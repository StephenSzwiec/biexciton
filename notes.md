# Notes from study of the codebase 

Personal notes for this codebase study. 

## Input Files Audit 

The following is an explanation and 'audit' of the current input files for this calculation. 

- `WAVECAR`: VASP output, a binary file containing wavefunction information. Essential single-source of truth for electronic structure data. 
- `energy_pop`: REDUNDANT file, contains [band energy population] data otherwise found in the WAVECAR's `eig()` array. Code should be refactored to avoid reliance on this file. 
- `input_overlap`: REDUNDANT file, contains three integers which should either be CLI arguments or read from a config file. Code should be refactored to avoid reliance on this file.
- `active_window`: REDUNDANT file. Consider the following:
```
136     # ddo - hole window size (bands below HO; free parameter)
90      # du  - electron window size (bands above LU; free parameter)
1112    # Nexcc - number of exciton states (redundant; should match excE count)
0.025   # ??? - (duplicate of gamma_b? what is this?)
0.025   # gamma_b - broadening parameter in eV; free parameter 
300     # ss - threshold scaling factor; free parameter 
```
- `excE` - Essential file, legitimate pre-computed data from BSE/TDDFT (?) Verify that these aren't just difference of eigenvalues from WAVECAR 
- `excwf_trunc` Essential file. Consider the structure:
```
1                    # exciton index
504                  # electron band index (j)
352                  # hole band index (i)  
-0.005556829...      # Re(\phi^k_ji)
0.005112732...       # Im(\phi^k_ji)
508                  # next e band
353                  # next h band
0.0                  # Re(...)
0.005212312...       # Im(...)
[...]
```

    - This is a _sparse representation_ of exciton wavefunction \phi^k_{ji}, where: k is the exciton state, j = exciton electron band index, i = exciton hole band index. The wavefunction is complex-valued, hence the real and imaginary parts. The sparse representation only lists non-zero coefficients, with the last two lines indicating the next (j,i) pair.
- `excwf_length`: REDUNDANT file. The structure:
```
1        # exciton index
1192     # number of non-zero entries
2        # next exciton
1175     # entries
[...]
```
    - contains data already inferable from `excwf_trunc`. Code should be refactored to avoid reliance on this file.
- `VCe1ehe2 & VCh1ehh2`: Essential files. These are massive (~103M and ~234M lines) pre-computed Coulomb matrix elements. These represent:
    - **VCe1ehe2**: V(e1, e, h, e2) - electron-electron interaction with hole spectator
    - **VCh1ehh2**: V(h1, e, h, h2) - hole-hole interaction with electron spectator
    - _Format_: `j l n k (Re, Im)` → Complex matrix element
    - These are expensive to compute (O(N^4) with integrals) and must be pre-computed.

### New idea for config management

- Current code relies on multiple redundant input files for configuration parameters. This is error-prone and hard to maintain.
- Propose a single `config` file to consolidate all parameters:
```
qq = 86             # initial exciton index 
gamma_b = 0.025     # broadening parameter in eV 
ss = 300            # threshold scaling factor 
ddo = 136           # bands above HO 
du = 90             # bands below LU
# consider that these 'free' parameters are redundant with respect to nbandmin/max 
# nbandmin = HO - N , so no need to input 
# nHO => inferrable from highest occ. 
# nbandmax = LU + M 
```

Other parameters, such as `OMP_NUM_THREADS` are system environment variables and could be read from the system environment during runtime. 

We could then take each parameter as command line arguments:
```
$ ./biexciton_tile  --wavecar /path/to/WAVECAR \
                    --excitons /path/to/exeC \
                    --coulomb-ee /path/to/VCe1ehe2 \
                    --coulomb-hh /path/to/VCh1ehh2 \
                    --config /path/to/config_file
``` 

So, the code would read the config file once at startup, parse the parameters, and use them throughout. This would eliminate the need for multiple redundant input files and make configuration much clearer.

## Code Structure Analysis 

The following is an analysis of the code at `tile_OMP_enRead.f` which is the main parallel drive for the calculation.

### Workflow Overview 

1. **Input Preparation** the system requires several input files:
    1. `WAVECAR` - VASP output file with wavefunction data
    2. `energy_pop` - band energy populations (redundant)
    3. `input_overlap` - overlap parameters (redundant)
    4. `active_window` - exciton window parameters (redundant)
    5. `excE` - pre-computed exciton energies
    6. `excwf_trunc` - sparse exciton wavefunctions
    7. `VCe1ehe2` & `VCh1ehh2` - Coulomb matrix elements 
2. **Parallel Tiling Setup**: The `tile_set.sh` script generates parallelized FORTRAN code:
    - Takes user input `max` (qq value) and `threads` (number of OpenMP threads, which is not good compared to `OMP_NUM_THREADS` env var) 
    - Divides the exciton loop range into two tiles for parallel execution:
        - Tile 1: upper 75% of exciton indices gets one thread 
        - Tile 2: lower 25% of exciton indices gets remaining threads (heavier workload) 
    - Concatenates code sections: `tile_script_top`, `tile_script_mid`, and `tile_script_bottom` to produce `tile_OMP_enRead.f` (the main parallel driver (barf inducing workflow))
3. **Main Computation Phase**: performed by the `tile_OMP_enRead.f` code:
    1. **Initialization and File Reads**: 
        - Read band window params from `input_overlap`
        - Read active window params from `active_window` 
        - Allocate energy arrays based on band range 
        - Read band energies from `energy_pop`
    2. **WAVECAR Processing**:
        - Open and read `WAVECAR` 
        - Read record length `irecl` and metadata 
        - Reopen and correct record length (vomit)
        - Read k-point, band, and unit-cell info 
        - Extract plane wave coefficients for Gamma-point only 
        - Build operator \nabla(r) via `find_r_dipmom` subroutine:
            - Calculate reciprocal lattice vectors 
            - Find cutoff sphere in k-space 
            - Generate (x,y,z) indices for plane waves within cutoff 
        - Store coefficients in array like `ac(npw, nbandmax-nbandmin+1)` 
        - Calculate and writes overlap integrals 
    3. **Exciton Data Loading**: 
        - Read exciton energies from `excE` 
        - Read truncated exciton wavefunctions from `excwf_trunc` 
        - Read wavefunction lengths from `excwf_length` 
        - Build index arrays:
            - `ccc(Nexcc, nLu+du+1, nLU+du+1)` maps (exciton, e, h) -> index 
            - `ce(Nexcc, (ddo+1)*(du+1))` electron band indices 
            - `ch(Nexcc, (ddo+1)*(du+1))` hole band indices
        - Store wavefunction coefficients in `psi(Nexcc, (du+1)*(ddo+1))` 
    4. **Coulomb Matrix Element Loading**:
        1. Read V(e1,e,h,e2) from `VCe1ehe2`:
            - Store in `VCeeh(nijp)` where `nijp = (du+1)^2 * (ddo+1) * (du+1)`
            - Build index array `cveeh(du+1, du+1, ddo+1, du+1)` mapping (e1,e2,h1,h2) -> index 
        2. Read V(h1,e,h,h2) from `VCh1ehh2`:
            - Store in `VCheh(nijp)` where `nijp = (ddo+1)^2 * (ddo+1) * (ddo+1)`
            - Build index array `cvheh(ddo+1, du+1, ddo+1, ddo+1)` mapping (h1,e,h2,h2) -> index
    5. **Biexciton Energy Calculation**:
        - For initial exciton state `qq` with energy `w = Eexc(qq)`:
            - OMP Parallel region over range of `ae` values:
                1. Energy conservation check:
                    - Calculate `dE = w - Eexc(ae) - Eexc(be)` 
                    - Calculate delta function broadening `deltaexc = gamma_b / (dE + gamma_b)^2`
                    - Skip if `deltaexc < 1/(gamma_b * ss)` 
                2. Double loop over exciton components:
                    - Loop over components of exciton `qq` indexed by `ii`
                    - Loop over components of exciton `be` indexed by `kk`
                3. Compute transition amplitude `Rh` \Theta_j part (electron spectator):
                    - Extract indices: `j=ce(qq,ii), i=ch(qq,ii), l=ce(be,kk), n=ch(be,kk)`
                    - Loop over intermediate states `k` (electron bands nLU to nLU+du)
                    - If state `(ae, k, i)` exists:
                        - calculate energy denominator and principal values 
                        - \Theta_l \Theta_{-n} part from `VCeeh` matrix elements: `Rh += VCeeh * PV * f * psi_a * conj(psi_c) * conj(psi_b)` 
                        - \Theta_{-l} \Theta_n part from `VCeeh` matrix elements: `Rh += conj(VCeeh) * PV * f * conj(psi_c) * psi_b` 
                    - \Theta_{-j} part (hole spectator):
                        - Swap indices: `i=ce(qq,ii), j=ch(qq,ii)` 
                        - Loop over intermediate states `k` (hole bands nHO-ddo to nHO)
                        - If state `(ae, j, k)` exists:
                            - calculate energy denominator and principal values 
                            - \Theta_l \Theta_{-n} part from `VCheh` matrix elements: `Rh += VCheh * PV * f * conj(psi_a) * psi_c * conj(psi_b)`
                            - \Theta_{-l} \Theta_n part from `VCheh` matrix elements: `Rh += conj(VCheh) * PV * f * psi_a * conj(psi_c) * psi_b`
                4. Accumulate contribution to biexciton energy:
                    - `RRh += |Rh|^2 * deltaexc`
                    - Write matrix element `M2cab_h_spec_w86` containing new line: `(qq, ae, be, |Rh|^2)` 
    6. **Output Generation**:
        - `M2cab_h_spec_w86`: biexciton transition matrix elements for all `(qq,ae,be)` combinations 
        - `R_h_spec_nb_w86`: Total rate `RRh` for the given `qq` exciton state and `w` energy 
        - `overlap` wavefunction normalization checks 
    
### Physical Interpretation 

This code performs the **Auger recombination rate calculation** for biexcitons from exciton->biexciton transitions mediated by Coulomb interactions, using 2nd order perturbation theory with this formula: 
R_h_spec = (2pi/hbar) * sum_{ae,be} |Rh|^2 * delta(E_exc(qq) - E_exc(ae) - E_exc(be))

Where:
- Rh is the transition amplitude for exciton qq to biexciton (ae, be) mediated by Coulomb interactions with hole spectator
- delta(...) is a broadened delta function enforcing energy conservation

This includes both the electron-spectator and hole-spectator channels for Rh, summing over all possible intermediate states and exciton components, and applies a **Lorentzian broadening scheme** for the delta function to account for finite lifetimes using `gamma_b` and `ss` parameters for energy conservation. 

### Weirdness Notes 

1. The tiling approach in `tile_set.sh` is convoluted and should be eliminated entirely. A more straightforward parallelization strategy using OpenMP pragmas directly in the main code would be preferable. 
2. The reliance on multiple redundant input files for configuration parameters is error-prone. Consolidating these into a single config file would improve usability and maintainability.
3. `WAVECAR` file reading is complex and clunky due to VASP binary format; "open twice" pattern is ugly.
4. Plane wave basis operations are currently buried in `find_r_dipmom`
    - `find_r_dipmom` is sort of a bad name; it builds PW basis, not dipole moments. 
    - Separate type for basis set management probably useful here 
    - Could be extended for actual dipole calculations later
    - On second look, this thing is complete cruft: it builds `npw` but that can be determined from WAVECAR header info directly. Please eliminate. 
5. You can eliminate `excwf_length` by counting during parse of `excwf_trunc`. 
6. Coulomb matrix element storage and lookup is inefficient; a type-safe access method could prevent index errors within nested loops. 
7. Core 200 line calculation loop is monolithic and could be refactored as:
```
!$OMP PARALLEL DO 
do ae = ae_start, ae_end 
    do be = 1, qq 
        call compute_trasition_amplitude(qq, ae, be, Rh, ...)
        if (deltaexc > threshold) then 
            RRh = RRh + abs(Rh)**2 * deltaexc
            call write_matrix_element(qq, ae, be, abs(Rh)**2) 
        end if 
    end do 
end do 
!$OMP END PARALLEL DO
```
8. OpenMP thread management could be abstracted and load balanced at runtime rather than code-gen time.

#### Current code issues

1. Double file I/O:  
```
! Wrong record length first read this is ugly  
 open(12,file="WAVECAR",status="old",form="unformatted", 
     c          access='direct',recl=1500000) 
 open(12,rec=1) rdum, rispin, rtag 
 irecl = rdum 
 close(12)
 open(12, ..., recl=rdum)
! instead, consider a single read with a helper function to parse the header:
call open_wavecar(filename, unit, header)
```
2. Array index gymnastics: current code has `j-nLU-1`,`n-(nHO-ddo)+1` everywhere. Store the offsets instead:
```
type :: index_helper
    integer :: nLU, nHO, ddo, du 
contains 
    procedure :: e_to_array => map_electron_index 
    procedure :: h_to_array => map_hole_index
end type index_helper 
! usage 
idx = helper%e_to_array(j) ! maps band index j to array index 
```
3. Repeated allocations without deallocation: 
``` 
! appears 5+ times: 
allocate(some_array(size1, size2))
do i = 1, N 
    read(unit, ...) some_array(i, :)
end do 
! make this a generic 
subroutine read_2d_array(unit, array, n1, n2) 
```
4. Magic numbers: 
```
! current 
if (deltaexc > 1/(gamma_b * ss)) then 
! better 
real, parameter :: delta_threshold = 1.0 / (gamma_b * ss) 
if (deltaexc > delta_threshold) then 
```
5. Duplicate Energy Check Logic: 
```
! Lines 750, 780, 810, all do: 
denom = Eaeki * Eaeki + gamma_b * gamma_b 
! Make a function 
function lorentzian_denominator(Eaeki, gamma_b) result(denom)
```

### Potential Improvements 

#### Key Areas for 1st-Pass Refactoring

1. **Command Line Argument Parsing**: for all parameters. Compile once, run many times. 
2. **Dynamic Memory Management:** with proper bounds checking to avoid overflows and allocatable arrays. 
3. **Error handling for file I/O:** to ensure robustness against missing/corrupted files. 
4. **Config File Consolidation:** to eliminate redundant input files, improving operation clarity. 
5. **Progress Reporting and Checkpointing:** to monitor long calculations and allow restarts. 
6. **Modular Structure:** break monolithic code into smaller subroutines and functions. 
7. **OpenMP Improvements:** use pragmas directly in code rather than code generation for better maintainability.
8. **Type Safety:** encapsulate complex data structures (e.g., Coulomb matrix elements) in derived types with accessor methods to prevent index errors. 
   
#### 2nd-Pass Enhancements
 
1. **Matrix Storage Optimization:** consider C-S-R storage for sparse matrices to reduce memory footprint, then pivot to HDF5 for matrix storage. 
