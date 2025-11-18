!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! biexciton_tile_read.f90 
!
! This program performs the Auger recombination rate calculation for biexcitons 
! from exciton->biexciton transitions mediated by Coulomb interactions, using 
! 2nd order perturbation theory, with this formula:
! 
! R_h_{spec} = (2pi / hbar ) 
!              * sum_{ae,be} [|Rh|^2 * delta(E_exc(qq) - E_exc(ae) - E_exc(be))]
! 
! where Rh is the transition amplitude for exciton qq to biexciton (ae, be)
! mediated by Coulomb interactions with hole spectator
! delta(...) is a broadened delta function enforcing energy conservation 
! 
! This includes both the electron and hole channel contributions to Rh using 
! spectator hole states (Tamm-Dancoff approximation ?) and applies a 
! Lorentzian broadening scheme to the delta function to account for 
! finite lifetimes of the states using two broadening parameters:
!   - gamma_b : broadening of exciton states
!   - ss : time threshold factor 
! 
!!! Input files: 
! - "WAVECAR" : contains KS wavefunctions, VASP binary, in plane-wave basis
! - "exeC" : contains pre-computed (BSE/TDDFT??) exciton energies, one per line
! - "excwf_trunc" : sparse representation of exciton wavefunctions(k,j,i):
! ```
! # exciton_index 
! # electron band_index j 
! # hold band index i 
! # Re(\phi^k_{ji})
! # Im(\phi^k_{ji})
! # next e band 
! # next h band 
! # Re(...)
! # Im(...) 
! [...]
! ```
! - "VCe1ehe2" and "VCh1ehh2" : electron-electron and hole-hole Coulomb matrices 
!   in sparse representation:
! ```
! # j l n k (Re, Im) 
! [...]
! ```
! - Config file "config.inp" or arguments:
!   - qq : exciton index for which to compute biexciton recombination rate
!   - gamma_b : broadening of exciton states (eV) 
!   - ss : time threshold scaling factor 
!   - ddo : bands above HO 
!   - du : bands below LU 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Module for configuration parsing and storage 
module mod_config
    implicit none 

    type :: biexciton_params 
        ! band structure parameters
        integer :: ddo      ! number of bands above HO to include 
        integer :: du       ! number of bands below LU to include 
        integer :: nHO      ! highest occupied band index 
        integer :: nLU      ! lowest unoccupied band index
        integer :: nbandmin ! minimum band index to consider 
        integer :: nbandmax ! maximum band index to consider 
        ! exciton parameters
        integer :: qq       ! exciton index for which to compute biexciton rate 
        real(8) :: gamma_b  ! exciton state broadening (eV) 
        integer :: ss       ! time threshold scaling factor 
        ! omp parameter 
        integer :: omp_threads ! number of omp threads
        ! file paths 
        character(len=*) :: wavecar_file 
        character(len=*) :: exciton_energy_file 
        character(len=*) :: exciton_wf_file 
        character(len=*) :: V_ee_file
        character(len=*) :: V_hh_file
        character(len=256) :: output_prefix
    end type biexciton_params

contains 
    subroutine read_config(filename, params)
        implicit none 
        character(len=*), intent(in) :: filename
        type(biexciton_params), intent(inout) :: params 
        integer :: ios 
        character(len=256) :: line, key, val 
        open(unit=101, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
    print *, 'Error: Unable to open config file: ', trim(filename)
            stop 
        end if
        do 
            read(101, '(A)', iostat=ios) line
            if (ios /= 0) exit 
            if (trim(line) == '' .or. line(1:1) == '#') cycle ! skip blank lines/comments
            read(line, '(A,A)', iostat=ios) key, val
            if (ios /= 0) cycle 
            select case (trim(adjustl(key)))
            case ('ddo')
                read(trim(adjustl(val)), *) params%ddo
            case ('du')
                read(trim(adjustl(val)), *) params%du
            case ('qq')
                read(trim(adjustl(val)), *) params%qq
            case ('gamma_b')
                read(trim(adjustl(val)), *) params%gamma_b
            case ('ss')
                read(trim(adjustl(val)), *) params%ss
            case ('omp_threads')
                read(trim(adjustl(val)), *) params%omp_threads
            end select
        end do
        close(101)
    end subroutine read_config 

    subroutine parse_commandline(params)
        implicit none 
        type(biexciton_params), intent(out) :: params 
        character(len=*) :: wavecar_file, exciton_energy_file 
        character(len=*) :: exciton_wf_file, V_ee_file, V_hh_file
        integer :: iarg, nargs 
        character(len=256) :: arg, command_name, path_buf
        logical :: exists

        nargs = command_argument_count()
        call get_command_argument(0, command_name)

        ! read flagged arguments 
        do iarg = 1, nargs 
            call get_command_argument(iarg, arg)
            select case (trim(arg))
            case ('--help', '-h')
                ! implement help message later 
            case ('--wavecar', '-w', '--WAVECAR')
                call get_command_argument(iarg+1, wavecar_file)
                params%wavecar_file = trim(wavecar_file)
            case ('--exciton_energy', '-e', '--exeC')
                call get_command_argument(iarg+1, exciton_energy_file)
                params%exciton_energy_file = trim(exciton_energy_file)
            case ('--exciton_wf', '-f', '--wf', '--excwf_trunc')
                call get_command_argument(iarg+1, exciton_wf_file)
                params%exciton_wf_file = trim(exciton_wf_file)
            case ('--V_ee', '--coulomb_ee', '--ee', '--VCe1ehe2')
                call get_command_argument(iarg+1, V_ee_file)
                params%V_ee_file = trim(V_ee_file)
            case ('--V_hh', '--coulomb_hh', '--hh', '--VCh1ehh2')
                call get_command_argument(iarg+1, V_hh_file)
                params%V_hh_file = trim(V_hh_file)
            case ('--output_prefix', '-o')
                call get_command_argument(iarg+1, path_buf)
                params%output_prefix = trim(path_buf)
            case ('-c', '--config')
                call get_command_argument(iarg+1, path_buf)
                call read_config(trim(path_buf), params) 
            case ('--ddo')
                call get_command_argument(iarg+1, arg)
                read(arg, *) params%ddo
            case ('--du')
                call get_command_argument(iarg+1, arg)
                read(arg, *) params%du
            case ('--qq')
                call get_command_argument(iarg+1, arg)
                read(arg, *) params%qq
            case ('--gamma_b')
                call get_command_argument(iarg+1, arg)
                read(arg, *) params%gamma_b
            case ('--ss')
                call get_command_argument(iarg+1, arg)
                read(arg, *) params%ss 
            case ('--omp_threads') 
                call get_command_argument(iarg+1, arg)
                read(arg, *) params%omp_threads
            end select
        end do 

        ! evaluate what we now have, and if that makes sense
        if (.not.allocated(params%wavecar_file)) params%wavecar_file = 'WAVECAR' 
        if (.not.allocated(params%exciton_energy_file)) params%exciton_energy_file = 'exeC'
        if (.not.allocated(params%exciton_wf_file)) params%exciton_wf_file = 'excwf_trunc'
        if (.not.allocated(params%V_ee_file)) params%V_ee_file = 'VCe1ehe2'
        if (.not.allocated(params%V_hh_file)) params%V_hh_file = 'VCh1ehh2'
        if (.not.allocated(params%output_prefix)) params%output_prefix = 'biexciton_output_'
        if (params%ddo <= 0) then
            print *, 'Error: ddo must be positive integer.'
            stop 
        end if
        if (params%du <= 0) then
            print *, 'Error: du must be positive integer.'
            stop 
        end if 
        if (params%qq < 1) then
            print *, 'Error: qq must be positive integer.'
            stop 
        end if
        if (params%gamma_b <= 0.0d0) then
            print *, 'Error: gamma_b must be positive real number.'
            stop 
        end if
        if (params%ss <= 0) then
            print *, 'Error: ss must be positive integer.'
            stop 
        end if
        if (.not.allocated(params%omp_threads)) then 
            ! try OMP_NUM_THREADS env variable 
            call get_environment_variable('OMP_NUM_THREADS', arg, status=iarg)
            ! if neither set, default to 1 
            if (iarg /= 0) then
                params%omp_threads = 1
            else
                read(arg, *) params%omp_threads
            end if 
        end if
        ! make sure all of these file paths actually exist 
        inquire(file=trim(params%wavecar_file), exist=exists)
        if (.not. exists) then
            print *, 'Error: WAVECAR file not found: ', trim(params%wavecar_file)
            stop 
        end if
        inquire(file=trim(params%exciton_energy_file), exist=exists)
        if (.not. exists) then
            print *, 'Error: Exciton energy file not found: ', trim(params%exciton_energy_file)
            stop 
        end if
        inquire(file=trim(params%exciton_wf_file), exist=exists)
        if (.not. exists) then
            print *, 'Error: Exciton wavefunction file not found: ', trim(params%exciton_wf_file)
            stop 
        end if
        inquire(file=trim(params%V_ee_file), exist=exists)
        if (.not. exists) then
            print *, 'Error: V_ee Coulomb matrix file not found: ', trim(params%V_ee_file)
            stop 
        end if
        inquire(file=trim(params%V_hh_file), exist=exists)
        if (.not. exists) then
            print *, 'Error: V_hh Coulomb matrix file not found: ', trim(params%V_hh_file)
            stop 
        end if
        if (params%omp_threads < 1) then
            print *, 'Error: omp_threads must be positive integer.'
            stop 
        end if
       
        ! cool, we should have everything now 
        ! because nbandmax = HO + ddo, nbandmin = LU - du, we set these later 
        ! final check: print out the parameters
        print *, 'Biexciton calculation parameters:'
        print *, '  WAVECAR file: ', trim(params%wavecar_file)
        print *, '  Exciton energy file: ', trim(params%exciton_energy_file)
        print *, '  Exciton wavefunction file: ', trim(params%exciton_wf_file)
        print *, '  V_ee file: ', trim(params%V_ee_file)
        print *, '  V_hh file: ', trim(params%V_hh_file)
        print *, '  Output prefix: ', trim(params%output_prefix)
        print *, '  ddo (bands above HO): ', params%ddo
        print *, '  du (bands below LU): ', params%du 
        print *, '  qq (exciton index): ', params%qq 
        print *, '  gamma_b (exciton broadening, eV): ', params%gamma_b 
        print *, '  ss (time threshold factor): ', params%ss 
        print *, '  omp_threads: ', params%omp_threads 
    end subroutine parse_commandline
end module mod_config

module mod_wavecar 
    implicit none 

    type :: wavecar_header 
        integer :: nkpt, nband, npw
        real(8) :: emax, efermi
        real(8) :: a1, b1, c1 ! lattice vectors
        integer :: rl         ! record length
        integer :: tag        ! precision tag 
    end type wavecar_header

    type :: kpoint_data
        real(8) :: kpt(3)                    ! k-point vector  
        real(8), allocatable :: eig(:)       ! eigenvalues
        real(8), allocatable :: occ(:)       ! occupations 
        complex(4), allocatable :: coef(:,:) ! psi coefficients (npw, nband) 
    end type kpoint_data

contains 
    subroutine read_wavecar(filename, header, kdata)
        implicit none
        character(len=*), intent(in) :: filename
        type(wavecar_header), intent(out) :: header
        type(kpoint_data), intent(out) :: kdata 
        integer :: ios, unit
        real(8) :: recordlen, rispin, rtag, rnpw
        real(8) :: inkpt, inband, lattice(3,3)
        integer :: i, j 
        real(8), allocatable :: gweight(:) ! not kept but read 
        real(8), parameter :: occ_threshold = 0.5d0

        ! Read record length 
        open(unit=102, file=trim(filename), form='unformatted', &
             access='direct', recl=8, status='old', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR: Unable to open WAVECAR file: ', trim(filename)
            stop 1
        end if
        read(102, rec=1) recordlen, rispin, rtag
        close(102)
        header%irecl = int(recordlen)
        header%tag = int(rtag)
        ! Warn about precision
        !if (header%rtag /= 45200 .and. header%rtag /= 45210) then
        !    print *, 'WARNING: Unexpected rtag value: ', header%rtag
        !    print *, 'Expected 45200 (single) or 45210 (double)'
        !end if
        open(unit=102, file=trim(filename), form='unformatted', &
             access='direct', recl=header%irecl, status='old', iostat=ios)
        read(102, rec=2) inkpt, inband, header%emax, &
                          ((lattice(i,j), i=1,3), j=1,3), header%efermi
        ! Store metadata
        header%nkpt = int(inkpt)
        header%nband = int(inband)
        ! Extract lattice parameters (assuming orthorhombic cell)
        header%a1 = lattice(1,1)
        header%b1 = lattice(2,2)
        header%c1 = lattice(3,3)
        ! allocate data structures 
        allocate(kdata%eig(header%nband))
        allocate(kdata%occ(header%nband))
        allocate(gweight(header%nkpt, header%nband)) ! not used 
        ! read k-point data 
        read(102, rec=3) rnpw, (kdata%kpt(i), i=1,3), &
                                  (kdata%eig(j),gweight(j),kdata%occ(j),&
                                  j=1,header%nband)
        header%npw = int(rnpw)
        allocate(kdata%coef(header%nband, header%npw))
        do i = 4, header%nband + 3
            read(102, rec=i) (kdata%coef(i-3,j), j=1, header%npw)
        end do
        close(102)
        ! deallocate gweight
        deallocate(gweight)
        ! report findings: 
        print *, 'WAVECAR read successfully:'
        print *, '  Number of k-points: ', header%nkpt
        print *, '  Number of bands: ', header%nband 
        print *, '  Max energy (eV): ', header%emax
        print *, '  Fermi energy (eV): ', header%efermi
        print *, '  Number of plane-waves: ', header%npw
        print *, '  Lattice parameters (Angstrom):'
        print *, '    a1: ', trim(header%a1), ' Å '
        print *, '    b1: ', trim(header%b1), ' Å '
        print *, '    c1: ', trim(header%c1), ' Å '
    end subroutine read_wavecar

    subroutine extract_active_bands(kdata, nbandmin, nbandmax, ac, overlaps)
        implicit none
        type(kpoint_data), intent(in) :: kdata 
        integer, intent(in) :: nbandmin, nbandmax 
        complex(8), allocatable, intent(out) :: ac(:,:) ! (nbandmax-nbandmin+1, npw)
        real(8), optional, intent(out) :: overlaps(:) ! (nbandmax-nbandmin+1) 
        integer :: npw, nactive, iband, ic, i 
        real(8) :: overlap 

        npw = size(kdata%coef, 2)
        nactive = nbandmax - nbandmin + 1 
        allocate(ac(nactive, npw))
        ic = 1 
        do iband = nbandmin, nbandmax 
            ac(ic, :) = cmplx(kdata%coef(iband, :), kind=8)
            ic = ic + 1 
        end do 
        if (present(overlaps)) then 
            allocate(overlaps(nactive))
            do ic = 1, nactive 
                overlap = 0.0d0 
                do i = 1, npw 
                    overlap = overlap + real(conjg(ac(ic,i)) * ac(ic,i))
                end do
                overlaps(ic) = overlap 
            end do 
        end if
    end subroutine extract_active_bands

    
    subroutine compute_ho_lu_indices(kdata, efermi, nho, nlu, egap, tolerance)
        implicit none
        type(kpoint_data), intent(in) :: kdata 
        real(8), intent(in) :: efermi      ! Fermi energy 
        integer, intent(out) :: nho        ! highest occupied band 
        integer, intent(out) :: nlu        ! lowest unoccupied band 
        real(8), intent(out) :: egap ! band gap 
        real(8), optional, intent(in) :: tolerance ! energy tolerance
        integer :: nband, iband 
        real(8) :: tol, E_ho, E_lu 
        integer, parameter :: default_tol = 0.01d0 
        
        if (present(tolerance)) then
            tol = tolerance
        else
            tol = default_tol
        end if

        nband = size(kdata%eig) 
        nho = -1
        nlu = -1 
        do iband = nband, 1, -1
            if (kdata%eig(iband) < efermi - tol) then 
                nho = iband
                E_ho = kdata%eig(iband)
                exit
            end if
        end do 

        do iband = 1, nband 
            if (kdata%eig(iband) > efermi + tol) then 
                nlu = iband 
                E_lu = kdata%eig(iband) 
                exit
            end if
        end do

        ! validate against errors 
        if(nho == -1 .or. nlu == -1) then 
            print *, 'Error: Unable to determine HO and LU band indices.'
            print *, 'E_Fermi (eV): ', efermi, ' eV'
            print *, 'Energy range: [', minval(kdata%eig), ', ', maxval(kdata%eig), '] eV'
            print *, 'Tolerance used (eV): ', tol, ' eV'
            stop 1 
        else if (nlu <= nho) then 
            print *, 'Error: Invalid band indices found: nHO=', nho, ', nLU=', nlu
            print *, 'Check Fermi energy and band energies (metallic system?)'
            print *, 'nHO = ', nho, ', E_HO = ', E_ho, ' eV'
            print *, 'nLU = ', nlu, ', E_LU = ', E_lu, ' eV'
            stop 1 
        end if

        ! double check versus occupations
        ho_check = -1 
        do i = nband, 1, -1 
            if (kdata%occ(i) > 0.5d0) then 
                ho_check = i 
                exit 
            end if
        end do
        if (ho_check /= nho) then 
            print *, 'Warning: Discrepancy in HO band index determination:'
            print *, '  From energies: nHO = ', nho
            print *, '  From occupations: nHO = ', ho_check
        end if

        ! degeneracy check ? I don't know if this is necessary but I guess it can't hurt
        n_degen = count(abs(kdata%eig - kdata%eig(nho)) < 1.0d-6)
        if (n_degen > 1) then 
            print *, 'Note: Highest occupied band is degenerate with ', n_degen, ' bands.'
        end if

        ! report findings 
        egap = E_lu - E_ho 
        print *, 'Determined band indices:'
        print *, '  Highest occupied band (nHO): ', nho, ' with energy (eV): ', E_ho 
        print *, '  Lowest unoccupied band (nLU): ', nlu, ' with energy (eV): ', E_lu
        print *, '  Band gap (eV): ', egap, ' eV '
        print *, '  E_Fermi (eV): ', efermi, ' eV (', efermi - E_ho, ' eV above HO)'
    end subroutine compute_ho_lu_indices


    subroutine free_kpoint_data(kdata)
        implicit none
        type(kpoint_data), intent(inout) :: kdata 
        if (allocated(kdata%eig)) deallocate(kdata%eig)
        if (allocated(kdata%occ)) deallocate(kdata%occ)
        if (allocated(kdata%coef)) deallocate(kdata%coef)
    end subroutine free_kpoint_data
end module mod_wavecar 

module mod_exciton
    implicit none 
    
    type :: exciton_basis 
        integer :: nexc ! number of exciton states
        real(8), allocatable :: exc(:) ! exciton energies (eV) exc(nexc)
        ! sparse representation of exciton wavefunctions 
        integer, allocatable :: n_components(:) ! nexc(k) number of (e,h) components for exciton k 
        integer, allocatable :: e_idx(:) ! flattened array of electron band indices ce(Nexc, max)
        integer, allocatable :: h_idx(:) ! flattened array of hole band indices ch(Nexc, max)
        complex(4), allocatable :: psi(:,:) ! wavefunctions 
        ! lookup map from ccc(exc, e, h) -> component index 
        integer, allocatable :: index_map(:,:,:) ! (nexc, nband_e, nband_h)
    end type exciton_basis

contains
    subroutine read_exciton_energies(filename, basis) 
        ! reads exciton energies from file, one per line 
        implicit none 
        character(len=*), intent(in) :: filename 
        type(exciton_basis), intent(inout) :: basis 
        integer :: ios, i 
        real(8) :: energy 
        ! how long is this file? 
        open(unit=103, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error: Unable to open exciton energy file: ', trim(filename)
            stop 
        end if
        basis%nexc = 0 
        do 
            read(103, *, iostat=ios) energy
            if (ios /= 0) exit 
            basis%nexc = basis%nexc + 1 
        end do 
        rewind(103) 
        allocate(basis%exc(basis%nexc))
        do i = 1, basis%nexc 
            read(103, *) basis%exc(i)
        end do 
        close(103)
        print *, 'Read ', basis%nexc, ' exciton energies from ', trim(filename)
    end subroutine read_exciton_energies
    
    
    subroutine read_exciton_wavefunctions(filename, basis, Nexcc)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: Nexcc
        type(exciton_basis), intent(inout) :: basis
        
        integer :: ios
        integer :: i, j, k, n_comp, current_k
        real(8) :: real_part, imag_part
        logical :: is_reading 
        
        ! Use linked-list style storage, then consolidate
        type :: comp_node
            integer :: j, i
            complex(4) :: psi
        end type comp_node
        
        type(comp_node), allocatable :: temp_data(:)
        integer :: temp_size, temp_capacity
        integer, allocatable :: exciton_start(:), exciton_count(:)
        
        temp_capacity = 10000  ! Initial capacity
        temp_size = 0

        allocate(temp_data(temp_capacity))
        allocate(exciton_start(Nexcc))
        allocate(exciton_count(Nexcc))
        
        exciton_count = 0
        is_reading = .false.
        current_k = 0
        
        open(unit=104, file=trim(filename), status='old', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR: Cannot open ', trim(filename)
            stop 1
        end if
        do
            read(104, *, iostat=ios) k
            if (ios /= 0) exit
            if (k >= 1 .and. k <= Nexcc) then
                ! New exciton
                current_k = k
                exciton_start(k) = temp_size + 1
                is_reading = .true.
            else if (is_reading) then
                ! Read quartet
                j = k
                read(104, *, iostat=ios) i, re_psi, im_psi
                if (ios /= 0) then
                    print *, 'ERROR: Incomplete quartet'
                    stop 1
                end if
                ! Grow temp_data if needed
                temp_size = temp_size + 1
                if (temp_size > temp_capacity) then
                    call grow_temp_data(temp_data, temp_capacity)
                end if
                temp_data(temp_size)%j = j
                temp_data(temp_size)%i = i
                temp_data(temp_size)%psi = cmplx(real_part, imag_part, kind=4)
                exciton_count(current_k) = exciton_count(current_k) + 1
            end if
        end do
        close(unit)
        
        ! Now consolidate into basis structure
        allocate(basis%n_components(Nexcc))
        basis%n_components = exciton_count
        
        ! Find max components
        max_comp = maxval(exciton_count)
        
        allocate(basis%e_idx(Nexcc, max_comp))
        allocate(basis%h_idx(Nexcc, max_comp))
        allocate(basis%psi(Nexcc, max_comp))
        
        basis%e_idx = -1
        basis%h_idx = -1
        basis%psi = cmplx(0.0, 0.0, kind=4)
        
        ! Copy data
        do k = 1, Nexcc
            do n_comp = 1, exciton_count(k)
                idx = exciton_start(k) + n_comp - 1
                basis%e_idx(k, n_comp) = temp_data(idx)%j
                basis%h_idx(k, n_comp) = temp_data(idx)%i
                basis%psi(k, n_comp) = temp_data(idx)%psi
            end do
        end do
        
        deallocate(temp_data, exciton_start, exciton_count)
        
        print *, 'Exciton wavefunctions loaded'

    contains
        subroutine grow_temp_data(arr, capacity)
            type(comp_node), allocatable, intent(inout) :: arr(:)
            integer, intent(inout) :: capacity
            type(comp_node), allocatable :: temp(:)
            
            allocate(temp(capacity * 2))
            temp(1:capacity) = arr
            call move_alloc(temp, arr)
            capacity = capacity * 2
        end subroutine grow_temp_data
    
    end subroutine read_exciton_wavefunctions 

    subroutine build_index_map(basis, nLU, du, nHO, ddo)
        implicit none 
        type(exciton_basis), intent(inout) :: basis 
        integer, intent(in) :: nLU, du, nHO, ddo 
        
        integer :: i, j, k, n 

        allocate(basis%index_map(basis%nexc, nLU:nLU+du, nHO-ddo:nHO))
        basis%index_map = 0 

        ! fill map 
        do k = 1, basis%nexc 
            do n = 1, basis%n_components(k) 
                i = basis%e_idx(k, n) ! electron band 
                j = basis%h_idx(k, n) ! hole band 
                basis%index_map(k, i, j) = n 
            end do 
        end do

        print *, 'Exciton index map built: ', &
            size(basis%index_map), ' entries (', &
            count(basis%index_map > 0), ' non-zero)'
    end subroutine build_index_map
end module mod_exciton


module mod_coulomb 
    implicit none 

    type :: columb_matrix 
        integer :: n_entries 
        complex(4), allocatable :: V(:)   ! matrix elements V(jlnk)
        integer, allocatable :: index_map(:,:,:,:) ! (j,l,n,k) -> index in V, 4D lookup 
        character(4) :: type ! 'eeeh' or 'hheh' 
    end type columb_matrix

contains 
    subroutine read_coulomb_matrix(filename, mat, nLU, nHO, du, ddo, mat_type)
        ! read Coulomb matrix from file in sparse format
        ! build index map for fast lookup 

        implicit none
        character(len=*), intent(in) :: filename
        type(columb_matrix), intent(out) :: mat 
        integer, intent(in) :: nLU, nHO, du, ddo 
        character(len=*), intent(in) :: mat_type ! 'eeeh' or 'hheh' 
        integer :: ios 
        integer :: j, l, n, k, idx 
        real(8) :: reV, imV

        ! input style per line 
        ! j l n k (Re,Im)  
        open(unit=105, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error: Unable to open Coulomb matrix file: ', trim(filename)
            stop 
        end if
        mat%n_entries = 0 
        ! First pass: count entries 
        do 
            read(105, *, iostat=ios) j, l, n, k, reV, imV
            if (ios /= 0) exit 
            mat%n_entries = mat%n_entries + 1 
        end do 
        rewind(105) 
        allocate(mat%V(mat%n_entries)) 
        allocate(mat%index_map(nLU:nLU+du, nLU:nLU+du, nHO-ddo:nHO, nHO-ddo:nHO))
        mat%index_map = -1 ! initialize to -1 (invalid index) 
        mat%type = mat_type
        idx = 1 
        do 
            read(105, *, iostat=ios) j, l, n, k, reV, imV
            if (ios /= 0) exit 
            mat%V(idx) = cmplx(reV, imV, kind=4)
            mat%index_map(j, l, n, k) = idx 
            idx = idx + 1 
        end do 
        close(105) 
        print *, 'Coulomb matrix ', trim(filename), ' read successfully:'
        print *, '  Number of entries: ', mat%n_entries 
        call check_hermiticity(mat)
    end subroutine read_coulomb_matrix

    function get_coulomb_element(mat, j, l, n, k) result(V) 
        ! safe indexed lookup with bounds checking
        implicit none 
        type(columb_matrix), intent(in) :: mat 
        integer, intent(in) :: j, l, n, k 
        complex(4) :: V 
        integer :: idx 

        if (j < lbound(mat%index_map,1) .or. j > ubound(mat%index_map,1) .or. &
            l < lbound(mat%index_map,2) .or. l > ubound(mat%index_map,2) .or. &
            n < lbound(mat%index_map,3) .or. n > ubound(mat%index_map,3) .or. &
            k < lbound(mat%index_map,4) .or. k > ubound(mat%index_map,4)) then 
            print *, 'Error: Coulomb matrix indices out of bounds: ', j, l, n, k
            stop 
        end if
        idx = mat%index_map(j, l, n, k) 
        if (idx == -1) then 
            V = cmplx(0.0d0, 0.0d0, kind=4) ! zero if not found 
        else 
            V = mat%V(idx)
        end if
    end function get_coulomb_element

    subroutine check_hermiticity(mat)
        ! check if Coulomb matrix is Hermitian [validate V(e2,h,e,e1) = conj(V(e1,h,e,e2))]

        implicit none 
        type(columb_matrix), intent(in) :: mat 
        integer :: j, l, n, k, idx1, idx2 
        complex(4) :: V1, V2 
        integer :: n_errors 
        n_errors = 0 

        do j = lbound(mat%index_map,1), ubound(mat%index_map,1)
            do l = lbound(mat%index_map,2), ubound(mat%index_map,2)
                do n = lbound(mat%index_map,3), ubound(mat%index_map,3)
                    do k = lbound(mat%index_map,4), ubound(mat%index_map,4)
                        idx1 = mat%index_map(j, l, n, k)
                        idx2 = mat%index_map(l, j, n, k)
                        if (idx1 /= -1 .and. idx2 /= -1) then 
                            V1 = mat%V(idx1)
                            V2 = mat%V(idx2)
                            if (abs(V1 - conjg(V2)) > 1.0d-6) then 
                                n_errors = n_errors + 1 
                                if (n_errors <= 10) then 
                                    print *, 'Hermiticity violation at indices: ', j, l, n, k
                                    print *, '  V(jlnk) = ', V1, ', conj(V(ljnk)) = ', conjg(V2)
                                end if
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end subroutine check_hermiticity
end module mod_coulomb

module mod_biexciton 
    use mod_exciton 
    use mod_coulomb 
    implicit none 

contains 
    subroutine compute_transition_amplitude(qq, ae, be, & 
            exc_basis, V_ee, V_hh, &
            energy_bands, gamma_b, ss, & 
            Rh, deltaexc) 
    end subroutine compute_transition_amplitude

    subroutine theta_j()
        ! Heaviside step function theta_j for electron spectator term
    end subroutine theta_j

    subroutine theta_minus_j()
        ! Heaviside step function theta_-j for hole spectator term 
    end subroutine theta_minus_j

    function energy_denominator(E1, E2, gamma_b) result(denom)
        ! Lorentzian energy denominator: gamma/(dE^2 + gamma^2) 
    end function energy_denominator 

end module mod_biexciton

module mod_parallel 
    use omp_lib 
    implicit none 

contains 
    subroutine setup_omp(nthreads)
        implicit none 
        integer, intent(in) :: nthreads 
        call omp_set_num_threads(nthreads)
        print *, 'OpenMP environment set up with ', nthreads, ' threads.'
    end subroutine setup_omp

    subroutine compute_load_balance(qq, nthreads, ranges)
        ! The "tiling" logic: upper 75% vs lower 25% of exciton states 
        ! where to split the work among threads 
        implicit none 
        integer, intent(in) :: qq, nthreads 
        integer, allocatable, intent(out) :: ranges(:,:) ! (nthreads, 2) start/end indices 
        integer :: i, total_work, work_per_thread, extra_work 
        integer :: start_idx, end_idx 
        ! we want to divide work, such that the first 75% of exciton states get 1 thread 
        ! and the remaining 25% get the rest of the threads 
        ! why? because lower exciton states are more computationally intensive 
        total_work = qq 
        ! skip the rest of nthreads is 1 or 2 
        if (nthreads == 1 ) then 
            allocate(ranges(1,2))
            ranges(1,1) = 1 
            ranges(1,2) = total_work 
            return 
        else if (nthreads == 2) then 
            allocate(ranges(2,2))
            ranges(1,1) = 1 
            ranges(1,2) = int(0.75d0 * total_work) 
            ranges(2,1) = ranges(1,2) + 1 
            ranges(2,2) = total_work 
            return 
        end if
        allocate(ranges(nthreads, 2)) 
        work_per_thread = total_work / nthreads 
        extra_work = mod(total_work, nthreads) 
        start_idx = 1 
        do i = 1, nthreads 
            end_idx = start_idx + work_per_thread - 1 
            if (i <= extra_work) end_idx = end_idx + 1 
            ranges(i, 1) = start_idx 
            ranges(i, 2) = end_idx 
            start_idx = end_idx + 1 
        end do 
    end subroutine compute_load_balance 

    subroutine parallel_biexciton_loop()
        ! wrapper for OMP section
    end subroutine parallel_biexciton_loop

end module mod_parallel

module mod_output
end module mod_output

program biexciton_calculation
    use mod_config
    use mod_wavecar
    use mod_exciton
    use mod_coulomb 
    use mod_biexciton
    use mod_parallel
    use mod_output
    implicit none 

    type(biexciton_params) :: params 
    type(wavecar_header) :: wavecar 
    type(exciton_basis) :: excitons 
    type(columb_matrix) :: V_ee, V_hh 
    type(kpoint_data) :: kdata
    complex(8), allocatable :: active_coefs(:,:) ! (nband_active, npw) 
    real(8), allocatable :: overlaps(:) ! (nband_active) 
    real(8) :: egap ! band gap (eV) 


    ! 1. Configuration and initialization 
    call parse_commandline(params) 

    ! 2. Read WAVECAR 
    call read_wavecar(params%wavecar_file, wavecar, kdata) 
    ! gamma point only for now 
    call compute_ho_lu_indices(kdata(1), wavecar%efermi, params%nHO, params%nLU, egap) 
    params%nbandmin = params%nLU - params%du 
    params%nbandmax = params%nHO + params%ddo 
    call extract_active_bands(kdata(1), params%nbandmin, params%nbandmax, &
                              active_coefs, overlaps)
    call free_kpoint_data(kdata)

    ! 3. Read exciton data 
    call read_exciton_energies(params%exciton_energy_file, excitons)
    call read_exciton_wavefunctions(params%exciton_wf_file, excitons, excitons%nexc)
    call build_index_map(excitions, params%nLU, params%du, params%nHO, params%ddo)
    
    ! 4. Read Coulomb matrices 
    call read_coulomb_matrix(params%V_ee_file, V_ee)
    call read_coulomb_matrix(params%V_hh_file, V_hh)

    ! 5. Setup parallel environment  
    call setup_omp(params%omp_threads) 
    call compute_load_balance(params%qq, params%omp_threads, ranges)

    ! 6. Compute biexciton recombination rate 
    call parallel_biexciton_loop(ranges, params, excitons, V_ee, V_hh, results)

    ! 7. Output results 
    call write_results(params, results)
    call write_summary(params, timing)

    ! Finalize
    call finalize_parallel()
end program biexciton_calculation
