!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Biexciton Tile Read Program 
! 
! This program reads data from several data sources to construct 
! from various criteria and physical parameters a series of output
! which uses the lower 25% of the `qq` value selected
! for each tile in the biexciton calculation.
!
! As a refresher, consider the nuclear-electronic Hamiltonian:
! H = T_n + T_e + V_ne + V_ee + V_nn 
! 
! We can discount the ionic kinetic energy T_n and the nuclear-nuclear repulsion V_nn
! as they are constants for a given nuclear configuration. 
! The electronic kinetic energy T_e and the electron-electron repulsion V_ee
! are also not directly relevant to the biexciton states we are interested in.
! Thus, we focus on the electron-nuclear attraction V_ne which primarily determines
! the biexciton binding energies and states.
!
! So, for any given two electron states |i> and |j>, the biexciton state can be approximated as: 
! |B_ij> = |i> ⊗ |j>
! 
! The biexciton binding energy E_B can be estimated as:
! 
! E_B = E_i + E_j - V_ne_ij
! 
! where E_i and E_j are the single exciton energies 
! and V_ne_ij is the electron-nuclear attraction
! which is a binding energy term, between the two states. 
!
! H_{xx} \phi = E_{xx} \phi -> E_b = 2E_x - E_{xx}
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


