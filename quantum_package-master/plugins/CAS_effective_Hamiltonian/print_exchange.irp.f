program print
 read_wf = .True.
 touch read_wf
 provide  mo_bielec_integrals_in_map
 call provide_all_stuffs
end
subroutine provide_all_stuffs
 implicit none
!provide ref_hamiltonian_matrix dressing_ref_hamiltonian
 integer :: i,j
 double precision :: hij,hmono,hdouble
 double precision               :: get_mo_bielec_integral     
 i = 9 
 j = 10
 
 print*,'H mono  = ',mo_mono_elec_integral(i,j) 
 print*,'kin     = ',mo_kinetic_integral(i,j)
 print*,'pot     = ',mo_nucl_elec_integral(i,j)
 print*,'pot 1   = ',mo_nucl_elec_integral_per_atom(i,j,1)
 print*,'pot 2   = ',mo_nucl_elec_integral_per_atom(i,j,2)
 print*,'J       = ',get_mo_bielec_integral(i,j,i,i,mo_integrals_map)
 call i_H_j_verbose(psi_ref(1,1,1),psi_ref(1,1,3),N_int,hij,hmono,hdouble)
 print*,'hij     = ',hij
 print*,'hmono   = ',hmono
 print*,'hdouble = ',hdouble


end
