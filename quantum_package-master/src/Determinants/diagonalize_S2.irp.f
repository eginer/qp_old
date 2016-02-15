program diag_s2
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
 implicit none
 use bitmasks
 integer :: i,j,k,l
 
 call diagonalize_s2_betweenstates(psi_det,psi_coef,N_det,size(psi_coef,1),n_states_diag)
 touch psi_coef
 call save_wavefunction
 




end

