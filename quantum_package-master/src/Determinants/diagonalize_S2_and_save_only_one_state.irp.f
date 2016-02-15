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
 
 print*,'PROGRAM TO EXTRACT ONLY THE VARIOUS STATES CORRESPONDING TO '
 print*,'expected_s2 = ',expected_s2
 call diagonalize_s2_betweenstates(psi_det,psi_coef,N_det,size(psi_coef,1),n_states_diag)

 integer :: i_good_state(1000),n_good_state
 n_good_state = 0
 double precision :: s2_scalar
 do i = 1, n_states_diag
  call get_s2_u0(psi_det,psi_coef(1,i),n_det,size(psi_coef,1),s2_scalar)
  if(dabs(s2_scalar-expected_s2).le.0.3d0)then
   n_good_state +=1
   i_good_state(n_good_state) = i
  endif
 enddo
 print*,'n_good_state = ',n_good_state
 double precision, allocatable :: psi_coefs_tmp(:,:),s2(:,:)
 allocate(psi_coefs_tmp(n_det,n_good_state),s2(n_good_state,n_good_state))
 do j = 1, n_good_state
  do i = 1, n_det
   psi_coefs_tmp(i,j) = psi_coef(i,i_good_state(j))
  enddo
 enddo

 call get_uJ_s2_uI(psi_det,psi_coefs_tmp,N_det,size(psi_coefs_tmp,1),size(psi_det,1),s2,n_good_state)
 print*,'S^2 matrix in the basis of the states considered'
 do i = 1, n_good_state
  write(*,'(10(F16.10,X))')s2(i,:)
 enddo
 
 call save_wavefunction_general(n_det,n_good_state,psi_det,n_det,psi_coefs_tmp)
 




end

