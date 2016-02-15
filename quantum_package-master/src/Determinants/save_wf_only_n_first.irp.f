program save_wf_only_n_first
 read_wf = .True.
 touch read_wf
 call routine
end

subroutine routine
 implicit none
 use bitmasks
 integer ::i,j
 integer :: n_det_restart
 integer(bit_kind),allocatable :: psi_det_tmp(:,:,:)
 double precision ,allocatable :: psi_coef_tmp(:,:),accu(:)
 write(*,*) 'How many determinants would you like ?'
 read(*,*) n_det_restart
 write(*,*) 'Excellent choice !'
 allocate (psi_det_tmp(N_int,2,N_det_restart),psi_coef_tmp(N_det_restart,N_states),accu(N_states))
 accu = 0.d0
 do i = 1, N_det_restart
  do j = 1, N_int
   psi_det_tmp(j,1,i) = psi_det(j,1,i)
   psi_det_tmp(j,2,i) = psi_det(j,2,i)
  enddo
  do j = 1,N_states
   psi_coef_tmp(i,j) = psi_coef(i,j)
   accu(j) += psi_coef_tmp(i,j) * psi_coef_tmp(i,j)
  enddo
 enddo
 do j = 1, N_states
  accu(j) = 1.d0/dsqrt(accu(j))
 enddo
 do j = 1,N_states
  do i = 1, N_det_restart
   psi_coef_tmp(i,j) = psi_coef_tmp(i,j) * accu(j)
  enddo
 enddo
 call save_wavefunction_general(N_det_restart,N_states,psi_det_tmp,N_det_restart,psi_coef_tmp)

 deallocate (psi_det_tmp,psi_coef_tmp,accu)


end
