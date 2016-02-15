program put_gess_always

 N_states = 2
 N_det= 4
 touch N_det N_states 
 call routine

end
subroutine routine
 use bitmasks
 implicit none
 integer :: i,j,N_det_tmp,N_states_tmp
 integer :: list(N_int*bit_kind_size,2)
 integer(bit_kind) :: string(N_int,2)
 integer(bit_kind) :: psi_det_tmp(N_int,2,4)
 double precision :: psi_coef_tmp(4,2)
 print*,'---------------'
 print*,'---------------'
 print*,'---------------'
 print*,'---------------'
 print*,'N_states = ',N_states

 integer :: iorb,jorb
 iorb = 6
 jorb = 7


 list = 0
 do i = 1, elec_alpha_num -1 
  list(i,1) = i
 enddo
 do i = 1, elec_beta_num -1 
  list(i,2) = i
 enddo
 list(elec_alpha_num,1) = iorb 
 list(elec_alpha_num,2) = jorb
 call list_to_bitstring( string(1,1), list(1,1), elec_alpha_num, N_int)
 call list_to_bitstring( string(1,2), list(1,2), elec_beta_num, N_int)
 call print_det(string,N_int)
 do j = 1,2
  do i = 1,  N_int 
   psi_det_tmp(i,j,1) = string(i,j)
  enddo
 enddo

 list = 0
 do i = 1, elec_alpha_num -1 
  list(i,1) = i
 enddo
 do i = 1, elec_beta_num -1 
  list(i,2) = i
 enddo
 list(elec_alpha_num,1) = jorb 
 list(elec_alpha_num,2) = iorb
 call list_to_bitstring( string(1,1), list(1,1), elec_alpha_num, N_int)
 call list_to_bitstring( string(1,2), list(1,2), elec_beta_num, N_int)
 call print_det(string,N_int)
 do j = 1,2
  do i = 1,  N_int 
   psi_det_tmp(i,j,2) = string(i,j)
  enddo
 enddo

 list = 0
 do i = 1, elec_alpha_num -1 
  list(i,1) = i
 enddo
 do i = 1, elec_beta_num -1 
  list(i,2) = i
 enddo
 list(elec_alpha_num,1) = jorb 
 list(elec_alpha_num,2) = jorb
 call list_to_bitstring( string(1,1), list(1,1), elec_alpha_num, N_int)
 call list_to_bitstring( string(1,2), list(1,2), elec_beta_num, N_int)
 call print_det(string,N_int)
 do j = 1,2
  do i = 1,  N_int 
   psi_det_tmp(i,j,3) = string(i,j)
  enddo
 enddo


 list = 0
 do i = 1, elec_alpha_num -1 
  list(i,1) = i
 enddo
 do i = 1, elec_beta_num -1 
  list(i,2) = i
 enddo
 list(elec_alpha_num,1) = iorb 
 list(elec_alpha_num,2) = iorb
 call list_to_bitstring( string(1,1), list(1,1), elec_alpha_num, N_int)
 call list_to_bitstring( string(1,2), list(1,2), elec_beta_num, N_int)
 call print_det(string,N_int)
 do j = 1,2
  do i = 1,  N_int 
   psi_det_tmp(i,j,4) = string(i,j)
  enddo
 enddo
 
 psi_coef_tmp(1,1) = 1.d0/dsqrt(2.d0)
 psi_coef_tmp(2,1) = 1.d0/dsqrt(2.d0)
 psi_coef_tmp(3,1) = 0.d0
 psi_coef_tmp(4,1) = 0.d0
 psi_coef_tmp(1,1) = 1.d0/dsqrt(2.d0)
 psi_coef_tmp(2,1) =-1.d0/dsqrt(2.d0)
 psi_coef_tmp(3,1) = 0.d0
 psi_coef_tmp(4,1) = 0.d0

 call save_wavefunction_general(n_det,n_states,psi_det_tmp,n_det,psi_coef_tmp)

end
