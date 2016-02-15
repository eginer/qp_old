BEGIN_PROVIDER [ double precision, H_matrix_all_dets,(N_det,N_det) ]
 implicit none
 BEGIN_DOC
 ! H matrix on the basis of the slater determinants defined by psi_det
 END_DOC
 integer :: i,j,degree
 double precision :: hij
 print*,'Providing H_matrix_all_dets'
 call  i_H_j(psi_det(1,1,1),psi_det(1,1,1),N_int,hij)
 !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,j,hij) &
 !$OMP SHARED (N_det, psi_det, N_int,H_matrix_all_dets)
 do i =1,N_det
  do j =1,N_det
   call  i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
   H_matrix_all_dets(i,j) = hij
  enddo
 enddo
 !$OMP END PARALLEL DO
END_PROVIDER

