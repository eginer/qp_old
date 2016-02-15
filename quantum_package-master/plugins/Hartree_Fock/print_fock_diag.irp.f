program print_fock_diag
 implicit none
 integer :: i,j
 do i = 1, mo_tot_num
  print*,'<F> = ',Fock_matrix_diag_mo(i),i
 enddo

end

