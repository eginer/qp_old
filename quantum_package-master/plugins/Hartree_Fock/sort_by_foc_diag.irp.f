program sort_by_fock_diag
 implicit none
 integer :: i,j
 character*(64) :: label
 label = "Orthonormalized"
 do i = 1, mo_tot_num
  print*,'<F> = ',Fock_matrix_diag_mo(i),i
 enddo
 call mo_sort_by_observable(Fock_matrix_diag_mo,mo_tot_num,mo_tot_num,label)
 mo_label = 'Orthonormalized'
 SOFT_TOUCH mo_coef mo_label

 call save_mos

end

