program reorder_mos
 implicit none
 integer :: i,j,k,l
 print*,'Enter the number of couple of orbital changes'
 integer :: n_change
 read(5,*) n_change
 print*,'Enter all the couple of orbital changes '
 integer,allocatable :: n_index_couple(:,:)
 allocate (n_index_couple(n_change,2))
 do i = 1, n_change
  read(5,*)n_index_couple(i,1),n_index_couple(i,2)
 enddo
 call change_order_mos(n_change,n_index_couple)
 soft_touch mo_coef
 call save_mos

end

