program print_OVB
 implicit none
 read_wf = .True.
 call provide_all

end
 
subroutine provide_all
 implicit none
 integer :: i,j,k,l,istate
 double precision, allocatable :: H_print(:,:)
 double precision, allocatable :: H_dressed(:,:)
 print*,'# nuclear_repulsion = ',nuclear_repulsion
 allocate (H_print(min_number_ionic:max_number_ionic,min_number_ionic:max_number_ionic))
 allocate (H_dressed(min_number_ionic:max_number_ionic,min_number_ionic:max_number_ionic))
 do istate = 1, N_states
   print*,'ISTATE = ',istate
   do i = min_number_ionic,max_number_ionic
    do j = min_number_ionic,max_number_ionic
     H_print(i,j) = H_OVB_naked(j,i,istate)
    enddo
   enddo
   do i = min_number_ionic,max_number_ionic
    H_print(i,i) -= H_OVB_naked(min_number_ionic,min_number_ionic,istate)
   enddo
  
   print*,'Ref Hamiltonian matrix emelent = ',H_OVB_naked(min_number_ionic,min_number_ionic,istate)
   print*,'-------------------'
   print*,'-------------------'
   print*,'CAS MATRIX         '
   print*,''
   do i = min_number_ionic,max_number_ionic
    write(*,'(I4,X,10(F8.5 ,4X))')i, H_print(i,:)
   enddo
   print*,''
   print*,'-------------------'
   print*,'-------------------'
   print*,'CAS MATRIX DRESSING '
   print*,''
   do i = min_number_ionic,max_number_ionic
    write(*,'(I4,X,10(F8.5 ,4X))')i, H_OVB_dressing(i,:,istate)
   enddo
   print*,''
   print*,'-------------------'
   print*,'-------------------'
   print*,'CAS MATRIX DRESSED  '
   print*,''
   do i = min_number_ionic,max_number_ionic
    do j = min_number_ionic,max_number_ionic
     H_dressed(i,j) = H_OVB_total_dressed(j,i,istate)
    enddo
   enddo
   do i = min_number_ionic,max_number_ionic
    H_dressed(i,i) -= H_OVB_total_dressed(min_number_ionic,min_number_ionic,istate)
   enddo
   do i = min_number_ionic,max_number_ionic
    write(*,'(I4,X,10(F8.5 ,4X))')i, H_dressed(i,:)
   enddo
   print*,''
!  print*,'-------------------'
 enddo


 deallocate (H_dressed)
 deallocate (H_print)
end
