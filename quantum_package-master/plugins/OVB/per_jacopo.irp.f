program print_OVB
 implicit none
 read_wf = .True.
 call provide_all

end

subroutine provide_all
 implicit none
 integer :: i,j,k,l,istate
 istate = 1
 print*,'H matrix in the contracted forms'
 print*,'Ionicity level       H matrix'
 do i = min_number_ionic, max_number_ionic 
  write(*,'(I3,X,10(F12.7,X))')i,H_OVB_naked(i,:,istate)
 enddo
 double precision, allocatable :: Hmatrix_tmp(:,:),eigvectors(:,:),eigvalues(:)
 allocate(Hmatrix_tmp(max_number_ionic+1,max_number_ionic+1),eigvectors(max_number_ionic+1,max_number_ionic+1),eigvalues(max_number_ionic+1))
 do i = min_number_ionic, max_number_ionic 
  do j = min_number_ionic, max_number_ionic 
   Hmatrix_tmp(i+1,j+1) = H_OVB_naked(i,j,istate)
  enddo
 enddo
 
 call lapack_diagd(eigvalues,eigvectors,Hmatrix_tmp,max_number_ionic+1,max_number_ionic+1)
 print*,'E = ',eigvalues(1) + nuclear_repulsion

end
