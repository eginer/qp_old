program print_OVB
 implicit none
 read_wf = .True.
 call provide_all

end
 
subroutine provide_all
 implicit none
 integer :: i,j,k,l,istate

 double precision, allocatable :: eigvalues(:),eigvectors(:,:)
 double precision, allocatable :: H_naked(:,:)
 double precision, allocatable :: H_dressed(:,:)
 allocate (H_dressed(max_number_ionic+1,max_number_ionic+1))
 allocate (H_naked(max_number_ionic+1,max_number_ionic+1))
 allocate (eigvalues(max_number_ionic+1))
 allocate (eigvectors(max_number_ionic+1,max_number_ionic+1))
 call diagonalize_CI


 istate = 1
 do istate = 1, N_states
 print*,'--------------'
 print*,'ISTATE = ',istate
 print*,'# nuclear_repulsion = ',nuclear_repulsion

  write(*,'(A27,F13.6)')'Naked Neutral energy      = ',H_OVB_naked(0,0,istate) + nuclear_repulsion
  write(*,'(A27,F13.6)')'Naked Ionic energy        = ',H_OVB_naked(1,1,istate) + nuclear_repulsion
  write(*,'(A27,F13.6)')'Naked Interaction         = ',H_OVB_naked(0,1,istate)
  write(*,'(A27,F13.6)')''
  do i = min_number_ionic,max_number_ionic
   do j = min_number_ionic,max_number_ionic
    H_dressed(j+1,i+1) = H_OVB_total_dressed(i,j,istate)
    H_naked(j+1,i+1) = H_OVB_naked(i,j,istate)
   enddo
  enddo
  call lapack_diagd(eigvalues,eigvectors,H_naked,max_number_ionic+1,max_number_ionic+1)
  write(*,'(A27,F13.6)')'Total energy naked        = ',eigvalues(istate) + nuclear_repulsion
  print*,''                                        
  write(*,'(A27,F13.6)')'dressing Neutral energy   = ',H_OVB_dressing(0,0,istate)
  write(*,'(A27,F13.6)')'dressing Ionic energy     = ',H_OVB_dressing(1,1,istate)
  write(*,'(A27,F13.6)')'dressing Interaction      = ',H_OVB_dressing(0,1,istate)
  write(*,'(A27,F13.6)')''                         
  write(*,'(A27,F13.6)')'dressed  Neutral energy   = ',H_OVB_total_dressed(0,0,istate)+ nuclear_repulsion
  write(*,'(A27,F13.6)')'dressed  Ionic energy     = ',H_OVB_total_dressed(1,1,istate)+ nuclear_repulsion
  write(*,'(A27,F13.6)')'dressed  Interaction      = ',H_OVB_total_dressed(0,1,istate)
 
  call lapack_diagd(eigvalues,eigvectors,H_dressed,max_number_ionic+1,max_number_ionic+1)
  write(*,'(A27,F13.6)')'Total energy dressed      = ',eigvalues(istate) + nuclear_repulsion
  write(*,'(A27,F13.6)')''
  write(*,'(A27,F13.6)')'Total energy diagonalized = ',CI_energy(istate)
  
 enddo


end
