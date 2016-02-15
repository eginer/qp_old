program print
 read_wf = .True.
 touch read_wf
 call provide_all_stuffs
end
subroutine provide_all_stuffs
 implicit none
 provide ref_hamiltonian_matrix 
 integer :: i,j,istate,k
 double precision, allocatable :: psi_restart_ref_normalized(:),psi_ref_zeroth_order(:),psi_ref_dressed(:)
 double precision, allocatable :: eigvalues(:),eigvectors(:,:)
 double precision, allocatable :: H_naked(:,:)
 double precision, allocatable :: H_dressed(:,:)
 double precision, allocatable :: H_print(:,:)
 double precision :: accu_norm
 allocate (H_dressed(max_number_ionic+1,max_number_ionic+1))
 allocate (H_print(min_number_ionic:max_number_ionic,min_number_ionic:max_number_ionic))
 allocate (H_naked(max_number_ionic+1,max_number_ionic+1))
 allocate (psi_restart_ref_normalized(min_number_ionic:max_number_ionic))
 allocate (psi_ref_zeroth_order(min_number_ionic:max_number_ionic))
 print*,'# nuclear_repulsion = ',nuclear_repulsion 
 allocate (psi_ref_dressed(min_number_ionic:max_number_ionic))
 allocate (eigvalues(max_number_ionic+1))
 allocate (eigvectors(max_number_ionic+1,max_number_ionic+1))
 double precision :: convert_hartree_ev
 convert_hartree_ev = 27.211399d0
 
 do istate= 1, N_states 
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
    write(*,'(I4,X,10(F16.10 ,4X))')i, H_print(i,:)*convert_hartree_ev
   enddo
  print*,''
  do i = min_number_ionic,max_number_ionic
   do j = min_number_ionic,max_number_ionic
    H_naked(j+1,i+1) = H_OVB_naked(i,j,istate)
   enddo
  enddo
 
  call lapack_diagd(eigvalues,eigvectors,H_naked,max_number_ionic+1,max_number_ionic+1)
  print*,'Energy = ',eigvalues(istate) + nuclear_repulsion
  do i = min_number_ionic,max_number_ionic
   psi_ref_zeroth_order(i) = eigvectors(i+1,istate)
  enddo
  print*,'Amplitudes of the various OVB components'
  print*,'Taking the neutral as reference '
  do i = min_number_ionic,max_number_ionic
   write(*,'(I4,X,10(F10.7 ,4X))') i,psi_ref_zeroth_order(i)/psi_ref_zeroth_order(min_number_ionic)
  enddo
  double precision,allocatable :: amplitudes(:)
  double precision, allocatable :: H_ovb_model_space(:,:)
  double precision, allocatable :: eigvalues_model_space(:),eigvectors_model_space(:,:)
  double precision, allocatable :: psi_coef_model_space(:)

  print*,'Choose the maximum level of ionicity you want in the MODEL space ...'
  integer :: max_level_ionic_model_space
  read(5,*) max_level_ionic_model_space

  print*,'Choose the maximum level of ionicity you want in the OUTER space ...'
  integer :: max_level_ionic_outer_space
  read(5,*) max_level_ionic_outer_space


  print*,'dressing the model space by the outer space !'
  print*,'Calculating the amplitudes'

  double precision :: accu
  allocate (amplitudes(1000))
  allocate (H_ovb_model_space(max_level_ionic_model_space+1,max_level_ionic_model_space+1))
  allocate (eigvalues_model_space(max_level_ionic_model_space+1))
  allocate (psi_coef_model_space(max_level_ionic_model_space+1))
  allocate (eigvectors_model_space(max_level_ionic_model_space+1,max_level_ionic_model_space+1))


  do i = max_level_ionic_model_space +1, max_level_ionic_outer_space
   accu = 0.d0
   do j = min_number_ionic,max_level_ionic_model_space
    accu += psi_ref_zeroth_order(j) * H_OVB_naked(i,j,istate)
   enddo
   amplitudes(i) = psi_ref_zeroth_order(i)/accu
   print*,'i = ',i
   print*,'amplitude = ',amplitudes(i)
  enddo
 
  do i = min_number_ionic,max_level_ionic_model_space 
   do j = min_number_ionic,max_level_ionic_model_space 
    H_ovb_model_space(i+1,j+1) = H_print(i,j)
    do k = max_level_ionic_model_space +1, max_level_ionic_outer_space
     H_ovb_model_space(i+1,j+1) += H_OVB_naked(i,k,istate) * H_OVB_naked(j,k,istate) * amplitudes(k)
    enddo
   enddo
  enddo
  
  print*,'Model space dressed by the outer space '
  do i = min_number_ionic+1,max_level_ionic_model_space+1
   write(*,'(I4,X,10(F16.10 ,4X))')i-1, H_ovb_model_space(i,:)*convert_hartree_ev
  enddo
 
  
  call lapack_diagd(eigvalues_model_space,eigvectors_model_space,H_ovb_model_space,max_level_ionic_model_space+1,max_level_ionic_model_space+1)
  do i = min_number_ionic,max_level_ionic_model_space 
   psi_coef_model_space(i+1) = eigvectors_model_space(i+1,istate)
  enddo
  
  print*,'Energy = ',eigvalues_model_space(istate) + nuclear_repulsion + H_OVB_naked(min_number_ionic,min_number_ionic,istate)
  print*,'Amplitudes of the various OVB components in the MODEL space '
  print*,'Taking the neutral as reference '
  do i = min_number_ionic,max_level_ionic_model_space
   write(*,'(I4,X,10(F16.10 ,4X))') i,psi_coef_model_space(i+1)/psi_coef_model_space(1)
  enddo
 enddo
 stop
 
 deallocate (H_dressed)
 deallocate (H_naked)
 deallocate (psi_restart_ref_normalized)
 deallocate (psi_ref_zeroth_order)
 deallocate (psi_ref_dressed)

 deallocate (eigvalues)
 deallocate (eigvectors)

end
