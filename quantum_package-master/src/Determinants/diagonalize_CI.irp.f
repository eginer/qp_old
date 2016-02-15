  use bitmasks
BEGIN_PROVIDER [ character*(64), diag_algorithm ]
  implicit none
  BEGIN_DOC
  ! Diagonalization algorithm (Davidson or Lapack)
  END_DOC
  if (N_det > N_det_max_jacobi) then
    diag_algorithm = "Davidson"
  else
    diag_algorithm = "Lapack"
  endif

  if (N_det < N_states_diag) then
    diag_algorithm = "Lapack"
  endif
  print*,'diag_algorithm = ',diag_algorithm
  print*,'N_det = ',N_det
  print*,'N_det_max_jacobi = ',N_det_max_jacobi
  print*,'n_states_diag    = ',n_states_diag
  
END_PROVIDER

BEGIN_PROVIDER [ double precision, CI_energy, (N_states_diag) ]
  implicit none
  BEGIN_DOC
  ! N_states lowest eigenvalues of the CI matrix
  END_DOC
  
  integer                        :: j
  character*(8)                  :: st
  call write_time(output_determinants)
  do j=1,N_states_diag
    CI_energy(j) = CI_electronic_energy(j) + nuclear_repulsion
    write(st,'(I4)') j
    call write_double(output_determinants,CI_energy(j),'Energy of state '//trim(st))
    call write_double(output_determinants,CI_eigenvectors_s2(j),'S^2 of state '//trim(st))
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, CI_electronic_energy, (N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors, (N_det,N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors_s2, (N_states_diag) ]
  implicit none
  BEGIN_DOC
  ! Eigenvectors/values of the CI matrix
  END_DOC
  integer                        :: i,j,k

  integer,allocatable :: list_order_good_state(:)
  integer,allocatable :: list_order_not_good_state(:)
  double precision, allocatable :: e_good_state(:)
  integer :: n_good_state, n_not_good_state
  double precision :: e_0
  integer(bit_kind), allocatable :: keys_tmp(:,:,:)
  allocate(keys_tmp(N_int,2,N_det))
  do i = 1, N_det
   do j = 1, N_int
    keys_tmp(j,1,i) = psi_det(j,1,i)
    keys_tmp(j,2,i) = psi_det(j,2,i)
   enddo
  enddo
  
  do j=1,N_states_diag
    do i=1,N_det
      CI_eigenvectors(i,j) = psi_coef(i,j)
    enddo
  enddo
  
  if (diag_algorithm == "Davidson") then
    
    call davidson_diag(keys_tmp,CI_eigenvectors,CI_electronic_energy, &
        size(CI_eigenvectors,1),N_det,N_states_diag,N_int,output_determinants)
    do j=1,N_states_diag
      call get_s2_u0(keys_tmp,CI_eigenvectors(1,j),N_det,size(CI_eigenvectors,1),CI_eigenvectors_s2(j))
    enddo

    
  else if (diag_algorithm == "Lapack") then
    
    double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:)
    double precision, allocatable  :: s2_eigvalues(:)
    allocate (eigenvectors(size(H_matrix_all_dets,1),N_det))
    allocate (eigenvalues(N_det))
    call lapack_diag(eigenvalues,eigenvectors,                       &
        H_matrix_all_dets,size(H_matrix_all_dets,1),N_det)
    
    do j = 1, N_states_diag
     CI_electronic_energy(j) = eigenvalues(j)
     do i=1,N_det
        CI_eigenvectors(i,j) = eigenvectors(i,j)
     enddo
    enddo
    integer :: i_state
    double precision :: s2
    if (s2_eig) then
      i_state = 0
      print*,'expected_s2= ',expected_s2
      do j=1,N_det
        call get_s2_u0(keys_tmp,eigenvectors(1,j),N_det,size(eigenvectors,1),s2)
        print*,'s2 = ',s2
        if(dabs(s2-expected_s2).le.0.3d0)then
        i_state += 1
        print*,'i_state = ',i_state
        do i=1,N_det
          CI_eigenvectors(i,i_state) = eigenvectors(i,j)
        enddo
        CI_electronic_energy(i_state) = eigenvalues(j)
        CI_eigenvectors_s2(i_state) = s2
        endif
        if (i_state.ge.N_states_diag) then
          exit
        endif
      enddo
      if(i_state == 0)then
       print*,'We did not found any interesting vectors ...'
       if(.not.diagonalize_s2)then
        print*,'stopping, sorry ... '
       else
        print*,'Diagonalizing the S^2 matrix within the available states !'
        allocate(list_order_good_state(0:n_det),list_order_not_good_state(0:n_det),s2_eigvalues(n_det))
        list_order_good_state = 0
        list_order_not_good_state = 0
        n_good_state = 0
        n_not_good_state = 0
        call diagonalize_s2_betweenstates(keys_tmp,eigenvectors,N_det,N_det,size(eigenvectors,1),N_det,s2_eigvalues)
        do j=1,N_det
          call get_s2_u0(keys_tmp,eigenvectors(1,j),N_det,size(eigenvectors,1),s2)
          if(s2_eig)then
              if(dabs(s2-expected_s2).le.0.3d0)then
               n_good_state +=1
               list_order_good_state(n_good_state) = j
              else  
               n_not_good_state +=1
               list_order_not_good_state(n_not_good_state) = j
              endif
          endif
        enddo
        print*,'n_good_state = ',n_good_state
        print*,'n_not_good_state = ',n_not_good_state
        if(n_good_state== 0)then
         print*,'We did not found any interesting vectors ...'
         print*,'stopping, sorry ... '
         stop
        endif
        do j = 1, n_good_state
         do k = 1, n_det
          CI_eigenvectors(k,j) = eigenvectors(k,list_order_good_state(j))
         enddo
        enddo
        i = 0
        do j = n_good_state+1, n_states_diag
         i+=1
         do k = 1, n_det
          CI_eigenvectors(k,j) = eigenvectors(k,list_order_not_good_state(i))
         enddo
        enddo
        deallocate(list_order_good_state,list_order_not_good_state,s2_eigvalues)
       endif
      endif
    else
      do j=1,N_states_diag
        call get_s2_u0(keys_tmp,eigenvectors(1,j),N_det,N_det,s2)
        do i=1,N_det
          CI_eigenvectors(i,j) = eigenvectors(i,j)
        enddo
        CI_electronic_energy(j) = eigenvalues(j)
        CI_eigenvectors_s2(j) = s2
      enddo
    endif
    deallocate(eigenvectors,eigenvalues)
  endif
! print*,'n_det.ge.n_states_diag = ',n_det.ge.n_states_diag
  if(diagonalize_s2.and.n_states_diag>1.and.(n_det.ge.n_states_diag))then
   allocate(s2_eigvalues(n_states_diag))
   if(s2_eig)then
    allocate(list_order_good_state(0:n_states_diag),list_order_not_good_state(0:n_states_diag))
    list_order_good_state = 0
    list_order_not_good_state = 0
    n_good_state = 0
    n_not_good_state = 0
   endif
   call diagonalize_s2_betweenstates(keys_tmp,CI_eigenvectors,N_det,N_det,size(CI_eigenvectors,1),n_states_diag,s2_eigvalues)
   do j=1,N_states_diag
     call get_s2_u0(keys_tmp,CI_eigenvectors(1,j),N_det,size(CI_eigenvectors,1),CI_eigenvectors_s2(j))
     if(s2_eig)then
         if(dabs(CI_eigenvectors_s2(j)-expected_s2).le.0.3d0)then
          n_good_state +=1
          list_order_good_state(n_good_state) = j
         else  
          n_not_good_state +=1
          list_order_not_good_state(n_not_good_state) = j
         endif
     else 
      call u0_H_u0(CI_electronic_energy(j),CI_eigenvectors(1,j),N_det,size(CI_eigenvectors,1),keys_tmp,N_int)
     endif
   enddo
   if(s2_eig)then
    print*,'n_good_state = ',n_good_state
    do i = 1, n_good_state
!!   print*,'list_order_good_state = ',list_order_good_state(i)
    enddo
    print*,'n_not_good_state = ',n_not_good_state
    do i = 1, n_not_good_state
!!   print*,'list_order_not_good_state = ',list_order_not_good_state(i)
    enddo
   endif
   
   if(s2_eig)then
    double precision, allocatable :: psi_coefs_tmp(:,:)
    integer, allocatable :: iorder(:),index_good_state_energetic_ordered(:)
    allocate(e_good_state(n_good_state),index_good_state_energetic_ordered(n_good_state))
    allocate(psi_coefs_tmp(N_det,n_states_diag),iorder(n_good_state))
!   print*,'good_state'
    do i = 1, n_good_state
     iorder(i) = i
     call u0_H_u0(e_0,CI_eigenvectors(1,list_order_good_state(i)),N_det,size(CI_eigenvectors,1),keys_tmp,N_int)
!    print*,'e_0 = ',e_0
     e_good_state(i) = e_0
    enddo
    do i = 1, N_states_diag
     do j = 1, N_det
      psi_coefs_tmp(j,i) = CI_eigenvectors(j,i)
     enddo
    enddo
!   print*,''

!   print*,'Sorting the good states '
    call dsort(e_good_state,iorder,n_good_state)   ! sort the good states by energy
    do i = 1, n_good_state
     index_good_state_energetic_ordered(i) = list_order_good_state(iorder(i))
    enddo
    do i = 1, n_good_state
!    print*,'index_good_state_energetic_ordered(i) = ',index_good_state_energetic_ordered(i)
     CI_electronic_energy(i) = e_good_state(i)
     do j = 1, N_det
      CI_eigenvectors(j,i) = psi_coefs_tmp(j,index_good_state_energetic_ordered(i))
     enddo
     call u0_H_u0(e_0,CI_eigenvectors(1,i),N_det,size(CI_eigenvectors,1),keys_tmp,N_int)
     call get_s2_u0(keys_tmp,CI_eigenvectors(1,i),N_det,size(CI_eigenvectors,1),s2)
!    print*,'e_0 = ',e_0
!    print*,'s2 = ',s2
     CI_eigenvectors_s2(i) = s2
    enddo
!   print*,''
!   print*,'not_good_state'
    do i = 1, n_not_good_state
!    print*,'list_order_not_good_state(i) = ',list_order_not_good_state(i)
     do j = 1, N_det
      CI_eigenvectors(j,i+n_good_state) = psi_coefs_tmp(j,list_order_not_good_state(i))
     enddo
     call get_s2_u0(keys_tmp,CI_eigenvectors(1,i+n_good_state),N_det,size(CI_eigenvectors,1),s2)
     call u0_H_u0(e_0,CI_eigenvectors(1,i+n_good_state),N_det,size(CI_eigenvectors,1),keys_tmp,N_int)
     CI_electronic_energy(i+n_good_state) = e_0
     CI_eigenvectors_s2(i+n_good_state) = s2
    enddo
    deallocate(psi_coefs_tmp,iorder,index_good_state_energetic_ordered)
    deallocate(list_order_good_state,list_order_not_good_state)
   endif
   deallocate(s2_eigvalues)
  endif
  deallocate(keys_tmp)
!  print*,''
!  print*,''
!  print*,''

! print*,'Checking the states'
! do i = 1, N_states_diag
!  call u0_H_u0(e_0,CI_eigenvectors(1,i),N_det,size(CI_eigenvectors,1),psi_det,N_int)
!  print*,'e_0 = ',e_0
!  print*,'CI_electronic_energy = ',CI_electronic_energy(i)
!  call get_s2_u0(psi_det,CI_eigenvectors(1,i),N_det,size(CI_eigenvectors,1),s2)
!  print*,'s2 = ',s2
!  print*,'CI_eigenvectors_s2(i) = ',CI_eigenvectors_s2(i)
! enddo
  
END_PROVIDER

subroutine diagonalize_CI
  implicit none
  BEGIN_DOC
!  Replace the coefficients of the CI states by the coefficients of the 
!  eigenstates of the CI matrix
  END_DOC
  integer :: i,j
  do j=1,N_states_diag
    do i=1,N_det
      psi_coef(i,j) = CI_eigenvectors(i,j)
    enddo
  enddo
  SOFT_TOUCH psi_coef CI_electronic_energy CI_energy CI_eigenvectors CI_eigenvectors_s2
end
