 BEGIN_PROVIDER [integer, pert_determinants, (N_states, psi_det_size) ]
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, lambda_mrcc, (N_states,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, lambda_pert, (N_states,psi_det_size) ] 
 implicit none
 BEGIN_DOC
 ! cm/<Psi_0|H|D_m> or perturbative 1/Delta_E(m)
 END_DOC
 integer :: i,k,j
 double precision               :: ihpsi(N_states), hii,delta_e_eff,ihpsi_current(N_states),hij
 integer :: i_ok,i_pert,i_pert_count
 i_ok = 0

 double precision :: phase_restart(N_states),tmp
 do i = 1, N_det_ref
  print*,'psi_ref_coef(1,k) = ',psi_ref_coef(i,1)
 enddo
 do k = 1, N_states
  phase_restart(k) = dsign(1.d0,psi_ref_coef_restart(1,k)/psi_ref_coef(1,k))
 enddo
 i_pert_count = 0
 
 do i=1,N_det_non_ref
   call i_h_psi(psi_non_ref(1,1,i), psi_ref_restart, psi_ref_coef_restart, N_int, N_det_ref,&
       size(psi_ref_coef_restart,1), n_states, ihpsi)
   call i_H_j(psi_non_ref(1,1,i),psi_non_ref(1,1,i),N_int,hii)
   do k=1,N_states
     lambda_pert(k,i) = 1.d0 / (psi_ref_energy_diagonalized(k)-hii)
     call i_h_psi(psi_non_ref(1,1,i), psi_ref, psi_ref_coef, N_int, N_det_ref,size(psi_ref_coef,1), n_states, ihpsi_current)
     tmp = psi_non_ref_coef(i,k)/ihpsi_current(k)
!    if(dabs(ihpsi_current(k)).ge.1.d-10)then
!     if(dabs(psi_non_ref_coef(i,k)) .le. 1.d-12 .or. dabs(ihpsi_current(k)) .le. 1.d-12)cycle
!     lambda_mrcc(k,i) = lambda_pert(k,i) * erf(tmp/lambda_pert(k,i)) + tmp * (1.d0 - erf(tmp/lambda_pert(k,i)))
!     if((ihpsi(k) * lambda_pert(k,i))/psi_non_ref_coef_restart(i,k) < 0.d0)then
!      lambda_mrcc(k,i) = lambda_pert(k,i)
!     endif
      i_pert = 1
      if((ihpsi(k) * lambda_pert(k,i))/psi_non_ref_coef_restart(i,k) .ge. 0.5d0 & 
         .and. (ihpsi(k) * lambda_pert(k,i))/psi_non_ref_coef_restart(i,k) > 0.d0 )then  ! test on the first order coefficient
       i_pert = 0
      endif
      do j = 1, N_det_ref
       call i_H_j(psi_non_ref(1,1,i),psi_ref(1,1,j),N_int,hij)
       if(dabs(hij * tmp).ge.0.5d0)then
        i_pert_count +=1
        i_pert = 1
        exit
       endif
      enddo
     if( i_pert == 1)then
      pert_determinants(k,i) = i_pert
     endif
     if(pert_determinants(k,i) == 1)then
       i_ok +=1
       lambda_mrcc(k,i) = lambda_pert(k,i)
     else
       lambda_mrcc(k,i) = psi_non_ref_coef(i,k)/ihpsi_current(k)
     endif

!    if(i == 185)then
!     print*,'det ',i,'not taken in perturbation'
!     print*,lambda_mrcc(k,i),psi_non_ref_coef(i,k),ihpsi_current(k)
!    endif
!    if(i == 187)then
!     print*,'det ',i,'not taken in perturbation'
!     print*,lambda_mrcc(k,i),psi_non_ref_coef(i,k),ihpsi_current(k)
!    endif
   enddo
 enddo
!if(oscillations)then
! print*,'AVERAGING the lambda_mrcc with those of the previous iterations'
! do i = 1, N_det_non_ref
!  do k = 1, N_states

!   double precision :: tmp
!   tmp = lambda_mrcc(k,i)
!   lambda_mrcc(k,i) += lambda_mrcc_tmp(k,i)
!   lambda_mrcc(k,i) = lambda_mrcc(k,i) * 0.5d0
!   if(dabs(tmp - lambda_mrcc(k,i)).ge.1.d-9)then
!   print*,''
!   print*,'i = ',i
!   print*,'psi_non_ref_coef(i,k) = ',psi_non_ref_coef(i,k)
!   print*,'lambda_mrcc(k,i)     = ',lambda_mrcc(k,i)
!   print*,'                 tmp = ',tmp
!   endif
!  enddo
! enddo
!endif
 print*,'N_det_non_ref = ',N_det_non_ref
 print*,'Number of Perturbatively treated determinants = ',i_ok
 print*,'i_pert_count = ',i_pert_count
 print*,'psi_coef_ref_ratio = ',psi_ref_coef(2,1)/psi_ref_coef(1,1)

END_PROVIDER


BEGIN_PROVIDER [ double precision, lambda_mrcc_tmp, (N_states,psi_det_size) ]
 implicit none
 lambda_mrcc_tmp = 0.d0
END_PROVIDER 

BEGIN_PROVIDER [ logical, oscillations ]
 implicit none
 oscillations = .False.
END_PROVIDER 



!BEGIN_PROVIDER [ double precision, delta_ij_non_ref, (N_det_non_ref, N_det_non_ref,N_states) ]
!implicit none
!BEGIN_DOC
!! Dressing matrix in SD basis
!END_DOC
!delta_ij_non_ref = 0.d0
!call H_apply_mrcc_simple(delta_ij_non_ref,N_det_non_ref)
!END_PROVIDER

 BEGIN_PROVIDER [ double precision, delta_ij, (N_det_ref,N_det_non_ref,N_states) ]
&BEGIN_PROVIDER [ double precision, delta_ii, (N_det_ref,N_states) ]
 implicit none
 BEGIN_DOC
 ! Dressing matrix in N_det basis
 END_DOC
 integer :: i,j,m
 delta_ij = 0.d0
 delta_ii = 0.d0
 print*,'Applying the T operator'
 call H_apply_mrcc(delta_ij,delta_ii,N_det_ref,N_det_non_ref)
 double precision :: max_delta
 double precision :: accu
 integer :: imax,jmax
 max_delta = 0.d0
 accu = 0.d0
 jmax = -1
 imax = -1
 do i = 1, N_det_ref
  do j = 1, N_det_non_ref
   accu += psi_non_ref_coef(j,1) * psi_ref_coef(i,1) * delta_ij(i,j,1)
   if(dabs(delta_ij(i,j,1)).gt.max_delta)then
    max_delta = dabs(delta_ij(i,j,1))
    imax = i
    jmax = j
   endif
  enddo
 enddo
 print*,''
 print*,''
 print*,'<psi| Delta H |psi> = ',accu
 if(imax .gt. 0 .and. jmax .gt. 0)then
  print*,'MAX VAL OF DRESING = ',delta_ij(imax,jmax,1)
  print*,'imax,jmax = ',imax,jmax
  print*,'psi_ref_coef(imax,1)     = ',psi_ref_coef(imax,1)
  print*,'psi_non_ref_coef(jmax,1) = ',psi_non_ref_coef(jmax,1)
 endif
 do i = 1, N_det_ref
  print*,'delta_ii(i,1)     = ',delta_ii(i,1)
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, h_matrix_dressed, (N_det,N_det,N_states) ]
 implicit none
 BEGIN_DOC
 ! Dressed H with Delta_ij
 END_DOC
 integer                        :: i, j,istate,ii,jj
 do istate = 1,N_states
   do j=1,N_det
     do i=1,N_det
       h_matrix_dressed(i,j,istate) = h_matrix_all_dets(i,j) 
     enddo
   enddo
   do ii = 1, N_det_ref
     i =idx_ref(ii)
     h_matrix_dressed(i,i,istate) += delta_ii(ii,istate)
    do jj = 1, N_det_non_ref
     j =idx_non_ref(jj)
     h_matrix_dressed(i,j,istate) += delta_ij(ii,jj,istate)
     h_matrix_dressed(j,i,istate) += delta_ij(ii,jj,istate)
    enddo
   enddo 
 enddo
!print*,'h_matrix_dressed'
!istate=1
!do i = 1, N_det
! write(*,'(100(F16.8,X))')h_matrix_dressed(i,:,istate)
!enddo
END_PROVIDER


 BEGIN_PROVIDER [ double precision, CI_electronic_energy_dressed, (N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors_dressed, (N_det,N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors_s2_dressed, (N_states_diag) ]
   implicit none
   BEGIN_DOC
   ! Eigenvectors/values of the CI matrix
   END_DOC
   integer                        :: i,j,k
   
   do j=1,N_states_diag
     do i=1,N_det
       CI_eigenvectors_dressed(i,j) = psi_coef(i,j)
     enddo
   enddo
   
   if (diag_algorithm == "Davidson") then
     
     integer                        :: istate
     istate = 1
     call davidson_diag_mrcc(psi_det,CI_eigenvectors_dressed,CI_electronic_energy_dressed,&
         size(CI_eigenvectors_dressed,1),N_det,N_states_diag,N_int,output_determinants,istate)
     
   else if (diag_algorithm == "Lapack") then
     
     double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:),H_mat_tmp(:,:)
    double precision, allocatable  :: s2_eigvalues(:)
     allocate (eigenvectors(N_det,N_det),H_mat_tmp(N_det,N_det))
     allocate (eigenvalues(N_det))
     do i = 1, N_det
      do j = 1,N_det
       H_mat_tmp(i,j) = H_matrix_dressed(i,j,1)
      enddo
     enddo
     call lapack_diag(eigenvalues,eigenvectors,H_mat_tmp,N_det,N_det)
     do j = 1, n_states_diag
      CI_electronic_energy_dressed(:) = eigenvalues(j)
      do i=1,N_det
        CI_eigenvectors_dressed(i,j) = eigenvectors(i,j)
      enddo
     enddo
     integer                        :: i_state
     double precision               :: s2
     i_state = 0
     if (s2_eig) then
       do j=1,N_det
         call get_s2_u0(psi_det,eigenvectors(1,j),N_det,N_det,s2)
!        print*,'s2 = ',s2
         if(dabs(s2-expected_s2).le.0.3d0)then
           i_state += 1
           do i=1,N_det
             CI_eigenvectors_dressed(i,i_state) = eigenvectors(i,j)
           enddo
           CI_electronic_energy_dressed(i_state) = eigenvalues(j)
           CI_eigenvectors_s2_dressed(i_state) = s2
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
        call diagonalize_s2_betweenstates(psi_det,eigenvectors,N_det,size(psi_det,1),size(eigenvectors,1),N_det,s2_eigvalues)
        do j=1,N_det
          call get_s2_u0(psi_det,eigenvectors(1,j),N_det,size(eigenvectors,1),s2)
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
         call get_s2_u0(psi_det,eigenvectors(1,j),N_det,N_det,s2)
         i_state += 1
         do i=1,N_det
           CI_eigenvectors_dressed(i,i_state) = eigenvectors(i,j)
         enddo
         CI_electronic_energy_dressed(i_state) = eigenvalues(j)
         CI_eigenvectors_s2_dressed(i_state) = s2
       enddo
     endif
     deallocate(eigenvectors,eigenvalues,H_mat_tmp)
   endif
  if(diagonalize_s2.and.n_states_diag>1)then
   integer,allocatable :: list_order_good_state(:)
   integer,allocatable :: list_order_not_good_state(:)
   double precision, allocatable :: e_good_state(:)
   integer :: n_good_state, n_not_good_state
   double precision :: e_0
   if(s2_eig)then
    allocate(list_order_good_state(0:n_states_diag),list_order_not_good_state(0:n_states_diag),s2_eigvalues(n_states_diag))
    list_order_good_state = 0
    list_order_not_good_state = 0
    n_good_state = 0
    n_not_good_state = 0
   endif
   call diagonalize_s2_betweenstates(psi_det,CI_eigenvectors_dressed,N_det,size(psi_det,1),size(CI_eigenvectors_dressed,1),n_states_diag,s2_eigvalues)
   do j=1,N_states_diag
     call get_s2_u0(psi_det,CI_eigenvectors_dressed(1,j),N_det,size(CI_eigenvectors_dressed,1),CI_eigenvectors_s2_dressed(j))
     if(s2_eig)then
         if(dabs(CI_eigenvectors_s2_dressed(j)-expected_s2).le.0.3d0)then
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
   if(s2_eig)then
    double precision, allocatable :: psi_coefs_tmp(:,:)
    integer, allocatable :: iorder(:),index_good_state_energetic_ordered(:)
    allocate(e_good_state(n_good_state),index_good_state_energetic_ordered(n_good_state))
    allocate(psi_coefs_tmp(N_det,n_states_diag),iorder(n_good_state))
!   print*,'good states'
    do i = 1, n_good_state
     iorder(i) = i
!    call u0_H_u_0_mrcc(e_0,CI_eigenvectors_dressed(1,list_order_good_state(i)),N_det,size(CI_eigenvectors_dressed,1),psi_det,N_int)
     call get_s2_u0(psi_det,CI_eigenvectors_dressed(1,list_order_good_state(i)),N_det,size(CI_eigenvectors_dressed,1),s2)
!    print*,'s2 = ',s2
!    print*,'e_0 = ',e_0
     e_good_state(i) = e_0
    enddo
    do i = 1, N_states_diag
     do j = 1, N_det
      psi_coefs_tmp(j,i) = CI_eigenvectors_dressed(j,i)
     enddo
    enddo

!   print*,'Sorting the good states '
    call dsort(e_good_state,iorder,n_good_state)   ! sort the good states by energy
    do i = 1, n_good_state
     index_good_state_energetic_ordered(i) = list_order_good_state(iorder(i))
    enddo
    do i = 1, n_good_state
     CI_electronic_energy_dressed(i) = e_good_state(i)
     do j = 1, N_det
      CI_eigenvectors_dressed(j,i) = psi_coefs_tmp(j,index_good_state_energetic_ordered(i))
     enddo
     call get_s2_u0(psi_det,CI_eigenvectors_dressed(1,i),N_det,size(CI_eigenvectors_dressed,1),s2)
     CI_eigenvectors_s2_dressed(i) = s2
    enddo
    do i = 1, n_not_good_state
     do j = 1, N_det
      CI_eigenvectors_dressed(j,i+n_good_state) = psi_coefs_tmp(j,list_order_not_good_state(i))
     enddo
     call get_s2_u0(psi_det,CI_eigenvectors_dressed(1,list_order_not_good_state(i)),N_det,size(CI_eigenvectors_dressed,1),s2)
!    call u0_H_u_0_mrcc(e_0,CI_eigenvectors_dressed(1,list_order_not_good_state(i)),N_det,size(CI_eigenvectors_dressed,1),psi_det,N_int)
     CI_electronic_energy_dressed(i+n_good_state) = e_0
     CI_eigenvectors_s2_dressed(i+n_good_state) = s2
    enddo
    deallocate(psi_coefs_tmp,iorder,index_good_state_energetic_ordered)
    deallocate(list_order_good_state,list_order_not_good_state,s2_eigvalues)
   endif
  endif
! print*,',CI_eigenvectors_dressed'
! do j = 1, N_states_diag 
!  print*,'State ',j
!  print*,'coef'
!  do i = 1, N_det
!   print*,CI_eigenvectors_dressed(i,j)
!  enddo
!   call get_s2_u0(psi_det,CI_eigenvectors_dressed(1,j),N_det,size(CI_eigenvectors_dressed,1),s2)
!   call u0_H_u_0_mrcc(e_0,CI_eigenvectors_dressed(1,j),N_det,size(CI_eigenvectors_dressed,1),psi_det,N_int)
!   print*,'s2                              = ',s2
!   print*,'CI_eigenvectors_s2_dressed(j)   = ',CI_eigenvectors_s2_dressed(j)
!   print*,'e_0                             = ',e_0
!   print*,'CI_electronic_energy_dressed(j) = ',CI_electronic_energy_dressed(j)
! enddo
   
END_PROVIDER

BEGIN_PROVIDER [ double precision, CI_energy_dressed, (N_states_diag) ]
  implicit none
  BEGIN_DOC
  ! N_states lowest eigenvalues of the dressed CI matrix
  END_DOC
  
  integer                        :: j
  character*(8)                  :: st
  call write_time(output_determinants)
! print*,'CI_energy_dressed'
  do j=1,N_states_diag
!   print*,'CI_electronic_energy_dressed(j) = ',CI_electronic_energy_dressed(j)
    CI_energy_dressed(j) = CI_electronic_energy_dressed(j) + nuclear_repulsion
  enddo

END_PROVIDER

subroutine diagonalize_CI_dressed
  implicit none
  BEGIN_DOC
!  Replace the coefficients of the CI states by the coefficients of the 
!  eigenstates of the CI matrix
  END_DOC
  integer :: i,j
! print*,'diagonalize_CI_dressed'
  do j=1,N_states_diag
!   print*,'j = ',j
    do i=1,N_det
      psi_coef(i,j) = CI_eigenvectors_dressed(i,j)
!     print*,'psi_coef(i,j) = ',psi_coef(i,j)
    enddo
  enddo
  SOFT_TOUCH psi_coef 

end
