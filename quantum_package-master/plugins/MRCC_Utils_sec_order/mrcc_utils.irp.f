
 BEGIN_PROVIDER [integer, pert_determinants, (N_states, psi_det_size) ]
 END_PROVIDER 


!BEGIN_PROVIDER [ double precision, lambda_mrcc, (N_states,psi_det_size) ]
!&BEGIN_PROVIDER [ double precision, lambda_pert, (N_states,psi_det_size) ] 
!implicit none
!BEGIN_DOC
!! cm/<Psi_0|H|D_m> or perturbative 1/Delta_E(m)
!END_DOC
!integer :: i,k
!double precision               :: ihpsi(N_states), hii,delta_e_eff
!integer :: i_ok
!i_ok = 0

!double precision :: phase_restart(N_states)
!do k = 1, N_states
! phase_restart(k) = dsign(1.d0,psi_ref_coef_restart(1,k)/psi_ref_coef(1,k))
!enddo
!
!do i=1,N_det_non_ref
!  call i_h_psi(psi_non_ref(1,1,i), psi_ref_restart, psi_ref_coef_restart, N_int, N_det_ref,&
!      size(psi_ref_coef_restart,1), n_states, ihpsi)
!  call i_H_j(psi_non_ref(1,1,i),psi_non_ref(1,1,i),N_int,hii)
!  do k=1,N_states
!    lambda_pert(k,i) = 1.d0 / (psi_ref_energy_diagonalized(k)-hii)
!    if((ihpsi(k) * lambda_pert(k,i))/psi_non_ref_coef_restart(i,k) .ge. 0.5d0 & 
!       .and. (ihpsi(k) * lambda_pert(k,i))/psi_non_ref_coef_restart(i,k) > 0.d0)then  ! test on the first order coefficient
!      call i_h_psi(psi_non_ref(1,1,i), psi_ref, psi_ref_coef, N_int, N_det_ref,size(psi_ref_coef,1), n_states, ihpsi)
!      lambda_mrcc(k,i) = psi_non_ref_coef(i,k)/ihpsi(k)
!    else
!      i_ok +=1
!      lambda_mrcc(k,i) = lambda_pert(k,i)
!    endif
!  enddo
!enddo
!print*,'N_det_non_ref = ',N_det_non_ref
!print*,'Number of Perturbatively treated determinants = ',i_ok
!print*,'psi_coef_ref_ratio = ',psi_ref_coef(2,1)/psi_ref_coef(1,1)

!END_PROVIDER

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
 END_PROVIDER 

 BEGIN_PROVIDER [ logical, oscillations ]
 implicit none
 oscillations = .False.
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, amplitudes_before, (N_det_non_ref,N_det_ref,N_states) ]

 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, amplitudes, (N_det_non_ref,N_det_ref,N_states) ]
 implicit none
 BEGIN_DOC
! Amplitudes t_Ii
 END_DOC
 integer :: i,j,k,l
 integer*2,allocatable  :: excitation_operators(:,:)
 double precision, allocatable :: amplitudes_phase_less(:)
 integer, allocatable :: index_connected(:)
 integer, save :: isave = 0
 if (isave == 0)then
  isave = 1
  allocate(excitation_operators(5,N_det_non_ref))
  allocate(amplitudes_phase_less(N_det_non_ref))
  allocate(index_connected(N_det_non_ref))
  amplitudes = 0.d0
  do i = 1, N_det_ref
   integer :: i_state,N_connect_ref
   do i_state = 1, N_states
    call get_excitation_operators_for_one_ref(psi_ref(1,1,i),i_state,N_det_non_ref,N_connect_ref,excitation_operators,amplitudes_phase_less,index_connected)
    do k = 1, N_connect_ref
     amplitudes( index_connected(k),i,i_state ) = amplitudes_phase_less(k)
    enddo
   enddo
  enddo
  print*,'Passed the AMPLITUDES !'
  deallocate(excitation_operators)
  deallocate(amplitudes_phase_less)
  deallocate(index_connected)
 else 
  integer :: ii
  double precision :: lambda_mrcc_new(N_det_non_ref),hijj,ihpsi
  double precision :: hij,sec_order,H_ref(N_det_ref),H_non_ref(N_det_non_ref)
  integer          :: idx(0:N_det_non_ref)
  do i_state = 1, N_states
   do i = 1, N_det_non_ref 
   lambda_mrcc_new(i) = 0.d0
   call filter_connected_i_H_psi0(psi_ref,psi_ref(1,1,i),N_int,N_det_ref,idx)  
   H_ref = 0.d0
   do ii=1,idx(0)
     k = idx(ii)
     !DEC$ FORCEINLINE
     call i_H_j(psi_ref(1,1,k),psi_ref(1,1,i),N_int,hij)
     lambda_mrcc_new(i) += hij * psi_ref_coef(k,i_state)
     H_ref(k)  = hij 
   enddo
    call filter_connected_i_H_psi0(psi_non_ref,psi_non_ref(1,1,i),N_int,N_det_non_ref,idx) 
    H_non_ref = 0.d0
    do ii=1,idx(0)
      k = idx(ii)
      if(k==i)cycle
      !DEC$ FORCEINLINE
      call i_H_j(psi_non_ref(1,1,k),psi_non_ref(1,1,i),N_int,hij)
      H_non_ref(k) = hij
    enddo
    double precision :: contrib_ref
    ! CALCULATION OF THE LAMBDA
    do ii = 1, N_det_ref
     contrib_ref = 0.d0
     do k = 1, N_det_non_ref
      contrib_ref +=  amplitudes_before(k,ii,i_state) * H_non_ref(k) 
     enddo
     lambda_mrcc_new(i) += contrib_ref * psi_ref_coef(ii,i_state)
    enddo
    lambda_mrcc_new(i) = psi_non_ref_coef(i,i_state) / lambda_mrcc_new(i)

    ! CALCULATION OF THE AMPLITUDES
    do ii = 1, N_det_ref
     contrib_ref = H_ref(ii)
     do k = 1, N_det_non_ref
      contrib_ref +=  amplitudes_before(k,ii,i_state) * H_non_ref(k) 
     enddo
     amplitudes(i,ii,i_state)  = contrib_ref * lambda_mrcc_new(i)
     if(dabs(amplitudes(i,ii,i_state)).ge.0.5d0)then
      amplitudes(i,ii,i_state) = 0.d0
     endif
    enddo

   enddo
  enddo

 endif

 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, delta_ij, (N_det_ref,N_det_non_ref,N_states) ]
&BEGIN_PROVIDER [ double precision, delta_ii, (N_det_ref,N_states) ]
 implicit none
 BEGIN_DOC
 ! Dressing matrix in N_det basis
 END_DOC
 integer :: i,j,m
 delta_ij = 0.d0
 delta_ii = 0.d0
 call H_apply_mrcc(delta_ij,delta_ii,N_det_ref,N_det_non_ref)
 double precision :: max_delta
 double precision :: accu
 integer :: imax,jmax
 max_delta = 0.d0
 accu = 0.d0
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
 print*,'MAX VAL OF DRESING = ',delta_ij(imax,jmax,1)
 print*,'imax,jmax = ',imax,jmax
 print*,'psi_ref_coef(imax,1)     = ',psi_ref_coef(imax,1)
 print*,'psi_non_ref_coef(jmax,1) = ',psi_non_ref_coef(jmax,1)
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
END_PROVIDER


 BEGIN_PROVIDER [ double precision, CI_electronic_energy_dressed, (N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors_dressed, (N_det,N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors_s2_dressed, (N_states_diag) ]
   implicit none
   BEGIN_DOC
   ! Eigenvectors/values of the CI matrix
   END_DOC
   integer                        :: i,j
   
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
     
     double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:)
     allocate (eigenvectors(size(H_matrix_dressed,1),N_det))
     allocate (eigenvalues(N_det))
     call lapack_diag(eigenvalues,eigenvectors,                      &
         H_matrix_dressed,size(H_matrix_dressed,1),N_det)
     CI_electronic_energy_dressed(:) = 0.d0
     do i=1,N_det
       CI_eigenvectors_dressed(i,1) = eigenvectors(i,1)
     enddo
     integer                        :: i_state
     double precision               :: s2
     i_state = 0
     if (s2_eig) then
       do j=1,N_det
         call get_s2_u0(psi_det,eigenvectors(1,j),N_det,N_det,s2)
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
     deallocate(eigenvectors,eigenvalues)
   endif
   
END_PROVIDER

BEGIN_PROVIDER [ double precision, CI_energy_dressed, (N_states_diag) ]
  implicit none
  BEGIN_DOC
  ! N_states lowest eigenvalues of the dressed CI matrix
  END_DOC
  
  integer                        :: j
  character*(8)                  :: st
  call write_time(output_determinants)
  do j=1,N_states_diag
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
  do j=1,N_states_diag
    do i=1,N_det
      psi_coef(i,j) = CI_eigenvectors_dressed(i,j)
    enddo
  enddo
  SOFT_TOUCH psi_coef 

end
