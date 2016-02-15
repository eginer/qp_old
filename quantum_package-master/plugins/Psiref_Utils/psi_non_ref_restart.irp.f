
 BEGIN_PROVIDER [ integer(bit_kind), psi_non_ref_restart,  (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_non_ref_coef_restart, (psi_det_size,n_states) ]
 implicit none
 BEGIN_DOC
  ! Set of determinants which are not part of the reference, defined from the application
  ! of the reference bitmask on the determinants. 
  ! idx_non_ref gives the indice of the determinant in psi_det.
  ! But this is with respect to the restart wave function. 
 END_DOC
 integer                        :: i_non_ref,j,k
 integer                        :: degree
 logical                        :: in_ref
 integer, save                  :: ifirst = 0 
 if(ifirst==0)then
  ifirst = 1
  i_non_ref =0
  do k=1,N_det
    in_ref = .False.
    do j=1,N_det_ref
      call get_excitation_degree(psi_ref(1,1,j), psi_det(1,1,k), degree, N_int)
      if (degree == 0) then
        in_ref = .True.
        exit
      endif
    enddo
    if (.not.in_ref) then
      double precision :: hij
      i_non_ref += 1
      do j=1,N_int
        psi_non_ref_restart(j,1,i_non_ref) = psi_det(j,1,k)
        psi_non_ref_restart(j,2,i_non_ref) = psi_det(j,2,k)
      enddo
      do j=1,N_states
        psi_non_ref_coef_restart(i_non_ref,j) = psi_coef(k,j)
      enddo
    endif
  enddo
 endif
END_PROVIDER
