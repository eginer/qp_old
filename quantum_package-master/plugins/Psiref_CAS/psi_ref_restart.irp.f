use bitmasks

 BEGIN_PROVIDER [ integer(bit_kind), psi_ref_restart, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_ref_coef_restart,  (psi_det_size,n_states) ]
  implicit none
  BEGIN_DOC
  ! Projection of the CAS wave function on the restart wave function. 
  END_DOC
  integer :: i,j,k
  integer, save                  :: ifirst

  if(ifirst == 0)then
   ifirst = 1
   do i=1,N_det_ref
     do k=1,N_int
       psi_ref_restart(k,1,i) = psi_cas(k,1,i)
       psi_ref_restart(k,2,i) = psi_cas(k,2,i)
     enddo
   enddo
   do k=1,N_states
     do i=1,N_det_ref
       psi_ref_coef_restart(i,k) = psi_cas_coef(i,k)
     enddo
   enddo
  endif

END_PROVIDER

