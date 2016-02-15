subroutine huckel_guess
  implicit none
  BEGIN_DOC
! Build the MOs using the extended Huckel model
  END_DOC
  integer                        :: i,j
  double precision               :: accu
  double precision               :: c
  character*(64)                 :: label

  integer :: iorb,jorb
  print*,'n_act_orb   = ' ,n_act_orb
  print*,'n_inact_orb = ',n_inact_orb
  print*,'n_virt_orb  = ',n_virt_orb
  do i = 1, n_act_orb
   iorb = list_act(i)
   do j = 1, n_inact_orb
    jorb = list_inact(j)
    Fock_matrix_mo(iorb,jorb) = 0.d0
    Fock_matrix_mo(jorb,iorb) = 0.d0
   enddo
   do j = 1, n_virt_orb
    jorb = list_virt(j)
    Fock_matrix_mo(iorb,jorb) = 0.d0
    Fock_matrix_mo(jorb,iorb) = 0.d0
   enddo
  enddo 
  touch Fock_matrix_mo fock_matrix_diag_mo
        
  call mo_as_eigvectors_of_mo_matrix(Fock_matrix_mo,          &
                                     size(Fock_matrix_mo,1),  &
                                     size(Fock_matrix_mo,2),label,1)
  mo_label = "Canonical"
  TOUCH mo_coef mo_label
  call save_mos

! c = 0.5d0 * 1.75d0

! do j=1,ao_num
!   !DIR$ VECTOR ALIGNED
!   do i=1,ao_num
!     Fock_matrix_ao(i,j) = c*ao_overlap(i,j)*(ao_mono_elec_integral_diag(i) + &
!                                                ao_mono_elec_integral_diag(j))
!   enddo
!   Fock_matrix_ao(j,j) = Fock_matrix_alpha_ao(j,j)
! enddo
! TOUCH Fock_matrix_ao
! mo_coef = eigenvectors_fock_matrix_mo
! SOFT_TOUCH mo_coef
! call save_mos

end
