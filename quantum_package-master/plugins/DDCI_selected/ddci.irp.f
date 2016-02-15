program ddci
  implicit none
  integer                        :: i,k

  
  double precision, allocatable  :: pt2(:), norm_pert(:), H_pert_diag(:)
  double precision,allocatable :: E_before(:)
  integer                        :: N_st, degree
  N_st = N_states
  allocate (pt2(N_st), norm_pert(N_st),H_pert_diag(N_st),E_before(N_st))
  character*(64)                 :: perturbation
  
  pt2 = 1.d0
! diag_algorithm = "Lapack"
  if (N_det > N_det_max) then
    print*,'Already too much determinant in the restart wf'
    call diagonalize_CI
    call save_wavefunction
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    N_det = N_det_max
    soft_touch N_det psi_det psi_coef
    call diagonalize_CI
    call save_wavefunction
    print *,  'N_det    = ', N_det
    print *,  'N_states = ', N_states
    print *,  'PT2      = ', pt2
    print *,  'E        = ', CI_energy
    print *,  'E+PT2    = ', CI_energy+pt2
    print *,  '-----'
  endif

  threshold_davidson = 1.d-6
  soft_touch threshold_davidson davidson_criterion
  double precision :: i_H_psi_array(N_states),diag_H_mat_elem,h,i_O1_psi_array(N_states)
  if(read_wf)then
   call i_H_psi(psi_det(1,1,N_det),psi_det,psi_coef,N_int,N_det,psi_det_size,N_states,i_H_psi_array)
   h = diag_H_mat_elem(psi_det(1,1,N_det),N_int)
   selection_criterion = dabs(psi_coef(N_det,1) *  (i_H_psi_array(1) - h * psi_coef(N_det,1))) * 0.1d0
   soft_touch selection_criterion
  endif
  do while (N_det < N_det_max.and.maxval(abs(pt2(1:N_st))) > pt2_max)
    call H_apply_DDCI_selection(pt2, norm_pert, H_pert_diag,  N_st)
    call save_wavefunction

    PROVIDE  psi_coef
    PROVIDE  psi_det
    PROVIDE  psi_det_sorted

    if (N_det > N_det_max) then
       psi_det = psi_det_sorted
       psi_coef = psi_coef_sorted
       N_det = N_det_max
       soft_touch N_det psi_det psi_coef
    endif
    call diagonalize_CI
    call save_wavefunction
    if(N_states_diag.gt.1)then
     print*,'Variational Energy difference'
     do i = 2, N_st
      print*,'Delta E = ',CI_energy(i) - CI_energy(1)
     enddo
    endif
    if(N_states.gt.1)then
     print*,'Variational + perturbative Energy difference'
     do i = 2, N_st
      print*,'Delta E = ',E_before(i)+ pt2(i) - (E_before(1) + pt2(1))
     enddo
    endif
    call ezfio_set_ddci_selected_energy(CI_energy)
    if (abort_all) then
      exit
    endif
  enddo
  threshold_davidson = 1.d-13
  soft_touch threshold_davidson davidson_criterion
  if(do_pt2_end)then
    call H_apply_DDCI_pt2(pt2, norm_pert, H_pert_diag,  N_st)
    print *,  'N_det    = ', N_det
    print *,  'N_states = ', N_states
    print *,  'PT2      = ', pt2
    print *,  'E        = ', CI_energy
    print *,  'E+PT2    = ', CI_energy+pt2
  endif
  deallocate(pt2,norm_pert,E_before)
end
