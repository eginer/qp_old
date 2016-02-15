program restart_more_singles
  ! Generates and select single excitations 
  ! on the top of a given restart wave function
  read_wf = .true.
  touch read_wf 
  print*,'ref_bitmask_energy = ',ref_bitmask_energy
  call routine

end 
subroutine routine
  implicit none
  integer                        :: i,k
  double precision, allocatable  :: pt2(:), norm_pert(:), H_pert_diag(:)
  integer                        :: N_st, degree
  double precision :: E_before
  integer :: n_det_before
  N_st = N_states
  allocate (pt2(N_st), norm_pert(N_st),H_pert_diag(N_st))
  i = 0
  print*,'N_det = ',N_det
  print*,'n_det_max = ',n_det_max
  print*,'pt2_max = ',pt2_max
  pt2=-1.d0
  E_before = ref_bitmask_energy
  do while (N_det < n_det_max.and.maxval(abs(pt2(1:N_st))) > pt2_max)
    n_det_before = N_det
    i += 1
    print*,'-----------------------'
    print*,'i = ',i
    call H_apply_just_1h_1p(pt2, norm_pert, H_pert_diag,  N_st)
    call diagonalize_CI
    print*,'N_det = ',N_det
    print*,'E        = ',CI_energy(1)
    print*,'pt2      = ',pt2(1)
    print*,'E+PT2    = ',E_before + pt2(1)
    E_before = CI_energy(1)
    call save_wavefunction
    if(n_det_before == N_det)then
     selection_criterion = selection_criterion * 0.5d0
    endif
    
  enddo

  call save_wavefunction
  deallocate(pt2,norm_pert)
end
