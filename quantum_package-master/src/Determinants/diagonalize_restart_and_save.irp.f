program diag_restart_and_save
 implicit none
 read_wf = .True.
 threshold_davidson = 1.d-10 
 touch threshold_davidson davidson_criterion
 touch read_wf
 print*,'----------'
 print*,'N_det = ',N_det
 call routine
end

subroutine routine
 implicit none
 call diagonalize_CI
 call save_wavefunction


end
