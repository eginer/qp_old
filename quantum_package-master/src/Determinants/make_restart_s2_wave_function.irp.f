program make_restart_s2_wave_function
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
 implicit none
 call make_s2_eigenfunction
 call diagonalize_CI
 call save_wavefunction
end
