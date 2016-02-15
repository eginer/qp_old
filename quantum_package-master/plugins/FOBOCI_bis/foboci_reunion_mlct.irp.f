program osoci_program
implicit none
   call FOBOCI_mlct
   call provide_all_the_rest
end
subroutine provide_all_the_rest
implicit none
integer :: i
   call update_one_body_dm_mo
   call provide_properties
   print*,''
   print*,''
   print*,''
   print*,'--------------------------'
   print*,''
   print*,'DOING THE REUNION OF ALL PERTINENT MLCT'
   print*,''
   print*,''
   print*,''
   print*,''
!  ! Update the generators 
   call set_mlct_to_generators_restart
   call set_psi_det_to_generators
   print*,'New generators !!'
   call save_wavefunction
   call set_bitmask_particl_as_input(reunion_of_bitmask)
   call set_bitmask_hole_as_input(reunion_of_bitmask)
   ! so all the mono excitation on the new generators 
   call all_single
   call provide_properties
   call save_wavefunction


 

end
