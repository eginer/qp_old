program osoci_program
call debug_det(ref_bitmask,N_int)

implicit none
   call FOBOCI_lmct_mlct
   call provide_all_the_rest
end
subroutine provide_all_the_rest
implicit none
integer :: i
   call update_one_body_dm_mo
   call provide_properties
   call save_osoci_natural_mos

 

end
