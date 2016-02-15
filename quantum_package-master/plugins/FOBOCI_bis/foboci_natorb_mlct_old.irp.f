program osoci_program

implicit none
   call FOBOCI_mlct_old
   call provide_all_the_rest
end
subroutine provide_all_the_rest
implicit none
integer :: i
   call update_one_body_dm_mo
   call provide_properties
   call save_osoci_natural_mos
   call save_mos

 

end
