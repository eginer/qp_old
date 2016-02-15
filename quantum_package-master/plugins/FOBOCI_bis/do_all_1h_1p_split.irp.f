program ohop_program
call debug_det(ref_bitmask,N_int)

implicit none
   call all_1h_1p_inact_split
   call provide_all_the_rest
end
subroutine provide_all_the_rest
implicit none
integer :: i
   call update_one_body_dm_mo_with_1h1p
   call provide_properties_1h1p

 

end
