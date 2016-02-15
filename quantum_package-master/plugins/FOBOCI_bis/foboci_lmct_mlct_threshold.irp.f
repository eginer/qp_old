
subroutine FOBOCI_lmct_mlct
  use bitmasks
 implicit none
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: unpaired_bitmask(:,:)
 integer, allocatable :: occ(:,:)
 integer       :: n_occ_alpha, n_occ_beta
 double precision :: norm_tmp,norm_total
 logical :: test_sym
 double precision :: thr
 double precision :: threshold
 logical :: verbose,is_ok
 verbose = .False.
 threshold = 1.d-2
 thr = 1.d-12
 allocate(unpaired_bitmask(N_int,2))
 allocate (occ(N_int*bit_kind_size,2))
 do i = 1, N_int
  unpaired_bitmask(i,1) = unpaired_alpha_electrons(i)
  unpaired_bitmask(i,2) = unpaired_alpha_electrons(i)
 enddo
 norm_total = 0.d0
 call initialize_density_matrix_osoci
 call bitstring_to_list(inact_bitmask(1,1), occ(1,1), n_occ_beta, N_int)
 print*,''
 print*,''
 print*,'DOING FIRST LMCT !!'
  do i = 1, n_inact_orb
   i_hole_osoci = list_inact(i)
   print*,'--------------------------'
   ! First set the current generators to the one of restart
   call set_generators_to_generators_restart
   call set_psi_det_to_generators
   call check_symetry(i_hole_osoci,thr,test_sym)
   if(.not.test_sym)cycle
   print*,'i_hole_osoci = ',i_hole_osoci
   ! Initialize the bitmask to the restart ones
   call initialize_bitmask_to_restart_ones
   ! Impose that only the hole i_hole_osoci can be done
   call modify_bitmasks_for_hole(i_hole_osoci)
   call print_generators_bitmasks_holes
   ! Impose that only the active part can be reached 
   call set_bitmask_particl_as_input(unpaired_bitmask)
   call all_single_h_core
!  ! Update the generators 
   call set_generators_to_psi_det
   call set_bitmask_particl_as_input(reunion_of_bitmask)
   call set_bitmask_hole_as_input(reunion_of_bitmask)
   call is_a_good_candidate(threshold,is_ok,verbose)
   print*,'is_ok = ',is_ok
   if(.not.is_ok)cycle
!  ! so all the mono excitation on the new generators 
   call all_single
   call set_intermediate_normalization_lmct(norm_tmp,i_hole_osoci)
   print*,'norm_tmp = ',norm_tmp
   norm_total += norm_tmp
   call update_density_matrix_osoci
 enddo

 print*,''
 print*,'DOING FIRST MLCT !!'
  do i = 1, n_virt_orb
   integer :: i_particl_osoci
   i_particl_osoci = list_virt(i)
   print*,'--------------------------'
   ! First set the current generators to the one of restart
   call set_generators_to_generators_restart
   call set_psi_det_to_generators
   call check_symetry(i_particl_osoci,thr,test_sym)
   if(.not.test_sym)cycle
   print*,'i_particl_osoci= ',i_particl_osoci
   ! Initialize the bitmask to the restart ones
   call initialize_bitmask_to_restart_ones
   ! Impose that only the hole i_hole_osoci can be done
   call modify_bitmasks_for_particl(i_particl_osoci)
   call print_generators_bitmasks_holes
   ! Impose that only the active part can be reached 
   call set_bitmask_hole_as_input(unpaired_bitmask)
   call all_single_h_core
!  ! Update the generators 
   call set_generators_to_psi_det
   call set_bitmask_particl_as_input(reunion_of_bitmask)
   call set_bitmask_hole_as_input(reunion_of_bitmask)
!  ! so all the mono excitation on the new generators 
   call is_a_good_candidate(threshold,is_ok,verbose)
   print*,'is_ok = ',is_ok
   if(.not.is_ok)cycle
   call all_single
   call set_intermediate_normalization_mlct(norm_tmp,i_particl_osoci)
   print*,'norm_tmp = ',norm_tmp
   norm_total += norm_tmp
   call update_density_matrix_osoci
 enddo

  if(.False.)then
   print*,'LAST loop for all the 1h-1p'
   print*,'--------------------------'
   ! First set the current generators to the one of restart
   call set_generators_to_generators_restart
   call set_psi_det_to_generators
   call initialize_bitmask_to_restart_ones
   ! Impose that only the hole i_hole_osoci can be done
   call set_bitmask_particl_as_input(inact_virt_bitmask)
   call set_bitmask_hole_as_input(inact_virt_bitmask)
!  call set_bitmask_particl_as_input(reunion_of_bitmask)
!  call set_bitmask_hole_as_input(reunion_of_bitmask)
   call all_single
   call set_intermediate_normalization_1h1p(norm_tmp)
   norm_total += norm_tmp
   call update_density_matrix_osoci
  endif


   print*,'norm_total = ',norm_total
   norm_total += 1.d0 
   norm_total = 1.d0/norm_total
   call rescale_density_matrix_osoci(norm_total)
   double precision :: accu
   accu = 0.d0
   do i = 1, mo_tot_num
    accu += one_body_dm_mo_alpha_osoci(i,i) + one_body_dm_mo_beta_osoci(i,i)
   enddo
   print*,'accu = ',accu
end


subroutine FOBOCI_mlct
  use bitmasks
 implicit none
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: unpaired_bitmask(:,:)
 integer, allocatable :: occ(:,:)
 integer       :: n_occ_alpha, n_occ_beta
 double precision :: norm_tmp,norm_total
 logical :: test_sym
 double precision :: thr
 double precision :: threshold
 logical :: verbose,is_ok
 verbose = .False.
 threshold = 1.d-2
 thr = 1.d-12
 allocate(unpaired_bitmask(N_int,2))
 allocate (occ(N_int*bit_kind_size,2))
 do i = 1, N_int
  unpaired_bitmask(i,1) = unpaired_alpha_electrons(i)
  unpaired_bitmask(i,2) = unpaired_alpha_electrons(i)
 enddo
 norm_total = 0.d0
 call initialize_density_matrix_osoci
 call bitstring_to_list(inact_bitmask(1,1), occ(1,1), n_occ_beta, N_int)
 print*,''
 print*,''
 print*,''
 print*,'DOING FIRST MLCT !!'
  do i = 1, n_virt_orb
   integer :: i_particl_osoci
   i_particl_osoci = list_virt(i)
   print*,'--------------------------'
   ! First set the current generators to the one of restart
   call set_generators_to_generators_restart
   call set_psi_det_to_generators
   call check_symetry(i_particl_osoci,thr,test_sym)
   if(.not.test_sym)cycle
   print*,'i_particl_osoci= ',i_particl_osoci
   ! Initialize the bitmask to the restart ones
   call initialize_bitmask_to_restart_ones
   ! Impose that only the hole i_hole_osoci can be done
   call modify_bitmasks_for_particl(i_particl_osoci)
   call print_generators_bitmasks_holes
   ! Impose that only the active part can be reached 
   call set_bitmask_hole_as_input(unpaired_bitmask)
   call all_single_h_core
!  ! Update the generators 
   call set_generators_to_psi_det
   call set_bitmask_particl_as_input(reunion_of_bitmask)
   call set_bitmask_hole_as_input(reunion_of_bitmask)
!  ! so all the mono excitation on the new generators 
   call is_a_good_candidate(threshold,is_ok,verbose)
   print*,'is_ok = ',is_ok
   if(.not.is_ok)cycle
   call all_single
   call set_intermediate_normalization_mlct(norm_tmp,i_particl_osoci)
   print*,'norm_tmp = ',norm_tmp
   norm_total += norm_tmp
   call update_density_matrix_osoci
 enddo

 print*,'norm_total = ',norm_total
 norm_total += 1.d0 
 norm_total = 1.d0/norm_total
 call rescale_density_matrix_osoci(norm_total)
 double precision :: accu
 accu = 0.d0
 do i = 1, mo_tot_num
  accu += one_body_dm_mo_alpha_osoci(i,i) + one_body_dm_mo_beta_osoci(i,i)
 enddo
 print*,'accu = ',accu
end


subroutine FOBOCI_lmct
  use bitmasks
 implicit none
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: unpaired_bitmask(:,:)
 integer, allocatable :: occ(:,:)
 integer       :: n_occ_alpha, n_occ_beta
 double precision :: norm_tmp,norm_total
 logical :: test_sym
 double precision :: thr
 double precision :: threshold
 logical :: verbose,is_ok
 verbose = .False.
 threshold = 1.d-2
 thr = 1.d-12
 allocate(unpaired_bitmask(N_int,2))
 allocate (occ(N_int*bit_kind_size,2))
 do i = 1, N_int
  unpaired_bitmask(i,1) = unpaired_alpha_electrons(i)
  unpaired_bitmask(i,2) = unpaired_alpha_electrons(i)
 enddo
 norm_total = 0.d0
 call initialize_density_matrix_osoci
 call bitstring_to_list(inact_bitmask(1,1), occ(1,1), n_occ_beta, N_int)
 print*,''
 print*,''
 print*,'DOING FIRST LMCT !!'
  do i = 1, n_inact_orb
   i_hole_osoci = list_inact(i)
   print*,'--------------------------'
   ! First set the current generators to the one of restart
   call set_generators_to_generators_restart
   call set_psi_det_to_generators
   call check_symetry(i_hole_osoci,thr,test_sym)
   if(.not.test_sym)cycle
   print*,'i_hole_osoci = ',i_hole_osoci
   ! Initialize the bitmask to the restart ones
   call initialize_bitmask_to_restart_ones
   ! Impose that only the hole i_hole_osoci can be done
   call modify_bitmasks_for_hole(i_hole_osoci)
   call print_generators_bitmasks_holes
   ! Impose that only the active part can be reached 
   call set_bitmask_particl_as_input(unpaired_bitmask)
   call all_single_h_core
!  ! Update the generators 
   call set_generators_to_psi_det
   call set_bitmask_particl_as_input(reunion_of_bitmask)
   call set_bitmask_hole_as_input(reunion_of_bitmask)
   call is_a_good_candidate(threshold,is_ok,verbose)
   print*,'is_ok = ',is_ok
   if(.not.is_ok)cycle
!  ! so all the mono excitation on the new generators 
   call all_single
!  call set_intermediate_normalization_lmct_bis(norm_tmp,i_hole_osoci)
   call set_intermediate_normalization_lmct(norm_tmp,i_hole_osoci)
   print*,'norm_tmp = ',norm_tmp
   norm_total += norm_tmp
   call update_density_matrix_osoci
 enddo

   print*,'norm_total = ',norm_total
   norm_total += 1.d0 
   norm_total = 1.d0/norm_total
   call rescale_density_matrix_osoci(norm_total)
   double precision :: accu
   accu = 0.d0
   do i = 1, mo_tot_num
    accu += one_body_dm_mo_alpha_osoci(i,i) + one_body_dm_mo_beta_osoci(i,i)
   enddo
   print*,'accu = ',accu
end
