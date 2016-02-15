
subroutine FOBOCI_spin_pol_thr
  use bitmasks
 implicit none
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: unpaired_bitmask(:,:)
 integer, allocatable :: occ(:,:)
 integer       :: n_occ_alpha, n_occ_beta
 double precision :: norm_tmp(N_states),norm_total(N_states)
 logical :: test_sym
 double precision :: thr
 double precision :: threshold
 logical :: verbose,is_ok
 verbose = .True.
 threshold = threshold_singles
 print*,'threshold = ',threshold
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
 print*,'Begining to look for the good 1h1p'
 integer :: i_particl_osoci
  do i = 1, n_inact_orb
   i_hole_osoci = list_inact(i)
   do j = 1, n_virt_orb
    i_particl_osoci = list_virt(j)
    print*,'--------------------------'
   ! First set the current generators to the one of restart
   call set_generators_to_generators_restart
   call set_psi_det_to_generators
   call check_symetry_1h1p(i_hole_osoci,i_particl_osoci,thr,test_sym)
   if(.not.test_sym)cycle
   print*,'couple hole, particle = ',i_hole_osoci,i_particl_osoci
   ! Initialize the bitmask to the restart ones
   call initialize_bitmask_to_restart_ones
   ! Impose that only the hole i_hole_osoci can be done
   call create_restart_1h_1p(i_hole_osoci,i_particl_osoci)
!  ! Update the generators 
   call set_generators_to_psi_det
   print*,'Passed set generators'
   call set_bitmask_particl_as_input(reunion_of_bitmask)
   call set_bitmask_hole_as_input(reunion_of_bitmask)
   call is_a_good_candidate(threshold,is_ok,verbose)
   print*,'is_ok = ',is_ok
   if(.not.is_ok)cycle
   ! so all the mono excitation on the new generators 
   if(.not.do_it_perturbative)then
    call all_single
   endif
   call set_intermediate_normalization_1h1p(norm_tmp,i_hole_osoci)
   do k = 1, N_states
    print*,'norm_tmp = ',norm_tmp(k)
    norm_total(k) += norm_tmp(k)
   enddo
   call update_density_matrix_osoci
  enddo
 enddo

end


