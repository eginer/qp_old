
subroutine all_1h_1p_inact_split
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
 call initialize_density_matrix_1h1p
 call bitstring_to_list(inact_bitmask(1,1), occ(1,1), n_occ_beta, N_int)
 print*,''
 print*,''
 print*,'mulliken spin population analysis'
 accu =0.d0
 do i = 1, nucl_num
  accu += mulliken_spin_densities(i)
  print*,i,nucl_charge(i),mulliken_spin_densities(i)
 enddo
  do i = 1, n_inact_orb
   i_hole_osoci = list_inact(i)
   print*,'--------------------------'
   call set_generators_to_generators_restart
   call set_psi_det_to_generators
   print*,'i_hole_osoci = ',i_hole_osoci
   call initialize_bitmask_to_restart_ones
   call modify_bitmasks_for_hole(i_hole_osoci)
   do k = 1, n_act_orb
    call modify_bitmasks_for_hole_in_out(list_act(k))
   enddo
   call print_generators_bitmasks_holes
   call set_bitmask_particl_as_input(reunion_of_bitmask)
   threshold_davidson = 1.d-10
   soft_touch threshold_davidson davidson_criterion
   call all_but_1h_1p_routine
   call set_intermediate_normalization_1h1p(norm_tmp)
   do k = 1, N_states
    print*,'norm_tmp = ',norm_tmp(k)
    norm_total(k) += norm_tmp(k)
   enddo
   call update_density_matrix_1h1p
 enddo

 print*,'norm_total = ',norm_total
 norm_total = 1.d0/(1.d0 + norm_total)
 call rescale_density_matrix_1h1p(norm_total)
 double precision :: accu
 accu = 0.d0
 do i = 1, mo_tot_num
  accu += one_body_dm_mo_alpha_1h1p(i,i) + one_body_dm_mo_beta_1h1p(i,i)
 enddo
 print*,'accu = ',accu
end


