
subroutine update_density_matrix_1h1p
 implicit none
 BEGIN_DOC
 ! one_body_dm_mo_alpha_osoci += Delta rho alpha
 ! one_body_dm_mo_beta_osoci  += Delta rho beta
 END_DOC
 integer :: i,j
 integer :: iorb,jorb
 do i = 1, mo_tot_num
  do j = 1, mo_tot_num
   one_body_dm_mo_alpha_1h1p(i,j) = one_body_dm_mo_alpha_1h1p(i,j) + (one_body_dm_mo_alpha(i,j) - one_body_dm_mo_alpha_generators_restart(i,j))
   one_body_dm_mo_beta_1h1p(i,j) = one_body_dm_mo_beta_1h1p(i,j) + (one_body_dm_mo_beta(i,j) - one_body_dm_mo_beta_generators_restart(i,j))
  enddo
 enddo

end


subroutine set_intermediate_normalization_1h1p(norm)
 implicit none
 double precision, intent(out) :: norm(N_states)
 integer :: i,j,degree,index_ref_generators_restart,k
 integer::  number_of_holes,n_h, number_of_particles,n_p
 logical :: is_the_hole_in_det,is_the_particl_in_det
 double precision :: inv_coef_ref_generators_restart(N_states),hij,hii,accu
 integer,allocatable :: index_good_hole_part(:)
 integer :: n_good_hole_part
 logical,allocatable :: is_a_ref_det(:)
 allocate(is_a_ref_det(N_det),index_good_hole_part(N_det))
 
 n_good_hole_part = 0
 ! Find the one holes and one hole one particle
 is_a_ref_det = .False.
 do i = 1, N_det
  ! Find the reference determinant for intermediate normalization
  call get_excitation_degree(ref_generators_restart,psi_det(1,1,i),degree,N_int)   
  if(degree == 0)then
   index_ref_generators_restart = i
   print*,'index_ref_generators_restart = ',index_ref_generators_restart
   do k = 1, N_states
    inv_coef_ref_generators_restart(k) = 1.d0/psi_coef(i,k)
    print*,'psi_coef(index_ref_generators_restart,k) = ',psi_coef(index_ref_generators_restart,k)
   enddo
  endif
  
  do j = 1, N_det_generators_restart
   call get_excitation_degree(psi_det(1,1,i),psi_det_generators_restart(1,1,j),degree,N_int)  
   if(degree == 0)then
    is_a_ref_det(i) = .True.
    exit
   endif
  enddo

  n_h = number_of_holes(psi_det(1,1,i))
  n_p = number_of_particles(psi_det(1,1,i))
  if(n_h == 1 .and. n_p == 1)then
   n_good_hole_part +=1
   index_good_hole_part(n_good_hole_part) = i
  else if(n_h == 2 .and. n_p == 2)then
   n_good_hole_part +=1
   index_good_hole_part(n_good_hole_part) = i
  else if(is_a_ref_det(i) == .False.)then
    do k = 1, N_states
     psi_coef(i,k) = 0.d0
    enddo
  endif
 enddo
 print*,''
 print*,'n_good_hole_part = ',n_good_hole_part
!do k = 1,N_states
! print*,'state ',k
! do i = 1, n_good_hole_part
!  print*,'index_good_hole_part(i)',index_good_hole_part(i)
!  print*,'psi_coef(index_good_hole_part) = ',psi_coef(index_good_hole_part(i),k)/psi_coef(index_ref_generators_restart,k)
! enddo
! print*,''
!enddo
 norm = 0.d0

 ! Set the wave function to the intermediate normalization
 do k = 1, N_states
  do i = 1, N_det
   psi_coef(i,k) = psi_coef(i,k) * inv_coef_ref_generators_restart(k)
  enddo
 enddo
 do k = 1,N_states
  print*,'state ',k
  do i = 1, N_det
   if (is_a_ref_det(i))then
    print*,'i,psi_coef_ref = ',psi_coef(i,k)
    cycle
   endif
   norm(k) += psi_coef(i,k) * psi_coef(i,k)
  enddo
  print*,'norm = ',norm(k)
 enddo
 deallocate(index_good_hole_part,is_a_ref_det)
 soft_touch psi_coef
end






subroutine initialize_density_matrix_1h1p
 implicit none
 one_body_dm_mo_alpha_1h1p = one_body_dm_mo_alpha_generators_restart
 one_body_dm_mo_beta_1h1p = one_body_dm_mo_beta_generators_restart
end

subroutine rescale_density_matrix_1h1p(norm)
 implicit none
 double precision, intent(in) :: norm(N_states)
 integer :: i,j
 double precision :: norm_tmp
 norm_tmp = 0.d0
 do i = 1, N_states
  norm_tmp += norm(i)
 enddo
 print*,'norm = ',norm_tmp
 
 do i = 1, mo_tot_num
  do j = 1,mo_tot_num
   one_body_dm_mo_alpha_1h1p(i,j) = one_body_dm_mo_alpha_1h1p(i,j) * norm_tmp
   one_body_dm_mo_beta_1h1p(j,i) = one_body_dm_mo_beta_1h1p(j,i) * norm_tmp
  enddo
 enddo
!soft_touch one_body_dm_mo_alpha_1h1p
end


 subroutine update_one_body_dm_mo_with_1h1p
   implicit none
   integer :: i
   double precision :: accu_tot,accu_sd
   print*,'touched the one_body_dm_mo_beta'
   one_body_dm_mo_alpha = one_body_dm_mo_alpha_1h1p
   one_body_dm_mo_beta = one_body_dm_mo_beta_1h1p
   touch one_body_dm_mo_alpha  one_body_dm_mo_beta 
   accu_tot = 0.d0
   accu_sd  = 0.d0
   do i = 1, mo_tot_num
    accu_tot += one_body_dm_mo_alpha(i,i) + one_body_dm_mo_beta(i,i)
    accu_sd  += one_body_dm_mo_alpha(i,i) - one_body_dm_mo_beta(i,i)
   enddo
   print*,'accu_tot = ',accu_tot
   print*,'accu_sdt = ',accu_sd 
 end
 
 subroutine provide_properties_1h1p
   implicit none
   integer :: i
   double precision :: accu
    accu= 0.d0
    do i = 1, nucl_num
     accu += mulliken_spin_densities(i)
     print*,i,nucl_charge(i),mulliken_spin_densities(i)
    enddo
    print*,'Sum of Mulliken SD = ',accu
 do i = 1, nucl_num
  write(*,'(I2,X,F3.1,X,4(F16.10,X))')i,nucl_charge(i),spin_density_at_nucleous(i),iso_hcc_gauss(i),iso_hcc_mhz(i),iso_hcc_cm_1(i)
 enddo
    
 end

