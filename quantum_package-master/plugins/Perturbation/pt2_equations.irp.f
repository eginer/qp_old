BEGIN_TEMPLATE

subroutine pt2_epstein_nesbet ($arguments)
  use bitmasks
  implicit none
  $declarations
  
  BEGIN_DOC
  ! compute the standard Epstein-Nesbet perturbative first order coefficient and second order energetic contribution
  !
  ! for the various N_st states.
  !
  ! c_pert(i) = <psi(i)|H|det_pert>/( E(i) - <det_pert|H|det_pert> )
  !
  ! e_2_pert(i) = <psi(i)|H|det_pert>^2/( E(i) - <det_pert|H|det_pert> )
  !
  END_DOC
  
  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock, h
  double precision               :: i_H_psi_array(N_st)
  PROVIDE  selection_criterion

  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  !call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array)
  call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)
  
  
  
  call debug_det(det_pert,N_int)
  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if(CI_electronic_energy(i)>h.and.CI_electronic_energy(i).ne.0.d0)then
      c_pert(i) = -1.d0
      e_2_pert(i) = selection_criterion*selection_criterion_factor*2.d0
    else if  (dabs(CI_electronic_energy(i) - h) > 1.d-6) then
        c_pert(i) = i_H_psi_array(i) / (CI_electronic_energy(i) - h)
        H_pert_diag(i) = h*c_pert(i)*c_pert(i)
        e_2_pert(i) = c_pert(i) * i_H_psi_array(i)
    else
      c_pert(i) = -1.d0
      e_2_pert(i) = -dabs(i_H_psi_array(i))
      H_pert_diag(i) = h
    endif
  enddo
! print*,'e_2_pert = ',e_2_pert
! print*,'c_pert = ',c_pert
! pause
  
end



subroutine pt2_hcc_epstein_nesbet ($arguments)
  use bitmasks
  implicit none
  $declarations
  
  BEGIN_DOC
  ! compute the first order contribution to the hcc 
  !
  ! for the various N_st states.
  !
  ! c_pert(i) = <psi(i)|H|det_pert>/( E(i) - <det_pert|H|det_pert> )
  !
  ! e_2_pert(i) = <psi(i)|hcc|det_pert> * c_pert(i)
  !
  END_DOC
  
  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock, h
  double precision               :: i_H_psi_array(N_st)
  double precision               :: i_O1_psi_array(N_st)
  PROVIDE  selection_criterion

  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  !call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array)
  call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)
  call i_O1_psi_alpha_beta(integral_hcc,det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_O1_psi_array)
  
  
  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if(CI_electronic_energy(i)>h.and.CI_electronic_energy(i).ne.0.d0)then
      c_pert(i) = -1.d0
      e_2_pert(i) = selection_criterion*selection_criterion_factor*2.d0
    else if  (dabs(CI_electronic_energy(i) - h) > 1.d-6) then
        c_pert(i) = i_H_psi_array(i) / (CI_electronic_energy(i) - h)
        H_pert_diag(i) = c_pert(i) * (i_O1_psi_array(i)+i_O1_psi_array(i) )
        e_2_pert(i) = -dabs(c_pert(i) * i_O1_psi_array(i))
    else
      c_pert(i) = -1.d0
      e_2_pert(i) = -dabs(i_O1_psi_array(i))
      H_pert_diag(i) = h
    endif
  enddo
  
end

subroutine pt2_epstein_nesbet_2x2 ($arguments)
  use bitmasks
  implicit none
  $declarations
  
  BEGIN_DOC
  ! compute the Epstein-Nesbet 2x2 diagonalization coefficient and energetic contribution
  !
  ! for the various N_st states.
  !
  ! e_2_pert(i) = 0.5 * (( <det_pert|H|det_pert> -  E(i) )  - sqrt( ( <det_pert|H|det_pert> -  E(i)) ^2 + 4 <psi(i)|H|det_pert>^2  )
  !
  ! c_pert(i) = e_2_pert(i)/ <psi(i)|H|det_pert>
  !
  END_DOC
  
  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock,delta_e, h
  double precision               :: i_H_psi_array(N_st)
  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  PROVIDE CI_electronic_energy

  !call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array)
  call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)
  
  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if (i_H_psi_array(i) /= 0.d0) then
      delta_e = h - CI_electronic_energy(i)
      if (delta_e > 0.d0) then
        e_2_pert(i) = 0.5d0 * (delta_e - dsqrt(delta_e * delta_e + 4.d0 * i_H_psi_array(i) * i_H_psi_array(i)))
      else
        e_2_pert(i) = 0.5d0 * (delta_e + dsqrt(delta_e * delta_e + 4.d0 * i_H_psi_array(i) * i_H_psi_array(i)))
      endif
      if (dabs(i_H_psi_array(i)) > 1.d-6) then
        c_pert(i) = e_2_pert(i)/i_H_psi_array(i)
      else
        c_pert(i) = 0.d0
      endif
      H_pert_diag(i) = h*c_pert(i)*c_pert(i)
    else
      e_2_pert(i) = 0.d0
      c_pert(i) = 0.d0
      H_pert_diag(i) = 0.d0
    endif
  enddo

end

subroutine pt2_moller_plesset ($arguments)
  use bitmasks
  implicit none
  $declarations
  
  BEGIN_DOC
  ! compute the standard Moller-Plesset perturbative first order coefficient and second order energetic contribution
  !
  ! for the various n_st states.
  !
  ! c_pert(i) = <psi(i)|H|det_pert>/(difference of orbital energies) 
  !
  ! e_2_pert(i) = <psi(i)|H|det_pert>^2/(difference of orbital energies) 
  !
  END_DOC
  
  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: phase,delta_e,h
  double precision               :: i_H_psi_array(N_st)
  integer                        :: h1,h2,p1,p2,s1,s2
  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  call get_excitation(ref_bitmask,det_pert,exc,degree,phase,Nint)
  if (degree == 2) then
    call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
    delta_e = Fock_matrix_diag_mo(h1) + Fock_matrix_diag_mo(h2) - &
            (Fock_matrix_diag_mo(p1) + Fock_matrix_diag_mo(p2))
    delta_e = 1.d0/delta_e
  else if (degree == 1) then
    call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
    delta_e = Fock_matrix_diag_mo(h1) - Fock_matrix_diag_mo(p1) 
    delta_e = 1.d0/delta_e
  else
    delta_e = 0.d0
  endif

  call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det,psi_selectors_size,n_st,i_H_psi_array)
  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,n_st
    H_pert_diag(i) = h
    c_pert(i) = i_H_psi_array(i) *delta_e
    e_2_pert(i) = c_pert(i) * i_H_psi_array(i)
  enddo
  
end


subroutine pt2_epstein_nesbet_SC2_projected ($arguments)
  use bitmasks
  implicit none
  $declarations
  BEGIN_DOC
  ! compute the Epstein-Nesbet perturbative first order coefficient and second order energetic contribution
  !
  ! for the various N_st states, 
  ! 
  ! but  with the correction in the denominator 
  ! 
  ! comming from the interaction of that determinant with all the others determinants 
  ! 
  ! that can be repeated by repeating all the double excitations
  !
  ! : you repeat all the correlation energy already taken into account in CI_electronic_energy(1)
  ! 
  ! that could be repeated to this determinant.
  !
  ! In addition, for the perturbative energetic contribution you have the standard second order
  !
  ! e_2_pert = <psi_i|H|det_pert>^2/(Delta_E)
  !
  ! and also the purely projected contribution 
  !
  ! H_pert_diag = <HF|H|det_pert> c_pert
  END_DOC
  
  double precision               :: i_H_psi_array(N_st)
  integer                        :: idx_repeat(0:ndet)
  integer                        :: i,j,degree,l
  double precision               :: diag_H_mat_elem_fock,accu_e_corr,hij,h0j,h,delta_E
  double precision               :: repeat_all_e_corr,accu_e_corr_tmp,e_2_pert_fonda

  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)

  call i_H_psi_SC2(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array,idx_repeat)
  accu_e_corr = 0.d0
  !$IVDEP
  do i = 1, idx_repeat(0)
   accu_e_corr = accu_e_corr + E_corr_per_selectors(idx_repeat(i))
  enddo
  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  h = h + accu_e_corr
  delta_E = 1.d0/(CI_SC2_electronic_energy(1) - h)


  c_pert(1) = i_H_psi_array(1) /(CI_SC2_electronic_energy(1) - h)
  e_2_pert(1) = i_H_psi_array(1) * c_pert(1)

  do i =2,N_st
    H_pert_diag(i) = h
    if  (dabs(CI_SC2_electronic_energy(i) - h) > 1.d-6) then
      c_pert(i) = i_H_psi_array(i) / (-dabs(CI_SC2_electronic_energy(i) - h))
      e_2_pert(i) = (c_pert(i) * i_H_psi_array(i))
    else
      c_pert(i) = i_H_psi_array(i)
      e_2_pert(i) = -dabs(i_H_psi_array(i))
    endif
  enddo

  degree = popcnt(xor( ref_bitmask(1,1), det_pert(1,1))) +                      &
           popcnt(xor( ref_bitmask(1,2), det_pert(1,2)))
  !DEC$ NOUNROLL
  do l=2,Nint
    degree = degree+ popcnt(xor( ref_bitmask(l,1), det_pert(l,1))) +            &
                     popcnt(xor( ref_bitmask(l,2), det_pert(l,2)))
  enddo
  if(degree==4)then
  ! <psi|delta_H|psi>
   e_2_pert_fonda = e_2_pert(1) 
   H_pert_diag(1) = e_2_pert(1) * c_pert(1) * c_pert(1)
   do i = 1, N_st
    do j = 1, idx_repeat(0)
     e_2_pert(i) += e_2_pert_fonda * psi_selectors_coef(idx_repeat(j),i) * psi_selectors_coef(idx_repeat(j),i)
    enddo
   enddo
  endif
  
end


subroutine pt2_epstein_nesbet_SC2_no_projected ($arguments)
  use bitmasks
  implicit none
  $declarations
  BEGIN_DOC
  ! compute the Epstein-Nesbet perturbative first order coefficient and second order energetic contribution
  !
  ! for the various N_st states, 
  ! 
  ! but  with the correction in the denominator 
  ! 
  ! comming from the interaction of that determinant with all the others determinants 
  ! 
  ! that can be repeated by repeating all the double excitations
  !
  ! : you repeat all the correlation energy already taken into account in CI_electronic_energy(1)
  ! 
  ! that could be repeated to this determinant.
  !
  ! In addition, for the perturbative energetic contribution you have the standard second order
  !
  ! e_2_pert = <psi_i|H|det_pert>^2/(Delta_E)
  !
  ! and also the purely projected contribution 
  !
  ! H_pert_diag = <HF|H|det_pert> c_pert
  END_DOC
  
  double precision               :: i_H_psi_array(N_st)
  integer                        :: idx_repeat(0:ndet)
  integer                        :: i,j,degree,l
  double precision               :: repeat_all_e_corr,accu_e_corr_tmp,e_2_pert_fonda
  double precision               :: diag_H_mat_elem_fock,accu_e_corr,hij,h0j,h,delta_E

  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)

  call i_H_psi_SC2(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array,idx_repeat)
  accu_e_corr = 0.d0
  !$IVDEP
  do i = 1, idx_repeat(0)
   accu_e_corr = accu_e_corr + E_corr_per_selectors(idx_repeat(i))
  enddo
  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  h = h + accu_e_corr
  delta_E = 1.d0/(CI_SC2_electronic_energy(1) - h)


  c_pert(1) = i_H_psi_array(1) /(CI_SC2_electronic_energy(1) - h)
  e_2_pert(1) = i_H_psi_array(1) * c_pert(1)

  do i =2,N_st
    H_pert_diag(i) = h
    if  (dabs(CI_SC2_electronic_energy(i) - h) > 1.d-6) then
      c_pert(i) = i_H_psi_array(i) / (-dabs(CI_SC2_electronic_energy(i) - h))
      e_2_pert(i) = (c_pert(i) * i_H_psi_array(i))
    else
      c_pert(i) = i_H_psi_array(i)
      e_2_pert(i) = -dabs(i_H_psi_array(i))
    endif
  enddo
end





subroutine pt2_epstein_nesbet_sc2 ($arguments)
  use bitmasks
  implicit none
  $declarations
  BEGIN_DOC
  ! compute the standard Epstein-Nesbet perturbative first order coefficient and second order energetic contribution
  !
  ! for the various N_st states, but with the CISD_SC2 energies and coefficients
  !
  ! c_pert(i) = <psi(i)|H|det_pert>/( E(i) - <det_pert|H|det_pert> )
  !
  ! e_2_pert(i) = <psi(i)|H|det_pert>^2/( E(i) - <det_pert|H|det_pert> )
  !
  END_DOC
  
  integer                        :: i,j
  double precision               :: i_H_psi_array(N_st)
  double precision               :: diag_H_mat_elem_fock, h
  PROVIDE  selection_criterion

  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  !call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array)
  call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)

  
  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if(CI_SC2_electronic_energy(i)>h.and.CI_SC2_electronic_energy(i).ne.0.d0)then
      c_pert(i) = -1.d0
      e_2_pert(i) = selection_criterion*selection_criterion_factor*2.d0
    else if  (dabs(CI_SC2_electronic_energy(i) - h) > 1.d-6) then
        c_pert(i) = i_H_psi_array(i) / (CI_SC2_electronic_energy(i) - h)
        H_pert_diag(i) = h*c_pert(i)*c_pert(i)
        e_2_pert(i) = c_pert(i) * i_H_psi_array(i)
    else
      c_pert(i) = -1.d0
      e_2_pert(i) = -dabs(i_H_psi_array(i))
      H_pert_diag(i) = h
    endif
  enddo
  
end



SUBST [ arguments, declarations ]

det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist ;

    integer, intent(in)             :: Nint
    integer, intent(in)             :: ndet
    integer, intent(in)             :: N_st
    integer, intent(in)             :: N_minilist
    integer(bit_kind), intent(in)   :: det_ref (Nint,2)
    integer(bit_kind), intent(in)   :: det_pert(Nint,2)
    double precision , intent(in)   :: fock_diag_tmp(2,mo_tot_num+1)
    double precision , intent(out)  :: c_pert(N_st)
    double precision , intent(out)  :: e_2_pert(N_st)
    double precision, intent(out)   :: H_pert_diag(N_st)
    integer, intent(in)             :: idx_minilist(0:N_det_selectors)
    integer(bit_kind), intent(in)   :: minilist(Nint,2,N_det_selectors)
;;


END_TEMPLATE

! Note : If the arguments are changed here, they should also be changed accordingly in
! the perturbation.template.f file.

