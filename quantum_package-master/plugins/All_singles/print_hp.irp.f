program print_hp
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
 implicit none
 integer :: i,j
 integer :: number_of_holes,n_h
 integer :: number_of_particles,n_p,degree
 logical :: is_a_1h1p
 integer          :: exc(0:2,2,2)
 double precision :: phase
 integer :: h1,p1,h2,p2,s1,s2
 
 print*,'BITMASK !!!'
 print*,'core'
 call debug_det(core_bitmask, N_int)
 
 print*,'inact'
 call debug_det(inact_bitmask, N_int)
 print*,'active'
 call debug_det(cas_bitmask(1,1,1),N_int)
 print*,'virt'
 call debug_det(virt_bitmask, N_int)
 print*,''
 print*,''
 print*,''
 do i = 1, N_det
  n_h = number_of_holes(psi_det(1,1,i))
  n_p = number_of_particles(psi_det(1,1,i))
  call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,1),degree,N_int)
  call get_excitation(psi_det(1,1,1),psi_det(1,1,i),exc,degree,phase,N_int)
  call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  print*,''
  print*,' i =',i
  print*,'degree = ',degree
  print*,'amplitude =',psi_coef(i,1)/psi_coef(1,1)
  print*,'n_h,n_p = ',n_h,n_p
  if(degree == 1)then
   print*,'s1',s1
   print*,'h1,p1 = ',h1,p1
  elseif (degree ==2)then
   print*,'s1',s1
   print*,'h1,p1 = ',h1,p1
   print*,'s2',s2
   print*,'h2,p2 = ',h2,p2
  endif
  call debug_det(psi_det(1,1,i),N_int)
  print*,''
 enddo




end
