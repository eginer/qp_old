subroutine run_mrcc
  implicit none
  call set_generators_bitmasks_as_holes_and_particles
  call mrcc_iterations
end

subroutine mrcc_iterations
  implicit none
  
  integer :: i,j

  double precision :: E_new, E_old, delta_e
  integer :: iteration,i_oscillations
  double precision :: E_past(4)
  E_new = 0.d0
  delta_E = 1.d0
  iteration = 0
  j = 1
  i_oscillations = 0
  do while (delta_E > 1.d-10)
    iteration += 1
    print *,  '===========================' 
    print *,  'MRCC Iteration', iteration
    print *,  '===========================' 
    print *,  ''
    E_old = sum(ci_energy_dressed)
    call write_double(6,ci_energy_dressed(1),"MRCC energy")
    call diagonalize_ci_dressed
    E_new = sum(ci_energy_dressed)
    delta_E = dabs(E_new - E_old)

    E_past(j) = E_new
    amplitudes_before = amplitudes 
    j +=1
    if(j>4)then
     j=1
    endif
    if(iteration > 4) then 
     if(delta_E > 1.d-10)then
      if(dabs(E_past(1) - E_past(3)) .le. delta_E .and. dabs(E_past(2) - E_past(4)).le. delta_E)then 
       print*,'OSCILLATIONS !!!'
       oscillations = .True.
       i_oscillations +=1
       lambda_mrcc_tmp = lambda_mrcc
      endif
     endif
    endif
    if (i_oscillations > 5) then
     exit
    endif
    if (iteration > 200) then
      exit
    endif
  enddo
  call write_double(6,ci_energy_dressed(1),"Final MRCC energy")
  call ezfio_set_mrcc_cassd_energy(ci_energy_dressed(1))
  call save_wavefunction

end

subroutine set_generators_bitmasks_as_holes_and_particles
 implicit none
 integer :: i,k
 do k = 1, N_generators_bitmask
  do i = 1, N_int
   ! Pure single part 
   generators_bitmask(i,1,1,k) = holes_operators(i,1)   ! holes for pure single exc alpha 
   generators_bitmask(i,1,2,k) = particles_operators(i,1) ! particles for pure single exc alpha 
   generators_bitmask(i,2,1,k) = holes_operators(i,2)   ! holes for pure single exc beta 
   generators_bitmask(i,2,2,k) = particles_operators(i,2) ! particles for pure single exc beta 

   ! Double excitation 
   generators_bitmask(i,1,3,k) = holes_operators(i,1)   ! holes for first single exc alpha 
   generators_bitmask(i,1,4,k) = particles_operators(i,1) ! particles for first single exc alpha 
   generators_bitmask(i,2,3,k) = holes_operators(i,2)   ! holes for first single exc beta 
   generators_bitmask(i,2,4,k) = particles_operators(i,2) ! particles for first single exc beta 

   generators_bitmask(i,1,5,k) = holes_operators(i,1)   ! holes for second single exc alpha 
   generators_bitmask(i,1,6,k) = particles_operators(i,1) ! particles for second single exc alpha 
   generators_bitmask(i,2,5,k) = holes_operators(i,2)   ! holes for second single exc beta 
   generators_bitmask(i,2,6,k) = particles_operators(i,2) ! particles for second single exc beta 

  enddo
 enddo
 touch generators_bitmask



end
