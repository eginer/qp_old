program print_mulliken
 implicit none
 read_wf = .True.
 touch read_wf
 print*,'Mulliken spin densities'

 call test
end
subroutine test
 double precision :: accu
 integer :: i
 integer :: j
 accu= 0.d0
 do i = 1, nucl_num
  print*,i,nucl_charge(i),mulliken_spin_densities(i)
  accu += mulliken_spin_densities(i)
 enddo
 print*,'Sum of Mulliken SD = ',accu
!accu = 0.d0
!do i = 1, mo_tot_num
! accu += one_body_spin_density_mo(i,i)
!enddo
!print*,'Sum of the diagonal elements of the SD in MOs '
!print*,accu
!accu = 0.d0
!do i = 1, ao_num
! do j = 1, ao_num
!  accu += one_body_spin_density_ao(i,j) * ao_overlap(i,j)
! enddo
!enddo
!print*,'Sum of the SD in the AO basis'
!print*,accu
 print*,'AO SPIN POPULATIONS'
 accu = 0.d0
 do i = 1, ao_num
  accu += spin_gross_orbital_product(i)
  write(*,'(I3,X,A4,X,I2,X,A4,X,F10.7)')i,trim(element_name(int(nucl_charge(ao_nucl(i))))),ao_nucl(i),trim(l_to_charater(ao_l(i))),spin_gross_orbital_product(i)
 enddo
 print*,'sum = ',accu
 accu = 0.d0
 print*,'Angular momentum analysis'
 do i = 0,  ao_l_max
  accu += spin_population_angular_momentum(i)
  print*,trim(l_to_charater(i)),spin_population_angular_momentum(i)
 print*,'sum = ',accu
 enddo

end

