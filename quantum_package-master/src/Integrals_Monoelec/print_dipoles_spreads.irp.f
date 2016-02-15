program print_dipoles
 implicit  none
 integer :: i,j
 print*,'Orbital      <Z>    <Z^2> - <Z>2'
 do i = 1, 2
  print*,i,mo_dipole_z(i,i),mo_spread_z(i,i) - mo_dipole_z(i,i)**2
 enddo

 print*,'Orbital      <x>    <x^2> - <x>2'
 do i = 1, 2
  print*,i,mo_dipole_x(i,i),mo_spread_x(i,i) - mo_dipole_x(i,i)**2
 enddo

 print*,'Orbital      <y>    <y^2> - <y>2'
 do i = 1, 2
  print*,i,mo_dipole_y(i,i),mo_spread_y(i,i) - mo_dipole_y(i,i)**2
 enddo


end
