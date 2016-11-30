subroutine calc_system_properties
! Calculates the various properties of the system
! Includes energy, angular momentum and orbits

use nbodydata

implicit none

real :: totalmass

if(heliocentric=='y') then
   call calc_heliocentric_frame(totalmass)
else
   call calc_centre_of_mass(totalmass)
endif
call angular_momentum
call energy

call calc_orbit_from_vector(totalmass)

end subroutine calc_system_properties

subroutine calc_heliocentric_frame(totalmass)
! Reset the system's frame so that the particle 1 is at the centre

use nbodydata


! Calculate centre of heliocentric frame
rcom(:) = pos(:,1)
vcom(:) = vel(:,1)
acom(:) = acc(:,1)

totalmass = sum(mass)

! Shift system into this frame frame
do ibody=1,N
pos(:,ibody) = pos(:,ibody)-rcom(:)
vel(:,ibody) = vel(:,ibody)-vcom(:)
acc(:,ibody) = acc(:,ibody)-acom(:)
enddo

end subroutine calc_heliocentric_frame


subroutine calc_centre_of_mass(totalmass)
! Calculate the orbital parameters of all the bodies
! Also compute energy and angular momentum for error tracking

use nbodydata

implicit none

real :: totalmass

! Calculate centre of mass
rcom(:) = 0.0
vcom(:) = 0.0
acom(:) = 0.0

totalmass = sum(mass)

do ibody=1,N
do ix=1,3
rcom(ix) = rcom(ix)+ mass(ibody)*pos(ix,ibody)
vcom(ix) = vcom(ix) + mass(ibody)*vel(ix,ibody)
acom(ix) = acom(ix) + mass(ibody)*acc(ix,ibody)
enddo

rcom(:) = rcom(:)/totalmass
vcom(:) = vcom(:)/totalmass
acom(:) = acom(:)/totalmass

enddo

! Shift system into centre of mass frame
do ibody=1,N
pos(:,ibody) = pos(:,ibody)-rcom(:)
vel(:,ibody) = vel(:,ibody)-vcom(:)
acc(:,ibody) = acc(:,ibody)-acom(:)
enddo

end subroutine calc_centre_of_mass



subroutine angular_momentum
! Calculate the angular momentum of all the bodies

use nbodydata

implicit none

! Compute angular momentum = r x v

angmom(1,:) = pos(2,:)*vel(3,:) - pos(3,:)*vel(2,:)
angmom(2,:) = pos(3,:)*vel(1,:) - pos(1,:)*vel(3,:)
angmom(3,:) = pos(1,:)*vel(2,:) - pos(2,:)*vel(1,:)
angmag(:) = sqrt(angmom(1,:)*angmom(1,:) + angmom(2,:)*angmom(2,:) + angmom(3,:)*angmom(3,:))

do ix=1,3
system_angmom(ix) = sum(angmom(ix,:))
enddo

system_ang =sqrt(system_angmom(1)*system_angmom(1) + system_angmom(2)*system_angmom(2)+system_angmom(3)*system_angmom(3))

if(initial_system_ang>1.0e-30)then
dL = (system_ang-initial_system_ang)/initial_system_ang
else
dL=0.0
endif

end subroutine angular_momentum






subroutine energy
! Calculate the orbital parameters of all the bodies
! Also compute energy and angular momentum for error tracking

use nbodydata

implicit none

real :: relpos
real,dimension(N) :: vmag
real,dimension(3) :: sep

! Compute kinetic energy of the bodies

vmag(:) = sqrt(vel(1,:)*vel(1,:) + vel(2,:)*vel(2,:) + vel(3,:)*vel(3,:))
ekin(:) = 0.5*vmag(:)*vmag(:)

! Compute potential energy

do ibody=1,N

epot(ibody)=0.0

do jbody=1,N

if(ibody==jbody) cycle

do ix=1,3
sep(ix) = pos(ix,ibody) - pos(ix,jbody)
enddo

relpos = sqrt(sep(1)*sep(1) + sep(2)*sep(2)+sep(3)*sep(3) +rsoft*rsoft)

do ix=1,3
epot(ibody) = epot(ibody) - G*mass(jbody)/(relpos)
enddo
enddo

enddo

etot(:) = ekin(:) + epot(:)

system_energy = sum(etot)

if(initial_system_energy >1.0e-30) then
dE = (system_energy-initial_system_energy)/initial_system_energy
else
dE = 0.0
endif

end subroutine energy





