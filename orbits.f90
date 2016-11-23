subroutine orbits
! Calculate the orbital parameters of all the bodies
! Also compute energy and angular momentum for error tracking

use nbodydata

implicit none

real :: totalmass, relpos, gravparam
real,dimension(3) :: sep
real,dimension(3,N) :: eccvector
real,dimension(N) :: vdotr,rmag,vmag


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

gravparam = G*totalmass

! Compute angular momentum = r x v

angmom(1,:) = pos(2,:)*vel(3,:) - pos(3,:)*vel(2,:)
angmom(2,:) = pos(3,:)*vel(1,:) - pos(1,:)*vel(3,:)
angmom(3,:) = pos(1,:)*vel(2,:) - pos(2,:)*vel(1,:)

! Compute kinetic energy of the bodies

rmag(:) =sqrt(pos(1,:)*pos(1,:) + pos(2,:)*pos(2,:) + pos(3,:)*pos(3,:))
vmag(:) = sqrt(vel(1,:)*vel(1,:) + vel(2,:)*vel(2,:) + vel(3,:)*vel(3,:))
angmag(:) = sqrt(angmom(1,:)*angmom(1,:) + angmom(2,:)*angmom(2,:) + angmom(3,:)*angmom(3,:))

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

do ix=1,3
system_angmom(ix) = sum(angmom(ix,:))
enddo

system_ang =sqrt(system_angmom(1)*system_angmom(1) + system_angmom(2)*system_angmom(2)+system_angmom(3)*system_angmom(3))


if(initial_system_energy >1.0e-30) then
dE = (system_energy-initial_system_energy)/initial_system_energy
else
dE = 0.0
endif

if(initial_system_ang>1.0e-30)then
dL = (system_ang-initial_system_ang)/initial_system_ang
else
dL=0.0
endif

! Calculate orbital parameters - a,e,i

vdotr(:) = 0.0
do ix=1,3
vdotr(:) = vdotr(:) + pos(ix,:)*vel(ix,:)
enddo

do ix=1,3
eccvector(ix,:) = (vmag(:)*vmag(:)*pos(ix,:) -vdotr(:)*vel(ix,:))/gravparam - pos(ix,:)/rmag(:)
enddo

ecc(:) = sqrt(eccvector(1,:)*eccvector(1,:) + eccvector(2,:)*eccvector(2,:) + eccvector(3,:)*eccvector(3,:))

semimaj(:) = angmag(:)*angmag(:)/(gravparam*(1.0- ecc(:)*ecc(:)))

inc(:) = 0.0

do ibody=1,N

if(angmag(ibody)<1.0e-15) cycle
    inc(ibody) = acos(angmom(3,ibody)/ angmag(ibody))
enddo

end subroutine orbits




