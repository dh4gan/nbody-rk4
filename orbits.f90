subroutine orbits
! Calculate the orbital parameters of all the bodies
! Also compute energy and angular momentum for error tracking

use nbodydata

implicit none

real :: totalmass, relpos, gravparam
real :: rdotV, ndotR,ndotV,edotn,edotR, nmag
real,dimension(3) :: sep,nplane
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


! Eccentricity first - calculate eccentricity (Laplace-Runge-Lenz) Vector
vdotr(:) = 0.0
do ix=1,3
   vdotr(:) = vdotr(:) + pos(ix,:)*vel(ix,:)
enddo

do ix=1,3
   eccvector(ix,:) = (vmag(:)*vmag(:)*pos(ix,:) -vdotr(:)*vel(ix,:))/gravparam - pos(ix,:)/rmag(:)
enddo

ecc(:) = sqrt(eccvector(1,:)*eccvector(1,:) + eccvector(2,:)*eccvector(2,:) + eccvector(3,:)*eccvector(3,:))


! Semimajor axis
semimaj(:) = angmag(:)*angmag(:)/(gravparam*(1.0- ecc(:)*ecc(:)))

inc(:) = 0.0

! Calculate the orbit's angles

do ibody=1,N

   ! Inclination

   if(angmag(ibody)<small) cycle
   inc(ibody) = acos(angmom(3,ibody)/ angmag(ibody))

   ! Longitude of the Ascending Node

   if (inc(ibody) <small) then
      longascend(ibody) = 0.0

      nplane(1) = angmag(ibody)
      nplane(2) = 0.0
      nplane(3) = 0.0
      nmag = angmag(ibody)

   else

      nplane(1) = -angmom(2,ibody)
      nplane(2) = angmom(1,ibody);
      nplane(3) = 0.0;

      nmag = sqrt(nplane(1)*nplane(1) + nplane(2)*nplane(2) + nplane(3)*nplane(3));

      longascend(ibody) = acos(nplane(1) / nmag);

      if (nplane(2) < 0.0) longascend(ibody) = 2.0 * pi - longascend(ibody);

   endif


   ! Calculate true anomaly

   !If orbit circular, no inclination, then use the position vector itself

   if (ecc(ibody) < small .and. abs(inc(ibody)) < small) then

      trueanom(ibody) = acos(pos(1,ibody) / rmag(ibody));
      if (vel(1,ibody) < 0.0) trueanom(ibody) = twopi - trueanom(ibody);

      ! If orbit circular and inclination non-zero, then use the orbital plane vector
   else if (ecc(ibody) < small) then

      ndotR = nplane(1)*pos(1,ibody) + nplane(2)*pos(2,ibody) + nplane(3)*pos(3,ibody)
      ndotR = ndotR / (rmag(ibody) * nmag);

      ndotV = nplane(1)*vel(1,ibody) + nplane(2)*vel(2,ibody) + nplane(3)*vel(3,ibody)

      trueanom(ibody) = acos(ndotR);

      if (ndotV > 0.0) trueanom(ibody) = twopi - trueanom(ibody);

      ! For non-circular orbits use the eccentricity vector
   else

      edotR = eccvector(1,ibody)*pos(1,ibody) + eccvector(2,ibody)*pos(2,ibody) + eccvector(3,ibody)*pos(3,ibody)
      edotR = edotR / (rmag(ibody) * ecc(ibody));

      rdotV = vel(1,ibody)*pos(1,ibody) + vel(2,ibody)*pos(2,ibody) + vel(3,ibody)*pos(3,ibody)

      trueanom(ibody) = acos(edotR);

      if (rdotV < 0.0)trueanom(ibody) = twopi - trueanom(ibody);
   endif

   ! Finally, calculate the longitude of periapsis - first calculate the argument of periapsis

   if (ecc(ibody) > small) then

      edotn = eccvector(1,ibody)*nplane(1) + eccvector(2,ibody)*nplane(2) + eccvector(3,ibody)*nplane(3);
      edotn = edotn / (nmag * ecc(ibody));

      argper(ibody) = acos(edotn);
      if (eccvector(3,ibody) < 0.0) argper(ibody) = twopi - argper(ibody);

      longper(ibody) = argper(ibody) + longascend(ibody)
   else

      argper(ibody) = 0.0
      longper(ibody) = 0.0

   endif

enddo

 end subroutine orbits




