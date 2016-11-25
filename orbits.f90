subroutine calc_orbit_from_vector(totalmass)
! Calculate the orbital parameters of all the bodies
! Also compute energy and angular momentum for error tracking

use nbodydata

implicit none

real :: totalmass, gravparam
real :: rdotV, ndotR,ndotV,edotn,edotR, nmag
real,dimension(3) :: nplane
real,dimension(3,N) :: eccvector
real,dimension(N) :: vdotr,rmag,vmag

gravparam = G*totalmass

rmag(:) =sqrt(pos(1,:)*pos(1,:) + pos(2,:)*pos(2,:) + pos(3,:)*pos(3,:))
vmag(:) = sqrt(vel(1,:)*vel(1,:) + vel(2,:)*vel(2,:) + vel(3,:)*vel(3,:))

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

end subroutine calc_orbit_from_vector

subroutine calc_vector_from_orbit(totalmass)
! Calculates body's position and velocity from orbital data

use nbodydata

implicit none

real :: totalmass,gravparam
real,dimension(N) :: rmag,vmag,semilatusrectum

rmag(:) = semimaj(:) * (1.0 - ecc(:) * ecc(:)) / (1.0 &
+ ecc(:) * cos(trueanom(:)));

! 2. Calculate position vector in orbital plane */

pos(1,:) = rmag(:)*cos(trueanom(:));
pos(2,:) = rmag(:) * sin(trueanom(:));
pos(3,:) = 0.0;

! 3. Calculate velocity vector in orbital plane */
semilatusrectum(:) = abs(semimaj(:) * (1.0 - ecc(:) * ecc(:)));
gravparam = G * totalmass;

do ibody=1,N

if (semilatusrectum(ibody) > small) then

vmag(ibody) = sqrt(gravparam / semilatusrectum(ibody));

else

vmag(ibody) = 0.0;
endif
enddo

vel(1,:) = -vmag(:) * sin(trueanom(:));
vel(2,:) = vmag(:) * (cos(trueanom(:)) + ecc(:));
vel(3,:) = 0.0;

! 4. Begin rotations to correctly align the orbit
! Firstly, Rotation around z axis by -argument of Periapsis */

call rotate_Z(pos, N, -1 * argper);
call rotate_Z(vel, N, -1 * argper);

! Secondly, Rotate around x by -inclination */

call rotate_X(pos, N, -1 * inc);
call rotate_X(vel, N, -1 * inc);

! Lastly, Rotate around z by longitudeAscendingNode */

call rotate_Z(pos,N,-1 * longascend);
call rotate_Z(vel,N,-1 * longascend);

end subroutine calc_vector_from_orbit

subroutine rotate_X(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in),dimension(nrows) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(:,:) = vector(:,:)
where(abs(angle)>1.0e-20)
newvector(1,:) = vector(1,:)
newvector(2,:) = vector(2,:)*cos(angle(:)) - vector(3,:)*sin(angle(:));
newvector(3,:) = vector(2,:)*sin(angle(:)) + vector(3,:)*cos(angle(:));
endwhere

end subroutine rotate_X

subroutine rotate_Y(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in),dimension(nrows) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(:,:) = vector(:,:)
where(abs(angle)>1.0e-20)

newvector(1,:) = vector(1,:)*cos(angle(:)) + vector(3,:)*sin(angle(:));
newvector(2,:) = vector(2,:)
newvector(3,:) = -vector(1,:)*sin(angle(:)) + vector(3,:)*cos(angle(:));

endwhere

end subroutine rotate_Y

subroutine rotate_Z(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in),dimension(nrows) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(:,:) = vector(:,:)

where(abs(angle)>1.0e-20)
newvector(1,:) = vector(1,:)*cos(angle(:)) - vector(2,:)*sin(angle(:));
newvector(2,:) = vector(1,:)*sin(angle(:)) + vector(2,:)*cos(angle(:));
newvector(3,:) = vector(3,:)
endwhere

end subroutine rotate_Z





