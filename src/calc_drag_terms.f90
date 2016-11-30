subroutine calc_drag_terms(position,velocity,acceleration)
! Calculates approximation to 3 drag terms given a migration timescale
! Migration drag at timescale tmig
! Eccentricity damping on timescale tmig/10
! Inclincation damping on timescale tmig/10
! (cf Alibert et al 2013)

use nbodydata
implicit none

real,dimension(3,N),intent(in) :: position,velocity
real,dimension(3,N),intent(inout) :: acceleration
real,dimension(N) :: vdotr,rmag

! Migration drag first

do ix=1,3   
   where(tmig(:) > small)
   acceleration(ix,:) = acceleration(ix,:)-velocity(ix,:)/(2.0*tmig(:))
   endwhere
enddo

! Eccentricity damping
vdotr(:) =  0.0
do ix=1,3
   vdotr(:) = vdotr(:)+ position(ix,:)*velocity(ix,:)
enddo

rmag(:) = 0.0
rmag(:) = position(1,:)*position(1,:) + &
     position(2,:)*position(2,:) + &
     position(3,:)*position(3,:)

do ix=1,3
where(tmig(:)*rmag(:)>small)
acceleration(ix,:) = acceleration(ix,:) - 2.0*dampfac*vdotr(:)*position(ix,:)/(rmag(:)*tmig(:))
endwhere
enddo

! Inclination damping

where(tmig(:) > small)
acceleration(3,:) = acceleration(3,:) - 2.0*dampfac*vel(3,:)/tmig(:)
endwhere

end subroutine calc_drag_terms
