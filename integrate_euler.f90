subroutine integrate(deltat, position,velocity,newposition,newvelocity)
! test routine - Integrates via Euler

use nbodydata
implicit none

real, intent(in) :: deltat
real, dimension(3,N), intent(in) :: position,velocity
real,dimension(3,N) :: acceleration
real,dimension(3,N),intent(out) :: newposition,newvelocity

call calc_grav_acceleration(position,acceleration)

newvelocity(:,:) = velocity(:,:) + acceleration*deltat
newposition(:,:) = position(:,:) + newvelocity*deltat


end subroutine integrate
