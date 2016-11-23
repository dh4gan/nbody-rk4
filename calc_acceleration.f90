subroutine calc_acceleration(position,velocity,acceleration)
! Routine drives calculation of all different acceleration terms

use nbodydata,only: N
implicit none
real,dimension(3,N), intent(in) :: position,velocity
real,dimension(3,N), intent(out) :: acceleration

acceleration(:,:) = 0.0

call calc_grav_acceleration(position,acceleration)

end subroutine calc_acceleration
