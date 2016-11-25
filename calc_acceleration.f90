subroutine calc_acceleration(position,velocity,acceleration)
! Routine drives calculation of all different acceleration terms

use nbodydata,only: N
implicit none
real,dimension(3,N), intent(in) :: position,velocity
real,dimension(3,N), intent(out) :: acceleration

real,dimension(N) :: tmig

acceleration(:,:) = 0.0

tmig(:) = 1.0e4

call calc_grav_acceleration(position,acceleration)
call calc_migration_drag(velocity,acceleration,tmig)

end subroutine calc_acceleration
