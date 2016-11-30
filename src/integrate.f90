subroutine integrate(deltat, position,velocity,newposition,newvelocity)
! This subroutine drives the RK4 integration

use nbodydata
implicit none

real, intent(in) :: deltat
real, dimension(3,N), intent(in) :: position,velocity

real,dimension(3,N),intent(out) :: newposition,newvelocity
real,dimension(3,N) :: nextpos,nextvel
real,dimension(3,N) :: k1vel, k2vel, k3vel, k4vel
real,dimension(3,N) :: k1pos, k2pos, k3pos, k4pos

! Begin calculating k-coefficients (for velocity and position)

! First k-coeff for velocity = acceleration
call calc_acceleration(position,velocity,k1vel)

! First k-coeff for position = velocity
k1pos(:,:) = velocity(:,:)

! Second k-coeff for velocity = acceleration at r=pos + 0.5*dt*k1pos
nextpos = position(:,:) + 0.5*deltat*k1pos(:,:)
nextvel = velocity(:,:) + 0.5*deltat*k1vel(:,:)

call calc_acceleration(nextpos,nextvel,k2vel)

! Second k-coeff for position = vvelocity + 0.5*k1vel
k2pos(:,:) = velocity(:,:) + 0.5*deltat*k1vel(:,:)

! Third velocity k-coeff
nextpos = position(:,:) + 0.5*deltat*k2pos(:,:)
nextvel = velocity(:,:) + 0.5*deltat*k2vel(:,:)
call calc_acceleration(nextpos,nextvel,k3vel)

! Third position k-coeff
k3pos(:,:) = velocity(:,:) + 0.5*deltat*k2vel(:,:)

! Fourth velocity k-coeff
nextpos = position(:,:)+ deltat*k3pos(:,:)
nextvel = velocity(:,:) + deltat*k3vel(:,:)
call calc_acceleration(nextpos,nextvel,k4vel)

! Fourth position k-coeff
k4pos(:,:) = velocity(:,:) + deltat*k3vel(:,:)

! New position and velocity
newposition(:,:) = position(:,:) + (deltat/6.0)*(k1pos(:,:) + 2.0*k2pos(:,:) + 2.0*k3pos(:,:) + k4pos(:,:))
newvelocity(:,:) = velocity(:,:) + (deltat/6.0)*(k1vel(:,:) + 2.0*k2vel(:,:) + 2.0*k3vel(:,:) + k4vel(:,:))

end subroutine integrate
