subroutine integrate(deltat, position,velocity,newposition,newvelocity)
! This subroutine drives the RK4 integration

use nbodydata
implicit none

real, intent(in) :: deltat
real, dimension(3,N), intent(in) :: position,velocity

real,dimension(3,N),intent(out) :: newposition,newvelocity
real,dimension(3,N) :: nextpos
real,dimension(3,N) :: k1vel, k2vel, k3vel, k4vel
real,dimension(3,N) :: k1pos, k2pos, k3pos, k4pos
real :: rad

! Begin calculating k-coefficients (for velocity and position)

print*, 'Position: ',position(:,2)
print*, 'DELTAT: ', deltat
! First k-coeff for velocity = acceleration
call gravforce(position,k1vel)
print*, 'k1vel: ', k1vel(:,2)

! First k-coeff for position = velocity
k1pos(:,:) = velocity(:,:)
!print*, 'k1pos: ', k1pos(:,2)

! Second k-coeff for velocity = acceleration at r=pos + 0.5*dt*k1pos
nextpos = position(:,:) + 0.5*deltat*k1pos(:,:)
print*, 'nextpos: ',nextpos(:,2)
call gravforce(nextpos,k2vel)
print*, 'k2vel: ', k2vel(:,2)

! Second k-coeff for position = vvelocity + 0.5*k1vel
k2pos(:,:) = velocity(:,:) + 0.5*deltat*k1vel(:,:)
print*, 'k2pos: ', k2pos(:,2)
! Third velocity k-coeff
nextpos = position(:,:) + 0.5*deltat*k2pos(:,:)
print*, 'nextpos: ',nextpos(:,2)
call gravforce(nextpos,k3vel)
!print*, 'k3vel: ', k3vel(:,2)

! Third position k-coeff
k3pos(:,:) = velocity(:,:) + 0.5*deltat*k2vel(:,:)

! Fourth velocity k-coeff
nextpos = position(:,:)+ deltat*k3pos(:,:)
call gravforce(nextpos,k4vel)

!print*, 'k4vel: ', k4vel(:,2)

! Fourth position k-coeff
k4pos(:,:) = velocity(:,:) + deltat*k3vel(:,:)

! New position and velocity
newposition(:,:) = position(:,:) + (deltat/6.0)*(k1pos(:,:) + 2.0*k2pos(:,:) + 2.0*k3pos(:,:) + k4pos(:,:))
newvelocity(:,:) = velocity(:,:) + (deltat/6.0)*(k1vel(:,:) + 2.0*k2vel(:,:) + 2.0*k3vel(:,:) + k4vel(:,:))

!print*, 'new position: ', newposition(:,2)

rad = sqrt(newposition(1,2)*newposition(1,2) + newposition(2,2)*newposition(2,2) + newposition(3,2)*newposition(3,2))

!print*, rad

!STOP
end subroutine integrate
