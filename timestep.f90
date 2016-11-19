subroutine timestep
! Adjusts the timestep based on a standard step doubling algorithm
! Takes two half timesteps and compares to the current result

use nbodydata

real :: halfdt

halfdt = dt/2

real, dimension(3,N) :: testpos1,testvel1,testpos2,testvel2

call integrate(halfdt,pos,vel,testpos1,testvel1)
call integrate(halfdt,testpos1,testvel1,testpos2,testvel2)

! Compare double step result with single step result


end subroutine timestep
