module nbodydata

!
! Unit data
!

real, parameter :: G=1
real, parameter :: msol = 1.99e33
real, parameter :: AU = 1.496e13

!
! Integers
!

integer :: N,ibody,jbody,ix

! Reals

real :: t, dt, tsnap, tdump, tend

! Body data

real, allocatable, dimension(:,:) :: pos,vel,acc

real,allocatable,dimension(:) :: mass, r,semimaj,ecc,inc,longascend,argper,longper

! Parameter filename

character(16),parameter :: paramfile = 'nbody_rk4.params'
character(len=6) :: snapshotformat
character(len=100) :: outputfile
character(len=50) :: outputprefix

end module nbodydata
