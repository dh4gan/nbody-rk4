module nbodydata

!
! Unit data
!

real, parameter :: G=1
real, parameter :: msol = 1.99e33
real, parameter :: AU = 1.496e13
real, parameter :: pi = 3.14159265
real, parameter :: twopi = 2.0*pi

!
! Integers
!

integer,parameter :: isnap = 10
integer,parameter :: ilog = 2

integer :: N,ibody,jbody,ix,snapshotcounter

! Reals

real :: small = 1.0e-20
real :: t, dt, tsnap, tdump, tend
real :: maxerror, tolerance,rsoft
real :: system_ang, system_energy, initial_system_ang,initial_system_energy
real :: dE, dL
integer :: start_clock,end_clock,clock_rate,clock_max

! Body data

real, allocatable, dimension(:,:) :: pos,vel,acc
real, allocatable,dimension(:,:) :: newpos,newvel
real,allocatable,dimension(:,:) :: angmom

real,dimension(3) :: system_angmom,rcom,vcom,acom
real,allocatable,dimension(:) :: mass, ekin,epot,etot,angmag
real,allocatable,dimension(:) :: r,semimaj,ecc,inc,longascend,argper,longper,trueanom

! Parameter filename

character(16),parameter :: paramfile = 'nbody_rk4.params'
character(1) :: snapshots
character(len=6) :: filenumformat
character(len=100) :: outputfile
character(len=50) :: outputprefix
character(len=10) :: fileno

contains

subroutine rotate_X(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(1,:) = vector(1,:)
newvector(2,:) = vector(2,:)*cos(angle) - vector(3,:)*sin(angle);
newvector(3,:) = vector(2,:)*sin(angle) + vector(3,:)*cos(angle);

vector(:,:) = newvector(:,:)

end subroutine rotate_X

subroutine rotate_Y(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(1,:) = vector(1,:)*cos(angle) + vector(3,:)*sin(angle);
newvector(2,:) = vector(2,:)
newvector(3,:) = -vector(1,:)*sin(angle) + vector(3,:)*cos(angle);

vector(:,:) = newvector(:,:)

end subroutine rotate_Y

subroutine rotate_Z(vector,nrows,angle)
! Rotates 3 x nrows vector of particle data around x-axis by angle

implicit none
integer, intent(in)::nrows
real,intent(in) :: angle
real,intent(inout),dimension(3,nrows) :: vector

real, dimension(3,nrows) :: newvector

newvector(1,:) = vector(1,:)*cos(angle) - vector(2,:)*sin(angle);
newvector(2,:) = vector(1,:)*sin(angle) + vector(2,:)*cos(angle);
newvector(3,:) = vector(3,:)

vector(:,:) = newvector(:,:)

end subroutine rotate_Z

end module nbodydata
