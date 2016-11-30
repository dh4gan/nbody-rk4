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

real,parameter :: small = 1.0e-20
real,parameter :: dampfac = 10.0

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
real,allocatable,dimension(:) :: mass, ekin,epot,etot,angmag,tmig
real,allocatable,dimension(:) :: r,semimaj,ecc,inc,longascend,argper,longper,trueanom

! Parameter filename

character(100) :: paramfile = 'nbody_rk4.params'
character(1) :: snapshots,heliocentric
character(len=6) :: filenumformat
character(len=100) :: outputfile
character(len=50) :: outputprefix
character(len=10) :: fileno

contains


end module nbodydata
