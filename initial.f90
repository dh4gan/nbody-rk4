subroutine initial
! Initialises the N-Body system
! During test phase, will hardwire setup

use nbodydata

implicit none

integer :: nzeros
real :: nfiles
character(1) :: zerostring


! TODO - read setup from file 

print*, 'Reading inputs from ',trim(paramfile)

open(10,file=paramfile,form='formatted')
read(10,*) outputprefix
read(10,*) snapshots
read(10,*) tend
read(10,*) tsnap
read(10,*) tolerance
read(10,*) rsoft
read(10,*) N


print*, N, ' bodies to be integrated'
print*, 'Files will be written with prefix ', outputprefix
if(snapshots=='y') then
    print*, 'Output in the form of snapshots'
else
    print*, 'Output in the form of individual files for each particle'
endif

print*, 'Maximum Runtime: ',tend,' years'
print*, 'Output every ',tsnap,' years'

print*, 'RK4 tolerance: ',tolerance
print*, 'Softening Length: ',rsoft, ' AU'

tend = tend*twopi
tsnap = tsnap*twopi

!snapshots = 'n'
!outputprefix = 'test'
!N = 2
!tolerance = 1.0e-4
!tend = 100*6.283
!tsnap = 0.1
!rsoft = 1.0e-5


allocate(pos(3,N),vel(3,N),acc(3,N))
allocate(newpos(3,N),newvel(3,N))
allocate(angmom(3,N),angmag(N),ekin(N),epot(N),etot(N))
allocate(mass(N),r(N),semimaj(N),ecc(N),inc(N))
allocate(longascend(N),argper(N),longper(N), trueanom(N))

pos(:,:) = 0.0
vel(:,:) = 0.0
acc(:,:) = 0.0

newpos(:,:) = 0.0
newvel(:,:) = 0.0

mass(1) = 1.0
mass(2) = 3.0e-6

pos(:,1) = 0.0

pos(:,2) = 0.0
pos(1,2) = 1.0

vel(:,:) = 0.0
vel(2,2) = 1.0


! If the output format is individual bodies
if(snapshots=='y') then
   
   ! Predicted number of files, and resulting format

   nfiles = tend/tsnap
   nzeros = int(log10(nfiles)) +2
   write(zerostring, '(I1)')nzeros
   filenumformat = "(I"//TRIM(zerostring)//"."//TRIM(zerostring)//")"

else

   nfiles = N
   nzeros = int(log10(nfiles)) +2
   write(zerostring, '(I1)')nzeros
   filenumformat = "(I"//TRIM(zerostring)//"."//TRIM(zerostring)//")"

   ! Open output files
   do ibody=1,N
      write(fileno, filenumformat) ibody

      outputfile = TRIM(outputprefix)//"."//TRIM(fileno)
      
      open(ibody+ilog,file=outputfile, form="formatted")
   enddo
endif

! Open log file

open(ilog,file=TRIM(outputprefix)//'.log',form='formatted')

! Calculate initial energy, angular momentum and store

call calc_acceleration(pos,vel,acc)
call orbits

initial_system_ang = system_ang
initial_system_energy = system_energy

print*, 'Initial Total Energy = ', initial_system_energy
print*, 'Initial Total Angular Momentum (Magnitude): ',initial_system_ang

! Measure wall clock time
call system_clock(start_clock,clock_rate,clock_max)


end subroutine initial
