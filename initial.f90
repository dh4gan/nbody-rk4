subroutine initial
! Initialises the N-Body system
! During test phase, will hardwire setup

use nbodydata

implicit none

integer :: nzeros,k
real :: nfiles
character(1) :: zerostring
character(1) :: varformat


! Read the simulation setup from file

print*, 'Reading inputs from ',trim(paramfile)

open(10,file=paramfile,form='formatted')
read(10,*) outputprefix
read(10,*) snapshots
read(10,*) tend
read(10,*) tsnap
read(10,*) tolerance
read(10,*) rsoft
read(10,*) N
read(10,*) varformat


! Allocate arrays to hold particle data
allocate(pos(3,N),vel(3,N),acc(3,N))
allocate(newpos(3,N),newvel(3,N))
allocate(angmom(3,N),angmag(N),ekin(N),epot(N),etot(N))
allocate(mass(N),r(N),semimaj(N),ecc(N),inc(N),tmig(N))
allocate(longascend(N),argper(N),longper(N), trueanom(N))


pos(:,:) = 0.0
vel(:,:) = 0.0
acc(:,:) = 0.0

if(varformat=='o') then

else
    do ibody = 1,N
    read(10,*) mass(ibody), (pos(k,ibody),k=1,3), (vel(k,ibody),k=1,3), tmig(ibody)
    enddo
endif

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
call calc_system_properties

initial_system_ang = system_ang
initial_system_energy = system_energy

print*, 'Initial Total Energy = ', initial_system_energy
print*, 'Initial Total Angular Momentum (Magnitude): ',initial_system_ang

call sleep(1)

! Measure wall clock time
call system_clock(start_clock,clock_rate,clock_max)


end subroutine initial
