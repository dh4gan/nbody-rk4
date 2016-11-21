subroutine initial
! Initialises the N-Body system
! During test phase, will hardwire setup

use nbodydata

implicit none

integer :: nzeros
real :: nfiles
character(1) :: zerostring


! TODO - read setup from file 

snapshots = 'n'
outputprefix = 'test'
N = 2

allocate(pos(3,N),vel(3,N),acc(3,N))
allocate(newpos(3,N),newvel(3,N))
allocate(mass(N),r(N),semimaj(N),ecc(N),inc(N))

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

tend = 10.0
tsnap = 0.1

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
      
      open(ibody,file=outputfile, form="formatted")
   enddo
endif

end subroutine initial
