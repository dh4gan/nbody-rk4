PROGRAM nbody_rk4
! Written 18/11/2016 by dh4gan
! This code implements a 4th order Runge Kutta integration of an N-Body system

use nbodydata

! Print code header
 IMPLICIT NONE

  !             Display header
  print*, " "
  print*, "-----------------------------------------------"
  print*, "     N-BODY INTEGRATOR (RUNGE KUTTA 4TH ORDER) "
  print*, "     Created by D.Forgan, 18th November 2016   "
  print*, "-----------------------------------------------"
  print*, " "
  print*, "-----------------------------------------------"
  call getarg(1,paramfile)

  if(paramfile=='') then
     print*, 'Parameter file name not found from command line'
     paramfile = 'nbody_rk4.params'
     print*, 'Reverting to default'
  endif
 
  print*, " input parameters to be read from ",trim(paramfile)
  print*, "-----------------------------------------------"
  call sleep(1)


! Read in parameter file and setup bodies

call initial

! Begin loop

t = 0.0
dt = 1.0e-3
tdump = tsnap

! Output initial conditions
snapshotcounter = 0
call output(snapshotcounter)

! Begin integration
do while(t<tend)

   call integrate(dt,pos,vel,newpos,newvel)

   ! Do a timestep check before updating position,velocity
   call timestep(newpos,newvel)

   ! IF timestep passes tolerance, update positions
   ! Otherwise, repeat step with new dt

   if(maxerror > tolerance) cycle

   pos = newpos
   vel = newvel

   t = t + dt

   if (t>tdump) then
         write(*,'(A,1P,2E12.3,A)'), 't, dt=',t/twopi,dt/twopi, ' years'
      snapshotcounter = snapshotcounter + 1
      call output(snapshotcounter)
      tdump = tdump + tsnap
   endif
   
enddo

! Once integration complete, close all files and exit

call endrun

END PROGRAM nbody_rk4
