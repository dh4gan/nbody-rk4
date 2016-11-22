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
  !print*, " input parameters in ",paramfile


! Read in parameter file and setup bodies

call initial

! Begin loop

t = 0.0
dt = 1.0e-5
tdump = tsnap

! Output initial conditions
snapshotcounter = 0

call output(snapshotcounter)

do while(t<tend)

   print*, 't=',t
   call integrate(dt,pos,vel,newpos,newvel)

   ! Do a timestep check before updating position,velocity
   call timestep(newpos,newvel)

   ! IF timestep passes tolerance, update positions
   ! Otherwise, repeat step with new dt

   if(maxerror > tolerance ) cycle

   pos = newpos
   vel = newvel

   t = t + dt

   
   if (t>tdump) then
      snapshotcounter = snapshotcounter + 1
      call output(snapshotcounter)
      tdump = tdump + tsnap
   endif
   
enddo

! Once integration complete, close all files and exit
! TODO - make a separate subroutine for end of integrations
do ibody=1,N
    close(ibody)
enddo


END PROGRAM nbody_rk4
