subroutine output(counter)
! Outputs data to file
! Currently writes each particle to separate file

use nbodydata

implicit none
integer, intent(in) :: counter

102 format (1P,19E15.5)
103 format (1P, 7E15.5)

call calc_grav_acceleration(pos,acc)
call orbits

! Either output timestep as a single snapshot
if(snapshots=='y') then
   write(fileno,filenumformat) counter

   outputfile = TRIM(outputprefix)//"."//TRIM(fileno)
   open(isnap,file=outputfile,form='formatted')

   do ibody=1,N
write(isnap,102) t/twopi, pos(:,ibody), vel(:,ibody), acc(:,ibody),&
    semimaj(ibody),ecc(ibody),inc(ibody), &
    ekin(ibody),epot(ibody),etot(ibody),angmom(:,ibody)
    call flush(isnap)
 enddo

   
else

   ! Either output individual bodies to separate files

   do ibody=1,N
      write(ibody+ilog,102) t/twopi, pos(:,ibody), vel(:,ibody), acc(:,ibody),&
semimaj(ibody),ecc(ibody),inc(ibody), &
ekin(ibody),epot(ibody),etot(ibody),angmom(:,ibody)
      call flush(ibody)
   enddo

endif

! Write log file containing global simulation data

write(ilog,103) t, dt, maxerror/tolerance, system_energy, dE, system_ang, dL


end subroutine output
