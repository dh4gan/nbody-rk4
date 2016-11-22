subroutine output(counter)
! Outputs data to file
! Currently writes each particle to separate file
! TODO - implement option of body output to one file

use nbodydata

implicit none
integer, intent(in) :: counter


call gravforce(pos,acc)

! Either output timestep as a single snapshot
if(snapshots=='y') then
   write(fileno,filenumformat) counter

   outputfile = TRIM(outputprefix)//"."//TRIM(fileno)
   open(isnap,file=outputfile,form='formatted')

   do ibody=1,N
      write(isnap,*) t, pos(:,ibody), vel(:,ibody), acc(:,ibody), semimaj(ibody),ecc(ibody),inc(ibody)
    call flush(isnap)
 enddo

   
else

   ! Either output individual bodies to separate files

   do ibody=1,N
      write(ibody+ilog,*) t, pos(:,ibody), vel(:,ibody), acc(:,ibody), semimaj(ibody),ecc(ibody),inc(ibody)
      call flush(ibody)
   enddo

endif

! Write log file containing global simulation data

write(ilog,*) t, dt, maxerror/tolerance


end subroutine output
