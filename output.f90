subroutine output(snapshotcounter)
! Outputs data to file
! Currently writes each particle to separate file
! TODO - implement option of body output to one file

use nbodydata

implicit none
integer, intent(in) :: snapshotcounter

! Either output timestep as a single snapshot
if(snapshots=='y') then
   write(fileno,filenumformat) snapshotcounter

   outputfile = outputprefix//"."//fileno
   open(isnap,file=outputfile,form='formatted')

   do ibody=1,N
      write(isnap,*) t, pos(ibody,:), vel(ibody,:), acc(ibody,:), semimaj(ibody),ecc(ibody),inc(ibody)
    call flush(isnap)
 enddo

   
else

   ! Either output individual bodies to separate files

   do ibody=1,N
      write(ibody,*) t, pos(ibody,:), vel(ibody,:), acc(ibody,:), semimaj(ibody),ecc(ibody),inc(ibody)
      call flush(ibody)
   enddo

endif


end subroutine output
