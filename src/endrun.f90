subroutine endrun
! Subroutine handles the end of the run - closes files etc

use nbodydata
implicit none

real :: elapsed_time

if(snapshots=='y') then

else
    do ibody=1,N
    close(ibody)
    enddo
endif

close(ilog)

call system_clock(end_clock,clock_rate,clock_max)

elapsed_time = real(end_clock-start_clock)/real(clock_rate)

write(*,'(A,I5,A)') 'Run complete on ', N, ' bodies'
write(*,'(A,1P,1E10.3,A)'), 'Elapsed time: ', elapsed_time, ' seconds'

end subroutine endrun
