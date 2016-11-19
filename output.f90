subroutine output
! Outputs data to file
! Currently writes each particle to separate file
! TODO - implement option of body output to one file

use nbodydata

implicit none

do ibody=1,N
    write(ibody,*) t, pos(ibody,:), vel(ibody,:), acc(ibody,:), semimaj(ibody),ecc(ibody),inc(ibody)
    call flush(ibody)
enddo

end subroutine output
