subroutine gravforce(position,acceleration)
! Calculates gravitational force between all bodies (brute force)
! given input positions

use nbodydata
implicit none

real, dimension(3,N),intent(in) :: position
real, dimension(3,N), intent(out) :: acceleration

real :: relpos
real, dimension(3) :: sep

do ibody=1,N

    acceleration(ibody,:)=0.0

    do jbody=1,N

        do ix=1,3
            sep(ix) = position(ix,ibody) - position(ix,jbody)
        enddo

        relpos = sqrt(sep(1)*sep(1) + sep(2)*sep(2)+sep(3)*sep(3))

        do ix=1,3
            acceleration(ix,ibody) = acceleration(ix,ibody) - mass(ibody)*mass(jbody)*sep(ix)/(relpos*relpos*relpos)
        enddo
    enddo

enddo

acceleration(:,:) = acceleration(:,:)*G

return
end subroutine gravforce
