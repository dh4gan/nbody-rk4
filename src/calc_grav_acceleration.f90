subroutine calc_grav_acceleration(position,acceleration)
! Calculates gravitational force between all bodies (brute force)
! given input positions

use nbodydata
implicit none

real, dimension(3,N),intent(in) :: position
real, dimension(3,N), intent(out) :: acceleration

real :: relpos,magipos,magjpos
real, dimension(3) :: sep

do ibody=1,N

    acceleration(:,ibody)=0.0

    ! If we're in the heliocentric frame
    if(heliocentric=='y') then
       if(ibody==1) cycle  ! Skip forces on the central body


       magipos = sqrt(position(1,ibody)*position(1,ibody) + &
                position(2,ibody)*position(2,ibody)+ &
                position(3,ibody)*position(3,ibody))

       acceleration(:,ibody) = acceleration(:,ibody) - &
            G*(mass(1)+mass(ibody))*position(:,ibody)/(magipos*magipos*magipos)

    endif

    do jbody=1,N

       if(ibody==jbody) cycle

        do ix=1,3
            sep(ix) = position(ix,ibody) - position(ix,jbody)
        enddo

        relpos = sqrt(sep(1)*sep(1) + sep(2)*sep(2)+sep(3)*sep(3) +rsoft*rsoft)


        if(heliocentric=='y') then

           if(jbody==1) cycle ! Central star force already calculated above

           magjpos = sqrt(position(1,jbody)*position(1,jbody) + &
                position(2,jbody)*position(2,jbody)+ &
                position(3,jbody)*position(3,jbody))

           do ix=1,3
              acceleration(ix,ibody) = acceleration(ix,ibody) - &
                   G*mass(jbody)*sep(ix)/(relpos*relpos*relpos) - &
                   G*mass(jbody)*position(ix,jbody)/(magjpos*magjpos*magjpos)
           enddo

        else

        do ix=1,3
            acceleration(ix,ibody) = acceleration(ix,ibody) - &
                 G*mass(jbody)*sep(ix)/(relpos*relpos*relpos)
        enddo

        endif

    enddo

enddo

acceleration(:,:) = acceleration(:,:)*G

return
end subroutine calc_grav_acceleration
