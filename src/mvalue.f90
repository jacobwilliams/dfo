double precision function mvalue( n, s, g, h, wrk, lwrk )
! 
!  ***********************************************
!     Computes the value of the quadratic model
!
!          S' * G + 0.5 * S' * H * S
!  ***********************************************
!
integer          n, lwrk
double precision s( n ), g( n ), h( n, n ), wrk( lwrk )

!
!  LOCAL VARIABLES
!

double precision ddot 

external         ddot
! 
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!    BLAS            :       DSYMV, DDOT
!


call dsymv( 'U', n, 0.5d0, h, n, s, 1, 0.0d0, wrk, 1 )

mvalue = ddot( n, s, 1, g, 1 ) + ddot( n, s, 1, wrk, 1 )


return
end

