


!      SUBROUTINE FUN( N, M, X, F, C, IFERR )
!C
!C  THIS SUBROUTINE COMPUTES
!C  THE OBJECTSIVE  FUNCTION VALUE AND THE VALUES OF  
!C  THE CONSTRAINTSFOR A GIVEN POINT X BY CALLING THE 
!C  CUTE INTERFACE SUBROUTINE CFN. 
!C  IF CFN RETURNS A NAN OF INF VALUE
!C  THE ROUTINE REPORTS AN ERROR BY SETTING IFERR=.TRUE.
!C
!
!
!   
!       DOUBLE PRECISION X(N), F, C(M)
!       LOGICAL          IFERR
!       INTEGER          N, M
!
!       INTEGER           I
!
!C
!C  COMPUTES THE VALUES OF FUNCTION AND CONSTRAINTS AND CHECKS IS 
!C  THE VALUES ARE NOT INF OR NAN
!C
!       IFERR=.FALSE.
!
!       CALL CFN( N , M , X , F , M, C)
!
!       IF (.NOT.( F .LT. 10.0D20 .AND. F .GT. -10.0D20 ))
!     +            IFERR = .TRUE.
!
!       DO 10 I=1, M
!         IF (.NOT.( C(I) .LT. 10.0D20 .AND. C(I) .GT. -10.0D20 ))
!     +              IFERR = .TRUE.
! 10    CONTINUE
!
!       RETURN
!       END




subroutine fmerit(m, val, objval, c, cu, cl, pp, method ) 

 double precision val, objval, c(*), cl(*), cu(*), pp
 integer          m, method


 integer          i
!
!  COMMON VARIABLES
!

!
!  PROBLEM CONTROL PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol 
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /
!
!  PP IS THE PENALTY PARAMETER
!

 val=objval
 if ( method .ne. 4 ) then 
   do 10 i =1, m
     if (cl(i)-cnstol.gt.c(i)) then
        val = val + pp*(cl(i)-c(i))
     elseif  (cu(i)+cnstol.lt.c(i)) then
        val = val + pp*(c(i)-cu(i))
     endif
10    continue
 else
   do 20 i =1, m
     if (cl(i)-cnstol.gt.c(i)) then
        val = val + pp*(cl(i)-c(i))**2
     elseif  (cu(i)+cnstol.lt.c(i)) then
        val = val + pp*(c(i)-cu(i))**2
     endif
20    continue
 endif
 return
 end









