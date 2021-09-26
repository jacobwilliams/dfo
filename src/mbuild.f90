subroutine mdbld( kappa , g   , h   , n    , nind , pntint, &
                  valint, poly, base, varnt, lpoly, lptint, &
                  neqcon, wrk , lwrk, iwrk , liwrk )  


!    
!  ***********************************************************************
!  THIS FUNCTION BUILDS A QUADRATIC INTERPOLATION
!  FUNCTION FOR THE SUPPLIED 'NIND' DATA POINTS THAT ARE
!  GIVEN IN ARRAY 'PNTINT' WITH VALUES IN ARRAY 'VALINT'.
!
!                     T         T
!      M(X)= KAPPA + G X + 0.5*X H X
!
!  IF 'VARNT' EQUALS 1 AND THE MODEL IS INCOMPLETE QUADRATIC
!  THEN MINIMUM FROBENIUS MODEL IS BUILT, OTHERWISE THE MODEL
!  IS BUILT FROM A LINEAR COMBINATION OF CONSTRUCTED NEWTON 
!  FUNDAMENTAL POLYNOMIALS. 
!
!  PARAMETERS:
!  
!   N      (INPUT)  DIMENSION OF THE PROBLEM
! 
!   NIND   (INPUT)  NUMBER OF INTERPOLATION POINTS
!
!   PNTINT (INPUT)  LIST OF  'NIND' DATA POINTS. THE I-TH POINT OCCUPIES
!                   POSITIONS ( I - 1 ) * N + 1 TO I * N.
!   VALINT (INPUT)  VALUES OF THE OBJECTIVE FUNCTION AT THE 'NIND'
!                   DATA POINTS CONTAINED IN PNTINT.
!   POLY   (INPUT)  THE ARRAY CONTAINING COEFFICIENTS OF NEWTON FUNDAMENTAL
!                   POLYNOMIALS (AS COMPUTED BY 'NBUILD')
!   BASE   (INPUT)  THE INDEX (IN PNTINT) OF THE CURRENT BASE (BEST)  POINT.
!
!   VARNT  (INPUT)  VARNT=1 -- FROBENIUS MODEL IS BUILT, VARNT=2 -- SUBQUADRATIC
!                   MODEL IS BUILT                    
!   KAPPA  (OUTPUT) THE CONSTANT TERM OF THE INTERPOLATION MODEL.
!
!   G      (OUTPUT) VECTOR OF THE LINEAR TERMS OF THE  INTERPOLATION MODEL.
!
!   H      (OUTPUT) MATRIX OF QUADRATIC TERMS OF THE  INTERPOLATION MODEL.
!
!   WRK    (INPUT)  REAL WORKSPACE OF LENGTH LWRK.
!
!   IWRK   (INPUT)  INTEGER WORKSPACE OF LENGTH LIWRK.
!
!  **************************************************************************
!

integer          n,       nind, lwrk  , liwrk, base, neqcon, &
                 iwrk( liwrk ), lptint, lpoly, varnt



double precision pntint( lptint ), kappa, g( n ), h( n, n ), &
                 valint( nind )  , wrk( lwrk )  , poly(lpoly)



integer          iout, iprint

double precision mcheps, cnstol

common /dfocm/   iout, iprint, mcheps, cnstol
save /dfocm/
!
!  SUBROUTINES AND FUNCTIONS CALLED:  
!
!  APPLICATION:  MBLDMF, MBLDNP 
!


!
!  LOCAL VARIABLES
!

integer          ilafla, imatr, ix, ivx, lenw, icurw, dd

logical          froben, fail


dd=(n+1)*(n+2)/2


!
!  SEE IF FROBENIUS MODEL SHOULD BE COMPUTED
!

froben=.false.
if (varnt==1 .and. nind > n+1 .and. nind<dd) froben=.true.


!
!  PARTITION REAL WORKSPACE
!

if ( froben ) then
  imatr  = 1
  ix     = imatr+nind**2
  ivx    = ix+(nind-1)*n
  icurw  = ivx+nind-1
  lenw   = lwrk-icurw+1
else 
  ilafla = 1
  icurw  = ilafla+nind
  lenw   = lwrk-icurw+1
endif

!
!  CHECK IF THE WORKING SPACE IS SUFFICIENT
!

if ( lenw < 1 ) then
   if ( iprint >= 0 ) write( iout, 1100 ) -lenw+1
   stop
endif

if ( liwrk < 1 ) then
   if ( iprint >= 0 )  write( iout, 1200 ) -liwrk+1
   stop
endif

fail = .true.
!
!  IF APPLICABLE, BUILD MINIMUM FROBENIUS MODEL
!

if ( froben ) then
  call mbldmf( kappa  , g ,   h , pntint    , valint    , base, &
               wrk(ix), wrk(ivx), wrk(imatr), n         , nind, &       ! 
               neqcon , lptint  , fail      , wrk(icurw), lenw, &
               iwrk   , liwrk   )
endif



!
!  IF MINIMUM FROBENIUS MODEL WAS NOT BUILT, THEN BUILD MODEL BASED ON 
!  NEWTON FUNDAMENTAL POLYNOMIALS
!

if ( fail ) then

  ilafla = 1
  call mbldnp( kappa  ,  g,  h, pntint, valint, poly  , &
               wrk(ilafla),  n, nind  , neqcon, lptint, lpoly )

endif          
       
return
!
!  NON-EXECUTABLE STATEMENTS
!

1100  format( ' MDBLD: *** ERROR: LWRK TOO SMALL!' / &
        '             IT SHOULD BE AT LEAST ',i5 )
1200  format( ' MDBLD: *** ERROR: LIWRK TOO SMALL!' / &
        '             IT SHOULD BE AT LEAST ', i5 )
end
!


!*******************************************************************************

!    NEXT SUBROUTINE

!*******************************************************************************

subroutine mbldnp( kappa, g, h   , pntint, valint, poly, &
                   lafla, n, nind, neqcon, lptint, lpoly )



!    
!  ***********************************************************************
!  THIS FUNCTION BUILDS A QUADRATIC INTERPOLATION
!  FUNCTION FOR THE SUPPLIED 'NIND' DATA POINTS THAT ARE
!  GIVEN IN ARRAY 'PNTINT' WITH VALUES IN ARRAY 'VALINT'.
!
!                     T         T
!      M(X)= KAPPA + G X + 0.5*X H X
!
!  THE MODEL  IS BUILT AS A LINEAR COMBINATION OF 'NIND'  NEWTON 
!  FUNDAMENTAL POLYNOMIALS STORED IN ARRAY 'POLY'. 
!
!  PARAMETERS:
!  
!   N      (INPUT)  DIMENSION OF THE PROBLEM
!
!   NIND   (INPUT)  NUMBER OF INTERPOLATION POINTS
! 
!   PNTINT (INPUT)  LIST OF  'NIND' DATA POINTS. THE I-TH POINT OCCUPIES
!                   POSITIONS ( I - 1 ) * N + 1 TO I * N.
!   VALINT (INPUT)  VALUES OF THE OBJECTIVE FUNCTION AT THE 'NIND'
!                   DATA POINTS CONTAINED IN PNTINT.
!   POLY   (INPUT)  THE ARRAY CONTAINING COEFFICIENTS OF NEWTON FUNDAMENTAL
!                   POLYNOMIALS (AS COMPUTED BY 'NBUILD')
!   LAFLA           AN AUXILIARY ARRAY OF LENGTH 'N' FOR FINITE DIFFERENCES
!
!   NEQCON (INPUT)  NUMBER OF LINEARLY INDEP. EQUALITY CONSTRAINTS
! 
!   KAPPA  (OUTPUT) THE CONSTANT TERM OF THE INTERPOLATION MODEL.
!
!   G      (OUTPUT) VECTOR OF THE LINEAR TERMS OF THE  INTERPOLATION MODEL.
!
!   H      (OUTPUT) MATRIX OF QUADRATIC TERMS OF THE  INTERPOLATION MODEL.
!
!  **************************************************************************
!


 integer           nind, n, neqcon, lptint, lpoly

 double precision lafla(nind) , poly(lpoly), pntint(lptint), &
                  valint(nind), kappa, g(n), h(n,n)
                 

!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!    APPLICATION:       FD,  GTERMS,  
!    BLAS       :       DCOPY
!



!
!  COMPUTE FINITE DIFFERENCES FROM NEWTON POLYNOMIALS, POINTS AND VALUES
!

    call fd(lafla, pntint, valint, poly, n, nind, neqcon, lptint, &
            lpoly)

!
!  COMPUTE TERMS OF THE MODEL FROM FINITE DIFFERENCES, PNTINT AND NEWTON
!  POLYNOMIALS
!

    call mterms(kappa , g   , h, lafla , poly , &
                pntint, nind, n, neqcon, lpoly, lptint)

 

  return
  end






!*******************************************************************************

!    NEXT SUBROUTINE

!*******************************************************************************




 subroutine mbldmf( kappa, g, h, pntint, valint, base, y   , vy  , &
                    matr , n, q, neqcon, lpnt  , fail, wrk , lwrk, &    ! 
                    iwrk , liwrk )




!    
!  ***********************************************************************
!  THIS FUNCTION BUILDS A QUADRATIC INTERPOLATION
!  FUNCTION FOR THE SUPPLIED 'Q' DATA POINTS THAT ARE
!  GIVEN IN ARRAY 'PNTINT' WITH VALUES IN ARRAY 'VALINT'.
!
!                     T         T
!      M(X)= KAPPA + G X + 0.5*X H X
!
!  THE MODEL IS INCOMPLETE QUADRATIC ( 'Q' IS LESS THAN (N+1)(N+2)/2 )
!  AND 'H' HAS THE SMALLEST  FROBENIUS NORM AMOUNG ALL HESSIANS OF ALL 
!  POSSIBLE QUADRATIC MODELS   SATISFYING THE INTERPOLATION CONDITION
!  FOR 'PNTINT' AND 'VALINT'.
!
!  PARAMETERS:
!  
!   N      (INPUT)  DIMENSION OF THE PROBLEM
! 
!   Q      (INPUT)  NUMBER OF INTERPOLATION POINTS
!
!   PNTINT (INPUT)  LIST OF  'NIND' DATA POINTS. THE I-TH POINT OCCUPIES
!                   POSITIONS ( I - 1 ) * N + 1 TO I * N.
!   VALINT (INPUT)  VALUES OF THE OBJECTIVE FUNCTION AT THE 'NIND'
!                   DATA POINTS CONTAINED IN PNTINT.
!   BASE   (INPUT)  THE INDEX (IN PNTINT) OF THE CURRENT BASE (BEST) POINT.
! 
!   KAPPA  (OUTPUT) THE CONSTANT TERM OF THE INTERPOLATION MODEL.
!
!   G      (OUTPUT) VECTOR OF THE LINEAR TERMS OF THE  INTERPOLATION MODEL.
!
!   H      (OUTPUT) MATRIX OF QUADRATIC TERMS OF THE  INTERPOLATION  MODEL.
!
!  FAIL    (OUTPUT) INDICATES IF THE THE SUBROUTINE HAS FAILED.
!  
!   Y               AUXILIARY ARRAY OF SHIFTED POINTS
!
!   VY              AUXILIARY ARRAY OF SHIFTED VALUES
!
!   MATR            AUXILIARY MATRIX QxQ
!
!   WRK    (INPUT)  REAL WORKSPACE OF LENGTH LWRK.
!
!   IWRK   (INPUT)  INTEGER WORKSPACE OF LENGTH LIWRK.
!
!  **************************************************************************
!


integer          n, q, lpnt, base, lwrk, iwrk(lwrk), liwrk, neqcon

double precision g(n), h(n,n), kappa, pntint(lpnt), valint(q), &
                 matr(q,q)   , wrk(liwrk), y(lpnt-n), vy(q-1)
logical          fail

double precision ddot

external         ddot

intrinsic        abs

integer          iout, iprint

double precision mcheps, cnstol

common /dfocm/   iout, iprint, mcheps, cnstol


!
!  LOCAL VARIABLES
!

double precision one, half, zero

parameter      ( one = 1.0d0,  half=0.5d0, zero=0.0d0 )

double precision anorm, rcond

integer          i, j, kk, ny, info, jbase 

!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!    APPLICATION:       RZRMAT, RZRVEC 
!    LAPACK     :       DPOTRF, DPOCON, DTRSV,   
!    BLAS       :       DCOPY , DAXPY , DGER , DGEMV
!    FORTRAN    :       ABS
!



fail = .false.
!
!  CHECK SUFFICIENCY OF THE WORKSPACE
!    

if (lwrk < 3*q) then
  if ( iprint >= 0 ) write(iout, 1000) 3*q
  stop
endif

if (liwrk < q) then
  if ( iprint >= 0 ) write(iout, 2000) q
  stop
endif 


if ( iprint >= 3 )   write(iout, 8000)

!
!  SHIFT ALL THE POINT SO THAT THE BASE IS AT THE ORIGIN AND PUT THE
!  POINTS, EXCEPT THE BASE, IN Y, AND THEIR VALUES IN VY. IGNORE THE
!  DUMMY POINTS IF THERE ARE ANY
!

 jbase=(base-1)*n
 ny=0
 do 5 i = 1,q
   if (i /= base .and. (i <= n+1-neqcon .or. i>n+1) ) then
     do 4 j=1,n
       y(ny*n+j) = pntint((i-1)*n+j)-pntint(jbase+j) 
4      continue
     ny=ny+1
   endif 
5  continue

!
!  FORM THE MATRIX 'MATR' BY FORMULAS: IF P IS THE MATRIX WITH
!  THE INTERPOLATION POINT AS THE COLUMNS THEN
!
!  *******************************
!  *          T            2     *
!  *     M = P P,   M=0.5*M +M   *
!  *                             *
!  *******************************

 do 20 i=1,ny
   do 10 j=1,ny
    matr(i,j)=ddot( n, y((i-1)*n+1), 1, y((j-1)*n+1), 1 )
10    continue
20  continue

 do 40 i=1,ny
   do 30 j=1,ny
    matr(i,j)=half*matr(i,j)**2 + matr(i,j)
30   continue
40  continue


!
!  USING LAPACK ROUTINES FIND AN ESTIMATE OF THE CONDITION NUMBER OF 'MATR'
!  'MATR' IS SYMMETRIC, SO, FIRST, COMPUTE THE LL^T FACTORIZATION IN 'DPOTRF'
!  IF FACTORIZATION FAILS, THEN PRINT A MESSAGE AND QUIT
!   



 call dpotrf('U', ny, matr, q, info)


 if ( info /= 0 ) then
   if( iprint>2 ) write(iout, 4000) 
   fail = .true.
   return
 endif
 
 anorm = zero
 do 46 i = 1, ny
   do 45 j = i, ny
     anorm=anorm+abs(matr(i,j))
45    continue 
46  continue  

!
!  FIND THE ESTIMATE OF THE COND. NUMBER, IF FAIL, THEN QUIT
!
 call dpocon('U' , ny  , matr, q , anorm, rcond, &
              wrk, iwrk, info)

 if ( info /= 0 ) then
   if( iprint>2 ) write(iout, 5000)  
   fail = .true.
   return
 endif
 if ( iprint >= 3 )   write(iout, 8010) rcond
!
!  CHECK IF THE CONDITION NUMBER IS NOT TOO BIG
!  IF NOT, THEN FORM THE RIGHT HAND SIDE FOR THE LINEAR SYSTEM
!  AND SOLVE THE SYSTEM USING LAPACK ROUTINES AND FACTORIZATION
!  COMPUTED BY 'DPOTRF'. STORE THE SOLUTION IN 'VY'
!
 if (rcond >= mcheps*1.0d4) then
   kk=0
   do 50 i = 1,q
     if (i/=base .and. (i <= n+1-neqcon .or. i>n+1)) then
       kk=kk+1
       vy(kk) = valint(i)-valint(base) 
     endif
50    continue
   call dtrsv('U', 'T', 'N', ny, matr, q, vy, 1)
   call dtrsv('U', 'N', 'N', ny, matr, q, vy, 1)


!
!  FIND THE COEFFICIENT OF THE SHIFTED MODEL BY FORMULAS:
!
!  ************************************ 
!  *     G <- 0,  H <- 0                *
!  *     FOR ALL POINTS P_i IN 'Y'    *
!  *                              T   *
!  *       H <- H + VY(I)*[P_i(P_i) ]  *
!  *                                  *
!  *       G <- G + VY(I)*P_i          *
!  *                                  *
!  ************************************
!

   call rzrvec(g, n)
   call rzrmat(h, n, n)

   do 60 i=1, ny
     call daxpy(n,  vy(i), y((i-1)*n+1),1, g, 1)
     call dger(n, n, vy(i), y((i-1)*n+1), 1, y((i-1)*n+1), &
               1, h, n)
60    continue

!
!  SHIFT THE MODEL BACK TO INTERPOLATE THE ORIGINAL POINTS AND FIND  KAPPA
! 
!  ************************************
!  *  G <- G - H*P_base                *
!  *                               T  *
!  *  KAPPA <- VALINT(BASE)-(P_base) G *
!  ************************************
!
   call dgemv('N',n, n, -one, h, n, pntint(jbase+1), 1, one, g ,1)
   call dcopy( n, g, 1, wrk, 1 )
   call dgemv('N',n, n, half, h, n, pntint(jbase+1), 1, one,wrk,1)
   kappa=valint(base)- ddot(n, wrk, 1, pntint(jbase+1), 1)

 else 
!
!  IF THE CONDITION NUMBER IS TOO LARGE, THEN STOP, SINCE IT SHOULD NOT
!  BE HAPPENING
!
   
   if (iprint>2) write (iout,3000) rcond
   fail = .true.
 endif
 return

1000 format( ' MBLDMF: *** ERROR: LWRK TOO SMALL !' / &
        '             IT SHOULD BE AT LEAST ',i5 )
2000 format( ' MBLDMF: *** ERROR: LIWRK TOO SMALL !' / &
        '             IT SHOULD BE AT LEAST ', i5 )

3000 format(' MINIMUM FROBENIUS NORM MODEL FAILED!' / &
        ' TOO ILL-CONDITIONED',  d14.7, /)
4000 format(' MINIMUM FROBENIUS NORM MODEL FAILED!' / &
        ' FACTORIZATION CANNOT BE COMPUTED' )
5000 format(' MINIMUM FROBENIUS NORM MODEL FAILED!' / &
        ' CONDITION NUMBER CANNOT BE COMPUTED' )
8000 format( ' MBLDMF: Minimum frobenius norm model computed ',/)
8010 format( ' MBLDMF: Condition number: ', d14.7,/ )  
 end









