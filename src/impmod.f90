subroutine impmod( poly  , pntint, valint, points, values, obfval, &
                   conval, sp2in , in2sp , n, m  , &
                   nq    , nind  , base  , pivthr, xchthr, pivval, &
                   dist  , delta , x     , a     , lda   , nclin , &    !  
                   ncnln , lb    , ub    , scale , scal  , nf    , &    !   
                   maxnf , impr  , pp    , neqcon, wrk   , lwrk  , &    ! 
                   iwrk  , liwrk , delmin, method)


!
!  *********************************************************************
!  THIS SUBROUTINE ATTEMPTS TO IMPROVE THE MODEL BY:
!      1. ADDING A NEW POINT TO THE INTERPOLATION SET, POSSIBLY
!         DROPPING ANOTHER POINT FROM THERE, SO THAT THE SET
!         IS BETTER POISED
!      2. RECOMPUTING THE WHOLE INTERPOLATION SET FROM SCRATCH,
!         CHOOSING POSSIBLY LESS AND POSSIBLY DIFFERENT POINTS
!         THE SET. THE BASE POINT IS ALWAYS GUARANTEED TO BE INCLUDED
!         IN THE NEW INTERPOLATION SET
!  PARAMETERS
!
!    POLY    (INPUT/OUTPUT) THE ARRAY OF NEWTON POLYNOMIALS
!    
!    PNTINT  (INPUT/OUTPUT) THE INTERPOLATION SET
!
!    VALINT  (INPUT/OUTPUT) THE VALUES OF THE FUNCTION AT POINTS  IN 'PNTINT'
!
!    POINTS  (INPUT/OUTPUT) SET OF ALL 'POINTS' WITH FUNCTION VALUES
!
!    VALUES  (INPUT/OUTPUT) VALUES AT POINTS IN 'POINTS'
!
!    SP2IN   (INPUT/OUTPUT) ARRAY OF 0/1 INDICATING FOR EVERY POINT IN
!                          'POINTS' IF IT BELONGS TO 'PNTINT' OR NOT
!    IN2SP   (INPUT/OUTPUT) ARRAY OF INDICES INDICATING FOR EVERY POINT IN
!                          'PNTINT' ITS POSITION IN 'POINTS'
!    NIND    (INPUT/OUTPUT) CARDINALITY OF INTERPOLATION SET
!
!    DIST    (INPUT/OUTPUT) DISTANCE OF ALL POINT IN 'POINTS' TO THE  BASE
!
!    NQ      (INPUT/OUTPUT) NUMBER OF POINTS IN 'POINTS'
!
!    PIVVAL  (INPUT/OUTPUT) ARRAY OF THE PIVOT VALUES ASSOCIATED WITH  EVERY
!                           POINT IN 'PNTINT'
!    X       (INPUT/OUTPUT) CUMULATIVE SHIFT OF THE INTERPOLATION CENTER
!
!    NF      (INPUT/OUTPUT) TOTAL NUMBER OF FUNCTION CALLS
!
!    MAXNF   (INPUT)        MAXIMUM  NUMBER OF FUNCTION CALLS
!
!    N       (INPUT)        PROBLEM DIMENSION
!
!    BASE    (INPUT)        THE INDEX OF THE BASE POINT IN 'PNTINT'
!
!    DELTA   (INPUT)        TRUST REGION RADIUS
!  
!    A       (INPUT)        MATRIX OF LINEAR CONSTRAINTS OF THE PROBLEM
!
!    NCLIN   (INPUT)        NUMBER OF LINEAR CONSTRAINTS
!   
!    NCNLN   (INPUT)        NUMBER OF NONLINEAR CONSTRAINTS
!
!    SCALE   (INPUT)        FLAG, INDICATING IF THE PROBLEM IS SCALED
!
!    SCAL    (INPUT)        ARRAY OF N SCALING FACTORS
!
!    LB      (INPUT/OUTPUT) LOWER BOUNDS OF THE PROBLEM
!
!    UB      (INPUT/OUTPUT) UPPER BOUNDS OF THE PROBLEM
!
!    IMPR    (INPUT/OUTPUT) ON INPUT
!              <0          THEN THE POINT WITH SMALLEST VALUE IS NOT IN
!                           IN THE INTERPOLATION SET, IMPR=-NQ, WHERE 
!                           NQ IS THE INDEX OF THE GOOD POINT IN "POINTS"
!
!                           THE OUTPUT INFORMATION
!               1           A POINT WAS ADDED TO THE INTERPOLATION SET
!               2           A POINT WAS REPLACED  IN THE INTERP. SET
!               3           WHOLE INTERPOLATION WAS RECOMPUTED, SINCE TOO OLD
!               4           WHOLE INTERPOLATION WAS RECOMPUTED SINCE
!                           NOTHING ELSE WORKED
!               0           NO CHANGE/IMPROVEMENT CAN BE MADE
!
!    WRK                    WORKING REAL SPACE
!  **************************************************************************
!

!
!  SUBROUTINE PARAMETERS
!

double precision  poly(lpoly)      , pntint(lptint)  , &
                  valint(lvlint)   , points(lpnts+n ), &
                  values( nq+1)    , a(lda*n)        , pp    , &
                  pivval(lvlint)   , dist(nq+1)      , x(n)  , &
                  lb(*), wrk(lwrk) , delta , &
                  ub(*), pivthr    , xchthr, &
                  conval(lconvl+m) , obfval(nq+1)    , scal(n), &
                  delmin

integer           sp2in(nq+1), in2sp(lvlint), n, nq, nind , base , &    ! 
                  iwrk(liwrk), nclin        , ncnln, lda  , impr , &    !   
                  nf         , liwrk        , lwrk , scale, maxnf, &
                  m          , neqcon       , method


!
!  COMMON VARIABLES
!
            

!
!  PRINTOUT PARAMETERS
! 
integer           iout  , iprint
double precision  mcheps, cnstol
common / dfocm /  iout  , iprint, mcheps, cnstol
save / dfocm /

!
!  INTERPOLATION CONTROL PARAMETERS
!
integer           npmin, layer, effort 
common / opti  /  npmin, layer, effort
save / opti  /
!
!  LENGTHS OF ARRAYS
!

integer           lpoly, lpnts, lvalue, lptint, lvlint, lconvl
common / rpart /  lpoly, lpnts, lvalue, lptint, lvlint, lconvl
save / rpart /


!
!  LOCAL VARIABLES
!

integer          ixg, j, np1, dd, inform, jpoly, &
                 minmax, ixgnew , icurw , lenw 

logical          fail, iferr

double precision vnew, del


double precision dnrmnf
external         dnrmnf


!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!    APPLICATION:       PTNEW , PTREPL, NBUILD, FUN, GETDIS,
!                       UNSCL , UNSHFT, SHIFT , SCL, SWAPNP,
!                       DNRMNF
!    BLAS:              DCOPY
!


!
!  SET THE POINTER TO SPACE FOR POINT 'XGNEW' - NEW POINT COMPUTED
!  TO IMPROVE GEOMETRY
!
ixgnew = 1

!
!  CHECK IF REMAINING REAL SPACE IS SUFFICIENT
!
icurw  = ixgnew + n
lenw   = lwrk - icurw + 1
if ( lenw < 0 ) then
  write( iout, 1000 ) - lenw 
  stop
endif
!
!  INITIALIZATION
!
np1   = n + 1 
dd    = np1*(n+2)/2
fail  =.true.
iferr =.false.
del   = delta
jpoly = nind + 1
minmax= 0
!
!  IF NO MODEL IMPROVEMENT IS NEEDED BUT THE LEVEL OF EFFORT
!  REQUIRES THE INTERPOLATION SET TO BE RECOMPUTED, PROCEED
!  DIRECTLY TO THE NBUILD PROCEDURE
!

if ( impr == 10 ) goto 25

!
!  IF THE MODEL IS NOT TOO OLD AND IMPR IS >= -1, THEN WE FIRST TRY TO 
!  FIND A "GEOMETRY" POINT WHICH CAN BE INCLUDED IN THE INTERPOLATION
!  TO IMPROVE MODEL/GEOMETRY.
!  WE SAY THAT THE MODEL IS OLD IF THE FIRST INTERPOLATION POINT
!  IS FURTHER THAN LAYER*DELTA AWAY FROM THE BASE; I.E. SHOULD
!  BE EXCLUDED FROM THE INTERPOLATION
!  IF THERE ARE NOT MORE THAT "NPMIN" POINTS IN THE MODEL WE 
!  NEED TO INCLUDE ANOTHER  POINT, EVEN IF THE MODEL IS "OLD".
!

10 if ( ( dist(in2sp(1)) <= layer*delta .or. nind <= npmin  ) &
      .and. impr >= -1 ) then
  if ( iprint >= 2 ) write( iout, 8000 )
  call ptnew( wrk(ixgnew), ixg  , vnew , pntint, lptint, n     , &
              poly ,lpoly, nind , base , del   , lb    , ub    , &      ! 
              a    , lda , nclin, ncnln, pivval, pivthr, xchthr, &
              wrk(icurw) , lenw , iwrk , liwrk , fail  , inform, &
              minmax     , x    )
!
!  IF A NEW POINT WAS NOT FOUND AND THE INTERPOLATION IF SUBLINEAR,
!  THEN WE MUST HAVE LESS DEGREES OF FREEDOM THAN WE THINK,
!  WE RETURN AND STOP OPTIMIZATION, OTHERWISE WE MAY GO INTO
!  LOOP
!      
  if ( fail .and. nind < n + 1 - neqcon ) then
    impr = 0
    goto 25
  endif
endif 
!
!  IF IFERR=.TRUE. IS MEANS THAT A NEW GEOMETRY POINT WAS FOUND
!  BUT A FUNCTION COMPUTATION FAILED. NOW ANOTHER POINT WAS TRIED
!  FAIL=.TRUE. MEANS THAT THE NEW POINT DOES NOT SATISFY
!  PIVOT THRESHOLD
!  IF  MINMAX = 0, THIS MEANS THAT WE CALLED PTNEW SECOND TIME
!  FOR THE SAME POLYNOMIAL. THUS, FIRST TIME WAS SUCCESSFUL AND
!  FUNCTION VALUE FAILED. SO, IF SEARCH FAILED NOW, DO NOT GIVE
!  UP AND CALL PTNEW AGAIN WITH A NEW POLYNOMIAL OR WITH A SMALLER
!  RADIUS.
!

if ( fail .and. iferr .and. minmax == 0 ) then
  if ( jpoly < n+1 ) then
!
!  IF WE STILL HAVE NOT TRIED ALL POLYNOMIALS  IN THE LINEAR BLOCK, 
!  PICK THE NEXT POLYNOMIAL, PUT IT IN APPROPRIATE PLACE  AND GO BACK TO 'PTNEW'
!
    jpoly = jpoly+1
    call swapnp( n, jpoly, nind+1, poly, lpoly )
    if ( iprint >= 2 ) write( iout, 8020 ) jpoly
    minmax = 0
    goto 10
  else
!
!  IF ALL POLYNOMIALS IN THIS BLOCK WERE TRIED, REDUCE RADIUS OF SEARCH
!  BY HALF AND REPEAT
!
    del   = del/2.0d0
    if ( del > delmin ) then
      jpoly = nind + 1
      minmax = 0
      goto 10
    else
      goto 25
    endif
  endif
endif

!
!  IF A "GEOMETRY" POINT WAS FOUND SUCCESSFULLY THEN 'IXG' IS THE
!  POSITION WHERE THE POINT SHOULD BE PLACED IN THE INTERPOLATION.
!


if ( .not. fail ) then 


!
!  COMPUTE THE FUNCTION VALUE AT THE NEW "GEOMETRY" POINT AND
!  IF THE FUNCTION EVALUATION FAILS, REDUCE RADIUS TO HALF
!  OF THE DISTANCE BETWEEN "GEOMETRY" POINT AND THE BASE AND
!  TRY TO COMPUTE ANOTHER "GEOMETRY" POINT.
!
  fail = .true.
  nf = nf + 1
  if ( nf > maxnf ) return
  call unshft(n, x, wrk(ixgnew))
  if ( scale /= 0 ) call unscl( n, wrk(ixgnew), scal )
  call fun( n, m, wrk(ixgnew), obfval(nq+1), conval(nq*m+1), &
            iferr)
  if ( scale /= 0 ) call scl( n, wrk(ixgnew), scal )
  call shift(n, x, wrk(ixgnew))
  if (iferr) then
    if ( iprint >= 2 ) write( iout, 8010 )
!
!  IF FUNCTION VALUE COULD NOT BE COMPUTED AT THE NEW GEOMETRY POINT
!  WE TRY TO FIND ANOTHER POINT
!
    if ( nind < n + 1 - neqcon ) then
!
!  IF WE DO NOT HAVE A FULLY LINEAR MODEL, THEN WE DO A BIT MORE WORK
!  TO GET A GEOMETRY POINT: 
!     MINMAX = 1 (-1) INDICATES THAT THE GEOMETRY POINT WAS FOUND BY
!     MAXIMIZING (-MINIMIZING) THE NEXT PIVOT.
!     IF FUNCTION COMPUTATION FAILS AT 'XGNEW', WE RETURN TO 'PTNEW' TO GET 
!     ANOTHER POINT: MINIMIZER IF 'XGNEW' WAS A MAXIMIZER (MINMAX=1)
!                    MAXIMIZER OF 'XGNEW' WAS A MINIMIZER (MINMAX=-1) 
!     MINMAX = 0 INDICATES THAT BOTH MINIMIZER AND MAXIMIZER WERE TRIED
!     SO THE SEARCH SHOULD BE REPEATED FOR A DIFFERENT POLYNOMIAL

      if ( minmax /= 0 ) then
        minmax = - minmax
        goto 10
      else if ( jpoly < n + 1 - neqcon ) then
!
!  IF WE STILL HAVE NOT TRIED ALL POLYNOMIALS  IN THE LINEAR BLOCK, 
!  PICK THE NEXT POLYNOMIAL, PUT IT IN APPROPRIATE PLACE  AND GO BACK TO 'PTNEW'
!
        jpoly = jpoly+1
        call swapnp( n, jpoly, nind+1, poly, lpoly )
        if ( iprint >= 2 ) write( iout, 8020 ) jpoly
        minmax = 0
        goto 10
      else
        impr = 0
        goto 25
      endif 
!
!  IF ALL POLYNOMIALS IN THIS BLOCK WERE TRIED, REDUCE RADIUS OF SEARCH
!  BY HALF AND REPEAT
!
!              DEL   = DEL/2.0D0
!              IF ( DEL > DELMIN ) THEN
!                JPOLY = NIND + 1
!                MINMAX = 0
!                GOTO 10
!              ELSE
!                GOTO 25
!              ENDIF
!            ENDIF
!          ELSE IF ( NIND < DD ) THEN
    else
      impr = 4
      goto 25
!
!  IF ALL POLYNOMIALS IN THIS BLOCK WERE TRIED, REDUCE RADIUS OF SEARCH
!  BY HALF AND REPEAT
!
!              DEL    = DEL/2.0D0
!              IF ( DEL > DELMIN ) THEN
!                JPOLY = NIND + 1
!                MINMAX = 0
!                GOTO 10
!              ELSE
!                GOTO 25
!              ENDIF

!            ENDIF
!          ELSE  
!            
!            CALL DCOPY( N, WRK(IXGNEW), 1, POINTS(NQ*N+1), 1 )
!            CALL GETDIS( N, NQ+1, NQ+1, POINTS, IN2SP(BASE), DIST, 
!     +                   WRK(ICURW), LENW )
!            DEL         = DIST( NQ+1 )/2.0D0
!            IF ( IPRINT >= 2 ) WRITE( IOUT, 8030 )
!            GOTO 10
    endif
  endif


!
!  DO  BOOKKEEPING RELATED WITH ADDING "XGNEW" TO "POINTS"
!
  call fmerit( m, values(nq+1), obfval(nq+1), conval(nq*m+1), &
               ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1),pp,method)
  call dcopy( n, wrk(ixgnew), 1, points(nq*n+1), 1 )
  nq          = nq + 1
  lpnts       = lpnts + n
  lvalue      = lvalue + 1
  lconvl      = lconvl + m
  call getdis( n  , nq  , nq, points, in2sp(base), dist, &
               wrk(icurw), lenw )
!
!  IF THE LEVEL OF EFFORT DOES NOT REQUIRE RECOMPUTING THE NEWTON
!  POLYNOMIALS FROM SCRATCH THEN  UPDATE THE NEWTON POLYNOMIALS TO
!  ACCOUNT FOR INCLUDING THE NEW POINT. 
!  THERE ARE TWO CASES: IXG=NIND+1 - THE NEW POINT IS ADDED TO 
!                                    INCOMPLETE INTERPOLATION SET
!                       IXG<NIND   - THE NEW POINT IS REPLACING 
!                                    INTERPOLATION POINT WITH INDEX IXG. 
!
  if ( ixg == nind+1 ) then
    impr        = 1
  else
    impr        = 2
  endif
  if  ( effort<=1 .or. &
      ( effort<=2 .and. ixg==nind+1)) then

    call ptrepl( wrk(ixgnew), ixg   , vnew, pntint, poly, nind, n, &    ! 
                 lpoly      , lptint)
    if ( iprint >= 2 ) write( iout, 8040 ) ixg, vnew
!
!  DO  BOOKKEEPING RELATED WITH ADDING "XGNEW" TO THE INTERP. SET
!
    fail = .false.
    if ( impr == 1 ) then
      pivval(ixg) = vnew
      nind        = nind+1
    else
      sp2in(in2sp(ixg)) = 0
      pivval(ixg)       = pivval(ixg)*vnew
    endif
    call dcopy( n, wrk(ixgnew), 1, pntint((ixg-1)*n+1), 1 )
    valint(ixg)     = values(nq)
    sp2in( nq ) = 1
    in2sp(ixg)  = nq
  endif 
endif

!
!  IF THE MODEL IS OLD OR WE COULD NOT FIND A GOOD "GEOMETRY" POINT THEN
!  WE RECOMPUTE THE INTERPOLATION SET (AND BASIS OF NEWTON POLYNOMIALS)
!  FROM SCRATCH. 
!
25 if ( fail ) then
  if ( impr == 0 ) goto 28
!
!  IF THE POINT FROM TRUST REGION MINIMIZATION HAS A GOOD VALUE, BUT
!  WAS NOT ADDED TO INTERPOLATION YET, THEN COPY IT AS SECOND POINT IN
!  INTERPOLATION SET (NOT THE FIRST, SINCE WE WANT TO INDICATE THAT
!  THE BASE HAD CHANGED) AND SET BASE=2 (-IMPR IS THE INDEX OF THE POINT) 
!
  if ( impr < 0 ) then
    if ( iprint >= 2 ) write( iout, 8050 )
    call dcopy(n, points((-impr-1)*n+1), 1, pntint(n+1), 1)
    base     = 2
    in2sp(2) = -impr
    call getdis( n, nq, 0, points, in2sp(base), dist, &
                 wrk(icurw), lenw )
  endif
!
!  SET APPROPRIATE VALUE FOR IMPR
!
  if (dist(in2sp(1))> layer*delta ) then
    if ( iprint >= 2 ) write( iout, 8060 ) dist(in2sp(1))
    impr = 3
  elseif ( impr /= 10 .and. impr /= 1 .and. impr /= 2 ) then
    impr = 4
  endif

!
!  SHIFT THE CURRENT BASE TO THE ORIGIN (BASE=1 MEANS IT IS AT THE
!                                        ORIGIN ALREADY)
!
28   if ( base /= 1 ) then
    call unshft(n, pntint((base-1)*n+1), x)
    do 30 j=1, nq
      call shift(n, pntint((base-1)*n+1), points((j-1)*n+1))
30     continue
  endif
  if ( iprint >= 2 ) write( iout, 8070 )
  call nbuild( poly  , points, values, pntint, valint, sp2in, &
               in2sp , nind  , n     , base  , dist  , delta, &
               pivthr, pivval, neqcon)
  if (( impr == 1 .or. impr == 2 ) .and. ( sp2in(nq) == 0 )) &
        impr = 0
endif 

return
1000 format( ' IMPMOD: *** ERROR: LWRK IS TOO SMALL!' / &
        '                    IT SHOULD BE AT LEAST ',i10,/ )
8000 format( ' IMPMOD: Compute a new point which improves geometry ',/)
8010 format( ' IMPMOD: Function computation failed at the new point',/)
8020 format( ' IMPMOD: Try to get a new point by  maximizing ',i4, &
                  '-th polynomial',/ )
8030 format( ' IMPMOD: Try to get a new point by reducing the radius ')
8040 format( ' IMPMOD: The new point is added as ', i4, &
                  '-th interpolation point, pivot=', d14.7,/)
8050 format( ' IMPMOD: The best point was not added to the ', &
         'interpolation yet',/ , 'It should be done here',/ )
8060 format( ' IMPMOD: The model is old, the base had moved=', d14.7 /)
8070 format( ' IMPMOD: Recompute the model ' )
end





!*******************************************************************************

!    NEXT SUBROUTINE

!*******************************************************************************


subroutine getdis( n, nq, idx, points, base, dist, wrk, lwrk )

!
!  ****************************************************************
!  THIS SUBROUTINE COMPUTES THE DISTANCE FROM THE BASE POINT TO
!  THE POINT WITH INDEX IDX. IF IDX IS 0 THEN DISTANCE TO ALL
!  POINTS IN 'POINTS' IS COMPUTED. THE DISTANCES ARE STORED
!  IN ARRAY 'DIST'
!
!  PARAMETERS
!  
!    N      (INPUT)  PROBLEM DIMENSION
!
!    NQ     (INPUT)  NUMBER OF POINTS IN 'POINTS'
!
!    IDX    (INPUT)  INDEX OF THE POINT FOR WHICH THE DISTANCE IS COMPUTED
!                    IF IDX=0 THEN DISTANCE IS COMPUTED FOR ALL POINTS
!    POINTS (INPUT)  ARRAY OF POINTS
!
!    BASE   (INPUT)  INDEX OF THE BASE POINT (FROM WHICH WE COMPUTE THE
!                    DISTANCES)
!    DIST   (OUTPUT) THE ARRAY OF DISTANCES
!  *****************************************************************
!

!
!  SUBROUTINE PARAMETERS 
!

integer          n, nq, base, lwrk, idx 
double precision points( lpnts+n ), dist( lvalue+1 ), wrk( lwrk )



!
!  COMMON VARIABLES
!

!
!  ARRAY LENGTHS
!
integer          lpoly, lpnts, lvalue, lptint, lvlint, lconvl
common /rpart/   lpoly, lpnts, lvalue, lptint, lvlint, lconvl
save /rpart/

!
!  PRINTOUT PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /



!
!  LOCAL VARIABLES
!
double precision zero
parameter      ( zero  = 0.0d0 )

integer          i , j , lb, li

double precision dnrmnf
external         dnrmnf

!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!       APPLICATION:       DNRMNF
!

!
!  CHECK SUFFICIENCY OF THE REAL WORKSPACE
!

if ( lwrk < n ) then
   if( iprint >= 0 ) write( iout, 1000 ) n
   stop
endif


!
!  CHECK IF THE IDX HAS A LEA GAL VALUE
!
if ( idx<0 .or. idx>nq ) then
   if( iprint >= 0 ) write( iout, 1010 ) idx
   stop 
endif

!
!  IF IDX=0 COMPUTE DISTANCE TO ALL POINTS
!
lb = ( base - 1 ) * n 

if ( idx == 0 ) then
  do 20 i = 1, nq
    if ( i == base ) then
      dist( i ) = zero
    else
      li = ( i - 1 ) * n 
      do 10 j = 1, n
        wrk( j ) = points( li + j ) - points( lb + j )
10       continue
!
!  COMPUTE THE NORM  BY CALLING DNRMNF (NOW L-INFINITY NORM)
!
      dist( i ) = dnrmnf( n, wrk )
    endif
20   continue

!
!  IF IDX>0 COMPUTE DISTANCE ONLY TO THAT POINT
!
else 
  if ( idx == base ) then
    dist( idx ) = zero
  else
    li = ( idx - 1 ) * n
      do 30 j = 1, n
        wrk( j ) = points( li + j ) - points( lb + j )
30       continue
!
!  COMPUTE THE NORM  BY CALLING DNRMNF (NOW L-INFINITY NORM)
!
    dist( idx ) = dnrmnf( n, wrk )
  endif
endif
return




1010 format( ' GETDIS: *** ERROR: IDX HAS ILLEGAL VALUE', i10, &
        '                    MUST BE A BUG!' ,/ )
1000 format( ' GETDIS: *** ERROR: LWRK IS TOO SMALL!' / &
        '                    IT SHOULD BE AT LEAST ',i10 )
end




!*************************************************************

!    NEXT SUBROUTINE

!*************************************************************



subroutine shift(n, x, y)

!
!  ***********************************************************
!  THIS SUBROUTINE SUBTRACTS VECTOR X FROM VECTOR Y: Y=Y-X
!
!  PARAMETERS
!     
!    X  (INPUT)        THE 'SHIFT' 
!
!    Y  (INPUT/OUTPUT) THE VECTOR WHICH IS 'SHIFTED'
!
!    N  (INPUT)        PROBLEM DIMENSION
!  **********************************************************
!

integer          n
double precision x(n), y(n)

!
!  LOCAL VARIABLES
! 
integer          i


do 10 i=1,n
  y(i)=y(i)-x(i)
10 continue
return
end




!**************************************************************

!    NEXT SUBROUTINE

!**************************************************************


subroutine unshft(n, x, y)

!
!  ***********************************************************
!  THIS SUBROUTINE ADDS   VECTOR X TO VECTOR Y: Y = Y + X
!
!  PARAMETERS
!     
!    X  (INPUT)        THE 'SHIFT' 
!
!    Y  (INPUT/OUTPUT) THE VECTOR WHICH IS 'SHIFTED'
!
!    N  (INPUT)        PROBLEM DIMENSION
!  **********************************************************
!

integer          n
double precision x(n), y(n)

!
!  LOCAL VARIABLES
! 
integer          i


do 10 i=1,n
  y(i)=y(i)+x(i)
10 continue
return
end




!**************************************************************

!    NEXT SUBROUTINE

!**************************************************************


double precision function dnrmnf( n, x )
!
!  *****************************************************
!  THIS FUNCTION COMPUTES THE INFINITY NORM OF VECTOR X:
!
!           DNRMNF = MAX_I {|X_I|}
!  
!  X (INPUT)  ARRAY OF LENGTH AT LEAST 'N' CONTAINING THE VECTOR
!
!  N (INPUT)  DIMENSION OF VECTOR IN X
!
!  *****************************************************

!
!  FUNCTION PARAMETERS
!
double precision x( n )
integer          n

!
!  LOCAL VARIABLES
!
double precision val
integer          j  
intrinsic        abs
!
!  FUNCTION CALLED
! 
!  BLAS:    ABS
!

val=0.0d0
do 10 j = 1, n
  if ( abs( x( j ) )>val) val=abs( x( j ) )
10 continue
dnrmnf=val
return
end




