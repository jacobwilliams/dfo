subroutine ptinit(n     , m     , x     , ldx   , nx    , np0   , &
                  nf    , points, obfval, conval, &
                  dist  , maxnf, delta , &
                  delmax, pivthr, lb    , ub    , a     , &
                  lda   , nclin , ncnln , scale , scal  , wrk   , &
                  lwrk  , iwrk  , liwrk , inform)


!
!  ******************************************************************
!  THIS SUBROUTINE BUILDS THE INITIAL INTERPOLATION SET.
!  GIVEN ONE STARTING POINT AND A TOTAL NUMBER OF POINT
!  IN INITIAL MODEL (NP0 - SPECIFIED BY THE USER) ADDITIONAL
!  POINTS ARE CONSTRUCTED IN SOME NEIGHBORHOOD OF THE INITIAL POINT.
!  AT LEAST  2 POINTS ARE NEEDED TO BUILD INITIAL MODEL.
!  AFTER DETERMINING THE POINTS AND THEIR FUNCTION VALUES WE
!  BUILD THE BASIS OF NEWTON FUNDAMENTAL POLYNOMIALS FOR THE
!  INITIAL INTERPOLATION SET
!  
!  PARAMETERS
!
!  N      (INPUT)  PROBLEM DIMENSION
!
!  NP0    (INPUT)  NUMBER OF POINTS REQUIRED FOR INITIAL MODEL BY THE
!                  USER
!         (OUTPUT) NUMBER OF POINTS FOR THE INITIAL MODEL SUPPLIED BY
!                  THE PROGRAM
!
!  X      (INPUT)  ARRAY OF LENGTH LDX*NX CONTAINING THE 'NX 'STARTING
!                  POINTS, PROVIDED BY THE USER.
!         (OUTPUT) CURRENT BEST POINT, POSSIBLY AFTER PROJECTION ONTO 
!                  FEASIBLE SET, IN FIRST N ENTREES
!
!  LDX    (INPUT)  LEADING DIMENSION OF ARRAY X
!
!  NX     (INPUT)  NUMBER OF INITIAL POINTS PROVIDED
!
!  DELTA  (INPUT)  TRUST REGION RADIUS
!
!  DELMAX (INPUT)  THE MAXIMUM TRUST REGION RADIUS ALLOWED
!
!  PIVTHR (INPUT)  PIVOT THRESHOLD VALUE
!
!  LB     (INPUT)  ARRAY OF LENGTH N+NCLIN+NCNLN OF LOWER BOUNDS 
!
!  UB     (INPUT)     ''       ''         ''        UPPER   ''
!
!  NCLIN  (INPUT)  NUMBER OF LINEAR ANALYTIC CONSTRAINTS
!
!  A      (INPUT)  (LDA X N) MATRIX OF LINEAR ANALYTIC CONSTRAINTS
!  
!  NCNLN  (INPUT)  NUMBER OF NONLINEAR ANALYTIC CONSTRAINTS
!
!  POINTS (OUTPUT) ARRAY (N*NP0) WITH THE POOL OF POTENTIAL INTERPOLATION 
!                  (SAMPLE) POINTS
!
!  OBJVAL (OUTPUT)  ARRAY (NP0) OF VALUES AT THE POINTS
!
!  CONVAL (OUTPUT)  ARRAY (M*NP0) OF VALUES OF CONSTRAINTS  AT THE POINTS
!
!  WRK             REAL SPACE WORKING ARRAY
!
!  IWRK            INTEGER SPACE WORKING ARRAY
!
!  INFORM (OUTPUT) INFORMATION ON EXIT
!              0    SUCCESSFUL MINIMIZATION
!              1    THE DERIVATIVES OF THE CONSTRAINT OR SOME PARAMETER
!                   SET BY THE USER IS INCORRECT
!              2    PROBLEM IS PROBABLY INFEASIBLE
!             -1    CANNOT COMPUTE FUNCTION VALUE AT ONE OF THE
!                   GENERATED POINTS OR ITS PROJECTION
!             -2    CANNOT COMPUTE FUNCTION VALUE AT ONE OF THE
!                   GENERATED POINTS SINCE RAN OUT OF FUNCTION EVALUATIONS
!  **********************************************************************
!


double precision x(nx*ldx)         , poly(lpoly), points(lpnts) , &
                 delta             , &
                 ub(n+nclin +ncnln), a(lda*n)   , scal(n)       , &
                 lb(n+nclin +ncnln), delmax        , &
                 pivthr            , wrk(lwrk)  , conval(lconvl), &
                 obfval(lvalue)    , dist(np0)  
integer          n     , np0   , nf   , nind, base, m, &
                 lda   , nclin , ncnln, lwrk      , iwrk, liwrk, &
                 inform, scale , nx   , ldx, maxnf       



!
!  COMMON VARIABLES:
!


!
!  MODEL PARAMETERS
!
include 'dfo_model_inc.inc'

!
!  PRINTOUT PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol 
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /
!
!  LENGTH OF ARRAYS
!

integer          lpoly, lpnts, lvalue, lptint, lvlint, lconvl
common /rpart/   lpoly, lpnts, lvalue, lptint, lvlint, lconvl
save /rpart/
!
!  INTERPOLATION CONTROL PARAMETERS
!      
integer          npmin, layer, effort
common / opti /  npmin, layer, effort
save / opti /

!
!  EXTERNAL ROUTINES
!

double precision ddot

external         ddot

!
!  LOCAL VARIABLES
!

logical          iferr, bdvltd, infeas

integer          i  , j    , np , inp, ix, iix  
double precision val, distb, del

double precision zero, half, two 

parameter       ( zero=0.0d0, half=0.5d0, two=2.0d0 )

!
!     SUBROUTINES AND FUNCTIONS CALLED:
!
!       APPLICATION:       MINTR , FUN   , GETDIS, SCL   , UNSCL,
!                          SHIFT , RANLUX, FUNCON, NBUILD
!                           
!                          
!       FORTRAN SUPPLIED:  MIN   , ABS
!       BLAS:              DCOPY , DDOT
!



!  **************************************************************
!  PROCESS THE FIRST POINT, PROVIDED BY THE USER, TO MAKE SURE IT
!  IS INCLUDED IN THE SAMPLE SET
!  **************************************************************


iferr  = .false.
bdvltd = .false. 
infeas = .false.

!
!  COPY THE POINT IN THE SET OF CURRENT SAMPLE POINTS
!

call dcopy( n, x, 1, points, 1 )
val = zero

!  ****************************************************
!  CHECK FEASIBILITY OF PROVIDED POINT
!  ****************************************************



!
!  CHECK IF THE INITIAL POINT IS FEASIBLE FOR SIMPLE BOUNDS
!  IF IT IS NOT, THEN PROJECT IT ON THE BOUNDS
!

do 10 i=1,n
  if (points(i).lt.lb(i)) then
    points(i) = lb(i)
    bdvltd    = .true.
  elseif (points(i).gt.ub(i)) then
    points(i) = ub(i)
    bdvltd    = .true.
  endif
10 continue


!
!  IF PROJECTED, PRINT WARNING MESSAGE 
!
if (bdvltd) then
  if ( iprint .ge. 0 ) write(iout, 1000)
  call dcopy( n, points, 1, x, 1 ) 
endif

!
!  CHECK FEASIBILITY OF FIRST POINT WRT LINEAR CONSTRAINTS
!

do 20 i=1, nclin
  val=ddot(n, a(i), lda, points, 1)
  if (val .gt. ub(n+i) .or. val .lt. lb(n+i)) infeas=.true.
20 continue


!
!  ******************************************************************
!  IF THERE ARE NONLINEAR CONSTRAINTS OR IF THE POINT VIOLATES LINEAR 
!  CONSTRAINTS, THEN  PROJECT IT ONTO FEASIBLE SET 
!  (IF IT IS ALREADY FEASIBLE, IS REMAINS THE SAME).
!  IF F_del  DENOTES INTERSECTION OF THE FEASIBLE REGION
!  AND TRUST REGION WITH CENTER AT X_1 AND RADIUS DEL, THEN WE SOLVE
!
!               2
!  MIN ||X-X_1||   S.T. { X IN F_del }
!
! 
!  **********************************************************
!


del=delta
if ( ncnln .gt. 0 .or. infeas) then
  do 40 i = 1, n
    gmod(i)   = -two*x( i )
    hmod(i,i) = two
    do 30 j = i+1, n
      hmod(i,j) = zero
      hmod(j,i) = zero
30    continue
40  continue
  
!
!  THIS MINIMIZATION PRODUCES PROJECTION OF THE  POINT ONTO
!  FEASIBLE REGION, INTERSECTED WITH TRUST REGION WITH RADIUS DEL
!


!
!  FIND THE PROJECTION
!
50  call mintr( n   , points, val   , del  , lb , ub  , &
             a   , lda   , nclin , ncnln, wrk, lwrk, &
             iwrk, liwrk , inform, 1) 



  if ( inform .eq. 1 ) then
    if ( iprint .gt. 0 ) write(iout,2000)
    return
  elseif (inform .eq. 2) then
!
!  IF FEASIBLE SOLUTION WAS NOT FOUND, TRY TO INCREASE 
!  TRUST REGION RADIUS 
!
    if (del .lt. delmax) then
      del=del*2
      goto 50
    else
!
!  IF THE TRUST REGION RADIUS IS AT IT'S MAXIMUM VALUE
!  THEN ASSUME THE PROBLEM IS INFEASIBLE AND QUIT
!
      if( iprint.gt.0 )   write(iout,2010)
      return
    endif
  endif
 
!
!  IF THE POINT HAS CHANGED, THEN THE STARTING POINT IS NOT
!  FEASIBLE AND WE FOUND ITS PROJECTION. PRINT A WARNING.
!
  val = zero
  do 60 i = 1, n 
    val = val + abs(points(i)-x(i))
60   continue
  if ( iprint .ge. 0 ) then
    if (val .gt. 100*n*mcheps) write(iout, 1010)
  endif
endif

if (val .gt. 100*n*mcheps .or. bdvltd) then

!
!  COMPUTE FUNCTION VALUE FOR THE FIRST POINT IF IT
!  HAD TO BE PROJECTED. IF SUCH COMPUTATION
!  FAILS, THEN QUIT THE PROGRAM
!
  nf = nf + 1
  if ( scale .ne. 0 ) call unscl( n, points, scal )
  call fun( n, m, points, obfval, conval, iferr )
  if ( scale .ne. 0 ) call scl( n, points, scal )
  if ( iferr ) then
    if( iprint.gt.0 ) write(iout,2020)
    inform = -1
    return
  endif
  if ( nf .ge. maxnf ) then
    inform = -2
    np0=1
    return
  endif
endif

!
!  INITIALIZE POINTERS TO CURRENT SAMPLE POINT (IN 'POINTS')
!
np  = 1
inp = n+1

!
!  CYCLE OVER ALL OTHER POINTS PROVIDED BY USER
!

do 120 ix=2, nx
!
!  POINTER TO CURRENT POINT PROVIDED BY USER (IN 'X')
!

  iix     = (ix-1)*ldx + 1

!
!  COPY THE POINT IN THE SET OF CURRENT SAMPLE POINTS
!

  call dcopy( n, x( iix ), 1, points( inp ), 1 )


!  ****************************************************
!  CHECK FEASIBILITY OF PROVIDED POINT
!  ****************************************************



!
!  CHECK IF THE INITIAL POINT IS FEASIBLE FOR SIMPLE BOUNDS
!  IF IT IS NOT, SKIP THIS POINT, AND MOVE TO THE NEXT ONE
!  PRINT A WARNING
!

  do 70 i=1,n
    if ( points(inp+i-1).lt.lb(i)-cnstol .or. &
         points(inp+i-1).gt.ub(i)+cnstol ) then 
      if ( iprint .ge. 0 ) write(iout, 1020) ix
      go to 120
    endif
70   continue


!
!  CHECK IF THE CURRENT POINT IS FEASIBLE WRT LINEAR CONSTRAINTS
!  IF IT IS NOT, SKIP THIS POINT, AND MOVE TO THE NEXT ONE
!  PRINT A WARNING
!

  do 80 i=1, nclin
    val=ddot(n, a(i), lda, points(inp), 1)
    if ( val .gt. ub(n+i)+cnstol .or. &
         val .lt. lb(n+i)-cnstol ) then
      if ( iprint .ge. 0 ) write(iout, 1030) ix
      go to 120
    endif
80  continue


!
!  CHECK IF THE POINT IS FEASIBLE WRT NONLINEAR CONSTRAINTS
!  BY PROJECTING IT ONTO THE FEASIBLE REGION (SEE HOW IT IS
!  DONE FOR THE FIRST POINT) AND CHECKING IF THE PROJECTION
!  IS DIFFERENT FROM THE POINT
!

  val = zero
  if ( ncnln .gt. 0 ) then
    do 100 i = 1, n
      gmod(i)   = -two*x( iix + i - 1 )
      hmod(i,i) = two
      do 90 j = i+1, n
        hmod(i,j) = zero
        hmod(j,i) = zero
90     continue
100   continue
  
!
!  THIS MINIMIZATION PRODUCES PROJECTION OF THE  POINT ONTO
!  FEASIBLE REGION, INTERSECTED WITH TRUST REGION WITH RADIUS DEL
!


!
!  FIND THE PROJECTION
!
    call mintr( n   , points(inp), val   , del  , lb , ub  , &
                a   , lda        , nclin , ncnln, wrk, lwrk, &
                iwrk, liwrk      , inform, 1) 



    if ( inform .eq. 1 ) then
      if ( iprint .gt. 0 ) write(iout,2000)
      return
    elseif (inform .eq. 2) then
!
!  IF NO FEASIBLE SOLUTION WAS FOUND THEN SKIP THIS POINTS AND
!  MOVE TO THE NEXT ONE
!
      if ( iprint .gt. 0 ) write( iout, 1040 ) ix
      go to 120 
    endif
 
!
!  IF THE POINT HAS CHANGED, THEN THE STARTING POINT IS NOT
!  FEASIBLE, WE SKIP IT AND MOVE TO THE NEXT POINT
!
    val = zero
    do 110 i = 1, n 
      val = val + abs(points(inp+i-1)-x(iix+i-1))
110     continue
    if (val .gt. 100*n*mcheps) then
      if ( iprint .gt. 0 ) write(iout, 1040) ix
      go to 120
    endif
  endif
!
!  IF THE POINT PASSED ALL FEASIBILITY TESTS, THEN ACCEPT IT AS
!  A SAMPLE POINT AND RECORD ITS FUNCTION VALUES AND CONSTRAINT VALUES
!
  np           = np  + 1
  inp          = inp + n
  if ( np .ne. ix ) then
    obfval( np ) = obfval( ix )
    call dcopy( m, conval((ix-1)*m+1), 1, conval((np-1)*m+1), 1 )
  endif

120 continue  

np0 = np

!  --------------------------------------------------------------
!
!  IF THERE IS ONLY ONE POINT IN THE SAMPLE SET, TRY TO FIND ANOTHER
!
!  --------------------------------------------------------------


!
!  GENERATE A POINT RANDOMLY WITHIN DELTA DISTANCE FROM X
!  AND SATISFYING THE SIMPLE BOUNDS
!  WRITE  THE POINT IN ARRAY 'POINTS'
!

if ( np0 .eq. 1 ) then
  if ( iprint .ge. 2 ) write( iout, 8000 ) 
  call dcopy( n, points, 1, x, 1 )
!
!  GENERATE N RANDOM NUMBERS BETWEEN 0 AND 1 AND RECORD THEM TO
!  POINTS STARTING FROM N+1-ST ENTREE
!       

  call ranlux(points(n+1), n)
  do 130 j = 1, n
    distb=zero
    distb=min(delta, ub(j)-x(j))
    if (distb.gt.mcheps) then
      points( n + j ) = distb * points(n+j) + x( j )
    endif 
    if (distb.le.mcheps) then
      distb=min(delta, x(j)-lb(j))
      points( n + j ) = -distb * points(n+j) + x( j ) 
    endif
130  continue

!
!  CHECK FEASIBILITY OF AUXILIARY POINT WRT LINEAR CONSTRAINTS
!
  infeas = .false.
  do 140 i=1, nclin
    val=ddot(n, a(i), lda, points(n+1), 1)
    if (val .gt. ub(n+i) .or. val .lt. lb(n+i)) infeas = .true.
140   continue
  val = zero
!  -------------------------------------------------------
!  FIND THE SECOND POINT FOR INTERPOLATION
!  -------------------------------------------------------

  if ( ncnln .gt. 0 .or. infeas) then          
    do 160 i = 1, n
      gmod(i)   = -two*points(n+i)
      hmod(i,i) =  two
      do 150 j = i+1, n
        hmod(i,j) = zero
        hmod(j,i) = zero
150       continue
160     continue
  
  
    call mintr ( n   , x    , val   , delta, lb , ub  , &
                 a   , lda  , nclin , ncnln, wrk, lwrk, &
                 iwrk, liwrk, inform, 1) 

    if (inform .eq. 1) then
      if( iprint.gt.0 )   write(iout,2000)
      return
    elseif (inform .eq. 2) then
      if( iprint.gt.0 )   write(iout,2010)
      return
    endif


!  -------------------------------------------------------
!  IF FIRST AND SECOND POINTS COINCIDE, FIND A DIFFERENT SECOND POINT
!  -------------------------------------------------------

    val = zero
    do 170 i = 1, n 
      val = val + abs(x(i)-points(i))
170     continue

    if (val .lt. 10*n*pivthr) then
      do 190 i = 1, n
        gmod(i)   =  two*points(n+i)
        hmod(i,i) = -two
        do 180 j = i+1, n
          hmod(i,j) = zero
          hmod(j,i) = zero
180         continue
190       continue


      call mintr ( n   , x    , val   , del  , lb , ub  , &
                   a   , lda  , nclin , ncnln, wrk, lwrk, &
                   iwrk, liwrk, inform, 1 ) 

      if ( inform .eq. 1 ) then
        if( iprint.gt.3 )   write(iout,2000)
        return
      elseif ( inform .eq. 2 ) then
        if( iprint.gt.3 )   write(iout,2010)
        return
      endif 
    endif
    call dcopy ( n, x, 1, points(n+1), 1 )
  endif

!  --------------------------------------------------------------
!
!  INTERPOLATION POINTS ARE COMPUTED
!
!  --------------------------------------------------------------




!
!  COMPUTE FUNCTION VALUE FOR THE  AUXILIARY  POINT. 
!  IF THERE ARE NONLINEAR CONSTRAINS, AND  IF FUNCTION EVALUATION 
!  FAILS FOR THE AUXILIARY  POINT, WE QUIT.
!  IF THERE ARE NO NONLINEAR CONSTRAINTS, THEN IF FUNCTION EVALUATION
!  FAILS FOR AUXILIARY POINT, WE SUBSTITUTE IT BY ITS CONVEX
!  COMBINATION WITH THE FIRST POINT:
!
!     X_k <-- 1/2(X_k+X_1)
!


200   if ( scale .ne. 0 ) call unscl( n, points( n+1 ), scal )
  call fun( n, m, points( n+1 ),  obfval( 2 ), conval( m+1 ), &
            iferr )
  if ( scale .ne. 0 ) call scl( n, points( n+1 ), scal )
  nf = nf + 1
  if ( nf .ge. maxnf ) then
    inform = -2
    np0=2
    return
  endif 
  if ( iferr ) then
    do 210 j = 1, n
      points(n+j) = half*(points(j) + points(n+j))
210     continue
!
!  IF THE SECOND POINT GETS TOO CLOSE TO FIRST POINT, QUIT
!
    call getdis( n, 2, 2, points, 1, dist, wrk, lwrk ) 
    if ( dist(2) .lt. pivthr ) then  
      inform = -1
      return 
    endif
!
!  IF THE SECOND POINTS BECOMES NON-FEASIBLE, QUIT
!
    if ( ncnln.gt. 0 ) then 
!            CALL FUNCON(1  , NCNLN , N  , NCNLN , IWRK,
!     *                  POINTS(N+1), WRK, WRK(NCNLN+1), 1  )
       call easycon(n, points(n+1), ncnln, wrk)
      do 220 j=1, ncnln 
        if ( wrk(j) .lt. lb(n+nclin+j) - cnstol .or. &
             wrk(j) .gt. ub(n+nclin+j) + cnstol     ) then
          inform = -1
          return 
        endif
220       continue
    endif  
    goto 200
  endif
!
!  SET THE NUMBER OF SAMPLE POINTS TO EQUAL 2
!
  np0 = 2
endif



!
!  CHECK IF ANY OF THE POINTS ARE TOO FAR, AND WOULD NOT BE INCLUDED
!  IN INTERPOLATION SET
!
bdvltd = .false.
do 250 i = 1, np0
  if ( dist(i) .gt. layer*delta ) bdvltd = .true.
250 continue
if ( bdvltd ) then
  if ( iprint .ge. 0 ) write(iout,1050)
endif

inform = 0
return 

1000 format(' DFO: WARNING: THE FIRST INITIAL POINT IS OUT OF BOUNDS',/ &
       ' DFO:',  10x, 'IT WILL BE PROJECTED ON THESE BOUNDS',/)
1010 format(' DFO: WARNING: THE FIRST INITIAL POINT IS NOT FEASIBLE', / &
       ' DFO:', 10x, 'IT WILL BE PROJECTED ON THE FEASIBLE SET',/)
1020 format(' DFO: WARNING: THE ', i4,'-TH INITIAL POINT IS OUT OF &
         BOUNDS',/ 'DFO:',    10x, 'IT WILL BE IGNORED',/)
1030 format(' DFO: WARNING: THE ', i4,'-TH INITIAL POINT DOES &
         NOT SATISFY','/DFO:',  10x, 'LINEAR CONSTRAINTS, &
         IT WILL BE IGNORED',/)
1040 format(' DFO: WARNING: THE ', i4,'-TH INITIAL POINT IS NOT &
         FEASIBLE',/ 'DFO:',    10x, 'IT WILL BE IGNORED',/)
1050 format(' DFO: WARNING: AT LEAST ONE  INITIAL POINT IS TOO FAR &
       FROM',/ 'DFO: THE BASE, IT WILL NOT BE IN THE INTERPOLATION &
       SET ',/ 'DFO: TO INCLUDE IT, INCREASE PARAMETER DELTA &
       OR LAYER',/)      
2000 format('DFO: SOME PARAMETER IN THE PROBLEM FORMULATION ',/ &
       'DFO: HAS ILLEGAL VALUE OR DERIVATIVES OF THE CONSTRAINTS'/ &
       'DFO: ARE WRONG. THE PROGRAM WILL STOP',/)
2010 format( 'DFO: FEASIBLE SET SEEMS TO BE EMPTY ',/ &
        'DFO:  THE PROGRAM WILL STOP',/)
2020 format( 'DFO: FUNCTION VALUE WAS NOT FOUND FOR INITIAL POINT',/ &
        'DFO: OR ITS PROJECTION. THE PROGRAM WILL STOP',/)
2030 format( 'DFO: FUNCTION VALUE WAS NOT FOUND FOR AN AUXILIARY &
        POINT.',/ 'DFO:  THE PROGRAM WILL STOP',/)
8000 format( 'DFO: PTINIT: Getting an auxiliary point '/ )
end





! *****************************************
!
!  NEW SUBROUTINE
!
! *****************************************




subroutine vlinit(n , m , nq, x, base, points, values, obfval, &
                  conval, cu, cl, dist, penpar, method, wrk, lwrk)





double precision points(lpnts) , x(n), values(lvalue), &
                 conval(lconvl), obfval(lvalue), dist(nq), &
                 penpar,         wrk(lwrk), cu(*), cl(*)
integer          n     , nq    , m,  base, lwrk, method

!
!  PRINTOUT PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol 
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /

!
!  LENGTH OF ARRAYS
!

integer          lpoly, lpnts, lvalue, lptint, lvlint, lconvl
common /rpart/   lpoly, lpnts, lvalue, lptint, lvlint, lconvl
save /rpart/


integer          i

!
!  COMPUTE THE VALUE OF THE MERIT FUNCTION AND CHOSE THE BASE
!

base = 1
do 10 i=1, nq
  call fmerit(m, values(i), obfval(i), conval((i-1)*m+1), cu, cl, &
              penpar, method ) 
  if ( values(i) .lt. values(base) + 1.0d2*mcheps ) base=i
10 continue


!  
!  THE BASE IS CHOSEN, SHIFT THE POINTS SO THAT THE BASE IS AT THE ORIGIN
! 

call dcopy(n, points((base-1)*n+1), 1, x, 1)
do 20 i=1,nq
  call shift(n, x, points((i-1)*n+1))
20 continue



!
!  GET THE DISTANCES OF ALL POINTS IN 'POINTS' TO THE BASE
!

call getdis( n, nq, 0, points, base, dist, wrk, lwrk)


return
end



