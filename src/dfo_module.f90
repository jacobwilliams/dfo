module dfo_module

  implicit none

contains

subroutine dfo( n     , nx   , x     , ldx, fx, conx, &
                ifiniv, m    , c     , &
                nclin , ncnln, lb    , ub , a   , lda , &
                xnames, pname, cnames, it , nf  , info, &
                maxit,  maxnf, stpcrtr, delmin, stpthr, &
                cnstolp, delta, pp, scale, ioutp , iprintp)

!
!  ****************************************************************
!  THIS SUBROUTINE MINIMIZES A NONLINEAR OBJECTIVE FUNCTION
!  SUBJECT TO LINEAR AND (POSSIBLY) NONLINEAR CONSTRAINTS
!  AND SIMPLE BOUNDS, WITHOUT USING DERIVATIVES OF 
!  THE OBJECTIVE FUNCTION AND NONLINEAR CONSTRAINTS.
!
!  THE PROBLEM IS CONSIDERED TO BE OF THE FORM
!                        
!                           MIN F( X )
!         
!        S.T.
!                               / X  \
!                        LB <= ( AX   ) <= UB
!                               \C(X)/
!
!  THIS PROGRAM IS BASED ON USING QUADRATIC INTERPOLATION OF THE
!  OBJECTIVE FUNCTION AND CONSTRAINTS IN COMBINATION WITH TRUST 
!  REGION FRAMEWORK
!
!  PARAMETERS
!  ==========
!
!  (INPUT) 
!
!    X         : ARRAY OF NX  INITIAL POINTS, GIVEN BY THE USER
!                'X' HAS TO CONTAIN AT LEAST ONE STARTING POINT
!                IF NO FUNCTION VALUES ARE PROVIDED FOR THE POINTS
!                IN 'X', THEN ALL BUT THE FIRST POINT ARE IGNORED.
!
!    LDX       : LEADING DIMENSION OF ARRAY 'X' (LDX>= N)
!
!    FX        : ARRAY OF FUNCTION VALUES AT THE INITIAL POINTS
!                OF LENGTH AT LEAST 1 ( MAY BE DUMMY IF NO VALUES  
!                ARE PROVIDED )
!    CONX      : ARRAY OF INITIAL VALUES OF THE DERIVTIVE FREE CONSTRAINTS
! 
!
!    NX        : NUMBER OF INITIAL POINTS FOR WHICH THE VALUE IS
!                PROVIDED. IF NO SUCH POINTS PROVIDED, THEN NX=1
!
!    IFINIV    : LOGICAL VARIABLE, .TRUE. IS INITIAL VALUES FOR
!                NX INITIAL POINTS ARE PROVIDED, .FALSE. OTHERWISE
!
!    N         : PROBLEM DIMENSION
!
!    NCLIN     : NUMBER F LINEAR CONSTRAINTS
!
!    NCNLN     : NUMBER OF NONLINEAR CONSTRAINTS
!
!    M         : NUMBER OF DERIVATIVE FREE CONSTRAINTS
!
!    LB        : ARRAY OF LOWER BOUNDS OF LENGTH >= (N+NCLIN+NCNLN)
!
!    UB        : ARRAY OF UPPER BOUNDS OF LENGTH >= (N+NCLIN+NCNLN)
!
!    A         : MATRIX OF LINEAR CONSTRAINTS,( DIMENSIONS: LDA X N)
!
!    LDA       : LEADING DIMENSION OF MATRIX A, HAS TO BE >= MAX(1, NCLIN)
!
!
!    XNAMES    : ARRAY OF STRINGS OF CHARACTERS, CONTAINING
!                NAMES OF THE VARIABLES
!
!    PNAME     : STRING OF CHARACTERS CONTAINING THE NAME OF THE PROBLEMS
!
!    CNAMES    : ARRAY OF STRINGS OF CHARACTERS, CONTAINING NAMES
!                OF CONSTRAINTS
!    MAXIT     : MAXIMUM NUMBER OF ITERATIONS ALLOWED
!   
!    MAXNF     : MAXIMUM NUMBER OF FUNCTION EVALUATIONS
!
!    STPCRTR   : STOPPING CRITERIA (1-BY MIN DELTA, 2-BY SLOW PROGRESS)
!    
!    DELMIN    : MINIMUM TRUST REGION RADIUS
!
!    STPTHR    : THRESHOLD FOR SLOW PROGRESS
!
!    CNSTOLP   : FEASIBILITY TOLERANCE FOR "EASY" CONSTRAINTS
!
!    DELTA     : INITIAL TRUST REGION RADIUS
!
!    PP        : THE PENALTY PARAMETER
!    
!    SCALE     : SCALING, 0-NO SCALING, 1-SCALE INITIAL POINT TO UNIT VECTOR
!                2- SCALE THE BOX TO UNIT BOX
!    IOUTP     : OUTPUT DEVICE, 6=SCREEN
! 
!    IPRINTP   : LEVEL OF PRINT <0 - SILENT MODE, 0 - NO OUTPUT, 
!                BESIDES FINAL OUTPUT AND SOME WARNINGS
!                AND ERROR MESSAGES, 1- SUMMARY OF EACH ITERATION
!                AND SOME WARNINGS, 2,3 - PROGRESSIVELY MORE OUTPUT.
!
!    (OUTPUT)
!
!     X        : THE OPTIMAL (OR BEST FOUND)  POINT 
!
!     FX       : THE VALUE OF OBJECTIVE FUNCTION AT X
!  
!     C        : THE ARRAY OF VALUES OF THE DERIVATIVE FREE CONSTRAITNS
!                AT X
!
!     INFO     : THE EXIT INFORMATION:
!                0     SUCCESSFUL MINIMIZATION
!                1     TOO MANY FUNCTION EVALUATIONS
!                2     TOO MANY ITERATIONS
!               -1     ILLEGAL VALUE OF AN INPUT PARAMETER
!               -2     ILLEGAL VALUE OF LDA OR LDX
!               -3     REAL WORKSPACE IS TOO SMALL
!               -4     INTEGER WORKSPACE IS TOO SMALL
!               -5     INCONSISTENT BOUNDS
!               -6     FUNCTION VALUE CANNOT BE COMPUTED
!                      AT INITIAL POINT
!               -7     MINIMIZATION OVER TRUST REGION FAILED
!               -8     MAXNF IS TOO SMALL OR T0O LARGE
!               -9     CANNOT BUILD A 2 POINT MODEL
!     IT       : THE NUMBER OF ITERATIONS
!
!     NF       : THE NUMBER OF FUNCTION EVALUATIONS
!
!  ********************************************************************
!
!     Programming:  K. Scheinberg, 1998-2003.
!
!  ********************************************************************
!

!
!  SUBROUTINE PARAMETERS
!
integer          n , m , nx  , ldx  , nclin, ncnln, lda, &
                 it, nf, info

double precision x(ldx*nx), fx(nx) , lb(*), &
                 ub(*) , a(lda*n), &
                 c(*), conx(*)        

logical          ifiniv

character*256     xnames(n), pname, cnames(*)
!
!  COMMON VARIABLES
!
!
include 'dfo_model_inc.inc'
!
!  PROCESS CONTROL PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /

!
!  INTERPOLATION CONTROL PARAMETERS
!      
integer          npmin, layer, effort
common / opti /  npmin, layer, effort
save / opti /

!
!  LENGTH OF ARRAYS
!

integer          lpoly, lpnts, lvalue, lptint, lvlint, lconvl
common /rpart/   lpoly, lpnts, lvalue, lptint, lvlint, lconvl
save /rpart/

!
!  LOCAL VARIABLES AND SOME USER DEFINED PARAMETERS
!
integer          lwrk , liwrk, ddmax , lmax  , rsmax0, rsmax1, &
                 rsmax2, rsmax3, &
                 rsmax4 , ismax0, ismax1, ismax2, nfpmax

parameter      ( lmax   = nvarmx+nlinmx+nnlnmx+nconmx )

parameter      ( ddmax  = (nvarmx+1)*(nvarmx+2)/2)

parameter      ( nfpmax = ddmax*(nvarmx+1)*nvarmx/2+(nvarmx+1)* &
                          nvarmx+1 )

parameter      ( rsmax0 = nfunmx*(nvarmx+3+nconmx) + &
                          ddmax*(nvarmx+3) + 2*nvarmx + nvarmx**2 &
                        + nvarmx**2*nconmx+nconmx+nfpmax)
parameter      ( rsmax1 = 3*(lmax) + (nvarmx+1)*(nnlnmx+nconmx) &
                          + nvarmx*nvarmx + 2*nvarmx + 1 )
parameter      ( rsmax2 = 2*nvarmx*nvarmx + nvarmx*nlinmx + &
                        2*nvarmx*(nnlnmx+nconmx) &
                      + 20*nvarmx + 11*nlinmx+21*(nnlnmx+nconmx)) 
parameter      ( rsmax3 = ddmax*ddmax + (ddmax-1)*nvarmx + &
                          4*ddmax-1 )

parameter      ( rsmax4 = 3*nvarmx+nvarmx**2 )

parameter      ( ismax0 = nfunmx + ddmax )
parameter      ( ismax1 = 3*nvarmx + nlinmx+2*nnlnmx )
parameter      ( ismax2 = lmax )
parameter      ( lwrk   = rsmax0 + rsmax1 +  rsmax2 +  rsmax4+ &
                          rsmax3  )
parameter      ( liwrk  = ismax0 + ismax1 +  ismax2 + ddmax )
!
!  ALLOCATE THE WORKING SPACE
! 
!    WRK       : REAL WORKSPACE ARRAY OF SIZE LWRK
!                THE REQUIRED LENGTH OF LWRK IS ROUGHLY APPROXIMATED 
!                FROM ABOVE  BY:
!                (2N+26)(N+NCLIN+NCNLN+M)+(2NCNLN+2M+6)*N+(MAXNF+NX)(N+2)+
!                [(N+1)**2 (N+2)(N+14)]   
! 
! 
!    IWRK      : INTEGER WORKSPACE ARRAY OF SIZE LIWRK
!                THE REQUIRED LENGTH OF INTEGER SPACE IS
!                4(N+NCLIN+NCNLN) + (MAXNF+NX) + (N+1)(N+2)
integer          iwrk(liwrk)
double precision wrk(lwrk)


integer          maxnf  , maxit , varnt , iospec, iscal, &
                 ipar(8), badbnd, scale , i     , j    , &
                 icurw  , lenw  , lscal , mjit  , nmjit, &
                 ip     , iv    , iob   , ic    , id   , &
                 lp     , lv    , lob   , lc    , ld   , &
                 inform , base  , nq    , method, stpcrtr, &
                 iprintp, ioutp 

double precision delmin  , lowbnd, anoise , rnoise, pivthr, &
                 addthr  , xchthr, rpar(8), delta , delmax, &
                 noise(2), avgsc , val, pp, stpthr, cnstolp

logical          iferr 
!
!  PARAMETERS
!   
double precision one, zero
parameter      ( one = 1.0d0, zero = 0.0d0 )
!
!  SUBROUTINES CALLED: 
!
!  APPLICATIONS: DFOSLV, PTINIT, FUN, SCL, UNSCL
!  BLAS        : DDOT
!
!
!  SET PARAMETERS FOR /MERIT/ COMMON BLOCK
! 
ncon=m
nlin=nclin
nnln=ncnln
!
!  DEFINE THE MINIMIZATION PARAMETERS THAT ARE NOT DEFINED BY THE USER
!
!   MAXMIMUM TR RADIUS
delmax=10*delta
!   LOWER BOUND ON FUNCTION VALUE
lowbnd=-1.0d0
!   ABSOLUTE NOISE
anoise=0.0d0
!   RELATIVE NOISE
rnoise=0.0d0
!   THRESHHOLD FOR POISEDNESS
pivthr=0.001d0
!   THRESHOLD FOR IMPROVING POISEDNESS BY ADDING POINT
addthr=100d0
!   THRESHOLD FOR IMPROVING POISEDNESS BY REPLACING POINT
xchthr=1000d0
!   MINIMUM NUMBER OF POINTS IN A MODEL 
npmin = 2
!   LAYER*DELTA IS THE MAXIMUM DISTANCE TO THE BASE ALLOWED
!   FOR THE INTERPOLATION POINTS
layer = 6
!   EFFORT LEVEL. 1,2 AND 3 INCREASING LEVELS OF LIN. ALGEBRA
effort = 1
!   METHOD FOR BULSING AN INCOMPLETE INTERP. MODEL. 1- MIN. FROBENIUS
!                                    2- INCOMPLETE MONOMIALS
varnt = 1
!   NUMBER OF MAJOR ITERATIONS
nmjit = 1
!   METHOD FOR HANDLING CONSTRAINTS. 1 - MODEL CONSTR, 2-SAME WITH SAFEGUARD,
!          3, 4 PENALTY FUNCTIONS
!
!    IF USING IPOPT DO NOT USE OPTIONS 2 OR 3 FOR NOW!!!
!
method=1
!   RECORD SOME OTHER PARAMETERES INTO COMMON BLOCK VARIABLES
iout=ioutp
iprint=iprintp
cnstol=cnstolp
!
!  STORE THE MINIMIZATION PARAMETERS, OBTAINED FROM SPECIFICATION FILE
!

mcheps     = 1.5d-16
ipar( 1 )    = maxit
ipar( 2 )    = maxnf
ipar( 3 )    = npmin
ipar( 4 )    = layer
ipar( 5 )    = effort
ipar( 6 )    = varnt
ipar( 7 )    = method
ipar( 8 )    = stpcrtr
rpar( 1 )    = delta
rpar( 2 )    = delmax
rpar( 3 )    = delmin
rpar( 4 )    = pivthr
rpar( 5 )    = addthr
rpar( 6 )    = xchthr
rpar( 7 )    = lowbnd
rpar( 8 )    = stpthr
noise( 1 )   = anoise
noise( 2 )   = rnoise

!  ----------------------------------
!  CHECK IF THE INITIAL DATA IS VALID
!  ----------------------------------


!
!  CHECK IF THE DIMENSIONS ARE POSITIVE
!
if ( n <= 0 .or. n > nvarmx ) then
   if (iprint >= 0 ) write( iout, 1000 ) n
   info = -1
   return
endif

if ( m < 0 .or. m > nconmx ) then
   if (iprint >= 0 ) write( iout, 1005 ) m
   info = -1
   return
endif

if ( nclin < 0  .or. nclin > nlinmx) then
   if (iprint >= 0 ) write( iout, 1010 ) nclin
   info = -1
   return
endif

if (  ncnln < 0 .or. ncnln > nnlnmx) then
   if (iprint >= 0 ) write( iout, 1020 ) ncnln
   info = -1
   return
endif
if ( maxnf > nfunmx-nx ) then
   if (iprint >= 0 ) write( iout, 1110 ) nfunmx
   info = -8
   return
endif
if ( maxnf <= 0 ) then 
   if (iprint >= 0 ) write( iout, 1050 ) 2
   info=-8
   return
endif
!
!  CHECK IF THE  BOUNDS ARE CONSISTENT
!
badbnd=0
do 10 i=1,n+nclin+ncnln+m
  if (lb( i ) > ub( i )+mcheps) badbnd=1
10 continue
if (badbnd==1) then
  if (iprint >= 0 ) write(iout, 1030)
  info=-5
  return
endif



!
!  CHECK IF DIMENSION OF LINEAR CONSTRAINT MATRIX IS SATISFACTORY
!

if ( lda < max(nclin,1) ) then 
   if (iprint >= 0 ) write(iout, 1060) max(nclin,1)
   info = -2              
   return
endif

!
!  CHECK IF LEADING DIMENSION OF INITIAL SET 'X' IS SATISFACTORY
!

if ( ldx < n ) then 
   if (iprint >= 0 ) write(iout, 1080) n
   info = -2              
   return
endif

!
!  CHECK IF OTHER DIMENSION OF INITIAL SET 'X' IS POSITIVE
!

if ( nx <= 0 ) then 
   if (iprint >= 0 ) write(iout, 1090) 
   info = -2              
   return
endif
!
!  CHECK IF OTHER INPUT PARAMETERS ARE VALID
!
if ( stpcrtr > 2 .or. stpcrtr < 1 .or. delmin< 0 .or. &
     stpthr < 0 .or. cnstolp < 0 .or. delta < 0 ) then
   if (iprint >= 0 ) write(iout, 1100) 
   info = -1           
   return
endif

!
!  PRINT OUT SOME SPECIFICATIONS INFO
!


if ( iprint >= 2 ) then
  if ( varnt ==1 ) then
    write ( iout, 4000 )
  else
    write ( iout, 4020 )
  endif
endif
!
!  REAL SPACE ALLOCATION
!

!
!  SET THE POINTERS TO SPACE FOR SCALING COEFFICIENTS
!
iscal = 1
!
!  POINTER TO  ARRAY OF ALL SAMPLE POINTS
!
ip  = iscal + n
!
!  POINTER TO MERIT FUNCTION VALUES AT SAMPLE POINTS
!
iv  = ip + ( maxnf + nx ) * n
!
!  POINTER TO OBJECTIVE FUNCTION VALUES AT SAMPLE POINTS
!
iob = iv + ( maxnf + nx ) 
!
!  POINTER TO CONSTRAINT  FUNCTION VALUES AT SAMPLE POINTS
!
ic  = iob + ( maxnf + nx ) 
!
!  POINTER TO DISTANCES BETWEEN CURRENT BASE AND  SAMPLE POINTS
!
id  = ic + max(1,( maxnf + nx )*m)
!
!  CHECK IF REAL SPACE IS SUFFICIENT
!
icurw = id + ( maxnf + nx )
lenw  = lwrk - icurw + 1
if ( lenw < 0 ) then
  write( iout, 1070 ) lwrk
  info = -3
  return
endif
lscal = iscal - 1
!
! INITIALIZE NF TO ZERO.
! IF NO INITIAL POINT WITH FUNCTION VALUE IS PROVIDED, COMPUTE
! THE FUNCTION VALUE AT THE FIRST POINT AND SET NF = 1
!
nf = 0
if ( .not. ifiniv ) then
  call fun( n, m, x, wrk(iob), wrk(ic), iferr )
  if ( iferr ) then
    if (iprint >= 0 ) write(iout, 1040)
    info = -6
    return
  endif
  fx(1)=wrk(iob)
  call dcopy(m, wrk(ic), 1,c, 1)
  nf = 1
  nx = 1
endif
!
!  CHECK IF THE INITIAL NUMBER OF INTERPOLATION POINTS
!  IS NOT LARGER THAN MAXIMUM NUMBER OF FUNCTION EVALUATIONS
!  
if ( nf >= maxnf ) then 
   if (iprint >= 0 ) write( iout, 1120 )
   info=1
   it=1
   scale =0
   goto 48
endif

!
!  IF THERE IS SOME INTIAL INFORMATION, THEN WE RECORD IT
!  IN ARRAYS WRK(IOB) AND WRK(IC) FOR FUTURE CONVINIENT HANDLING.
!

if ( ifiniv ) then
  call dcopy( nx, fx, 1, wrk(iob), 1 )
  call dcopy( nx*m, conx, 1, wrk(ic), 1 )
endif

!       
!  CHECK IF SCALING AND 'EASY' NONLINEAR CONSTRAINTS ARE NOT USED AT
!  THE SAME TIME (NOT SUPPORTED BY CURRENT VERSION)
!
if ( scale /= 0 .and. ncnln > 0 ) then
   if (iprint >= 0 ) write(iout, 5000) 
   scale = 0
endif

!       
!  SCALE THE PROBLEM IF REQUIRED
!

if ( scale /= 0 ) then
  avgsc=one
  do 20 i=1, n
    if ( scale == 1 ) then
      wrk( lscal + i ) = one + abs(x(i))
      avgsc   = avgsc * wrk( lscal + i )
    else 
      wrk( lscal + i ) = ub(i) - lb(i)
      if ( wrk( lscal + i) < cnstol .or. ub(i) >= 1.0d+20 &
                           .or. lb(i)  <= -1.0d+20 ) then
        if (iprint >= 0 ) write( iout, 5200 )  i
        wrk( lscal + i ) = one
        avgsc   = avgsc * wrk( lscal + i )
      endif
    endif
20   continue
endif

!
!  SCALE THE POINT AND THE BOUNDS
!
if ( scale /= 0 ) then
  do 25 i=1,nx
    call scl( n, x((i-1)*n+1), wrk( iscal ) )
25   continue  
  call scl( n, lb, wrk( iscal ) )
  call scl( n, ub, wrk( iscal ) )
!
!  SHOULD SCALE THE TRUST REGION RADIUS AS WELL. 
!

  rpar( 1 ) = rpar( 1 )/(exp(log(avgsc)/n))
  rpar( 2 ) = rpar( 2 )/(exp(log(avgsc)/n))
  rpar( 3 ) = rpar( 3 )/(exp(log(avgsc)/n))
  do 40 i=1,n
    do 30 j=1,nclin
      a(lda*(i-1)+j)=a(lda*(i-1)+j)*wrk( lscal + i )
30     continue  
40   continue  
  if ( iprint >= 2 ) then
    write( iout, 8000 )
    do 45 i = 1, n
      write( iout, 8010 ) wrk( lscal + i )
45     continue
  endif  
endif
!
!  INITIALIZE: NUMBER OF SAMPLE POINTS
!              LENGTH OF ARRAY OF POINTS AND VALUES
!             

nq     = max( nx, 2 )
lvalue = nq
lpnts  = nq*n
lconvl = max(1,nq*m)

!
!  COMPUTE INITIAL SET OF SAMPLE POINTS POINTS
!

call ptinit(n         , m      , x      , ldx     , nx     , &
            nq        , nf     , wrk(ip), wrk(iob), wrk(ic), &
            wrk(id)   , maxnf, delta  , &
            delmax    , pivthr*delta    , lb    , ub    , a , &
            lda       , nclin  , ncnln  , scale , wrk(iscal), &
            wrk(icurw), lenw   , iwrk(1)   , liwrk , inform )


!
!  IF A SECOND POINT FOR WHICH FUNCTION VALUE  CAN BE COMPUTED
!  CANNOT BE FOUND, THEN PRINT A MESSAGE AND STOP
!
if ( inform == -1 ) then
  if ( iprint >= 0 ) write(iout, 2030)
  info = -9
  goto 48
endif
!
!  SOME PARAMETERS DATA PROBLEM IS SET INCORRECTLY
!  IT IS LIKELY THAT SOMETHING IS WRONG WITH CONSTRAINTS
!  PRINT A MESSAGE AND STOP
!

if ( inform == 1 ) then
  if ( iprint >= 0 ) write(iout, 2040)
  info = -7
  return
endif


!
!  IF WE COULD NOT FIND AT LEAST TWO FEASIBLE POINTS
!  TO INITIATE THE PROCESS WE STOP AND PRINT A MESSAGE      
!

if ( inform == 2 ) then
  if ( iprint >= 0 ) write(iout, 2050)
  info = -9
  return
endif
lp  = 1
lv  = iv  - ip + 1
ld  = id  - ip + 1
lob = iob - ip + 1
lc  = ic  - ip + 1

do 47 mjit=1, nmjit

!
!  FIND THE VALUE OF THE MERIT FUNCTION FOR EVERY POINT
!

  call vlinit(n, m, nq, x , base, wrk(ip), wrk(iv), wrk(iob)  , &
              wrk(ic) , ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1), &
              wrk(id), pp, method,  wrk(icurw), &
              lenw)

  fx(1)=wrk(iv+base-1)
!
!  CALL THE DFO SOLVER
!  
!
!  CHECK IF WE RAN OUT OF CUNTION EVALUATIONS
! 
if ( inform == -2 ) then
  if ( iprint >= 0 ) write(iout, 1120)
  info = 2
  call bestpt( n, m, x, fx(1), c, wrk(ip), wrk(iv), wrk(iob), &
             wrk(ic), nq )
  goto 48
endif
  call dfoslv( n  , m  , nclin, ncnln     , x     , fx(1) , &
               c  , nq , base , lp        , lv    , ld    , &
               lob, lc , lb    , ub   , &
               a  , lda, scale, wrk(iscal), pp, info , &
               it , nf , noise, ipar      , rpar  , icurw, &
               wrk(ip) , lwrk-ip+1 , iwrk , liwrk )

  if ( nf >= maxnf ) go to 48
  if ( it >= maxit ) go to 48
  pp = pp * 1.0d1

47 continue  
!
!  UNSCALE THE PROBLEM IF SCALED
!

48 continue

if (scale /= 0) then
  call unscl( n, x, wrk(iscal) )
  call unscl( n, lb, wrk(iscal) )
  call unscl( n, ub, wrk(iscal) )
  do 60 i=1,n
    do 50 j=1,nclin
      a(lda*(i-1)+j)=a(lda*(i-1)+j)/wrk( lscal + i )
50     continue  
60   continue  
endif

!
!  IF DEMANDED, PRINT OUT THE FINAL OUTPUT
!
if ( iprint >= 0 ) then
  write ( iout, 6000 ) pname
  do 70 i = 1, n
    write( iout, 6200 ) xnames(i), x( i )
70   continue
  do 80 i = 1, m
    write( iout, 6300 ) cnames(i), lb(n+nclin+ncnln+i), &
           c( i ), ub(n+nclin+ncnln+i)
80   continue
  write( iout, 6400 )  pname, n, nf, it, fx(1), info
!       WRITE( IOUT, 6600 )  PNAME, N, NF, IT, FX(1), INFO 
endif

return


!
!  NON-EXECUTABLE STATEMENTS
!



1000 format(/' DFO: *** ERROR: NUMBER OF VARIABLES',/ &
        ' DFO:            HAS ILLEGAL VALUE:', i6,/ )
1005 format(/' DFO: *** ERROR: NUMBER OF GENERAL CONSTRAINTS',/ &
        ' DFO:            HAS ILLEGAL VALUE:', i6,/)
1010 format(/' DFO: *** ERROR: NUMBER OF LINEAR CONSTRAINTS',/ &
        ' DFO:            HAS ILLEGAL VALUE:', i6,/)
1020 format(/' DFO: *** ERROR: NUMBER OF NONLINEAR CONSTRAINTS',/ &
        ' DFO:            HAS ILLEGAL VALUE:', i6,/ )
1030 format(/' DFO: *** ERROR: THE UPPER AND LOWER BOUNDS ARE',/, &
        ' DFO:            INCONSISTENT', /)
1040 format(/' DFO: *** ERROR: FUNCTION VALUE CANNOT BE OBTAINED',/, &
        ' DFO:            FOR INITIAL POINT!')
1050 format(/' DFO: *** ERROR: MAXNF IS TOO SMALL!' / &
        ' DFO:            IT SHOULD BE AT LEAST ', i5 )
1060 format(/' DFO: *** ERROR: THE LEADING DIMENSION OF ',/, &
        ' DFO:            LINEAR CONSTRAINTS SHOULD BE AT LEAST', &
                          i6,/)
1070 format('  DFO: *** ERROR: LWRK IS TOO SMALL EVEN TO START!' / &
        ' DFO:             IT IS: ',i10 )
1080 format(/' DFO: *** ERROR: THE LEADING DIMENSION OF ',/, &
        ' DFO:            INITIAL SET "X" SHOULD BE AT LEAST', &
                          i6,/)
1090 format(/' DFO: *** ERROR: THE NUMBER  OF INITIAL ',/, &
        ' DFO:            POINTS "X" SHOULD BE POSITIVE')
1100 format(/' DFO: *** ERROR: AN INPUT PARAMETER ',/ &
        ' DFO:            HAS AN ILLEGAL VALUE:', i6,/ )
1110 format(/' DFO: *** ERROR: MAXNF OR NX IS TOO BIG!' / &
        ' DFO:            MAXNF+NX IT SHOULD BE AT MOST ', i5 )
1120 format(/' DFO: *** WARNING: MAXNF IS TOO SMALL' / &
        ' DFO:              TO PERFORM OPTIMIZATION')
2030 format(/' DFO: CANNOT COMPUTE THE FUNCTION VALUE FOR INITIAL', / &
        ' DFO: INTERPOLATION SET OF AT LEAST TWO POINTS')
2040 format(/' DFO: SOME PARAMETER OR DATA  IS SET INCORRECTLY', / &
        ' DFO  CHECK CONSTRAINT COMPUTATIONS')
2050 format(/' DFO: CANNOT FIND TWO POINTS TO BUILD INITIAL MODEL!', / &
        ' DFO: CHECK FEASIBILITY AND ACCURACY REQUIREMENTS')
3000 format( 4( d12.4, /), 2( i12, /), 6( d12.4, /), 8(i12, /), d12.4, &
        /, i12 )
5000 format(/' DFO: *** WARNING: SCALING IS NOT SUPPORTED ',/, &
        '          FOR PROBLEMS WITH NONLINEAR CONSTRAINTS',/, &
        '          NO SCALING WILL BE PERFORMED',/)
5200 format( ' DFO: *** WARNING:',i4,'-th variable cannot be scaled,' &
                   ,/, &
        '          scaling coefficient is too small or',/, &
        '          too large',/)
4000 format( ' Minimum frobenius norm models will be used', / )
4020 format( ' Models based on incomplete Newton Polynomial basis',/ &
        ' will be used',/ )

6000 format( /, a10,'- OPTIMAL POINT:', / ) 
6200 format( a10, d13.6, / )
6300 format( a10, 3( d13.6), /) 
6400 format( /, 24('*'), 'DFO: FINAL OUTPUT ', 24('*') // &
    ' problem                 =      ', a8 ,   /, &
    ' # variables             =      ', i10,   /, &
    ' # objective functions   =      ', i10,   /, &
    ' # iterations            =      ', i10,   /, &
    ' Final f                 =      ', d15.7, /, &
    ' Exit code               =      ', i10,   /, &
     65('*') / )
6600 format( a10, i10, 1x, i10, 1x, i10, 1x,  d14.7, 1x, i10 )
8000 format( 'SCALING COEFFICIENTS:',/)
8010 format( d14.7, 1x )


end
!


subroutine dfoslv( n  , m   , nclin, ncnln, x   , fx   , &
                   c  , nq  , base , ip   , iv  , id   , &
                   iob, ic  , lb  , ub   , &
                   a  , lda , scale, scal , pp  , exit , &
                   it , nf  , noise, ipar , rpar, ibegw, &
                   wrk, lwrk, iwrk , liwrk)  
!
!  ******************************************************************
!  THIS SUBROUTINE MINIMIZES A NONLINEAR OBJECTIVE FUNCTION
!  SUBJECT TO LINEAR AND (POSSIBLY) NONLINEAR CONSTRAINTS
!  AND SIMPLE BOUNDS, WITHOUT USING DERIVATIVES OF 
!  THE OBJECTIVE FUNCTION.
!
!  THIS ALGORITHM  IS BASED ON USING QUADRATIC INTERPOLATION OF THE
!  OBJECTIVE FUNCTION IN COMBINATION WITH TRUST REGION FRAMEWORK
!
!
!  PARAMETERS
!  ==========
!
!  (INPUT) 
!
!    X         : ARRAY OF NX  INITIAL POINTS, GIVEN BY THE USER
!
!    LDX       : LEADING DIMENSION OF ARRAY 'X' (LDX>= N)
!
!    FX        : ARRAY OF FUNCTION VALUES AT THE INITIAL POINTS
!                OF LENGTH AT LEAST 1
!
!    NX        : NUMBER OF INITIAL POINTS FOR WHICH THE VALUE IS
!                PROVIDED.
!
!    N         : PROBLEM DIMENSION
!
!    NCLIN     : NUMBER OF LINEAR CONSTRAINTS
!
!    NCNLN     : NUMBER OF NONLINEAR CONSTRAINTS
!
!    NOISE     : NOISE(1) - ABSOLUTE NOISE LEVEL
!                NOISE(2) - RELATIVE NOISE LEVEL
!                ( IF KNOWN BY THE USER )
!
!    IPAR      : ARRAY OF INITIAL INTEGER PARAMETERS
!
!    RPAR      : ARRAY OF INITIAL REAL PARAMETERS
!
!    LB        : ARRAY OF LOWER BOUNDS OF LENGTH >= (N+NCLIN+NCNLN)
!
!    UB        : ARRAY OF UPPER BOUNDS OF LENGTH >= (N+NCLIN+NCNLN)
!
!    A         : MATRIX OF LINEAR CONSTRAINTS,( DIMENSIONS: LDA X N)
!
!    LDA       : LEADING DIMENSION OF MATRIX A, HAS TO BE >= MAX(1, NCLIN)
!
!    WRK       : REAL WORKSPACE ARRAY OF SIZE LWRK
! 
!    IWRK      : INTEGER WORKSPACE ARRAY OF SIZE LIWRK
!
!    (OUTPUT)
!
!     X        : THE OPTIMAL (OR BEST SO FAR)  POINT 
!
!     FX       : THE VALUE OF OBJECTIVE FUNCTION AT X
!
!     EXIT     : THE EXIT INFORMATION:
!                0     SUCCESSFUL MINIMIZATION
!                1     TOO MANY FUNCTION EVALUATIONS
!                2     TOO MANY ITERATIONS
!               -3     REAL WORKSPACE IS TOO SMALL
!               -4     INTEGER WORKSPACE IS TOO SMALL
!               -7     MINIMIZATION OVER TRUST REGION FAILED
!               -9     CANNOT FIND TWO POINTS TO BUILD INITIAL MODEL
!     IT       : THE NUMBER OF ITERATIONS
!
!     NF       : THE NUMBER OF FUNCTION EVALUATIONS
!  ********************************************************************
!

!
!  SUBROUTINE PARAMETERS
!
integer          n, m , nclin    , ncnln, it  , exit  , nf   , &
                 maxit, maxnf    , liwrk, lwrk, iwrk( liwrk ), &
                 varnt, ipar( 8 ), scale, lda , ic, iob, ibegw

double precision x(n)   , noise( 2 ), rpar( 8 )  , delmin , &
                 lowbnd , fx        , wrk( lwrk ), scal(n), &
                 lb( * )    , &
                 ub(*), a( lda*n ) , &
                 pp , c(*)



!
!  COMMON VARIABLES
!
include 'dfo_model_inc.inc'
!
!  LENGTH OF ARRAYS
!
integer          lpoly, lpnts, lvalue, lptint, lvlint, lconvl
common /rpart/   lpoly, lpnts, lvalue, lptint, lvlint, lconvl
save /rpart/

!
!  PROBLEM CONTROL PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /

!
!  INTERPOLATION CONTROL PARAMETERS
!      
integer          npmin, layer, effort
common / opti /  npmin, layer, effort
save / opti /


!
!
!  LOCAL VARIABLES
!  ---------------
!


integer          impr  , nstop , base  , i     , nq    , iter  , &
                 lp    , lv    , lk    , ik    , lg    , lh    , &
                 lcurw , ig    , ih    , ip    , iv    , icurw , &
                 ibeg  , lsp2in, next  , inq   , lnq   , ibase , &
                 j     , obase , ipend , lbase , lrs   , lis   , &
                 ld    , id    , ii    , icuriw, isp2in, lpiv  , &
                 inp   , nind  , ivi   , inext , ipoly , ifadd , &
                 ifxg  , ipi   , dd    , rsneed, lnp   , lpi   , &
                 mntris, npslis, impdrs, lvi   , ipiv  , lin2sp, &
                 iin2sp, inform, oldnq , lcuriw, ptnwrs, mntrrs, &
                 npslrs, mdblrs, gtdsrs, mdblis, neqcon, lob   , &
                 lc    , lcq   , lcl   , lcc   , lci   , k     , &
                 kk    , ici   , icq   , icl   , icc   , nindr , &
                 method, stopcnt, howstop

logical          fail  , okgeom, linear, iferr , usemer

double precision delta , delmax, anoise, rnoise, ratio , rho   , &
                 fbase , mval  , prered, noisef, fnq   , theta , &
                 ft    , stnorm, rhow  , delthr, val   , snorm , &
                 kappa , pivthr, addthr, pivt  , addt  , xchthr, &
                 kmod  , del   , odelta, whenstop

double precision zero  , one   , two   , ten   ,  half , thous
parameter      ( zero  = 0.0d0 , ten = 10.0d0  ,  two   =  2.0d0 )
parameter      ( half  = 0.5d0 , one = 1.0d0   ,  thous =  1.0d3 )

!
!     TRUST REGION METHOD PARAMETERS
!
double precision rhomin, rhoinc, rhojmp,  maxjmp

!
parameter      ( rhomin = 0.05d0, rhojmp =  0.9d0 ) 
parameter      ( rhoinc = 0.25d0, maxjmp =  1.0d2 )

!
integer          seed
parameter      ( seed = 926865394 )

!
!     SUBROUTINES AND FUNCTIONS CALLED:
!
!       APPLICATION:       PTINIT, MDBLD , MINTR , FUN   , GETDIS, 
!                          PTREPL, PTEXCH, NEXTNP, EVALX , IMPMOD
!                          SCL   , UNSCL , SHIFT , UNSHFT, DNRMNF,
!                          EVALNP, BESTPT, SWAPNP, MVALUE
!       FORTRAN SUPPLIED:  MAX   , MIN   , ABS
!       BLAS:              DCOPY 
!


!
!  SET VARIOUS PARAMETERS
!

!
!  MAXIMUM NUMBER OF ITERATIONS
!
maxit  = ipar( 1 )
!
!  MAXIMUM NUMBER OF FUNCTION EVALUATIONS
!
maxnf  = ipar( 2 )
!
!  THE DESIRED MINIMUM NUMBER OF POINT IN AN INTERPOLATION SET
!
npmin  = ipar( 3 )
!
!  THE CONTROL PARAMETER OF MAXIMUM DISTANCE OF POINTS TO THE BASE
!
layer  = ipar( 4 )
!
!  THE LEVEL OF EFFORT OF INTERPOLATION COMPUTATIONS
!  
effort = ipar( 5 )
!
!  FLAG RESPONSIBLE FOR THE CHOICE OF INCOMPLETE MODEL: 
!  1 - MINIMUM FROBENIUS NORM, 2 - DERIVED FROM NEWTON POLYNOMIAL BASIS
!
varnt  = ipar( 6 )
!
!  METHOD USED TO HANDLE CONSTRAINTS 3,4 - PENALTY FUNCTION OPTIMIZATION
!                                    1,2 - MODELING CONSTRAINTS
!
method  = ipar( 7 )
!
!  WHICH STOPPING CRITERIA TO USE   1 - SIZE OF THE TRUST REGION, TENDS TO
!                                       PRODUCE MORE ACCURATE RESULTS
!                                   2 - SLOW PROGRESS, TENDS TO STOP SOONER
!
howstop  = ipar( 8 )
!
!  THE INITIAL TRUST REGION RADIUS
!
delta  = rpar( 1 )
!
!  MAXIMUM TRUST REGION RADIUS
!
delmax = rpar( 2 )
!
!  MINIMUM TRUST REGION RADIUS (STOPPING CRITERIA)
!
delmin = rpar( 3 )
!
!  RELATIVE (W.R.T THRUST REGION RADIUS SIZE) PIVOT THRESHOLD
!
pivt   = rpar( 4 )
!
!  RELATIVE THRESHOLD OF ADDING A TRUST REGION POINT WITH BAD RHO 
!
addt   = rpar( 5 )
!
!  MINIMUM REQUIRE FACTOR OF PIVOT IMPROVEMENT WHILE REPLACING
!  AN INTERPOLATION POINT BY A  GEOMETRY POINT
!
xchthr = rpar( 6 )
!
!  LOWER BOUND ON DESIRED MINIMUM FUNCTION VALUE 
!
lowbnd = rpar( 7 )
!
!  THRESHHOLD FOR CONSIDERING INSUFFICIENT PROGRESS
!
whenstop = rpar( 8 )
!
!  ABSOLUTE NOISE LEVEL
!
anoise = noise( 1 )
!
!  RELATIVE NOISE LEVEL
!
rnoise = max( noise( 2 ), mcheps )




!  ------------------------------------------------------------
!  PARTITION THE REAL WORKSPACE
!  ------------------------------------------------------------

!
!  MAXIMUM NUMBER OF INTERPOLATION POINTS/POLYNOMIALS 
!  IN A QUADRATIC MODEL
!
dd     = ((n + 1)*(n + 2))/2
!
!  MAXIMUM LENGTH OF ARRAY OF  NEWTON FUNDAMENTAL POLYNOMIALS
! 
lpoly  = dd*(n+1)*n/2+(n+1)*n+1
!
!  MAXIMUM LENGTH OF ARRAY OF INTERPOLATION POINTS
!
lptint = dd*n
!
!  MAXIMUM LENGTH OF ARRAY OF FUNCTION  VALUES OF INTERP. POINTS
!
lvlint = dd



! 
!  POINTERS OF ARRAYS IN REAL WORKSPACE
!  -----------------------------------


!
!  POINTER TO NEWTON FUNDAMENTAL POLYNOMIALS
!
inp = ibegw
!
!  POINTER TO ARRAY OF INTERPOLATION POINTS
!
ipi = inp + lpoly
!
!  POINTER TO FUNCTION VALUES AT INTERPOLATION POINTS
!
ivi = ipi + dd*n
!
!  POINTS TO CONSTRAINT VALUES AT INTERPOLATION POINTS
!
ici = ivi + dd
!
!  POINTER TO PIVOT VALUES ASSOCIATED WITH INTERPOLATION POINTS
!
ipiv= ici + dd
!
!  POINTER TO LINEAR TERMS OF THE QUADRATIC MODEL
!
ig  = ipiv + dd
!
!  POINTER TO QUADRATIC TERMS OF THE QUADRATIC MODEL
!
ih  = ig  + n
!
!  POINTER TO CONSTANT TERMS OF THE MODEL OF CONSTRAINTS
!
icc = ih  + n * n
!
!  POINTER TO LINEAR TERMS OF THE MODEL OF CONSTRAINTS 
!
icl = icc + m
!
!  POINTER TO QUADRATIC TERMS OF THE MODEL OF CONSTRAINTS 
!
icq = icl + m * n
!
!  POINTER TO AUXILIARY SPACE FOR NEW POINTS
!
ik  = icq  + n * n * m
!
!  POINTER TO BEGINNING OF WORKING SPACE
!
icurw = ik + n
!
!  THE LENGTH OF THE REMAINING SPACE
!
lrs    = lwrk - icurw + 1
!
!  SPACE REQUIRED BY 'MINTR' SUBROUTINE
!
mntrrs = 3*(n+nclin+m+ncnln) + (n+1)*max(1,ncnln+m) + n*n + n
!
!  SPACE REQUIRED BY 'NPSOL' SUBROUTINE
!
npslrs = 2*n*n + n*nclin + 2*n*(ncnln+m) + 20*n + &
         11*nclin + 21*(ncnln+m)
!
!  SPACE REQUIRED BY 'MDLBLD' SUBROUTINE
!
mdblrs = dd*dd + (dd-1)*n + dd - 1
!
!  SPACE REQUIRED BY 'PTNEW' SUBROUTINE
!
ptnwrs = 2*n + n*n
!
!  SPACE REQUIRED BY 'GETDIS' SUBROUTINE 
!
gtdsrs = n
!
!  SPACE REQUIRED BY 'IMPMOD' SUBROUTINE
!
impdrs = n
!
!  CHECK IF THE REMAINING REAL SPACE IS SUFFICIENT  
!  BY COMPARING IT  WITH THE "CRITICAL PATH"
!  -----------------------------------------------

rsneed = ptnwrs + mntrrs + npslrs + impdrs
if ( lrs <  max( rsneed, mdblrs ) ) &
   then
   if ( iprint >= 0 ) write( iout, 2000 ) &
          max(rsneed, mdblrs) - lrs
   exit = -3
   return
endif

!
!  HELPFUL 'SHIFTED' POINTERS TO ARRAYS
!
ld    = id    - 1
lp    = ip    - 1
lv    = iv    - 1
lc    = ic    - 1
lob   = iob   - 1
lg    = ig    - 1
lh    = ih    - 1
lcc   = icc   - 1
lcl   = icl   - 1
lcq   = icq   - 1
lp    = ip    - 1
lv    = iv    - 1
lpi   = ipi   - 1
lvi   = ivi   - 1
lci   = ici   - 1
lpiv  = ipiv  - 1
lk    = ik    - 1
ld    = id    - 1
lnp   = inp   - 1
lcurw = icurw - 1



!  ---------------------------------------------------------------
!  PARTITION THE INTEGER WORKSPACE
!  ---------------------------------------------------------------

!
!  POINTER TO ARRAY OF 0/1 INDICATORS: 1 - IF A SAMPLE POINT 
!  IS INCLUDED IN THE CURRENT INTERPOLATION SET, 0 - OTHERWISE
!
isp2in = 1
!
!  POINTER TO ARRAY INDICATING FOR EVERY INTERPOLATION POINTS 
!  ITS LOCATION IN SAMPLE SET
!
iin2sp = isp2in+maxnf
!
!  POINTER TO THE INTEGER WORKING SPACE
!  ------------------------------------

icuriw = iin2sp  + dd

!
!  LENGTH OF THE REMAINING INTEGER SPACE
!
lis = liwrk - icuriw + 1
!
!  INTEGER SPACE REQUIRES BY 'MINTR' PROCEDURE
!
mntris = n + nclin + ncnln + m
!
!  INTEGER SPACE REQUIRES BY 'NPSOL' PROCEDURE
!
npslis = 3*n + nclin + 2*(ncnln + m)
!
!  INTEGER SPACE REQUIRES BY 'MDBLD' PROCEDURE
!
mdblis = dd
!
!  CHECK IF THE REMAINING INTEGER SPACE IS SUFFICIENT
!  BY COMPARING IT WITH THE "CRITICAL PATH"
!
if ( lis < max(mntris + npslis, mdblis) ) then
   if ( iprint >= 0 ) write( iout, 2020 ) &
        max(mntris + npslis, mdblis) - lis
   exit = -4
   return
endif
!
!  HELPFUL 'SHIFTED' POINTERS TO INTEGER ARRAYS
!
lsp2in = isp2in - 1
lin2sp = iin2sp - 1
lcuriw = icuriw - 1

!
!  SET PARAMETERS FOR /MERIT/ COMMON BLOCK
! 
ncon=m
nlin=nclin
nnln=ncnln
penpar=pp
do 5 i=1, m
  conl(i)=lb(n+nclin+ncnln+i)
  conu(i)= ub(n+nclin+ncnln+i)
5 continue  



!  ----------------------------------------------------
!  INITIAL INTERPOLATION SET IS BUILD HERE
!  ----------------------------------------------------





!
!  CALL SUBROUTINE 'NBUILD' TO COMPUTE INITIAL: 
!       INTERPOLATION SET AND 
!       NEWTON FUNDAMENTIAL POLYNOMIALS
!       DISTANCES OF SAMPLE POINTS TO BASE POINT
!       PIVOT VALUES ASSOCIATED WITH INTERPOLATION SET
!       ARRAYS IN2SP AND SP2IN CONNECTING SAMPLE SET WITH INTERP. SET
!


if ( iprint >= 2 ) write( iout, 8000 ) 
pivthr = pivt*delta
call dcopy (n, wrk(ip+(base-1)*n), 1, wrk(ipi), 1)
wrk(ivi) = wrk(iv + base - 1)
iwrk(iin2sp)=base
iwrk( isp2in + base - 1 ) = 1
base = 1
!
!  COMPUTE THE NUMBER OF LINEARLY INDEPENDENT INEQUALITY CONSTRAINTS
!           
call dcopy ( n, wrk(ip+(base-1)*n), 1, wrk(ik), 1 )
call unshft( n, x, wrk(ik) )
call intdim( wrk(ik)   , n   , nclin , ncnln , neqcon, &
             a, lda    , lb  , ub    , pivthr, &
             wrk(icurw), lrs , iwrk(icuriw)  , lis    )
!
!  COMPUTE THE BASIS OF POLYNOMIALS
!
call nbuild( wrk(inp)    , wrk(ip)     , wrk(iv)     , wrk(ipi), &
             wrk(ivi)    , iwrk(isp2in), iwrk(iin2sp), nind    , &
             n           , base        , wrk(id)     , delta   , &
             pivthr      , wrk(ipiv)   , neqcon )   



!
!  IF THE BASIS CONSISTS ONLY OF A SINGLE POLYNOMIAL (CONSTANT TERM)
!  THEN SOMETHING IS WRONG. WE STOP AND PRINT A MESSAGE      
!
if (nind<2) then
  if ( iprint >= 0 ) write(iout, 2060)
  exit=-9
  return
endif


!
!  INITIALIZE  POINTERS TO BASE POINT AND TO THE LAST SAMPLE POINT
!
lbase = lpi + ( base - 1 ) * n
ibase = lbase + 1
lnq   = lp + (nq-1)*n
inq   = lnq + 1

!
!  SOME MORE INITIALIZATION
!
oldnq = nq
obase = base
odelta=delta
nstop = 0
stopcnt = 0
ratio = two 
rho   = zero
impr  = 3



!
! ===========================================================================
!
!                    MAIN ITERATION
!
! ===========================================================================
!
do 200 iter = 1, maxit

   if ( iprint >= 2 ) write( iout, 8010 ) iter
!
!  SET THE CURRENT VALUES FOR PIVOT-RELATES THRESHOLDS
!       
   pivthr = pivt*min(one,delta)
   addthr = addt*min(one,delta)


!
!  SELECT AN INTERPOLATION POINT WITH THE SMALLEST FUNCTION VALUE AS BASE
! 
   fbase  = wrk( lv + iwrk( lin2sp + base ))
   do 10 i = 1, nind
     if (( i <= n+1-neqcon .or. i > n+1 ) .and. &
     ( wrk( lv + iwrk( lin2sp + i )) <= fbase - ten * mcheps)) &
          then
       fbase = wrk(lv + iwrk( lin2sp + i ))
       base  = i
     endif
10    continue
   if ( iprint >= 2 ) write( iout, 8020 ) base, nind

!
!  IF THE BASE POINT HAD CHANGED, WE NEED TO RECOMPUTE THE DISTANCES
!  FROM BASE TO ALL SAMPLE POINTS IN "POINTS" 
!
   if ( obase /= base ) then 
      lbase = lpi + ( base - 1 ) * n
      ibase = lbase + 1
      obase = base
      call getdis( n, nq, 0 , wrk( ip ) , iwrk(lin2sp+base), &
                   wrk( id ), wrk(icurw), lrs )
   endif

!
!  THE MAXIMUM POSSIBLE NUMBER OF INDEPENDENT INTERPOLATION POINTS
!  DEPENDS ON THE ACTUAL DIMENSION OF THE FEASIBLE SET, RATHER THEN
!  ON THE DIMENSION OF THE WHOLE SPACE. SUBROUTINE 'INTDIM' COMPUTES
!  THE DIMENSION OF THE SPACE SPANNED BY EQUALITY CONSTRAINTS AT
!  THE BASE. WE CONSIDER A CONSTRAINT AS  EQUALITY IF THE DIFFERENCE
!  BETWEEN ITS UPPER AND LOWER BOUNDS IS LESS THAN PIVOT THRESHOLD
!
!  NEQCON - NUMBER OF LINEARLY INDEPENDENT EQUALITY CONSTRAINTS
!
   if ( odelta /= delta ) then
     call dcopy ( n, wrk(ibase), 1, wrk(ik), 1 )
     call unshft( n, x, wrk(ik) )
     call intdim( wrk(ik)   , n   , nclin , ncnln , neqcon, &
                  a, lda    , lb  , ub    , pivthr, &
                  wrk(icurw), lrs , iwrk(icuriw)  , lis    )
     odelta = delta
   endif           

!
!  SET 'OKGEOM' TO FALSE BY DEFAULT
!  SET 'LINEAR' TO TRUE IF THE INTERPOLATION IS FULLY LINEAR
!  (I.E., THERE IS MAXIMUM POSSIBLE NUMBER OF LINEAR NEWTON
!   POLYNOMIALS THAT CAN FIT INTO THE FEASIBLE SET. IN CASE 
!   WHEN THERE ARE NO EQUALITY CONSTRAINTS NIND SHOULD BE AT
!   LEAST N+1, OTHERWISE - AT LEAST N+1 - NUMBER OF EUALITY
!   CONSTRAINTS
!

   okgeom= .false.
   linear=( nind > n - neqcon ) 
!
!  IF THERE ARE EQUALITY CONSTRAINTS, THEN UPON REACHING THE MAXIMUM
!  POSSIBLE NUMBER OF POINTS FOR LINEAR INTERPOLATION WE COMPLETE
!  THE LINEAR INTERPOLATION TO N+1 ELEMENTS WITH DUMMY POLYNOMIALS
!  AND POINTS. THIS IS DONE BECAUSE THE PROGRAM IS WRITTEN UNDER
!  ASSUMPTION THAT QUADRATIC BLOCK IS NOT CONSTRAUCTED UNTILL THE
!  LINEAR BLOCK IS FULL. BY COMPLETETING THE LINEAR BLOCK BY DUMMY
!  ELEMENTS WE MAKE SURE THAT WE CAN MOVE ON TO QUADRATIC ELEMENTS.  
!
   
   if ( neqcon > 0 .and. nind == n + 1 - neqcon ) then
     call complt( wrk(inp), wrk(ivi), wrk(ipiv), n, nind, lpoly )
   endif

!  --------------------------------------
!  PRINT OUT SUMMARY OF CURRENT ITERATION
!  --------------------------------------
!
!  NINDR - NUMBER OF ELEMENTS IN THE INTERPOLATION (WITHOUT THE DUMMY)  
!
   
   if ( iprint >= 1 ) then
     nindr = nind
     if ( nind > n + 1 - neqcon ) nindr = nind - neqcon
!           IF (ITER==( ITER/10 )*10 +1 ) WRITE( IOUT, 1000 )
     write( iout, 1000 )
     write( iout, 1010) iter, nf, fbase, rho, delta, impr, nindr, &
                        neqcon
   endif
!
!  STOPPING TEST: IF THE TRUST REGION RADIUS IS SMALL ENOUGH, WE STOP
!                 IF THE BEST VALUE HAS REACHED THE LOWER BOUND, WE STOP
!
   if ( iprint >= 2 ) write( iout, 8030 ) delta 
   if ( ( howstop==2 .and. stopcnt>=2) &
         .or. (delta < delmin) ) then 
     exit = 0
     call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
                 wrk(ic), nq)
     it   = iter
     return               
   endif
!
!  IN DEBUGGING MODE, CHECK THE ACCURACY OF NEWTON FUNDAMENTAL POLYNOMIALS
!
 if ( iprint >= 3 ) then
   do 30 i = 2, nind
     if( i <=  n + 1 - neqcon) then
       ii = min( n + 1 - neqcon, nind )
     else 
       ii = nind
     endif
     do 20 j = 1, ii
      if (j <= n + 1 - neqcon .or. j > n + 1 ) then
        call  evalnp( val, wrk(ipi), j, wrk(inp), i, n, &
                      lptint , lpoly )
        if ( (i==j .and. abs(val-one) > 1d-10) .or. &
             (i/=j .and. abs(val)     > 1d-10) ) then
          write(iout,1020) i, j, val
        endif
      endif
20      continue
30    continue
 endif

 if ( method <= 2 ) then
!  ---------------------------------------------------------------------
!  BUILD THE QUADRATIC MODEL OF THE OBJECTIVE FUNCTION
!  ---------------------------------------------------------------------

   if ( iprint >= 2 ) write( iout, 8040 ) 
   do 31 j=1, nind
     ii = iwrk(lin2sp + j)
     wrk(lvi + j) = wrk( lob + ii ) 
31    continue  
   call mdbld( kappa       , wrk(ig) , wrk(ih) , n   , nind      , &
               wrk(ipi)    , wrk(ivi), wrk(inp), base, varnt     , &
               lpoly       , lptint  , neqcon  , wrk(icurw) , lrs, &    !   
               iwrk(icuriw), lis     )

!  ---------------------------------------------------------------------
!  BUILD THE QUADRATIC MODEL OF THE CONSTRAINTS  
!  ---------------------------------------------------------------------

   do 33 i = 1, m
     do 32 j=1, nind
       ii = (iwrk(lin2sp + j)-1)*m
       wrk(lci + j) = wrk(lc + ii + i) 
32     continue  
     call mdbld(wrk(lcc+i), wrk(icl+(i-1)*n) , wrk(icq+(i-1)*n*n), &    ! 
                n         , nind             , wrk(ipi)          , &
                wrk(ici)  , wrk(inp)  , base , varnt             , &
                lpoly     , lptint           , neqcon            , &    ! 
                wrk(icurw), lrs              , iwrk(icuriw)      , &    ! 
                lis       )

33   continue

 else 
!  ---------------------------------------------------------------------
!  BUILD THE QUADRATIC MODEL OF THE MERRIT FUNCTION
!  ---------------------------------------------------------------------

   if ( iprint >= 2 ) write( iout, 8040 ) 
   do 34 j=1, nind
     ii = iwrk(lin2sp + j)
     wrk(lvi + j) = wrk( lv + ii ) 
34   continue  
   call mdbld( kappa       , wrk(ig) , wrk(ih) , n   , nind      , &
               wrk(ipi)    , wrk(ivi), wrk(inp), base, varnt     , &
               lpoly       , lptint  , neqcon  , wrk(icurw) , lrs, &    !   
               iwrk(icuriw), lis     )
 endif
!
!  IN DEBUGGING MODE CHECK THE ACCURACY OF QUADRATIC INTERPOLATION
!
 if ( iprint >= 3 ) then
   do 41 i=1, m
     do 40 j = 1, nind
       if ( j <= n + 1 - neqcon .or. j > n + 1 ) then
         ii = (iwrk(lin2sp + j)-1)*m
         wrk(lci + j) = wrk(lc + ii + i) 
         val = wrk(lci + j) - (wrk(lcc+i) + &
               mvalue( n, wrk( ipi+(j-1)*n), &
               wrk(icl+(i-1)*n), wrk( icq+(i-1)*n*n ), &
               wrk( icurw ), lrs))
         if ( abs(val) > 1.0d-6 ) write(iout,1030)  j, val 
       endif
40      continue
41    continue  
 endif     
 

!  -------------------------------------------------------------------------
!
!  FIND THE NEXT CANDIDATE POINT FOR SAMPLING 
!  BY MINIMIZING THE MODEL OVER TRUST REGION (INTERSECTED WITH FEASIBLE SET)
!
!  -------------------------------------------------------------------------


!
!  SET THE RADIUS TO DELTA
!        

 del = delta

!
!  -------------------  MODIFY THE MODEL  ----------------------------
!  THE MODEL INTERPOLATES A SHIFTED SET OF POINTS, SINCE THE FEASIBLE
!  SET IS NOT SHIFTED, WE MINIMIZE A MODIFIED MODEL, WHICH INTERPOLATES
!  NON-SHIFTED SET OF POINTS OVER NOT-SHIFTED REGION
!  THE PARAMETERS OF THE MODIFIED MODEL ARE STORED IN KMOD, GMOD, HMOD  
!  THE MODIFIED MODEL INTERPOLATES POINTS X_i+X, WHERE X_i ARE INTERP. POINTS,
!  THEN THE MODIFIED MODEL IS COMPUTED AS
!
!                HMOD <-- H
!                GMOD <-- G - H'*X
!                KMOD <-- KAPPA - G'*X - 0.5X*H'*X
!  -------------------------------------------------------------------
!


 kmod=kappa
 do 51 i=1,n
   gmod(i)= wrk(ig+i-1)
   kmod   = kmod-wrk(ig+i-1)*x(i)
   ii=(i-1)*n
   do 50 j=1,n
     hmod(i,j)= wrk(ih+ii+j-1)
     gmod(i)  = gmod(i)-hmod(i,j)*x(j)
     kmod     = kmod+half*hmod(i,j)*x(j)*x(i)
50    continue
51    continue   

 
 if ( method <= 2 ) then
   do 60 k=1, m
     ccon(k)=wrk(lcc+k)
     kk=(k-1)*n
     do 56 i=1,n
       lcon(kk+i)= wrk(lcl+kk+i)
       ccon(k)   = ccon(k)-wrk(lcl+kk+i)*x(i)
       ii=(i-1)*n
       do 55 j=1,n
         qcon(kk*n+ii+j)= wrk(lcq+kk*n+ii+j)
         lcon(kk+i)  = lcon(kk+i)-qcon(kk*n+ii+j)*x(j)
         ccon(k)     = ccon(k)+half*qcon(kk*n+ii+j)*x(j)*x(i)
55        continue
56      continue   
60    continue 
 endif 
!
!  PRINT OUT MODEL COEFFICIENT IF DEMANDED
!
 if ( iprint >= 3 ) then
   write( iout, 8050 )
   write( iout, 8051 )
   do 62  i = 1, n
     write(iout, 8052)
     do 61 j = 1, n
       write(iout, 8053) hmod(i,j)
61      continue
62    continue
   write( iout, 8054 )
   do 63 j = 1, n
     write(iout, 8055) gmod(j)
63    continue
   write( iout, 8056 ) kmod
   write( iout, 8060 ) 
 endif
!
!  SET STARTING POINT FOR MINIMIZATION TO THE SHIFTED BASE POINT
!  AND SHIFT IT BACK TO ITS ORIGINAL POSITION
!  ( WE HAVE TO DO IT SINCE THE FEASIBLE SET IS NOT SHIFTED )
!

70  call dcopy ( n, wrk(ibase), 1, wrk(ik), 1 )
 call unshft( n, x, wrk(ik) )



!
!  MINIMIZE THE MODIFIED MODEL 
!

 if ( iprint >= 2 ) write( iout, 8065 ) 
 if ( method <= 2 ) then
   call mintr (n         , wrk( ik ), mval        , del  , lb, &
               ub        , a        , lda         , nclin,m+ncnln, &
               wrk(icurw), lrs      , iwrk(icuriw), lis  , inform, &
               1    )
 else 
   call mintr (n         , wrk( ik ), mval        , del  , lb, &
               ub        , a        , lda         , nclin, ncnln , &
               wrk(icurw), lrs      , iwrk(icuriw), lis  , inform, &
               method    )
 endif
!
!  WARNING!! BE CAREFULL, DON'T USE WRK(ICURW) INTIL FMERIT IS CALLED
!  THE CURRENT VALUES OF THE CONSTRAINED ARE STORED THERE
!

!
!  IF NO FEASIBLE POINT CAN BE FOUND FOR THE MODELS OF  CONSTRAINTS
!  THEN RUN MINTR AGAIN MINIMIZING UNCONSTRAINED MERIT FUNCTION
!

 if ( method == 2 .and. inform== 3 ) then 
   usemer=.true.
   call mintr (n         , wrk( ik ), mval        , del  , lb, &
               ub        , a        , lda         , nclin, ncnln , &
               wrk(icurw), lrs      , iwrk(icuriw), lis  , inform, &
               2         )
 else
   usemer=.false.
   if ( inform == 3 ) inform = 0
 endif
!
!  IF NO FEASIBLE POINT FOUND, THEN CONSIDER REDUCTION TO BE ZERO
!

 if (inform == 2 ) then
   val = fbase
   goto 81
 endif 

!
!  IF THE MINIMIZATION FAILED, THEN RECORD THE BEST POINT AND EXIT
!
 if (inform/=0) then
   call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
               wrk(ic), nq)
   it    = iter
   exit=-7
   return
 endif
!
!  COMPUTE THE PREDICTED VALUE OF THE MERIT FUNCTION
!
 if ( method == 1 .or. (method == 2 .and. .not. usemer) ) then
   call fmerit(m, val, mval+kmod, wrk(icurw), &
               ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1), pp, &
               method )
 else
   val=mval+kmod
 endif
!
!  SHIFT THE NEW POINTS
! 
 call shift ( n, x, wrk(ik))

!
!  COMPUTE THE DISTANCE FROM THE NEW POINT TO THE BASE
!
 do 80 i=1,n
   wrk(lcurw+i)=wrk(lk+i)-wrk(lbase+i)
80  continue

 snorm  = dnrmnf( n, wrk( icurw ) )

!
!  COMPUTE PREDICTED REDUCTION
!
81  prered = fbase - val

!
!  IF THE MODEL REDUCTION IS SMALL, THEN WE DO NOT SAMPLE FUNCTION
!  AT THE NEW POINT. WE THEN WILL TRY TO IMPROVE THE MODEL.
!
 noisef = half * max( anoise, abs( fbase ) ) * rnoise

 if ( prered < max( noisef,delmin*1.0d-2) .or. &
      prered < delmin*abs(fbase) .and. snorm < delmin .or. &
      prered < cnstol*abs(fbase)*1.0d-3 ) then
   if ( iprint >= 2 ) write( iout, 8067 ) prered
   rho   = -thous
!
!  SET INDICATORS THAT NO NEW POINT WAS ADDED TO INTERPOLATION SET
!  AND NO OTHER POINT WAS COMPUTED TO IMPROVE GEOMETRY 
!
   ifadd = 0
   ifxg  = 0

!
!  IF THE MODEL REDUCTION IS SUFFICIENT, THEN SAMPLE THE FUNCTION
!  AT THE NEW POINT
!
 else
!
!  RECORD THE NEW POINT INTO THE SAMPLE SET
!
   call dcopy(n, wrk(ik),1,wrk(inq+n),1)
!
!  IF THE NUMBER OF FUNCTION EVALUATIONS EXCEEDED THE LIMIT, RECORD
!  THE BEST POINT AND EXIT
!
   nf = nf + 1
   if ( nf > maxnf ) then
     call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
                 wrk(ic), nq)
     it   = iter
     exit = 1
     return
   endif

!
!  COMPUTE THE FUNCTION VALUE AT THE NEW POINT
!
   call unshft(n, x, wrk(ik))
   if ( scale /= 0 ) call unscl( n, wrk(ik), scal )
   call fun( n, m, wrk( ik ), wrk(iob+nq), wrk(ic + nq*m), iferr )
   if ( scale /= 0 ) call scl( n, wrk(ik), scal )
!
!  IF FUNCTION COMPUTATION FAILS, THEN REDUCE THE RADIUS FOR 
!  MINIMIZATION TO HALF OF THE DISTANCE BETWEEN THE NEW POINT 
!  AND THE BASE. REPEAT THE MINIMIZATION
!
   if (iferr) then
     del  = snorm/2
     if ( del >= delmin ) then 
       goto 70
     else
       rho   = -thous
       ifadd = 0
       ifxg  = 0
       goto 170
     endif
   endif
!
!  IF FUNCTION VALUE IS COMPUTED, THEN RECORD THE VALUE AND DO
!  APPROPRIATE BOOKKEEPING AND UPDATES
!
   call fmerit(m, fnq, wrk(iob+nq), wrk(ic + nq*m), &
               ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1), pp, &
               method )
   nq     = nq + 1
   lnq    = lp + ( nq - 1 ) * n
   inq    = lnq + 1
   lpnts  = lpnts+n
   lvalue = lvalue+1
   lconvl = lconvl+m
   impr   = 10
   wrk ( lv + nq ) = fnq
   next            = nq
   wrk ( ld + nq ) = snorm

!  -----------------------------------------------------------
!  COMPUTE RHO = ACHIEVED REDUCTION VS. PREDICTED REDUCTION
!  -----------------------------------------------------------

   rho = ( fbase - fnq ) / prered
   if ( iprint >= 2 ) write( iout, 8080 ) rho

!  ------------------------------------------------------------
!  COMPUTE RELATIVE ACHIEVED REDUCTION AND CHECK IF IT IS SMALL
!  ------------------------------------------------------------
   val = ( fbase - fnq ) /(one+abs(fbase))
   if ( ( val > 0 ).and. ( val < whenstop).and. &
        ( rho > rhomin )) then
     stopcnt=stopcnt+1
   else
     stopcnt=0
   endif        

!  -----------------------------------------------------------
!  IF MODEL AGREEMENT IS VERY GOOD, ATTEMPT A JUMP
!  -----------------------------------------------------------

!
!  A JUMP IS DONE BY MINIMIZING THE SAME MODEL OVER A LARGER
!  TRUST REGION. IT IS ONLY DONE IF THE MINIMIZATION WITH RADIUS
!  DELTA LANDED ON THE TRUST REGION BOUND
!

   if ( rho >= rhojmp  .and. rho <= two - rhojmp &
            .and. abs(snorm-delta)<cnstol ) then
     if ( iprint >= 2 ) write( iout, 8085 )
     theta = zero
     do 90 i = 1, nq - 1
!
!  FOR ALL SAMPLE POINTS (EXCEPT LAST) WHICH ARE NOT IN INTERPOLATION
!  AND ARE FURTHER THAT DELTA AWAY  FROM THE BASE COMPUTE THE AGREEMENT
!  BETWEEN THE MODEL AND FUNCTION
!
       if ( iwrk( lsp2in+ i ) == 0   .and. &
            wrk(  ld + i ) >= delta ) then
         ibeg = lp + (i-1) * n 
!
!  COMPUTE THE PREDICTED VALUE OF THE MERIT FUNCTION
!
         mval= mvalue( n, wrk( ibeg+1 ), wrk( ig ), &
                       wrk( ih ), wrk( icurw), lrs )
         if ( method <= 2  ) then
           call fmerit(m, val, mval+kappa, wrk(ici+(i-1)*m), &
                       ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1), &
                       pp, method )
         else
           val = mval + kappa
         endif
         val=val-fbase
         if ( abs(val)>mcheps ) then
           rhow = ( wrk( lv + i ) - fbase ) / val
!
!  IF FOR A GIVEN POINT THE AGREEMENT IF GOOD, THEN WE TRY TO JUMP 
!  AT LEAST AS FAR AS THIS POINT FROM THE BASE
!
           if ( rhow >= one - rhomin .and. &
                rhow <= one + rhomin      ) then
             theta = max( theta, wrk( ld + i ) )
           endif 
         endif                  
       endif
90      continue

!
!  MAKE SURE THAT THE SIZE OF THE JUMP DOES NOT EXCEED THE LIMIT
!          
     theta = min( theta, maxjmp * delta )
     if ( iprint >= 2 ) write( iout, 8090 ) theta
       
     del = theta

!
!  IF THE POSSIBLE SIZE OF THE JUMP BIG ENOUGH (SO IT COVERS AREA LARGER
!  THAN THE AREA THAT WOULD BE COVERED BY THE NEXT ITERATION) THEN
!  COMPUTE THE JUMP
!
100      if ( del > snorm + ratio * delta) then

!
!  SET INITIAL POINT TO THE BASE POINT AND SHIFT IT TO ORIGINAL POSITION
!
       call dcopy(n, wrk(ibase), 1, wrk( ik), 1)        
       call unshft(n, x, wrk(ik) )
!
!  CALL MINIMIZATION
!
       if ( method <= 2 ) then 
         call mintr( n           , wrk(ik), mval      , del, &
                     lb          , ub     , a         , lda, &
                     nclin       , m      , wrk(icurw), lrs, &
                     iwrk(icuriw), lis    , inform    , 1  )   
       else 
         call mintr( n           , wrk(ik), mval      , del, &
                     lb          , ub     , a         , lda, &
                     nclin       , ncnln  , wrk(icurw), lrs, &
                     iwrk(icuriw), lis    , inform    , method ) 
       endif
       call shift(n, x, wrk(ik) )
       if ( inform == 3 ) inform=0
       if ( inform == 2 ) goto 115
!
!  COMPUTE THE PREDICTED VALUE OF THE MERIT FUNTION
!
       if ( method <= 2 ) then
         call fmerit(m, val, mval+kmod, wrk(icurw), &
          ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1), pp, &
          method )
       else
         val=mval+kmod
       endif
!
!  IF MINIMIZATION HAD FAILED, THEN RECORD THE BEST POINT AND EXIT
!
       if (inform/=0) then
         call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
                 wrk(ic), nq)
         it   = iter
         exit = -7
         return
       endif

!
!  COMPUTE THE DISTANCE FROM THE NEW POINT TO THE BASE
!
       do 110 i=1,n
         wrk(lcurw+i)=wrk(lk+i)-wrk(lbase+i)
110        continue
       stnorm  = dnrmnf( n, wrk( icurw ) )
       if ( iprint >= 2 ) write( iout, 8092 ) stnorm
!
!  IF THE ACTUAL JUMP IS LARGE ENOUGH THEN SAMPLE THE FUNCTION VALUE
!  AT THE NEW POINT
!
       if ( stnorm >  snorm + ratio * delta ) then
         call dcopy(n, wrk(ik), 1, wrk( inq + n ), 1)
         nf = nf + 1
!
!  IF THE NUMBER OF FUNCTION CALLS EXCEEDS MAXIMUM, THEN RECORD THE
!  BEST POINT AND EXIT
!
         if ( nf >= maxnf ) then
           call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
                       wrk(ic), nq)
           it   = iter
           exit = 1
           return
         endif

!
!  COMPUTE THE FUNCTION VALUE 
!
         call unshft(n, x, wrk(ik))
         if ( scale /= 0 ) call unscl( n, wrk(ik), scal )
         call fun( n, m, wrk( ik ), wrk(iob+nq), wrk(ic + nq*m), &
                   iferr )
         if ( scale /= 0 ) call scl( n, wrk(ik), scal )
!
!  IF THE FUNCTION COMPUTATION FAILS THEN REDUCE THE SIZE OF THE JUMP
!  AND TRY TO COMPUTE A NEW JUMP
!
         if ( iferr ) then 
           del=stnorm/2
           if ( del >= delta ) then 
             goto 100
           else
             goto 115
           endif
         endif

!
!  IF THE FUNCTION VALUE IS COMPUTED THEN RECORD IT AND DO APPROPRIATE
!  BOOKKEEPING AND UPDATES
!
         if ( iprint >= 2 ) write( iout, 8095 )
         call fmerit(m, ft, wrk(iob+nq), wrk(ic + nq*m), &
              ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1), &
              pp, method )
         nq     = nq + 1
         lnq    = lp + ( nq - 1 ) * n
         inq    = lnq + 1
         lpnts  = lpnts + n
         lvalue = lvalue + 1
         lconvl = lconvl + m
         wrk( lv + nq )     = ft
         iwrk( lsp2in+ nq ) = 0
         wrk ( ld + nq )    = stnorm

!
!  IF THE NEW POINT IS BETTER THAN HE PREVIOUS ONE
!  THE TAKE THE NEW POINTS AND NEXT TO BE INCLUDED IN
!  HE INTERPOLATION SET
!
            
         if ( ft <  fnq ) then
         
!
!  COMPUTE 'RHO' - ACHIEVED REDUCTION VS. PREDICTED REDUCTION
!
           rho   = ( ft - fbase ) / (val - fbase)
           next  = nq
           snorm = stnorm
           impr  = 20
           if ( iprint >= 2 ) write( iout, 8097 ) rho
         endif
       endif
     endif
   endif
!  ------------------------------------------------------------------
!         END OF JUMPING ATTEMPT
!  ------------------------------------------------------------------
  
!  ------------------------------------------------------------------
!
!  UPDATE THE INTERPOLATION SET
!
!  ------------------------------------------------------------------

!
!  SET FLAGS TO INDICATE THAT NOTHING HAS BEEN DONE TO THE MODEL YET
!     
115    ifadd = 0
   ifxg  = 0
   inext =(next-1)*n

!
!  IF THE IMPROVEMENT IN FUNCTION VALUE IS SUFFICIENT WE INCLUDE
!  THE NEXT POINT IN THE MODEL
!

   if ( rho > rhomin ) then
     if ( iprint >= 2 ) write( iout, 8100 ) 
!
!  IF THE MODEL IS INCOMPLETE, THEN WE ADD THE NEXT POINT TO THE
!  MODEL ( IF EFFORT >= 4, THEN IT WILL BE DONE LATER )
!
     if ( nind < dd .and. effort < 4 ) then
!
!  DETERMINE THE INDEX OF LAST POLYNOMIAL IN THE CURRENT BLOCK
!
       if ( nind <= n ) then
         ipend = n + 1
       else
         ipend = dd
       endif
!
!  TRY TO FIND A POLYNOMIAL WHICH PRODUCES ACCEPTABLE PIVOT 
!  TO INCLUDE THE NEW POINT IN INTERPOLATION.
!  WE MAY NEED TO CHECK ALL 'AVAILABLE' POLYNOMIALS IN CURRENT BLOCK
! 
       do 120 ipoly = nind+1, ipend
!
!  UPDATE THE NEXT NEWTON POLYNOMIAL, SO THAT IT IS 'ORTHOGONAL'
!  TO THE CURRENT INTERPOLATION SET
!

         call nextnp( ipoly, wrk(inp), wrk(ipi), nind, n, &
                      lpoly, lptint )
!
!  EVALUATE NEXT POLYNOMIAL AT THE NEW POINT
!

         call evalx( val, wrk(ip+inext), wrk(inp), ipoly, &
                     n  , lpoly )
!
!  IF THE VALUE (THE PIVOT) IS ACCEPTABLE, ACCEPT CURRENT POLYNOMIAL
!  AS THE NEXT AND DO NOT CHECK OTHER POLYNOMIALS.
!
          
         if (abs(val) > pivthr ) go to 130

120        continue
       go to 140
!
!  PLACE THE IPOLY-TH POLYNOMIAL IN THE PLACE OF NIND+1-ST
!  POLYNOMIAL, BY SWAPPING THEM IN ARRAY POLY
!
130        if ( ipoly /= nind+1 ) then
         call swapnp( n, nind+1, ipoly, wrk(inp), lpoly )
       endif
!
!  PERFORM THE APPROPRIATE UPDATES TO THE SET OF THE FIRST NIND
!  NEWTON POLYNOMIALS, WHICH ARE REQUIRED TO INCLUDE THE NEW
!  POINT-POLYNOMIAL PAIR IN THE INTERPOLATION
!
       call ptrepl( wrk(ip+inext), nind+1, val, wrk(ipi), &
                    wrk(inp)     , nind  , n  , lpoly   , lptint )
!
!  RECORD THE NEW POINT IN THE INTERPOLATION SET
!
       call dcopy(n, wrk(ip+inext), 1, wrk(ipi+nind*n), 1)
       nind              = nind+1
       wrk(lvi+nind)     = wrk(lv+next)
       iwrk(lsp2in+next) = 1
       iwrk(lin2sp+nind) = next
       wrk(lpiv+nind)    = val
       ifadd             = 1
!
!  IF EFFORT = 4 THEN WE CHANGE THE BASE TO THE NEW BEST POINT
!  ( OTHERWISE THE BASE WILL BE CHANGED LATER )
!
       if ( effort == 4 ) then
         base              = nind 
         lbase = lpi + ( base - 1 ) * n
         ibase = lbase + 1
         obase = base
         call getdis( n, nq, 0 , wrk( ip ) , iwrk(lin2sp+base), &
                      wrk( id ), wrk(icurw), lrs )
       endif
       if ( iprint >= 2 ) write( iout, 8110 ) val
     endif


140      oldnq = nq          

!
!  IF THE MODEL IS FULL OR THE PIVOT VALUE FOR ADDING THE POINT
!  IS TOO SMALL, THEN TRY TO INCLUDE THE NEW POINT BY REPLACING
!  ANOTHER POINT BY IT. IF SUCCEED, UPDATE THE NEWTON POLYNOMIALS
!  ( IF EFFORT >= 4, THEN IT WILL BE DONE LATER )
!

     if ( ifadd==0 ) then
!
!  IF THE PIVOT IS TOO SMALL AND THE MODEL IS NOT FULLY LINEAR
!  SET A FLAG TO MAKE SURE THAT WE TRY TO ADD A 'GEOMETRY' POINT
!  LATER IN 'IMPMOD' PROCEDURE
!
       if (.not. linear) ifxg=1
       ipoly=base
       call ptexch( wrk(inp), wrk(ip+inext), wrk(ipi), ipoly , &
                    nind    , n            , lpoly   , lptint, &
                    pivthr  , wrk(ipiv)    , val     , fail  )

!
!  IF WE FIND A POINT TO BE REPLACED BY THE NEW POINT, THEN
!  RECORD THE EXCHANGE IN INTERPOLATION SET
!
       if (.not.fail) then
         call dcopy(n, wrk(ip+inext), 1 , &
                    wrk(ipi+(ipoly-1)*n), 1)
         wrk(lvi+ipoly)     = wrk(lv+next)
         wrk(lpiv+ipoly)    = val*wrk(lpiv+ipoly)
         iwrk(lin2sp+ipoly) = next
         iwrk(lsp2in+next)  = 1
         iwrk(lsp2in+iwrk(lin2sp+ipoly))=0
         ifadd              = 1
!
!  IF EFFORT = 4 THEN WE CHANGE THE BASE TO THE NEW BEST POINT
!  ( OTHERWISE THE BASE WILL BE CHANGED LATER )
!
         if ( effort == 4 ) then
           base              = ipoly 
           lbase = lpi + ( base - 1 ) * n
           ibase = lbase + 1
           obase = base
           call getdis( n, nq, 0 , wrk( ip ) , iwrk(lin2sp+base), &
                        wrk( id ), wrk(icurw), lrs )
         endif
         if ( iprint >= 2 ) write( iout, 8120 ) ipoly,  val
       endif
     endif
!
!  IF THE MODEL REDUCTION IS NOT SUFFICIENTLY GOOD, WE STILL
!  TRY TO ADD THE NEW POINT (SINCE IT CONTAINS NEW INFORMATION)
!

   else
     if ( iprint >= 2 ) write( iout, 8105 )
!
!  IF THE MODEL IS INCOMPLETE TRY ADDING THE NEW POINT
!  IF THE EFFORT >= 3 TRY TO DO IT LATER
!

     if ( nind<dd .and. effort < 3 ) then
!
!  DETERMINE THE INDEX OF LAST POLYNOMIAL IN THE CURRENT BLOCK
!
       if ( nind <= n ) then
         ipend = n + 1
       else
         ipend = dd
       endif
!
!  TRY TO FIND A POLYNOMIAL WHICH PRODUCES ACCEPTABLE PIVOT 
!  TO INCLUDE THE NEW POINT IN INTERPOLATION.
!  WE MAY NEED TO CHECK ALL 'AVAILABLE' POLYNOMIALS IN CURRENT BLOCK
! 
       do 150 ipoly = nind+1, ipend
!
!  UPDATE THE NEXT NEWTON POLYNOMIAL, SO THAT IT IS 'ORTHOGONAL'
!  TO THE CURRENT INTERPOLATION SET
!

         call nextnp( ipoly, wrk(inp), wrk(ipi), nind, n, &
                      lpoly, lptint )
!
!  EVALUATE NEXT POLYNOMIAL AT THE NEW POINT
!

         call evalx( val, wrk(ip+inext), wrk(inp), ipoly, &
                     n  , lpoly )
!
!  IF THE VALUE (THE PIVOT) IS ABOVE A CERTAIN THRESHOLD
!  (WHICH MAY BE SET HIGHER THAN PIVOT THRESHOLD, TO IMPOSE
!  STRICTER GEOMETRY REQUIREMENTS ON  A POINT WHICH GIVES POOR
!  REDUCTION)  THEN ADD THE POINT TO THE INTERPOLATION SET,
!  UPDATING OTHER NEWTON POLYNOMIALS
!                
         if (abs(val) > addthr ) go to 160

150        continue
!
!  IF THERE IS NO POLYNOMIAL WHICH GIVES ACCEPTABLE PIVOT, THEN
!  MOVE ON TO MODEL IMPROVEMENT

       go to 170
!
!  PLACE THE IPOLY-TH POLYNOMIAL IN THE PLACE OF NIND+1-ST
!  POLYNOMIAL, BY SWAPPING THEM IN ARRAY POLY
!
160        if (ipoly /= nind+1) then
         call swapnp(n, nind+1, ipoly, wrk(inp), lpoly)
       endif


!
!  ADD THE POINT TO THE INTERPOLATION SET,
!  UPDATING OTHER NEWTON POLYNOMIALS
!


       call ptrepl( wrk(ip+inext), nind+1, val, wrk(ipi), &
                    wrk(inp)     , nind  , n  , lpoly   , lptint)

!
!  DO APPROPRIATE UPDATES
!
       call dcopy(n, wrk(ip+inext), 1, wrk(ipi+nind*n), 1)
       nind              = nind+1
       wrk(lvi+nind)     = wrk(lv+next)
       iwrk(lsp2in+next) = 1
       iwrk(lin2sp+nind) = next
       wrk(lpiv+nind)    = val
       ifadd             = 1
       if ( iprint >= 2 ) write( iout, 8110 ) val
     endif
 
!
!  END OF ATTEMPT OF INCLUDING A NEW POINT
!
   endif
 endif
170  continue

 if ( iprint >= 3 ) then
   write( iout, 8170 )
   do 175 i = 1, nind
     write( iout, 8055 ) wrk(ipiv + i - 1)
175    continue
 endif 
!
!  IF THE IMPROVEMENT OF THE MODEL IS REQUIRED, SET IMPR TO 0
! 
 if ( (ifadd == 0) .or. (ifxg == 1) ) impr = 0
!
!  IF THE NEW POINT HASN'T BEEN ADDED, BUT IT GIVES REDUCTION, THEN
!  SET A VALUE OF IMPR, SO THAT THIS POINT IS SET TO BASE IN 'IMPMOD' 
! 
 if ( ifadd==0 .and. rho> rhomin ) impr=-next
!
!  TRY TO IMPROVE THE MODEL (BY FINDING ANOTHER POINT OR  RECOMPUTING
!  THE BASIS OF NEWTON POLYNOMIALS AND THE INTERPOLATION SET IF:
!      1. NOTHING WAS ADDED TO THE POINT
!      2. IF A NEW 'GEOMETRY' POINT IS REQUIRED 
!      3. THE MODEL IS TOO OLD, THAT IS THE CURRENT BASE IS LAYER*DELTA
!         AWAY FROM THE FIRST INTERPOLATION POINT
!      4. IF THE EFFORT LEVEL REQUIRES INTERPOLATION TO BE RECOMPUTED      
!
 if (  impr <= 0                            .or. &
      (wrk(ld+iwrk(iin2sp)) > layer*delta) .or. &
      (effort == 3 .and. rho <= rhomin)   .or. &
       effort == 4 ) then

!
!  CHECK IF THE NUMBER IF FUNCTION CALLS REACHED MAXIMUM
!  IF IT DOES RECORD THE BEST POINT AND EXIT
!
   if ( nf >= maxnf ) then
     call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
                 wrk(ic), nq)
     it   = iter
     exit = 1
     return
   endif
!
!  COMPLETE THE INTERPOLATION WITH DUMMY ELEMENTS IF THE LINEAR
!  BLOCK REACHES THE MAXIMUM POSSIBLE SIZE (MORE DETAILS ABOVE)
!

   if ( neqcon > 0 .and. nind == n + 1 - neqcon ) then
     call complt( wrk(inp), wrk(ivi), wrk(ipiv), n, nind, lpoly )
   endif              
   if ( iprint >= 2 ) write( iout, 8140 )
   call impmod( wrk(inp)    , wrk(ipi)    , wrk(ivi)    , wrk(ip), &    ! 
                wrk(iv)     , wrk(iob)    , wrk(ic)     , &
                iwrk(isp2in), iwrk(iin2sp), n, m   , &
                nq          , nind        , base        , pivthr , &    ! 
                xchthr      , wrk(ipiv)   , wrk(id)     , delta  , &
                x           , a           , lda         , nclin  , &    ! 
                ncnln       , lb          , ub          , scale  , &
                scal        , nf          , maxnf       , impr   , &
                pp          , neqcon      , &
                wrk(icurw+1), lrs         , iwrk(icuriw), lis    , &
                delmin      , method      )

!
!  CHECK IF THE NUMBER IF FUNCTION CALLS EXCEEDS MAXIMUM
!  IF IT DOES RECORD THE BEST POINT AND EXIT
!
   if ( nf >= maxnf ) then
     call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
                 wrk(ic), nq)
     it   = iter
     it   = iter
     exit = 1
     return
   endif


!
!  IF THE BASIS CONSISTS ONLY OF A SINGLE POLYNOMIAL (CONSTANT TERM)
!  THEN SOMETHING IS WRONG. WE RETURN AND PRINT A MESSAGE      
!

   if (nind<2) then
     if ( iprint >= 0 ) write(iout, 2060)
     call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
                 wrk(ic), nq)
     it    = iter
     exit  = -9
     return
   endif
!
!  IF THE NUMBER OF SAMPLE POINTS INCREASED AFTER APPLYING IMPMOD
!  WE UPDATE THE POINTERS TO THE LAST POINT ACCORDINGLY
!

   if ( oldnq /= nq ) then
     lnq   = lp  + (nq-1)*n
     inq   = lnq + 1
     oldnq = nq
   endif
 endif

!
!  IF THE MODEL WAS CHANGED EITHER BY ADDING A NEW "GOOD" POINT
!  OR BY THROWING AWAY OLD POINTS OR IF THE MODEL WAS NOT LINEAR
!  WE KEEP 'OKGEOM' EQUAL TO FALSE TO INDICATE THAT REDUCING TRUST 
!  REGION IS NOT NECESSARY, SINCE THE MODEL HAS IMPROVED. OTHERWISE
!  SET 'OKGEOM' TO TRUE, SINCE THE MODEL WAS OK AND TR REDUCTION 
!  MAY BE NEEDED. IMPR=0 MEANS THAT THE MODEL HAS NOT IMPROVED EVEN THOUGH
!  THE MODEL WAS NOT FULLY LINEAR. THIS CAN HAPPEN IN "DEGENERATE"
!  CONSTRAINED  CASES
!
 if ((( impr >= 4  .or. rho == - thous) .and. linear) &
       .or. impr == 0 ) okgeom=.true.


!  --------------------------------------------------------------
!
!  UPDATING TRUST REGION RADIUS
!
!  --------------------------------------------------------------


!
!  IF THE REDUCTION IF GOOD, INCREASE DELTA ACCORDING TO 
!  THE STEP TAKEN BY TRUST REGION  MINIMIZATION  
!
 if ( rho >= rhoinc .or. impr == 20 ) then
   if ( iprint >= 2 ) write( iout, 8130 )
   delta = max( delta, min( delmax, ratio*snorm ))

!
!  IF THE REDUCTION WAS BAD, AND THE MODEL WAS ALREADY GOOD ENOUGH
!  THEN REDUCE DELTA
!
 elseif ( rho < rhomin .and. okgeom ) then
   if ( iprint >= 2 ) write( iout, 8150 )
   delthr =  ratio  * delmin
   delta = max( delta / ratio,  delthr )
 endif

!
!  IF DELTA IS JUST ABOVE THE MINIMUM THRESHOLD, MAKE SURE THAT IT
!  DOES NOT GET DECREASED UNTIL WE ARE SURE THAT WE CANNOT GET
!  FURTHER IMPROVEMENT
!

 if ( delta <= ratio * delmin + 10*mcheps .and. &
      rho  < rhomin  ) then
   if ( iprint >= 2 ) write( iout, 8160 )
   nstop = nstop + 1
   if (  nstop >= 5 .and. okgeom ) then
     delta = half * delmin
   endif
 else
   nstop = 0
 endif

200 continue

call bestpt( n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
             wrk(ic), nq )
it   = iter
exit = 2
return


2000 format(/, ' DFO: *** ERROR: LWRK IS TOO SMALL!' / &
         ' DFO:        ADDITIONAL SPACE OF ',i8 ,' REQUIRED')
2020 format(/, ' DFO: *** ERROR: LIWRK IS TOO SMALL!' / &
         ' DFO:        ADDITIONAL SPACE OF ',i8 ,' REQUIRED')
2060 format('DFO: CANNOT FIND TWO POINTS TO BUILD A MODEL!', / &
        '    THE PIVOT THRESHOLD MAY BE TOO HIGH')
1000 format('DFO: ', 1x,' it ', 1x, ' nf ', 1x, '    value     ', 1x, &
        '     rho      ', 1x, '     delta   ', 1x,'impr', &
        1x, 'nind', 1x, 'neqc'  )
1010 format( 'DFO:', 1x, i4, 1x, i4, 1x, d14.7, 1x, d14.7, 1x, d14.7, &
         1x, i4, 1x, i4, 1x, i4  )
1020 format('DFO: INEXACT NEWTON POLYNOMIAL', 1x, i4, 1x,'AT POINT:', &
              i4,1x, 'ERROR IS:', d13.6)
1030 format('DFO: INEXACT MODEL AT POINT:', 1x, i4,  1x,'ERROR IS:', &
             d13.6)
8000 format( /, ' DFO: computing initial interpolation set')
8010 format( /, ' DFO: ************starting iteration ', i4,'******')
8020 format( /, ' DFO: base point = ',i4,' out of ', i4,/ )
8030 format( /, ' DFO: convergence test: DELTA = ', d14.7 )
8040 format( /, ' DFO: building the new model' )
8050 format( /, ' DFO: **********the model coefficients:********* ' )
8051 format( /, ' quadratic terms:' )
8052 format( / )
8053 format( d14.7 )
8054 format( /, ' linear terms:',/ )
8055 format( d14.7, 1x )
8056 format( /, ' constant term: ', /, d14.7,/ )
8060 format( /, ' DFO: **********end of model coefficients:********')
8065 format( /, ' DFO: trust region minimization' )
8067 format( /, ' DFO: PRERED=', d14.7,' is too small',/ &
           ' DFO:       a new points is not produced' )
8070 format( /, ' DFO: sufficient reduction: take the step' )
8080 format( /, ' DFO: ared / prered = ', d14.7 )
8085 format( /, ' DFO: attempt a jump' )
8090 format( /, ' DFO: attempt jump distance = ', d14.7 )
8092 format( /, ' DFO: actual jump distance = ', d14.7 )
8095 format( /, ' DFO: a new point produced by jump' )
8097 format( /, ' DFO: jump accepted, new prered/ared=', d14.7 )
8100 format( /, ' DFO: the new point gives good reduction' )
8105 format( /, ' DFO: the new point does not give good reduction' )
8110 format( /, ' DFO: the new points is added to the model,',/ &
           ' DFO:       pivot=', d14.7 )
8120 format( /, ' DFO: the new point replaced ',i4,'-th point,',/ &
           ' DFO:       pivot=', d14.7  )
8140 format( /, ' DFO: improve the geometry' )
8130 format( /, ' DFO  : increase the trust region radius' )
8150 format( /, ' DFO  : decrease the trust region radius' )
8160 format( /, ' DFO  : trust region radius is at the lower bound',/ &
           ' DFO:       stopping count = ', i2 )
8170 format( /, ' DFO: relative pivots of the interpolation set:',/)
end








!*************************************************************

!    NEXT SUBROUTINE

!*************************************************************


subroutine bestpt( n, m, x, fx, c, points, values, objval, &
                   conval, nq )
!  ----------------------------------------------------------
!
!  THIS SUBROUTINE SELECTS THE POINT WITH THE BEST VALUE FROM
!  ARRAY 'POINTS' WITH VALUE 'VALUES'. THE POINT IS RECORDED 
!  IN X AND THE VALUE - IN FX.
!
!  ----------------------------------------------------------


!
!  SUBROUTINE PARAMETERS
!

integer          n, nq, m
double precision x(n), fx, points(nq*n), values(nq), &
                 conval(nq*m), c(m), objval(nq)

!
!  COMMON VARIABLES
!
integer          iout  , iprint
double precision mcheps, cnstol
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /

!
!  LOCAL VARIABLES
!  
integer          i, imin
double precision val

!
!  SUBROUTINE CALLED
!
!  BLAS:      DCOPY
!


val   = values(nq)
call unshft( n, x, points((nq-1)*n + 1) )
imin = nq
do 10 i = nq-1, 1, -1
  call unshft( n, x, points((i-1)*n + 1) )
  if ( values(i) < val - mcheps ) then
    val  = values(i)
    imin = i
  endif
10 continue 
fx   = objval(imin) 
call dcopy(  n, points((imin-1)*n+1), 1, x, 1 )
call dcopy(  m, conval((imin-1)*m+1), 1, c, 1)

return 
end


!*************************************************************

!    NEXT SUBROUTINE

!*************************************************************


subroutine scl(n, x, scal)
!  **********************************************************
!  THIS SUBROUTINE SCALES (BY DOING ENTREE-WISE MULTIPLICATION) 
!  VECTOR X BY VECTOR 1/SCAL
!  **********************************************************
integer          n
double precision x(n), scal(n)
integer          i

do 10 i=1,n
  x(i)=x(i)/scal(i)
10 continue

return 
end



!*************************************************************

!    NEXT SUBROUTINE

!*************************************************************


subroutine unscl(n, x, scal)
!  **********************************************************
!  THIS SUBROUTINE SCALES (BY DOING ENTREE-WISE MULTIPLICATION) 
!  VECTOR X BY VECTOR SCAL
!  **********************************************************
integer          n
double precision x(n), scal(n)
integer          i

do 10 i=1,n
  x(i)=x(i)*scal(i)
10 continue

return  
end



!*************************************************************

!    NEXT SUBROUTINE

!*************************************************************


subroutine complt( poly, valint, pivval, n, nind, lpoly)

!  *********************************************************** 
!  THIS SUBROUTINE COMPLETES THE LINEAR BLOCK OF THE INTERPOLATION
!  WITH DUMMY ELEMENTS, SO THAT THEY ARE NEVER CONSIDERED AND
!  NEVER HAVE EFFECT ON INTERPOLATION, AND SO THAT THE LINEAR
!  BLOCK HAS SIZE N + 1, I.E., IS COMPLETE.
!
!  PARAMETERS
!  
!    POLY   (INPUT/OUTPUT)  NEWTON FUNDAMENTAL POLYNOMAILS
!
!    VALINT (INPUT/OUTPUT)  VALUES OF INTERPOLATION POINTS
!
!    PIVVAL (INPUT/OUTPUT)  PIVOT VALUES OF INTERPOLATION POINTS
!
!     N     (INPUT)         PROBLEM DIMENSION
!
!    NIND   (INPUT/OUTPUT)  NUMBER OF INTERPOLATION POINTS/POLYNOMIALS
!  ***********************************************************


!
!  SUBROUTINE PARAMETERS
!

integer          n, nind, lpoly
double precision poly(lpoly), valint(n+1), pivval(n+1)

!
!  LOCAL VARIABLES
!
integer          i, ibeg, iend
double precision zero, huge
parameter       (zero = 0.0d0, huge = 1.0d20)

!
!  SET ALL POLYNOMIALS FROM NIND+1-ST TILL N+1-ST TO ZERO
!

ibeg = 2 + (n+1)*(nind-1) 
iend = 1 + (n+1)*n
do 10 i=ibeg, iend
  poly(i)=zero
10 continue 
!
!  SET THE CORRESPONDING FUNCTION VALUES TO HUGE (SO THEY WILL NOT
!  BE PICKED UP AS THE  SMALLEST CURRENT VALUE), AND PIVOT VALUE TO
!  HUGE, SO THAT THEY WILL NOT BE REPLACED LATER
! 
do 20 i = nind+1, n+1
  valint(i)=huge
  pivval(i)=huge
20 continue 
!
!  SET NING TO THE SIZE OF FULLY LINEAR MODEL
! 
nind = n+1
return
end




!*************************************************************

!    NEXT SUBROUTINE

!*************************************************************

subroutine intdim( x    , n , nclin , ncnln , neqcon, a   , lda , &     ! 
                   lb   , ub, pivthr, wrk   , lwrk  , iwrk, liwrk)

!  ******************************************************************
!  THIS SUBROUTINE IDENTIFIES THE NUMBER OF LINEARLY INDEPENDENT 
!  EQUALITY CONSTRAINTS OF THE FEASIBLE SET. A CONSTRAINT IS CONSIDERED
!  EQUALITY IF THE UPPER AND LOWER BOUND DIFFER BY LESS THAT PIVTHR.
!  TO FIND THIS NUMBER WE CONSTRUCT THE JACOBIAN MATRIX OF EQUALITY
!  CONSTRAINTS AND THEN COMPUTE ITS QR DECOMPOSITION AND COUND THE
!  NUMBER OF ZEROS IN THE DIAGONAL OF R.
!  THIS NUMBER OF LIN. IND. EQUAL. IS COMPUTED AT A GIVEN POINT, SINCE
!  IN THE PRESENCE OF NONLINEAR CONSTRAINTS IT CAN BE DIFFERENT 
!  (IN THEORY)
!
!   PARAMETERS
!
!    X      (INPUT)  CURRENT POINT
!
!    N      (INPUT)  DIMENSION OF THE PROBLEM
!
!    NCLIN  (INPUT)  NUMBER OF LINEAR CONSTRAINTS
!    
!    NCNLN  (INPUT)  NUMBER OF NONLINEAR CONSTRAINTS
!
!    NEQCON (OUTPUT) THE NUMBER OF LINEAR INDEP. EQUALITY CONSTR.
!
!    A      (INPUT)  MATRIX OF LINEAR CONSTRAINTS
!
!    LDA    (INPUT)  LEADING DIMENSION OF MATRIX A 
!                   
!    LB     (INPUT)  ARRAY OF LOWER BOUNDS
!
!    UB     (INPUT)  ARRAY OF UPPER BOUNDS
!
!    PIVTHR (INPUT)  PIVOT THRESHOLD
!
!    WRK    (INPUT)  WORKING REAL SPACE ARRAY
!
!    IWRK   (INPUT)  WORKING INTEGER SPACE ARRAY
!  *******************************************************************

!
!  SUBROUTINE PARAMETERS
!
double precision a(lda*n) , lb(n + nclin + ncnln), &
                 pivthr   , ub(n + nclin + ncnln), &
                 wrk(lwrk), x(n) 

integer          n, nclin, ncnln, neqcon, lda,  lwrk, liwrk, &
                 iwrk(liwrk)


!
!  COMMON VARIABLES
!
include 'dfo_model_inc.inc'

integer          iout  , iprint
double precision mcheps, cnstol
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /
!
!  LOCAL VARIABLES
!

integer          info, i, j, m, ic, icjac, lcjac, ldmat, imat, &
                 lmat, icurw, lrs, ldcj


double precision zero, one
parameter      ( zero = 0.0d0, one =1.0d0 ) 

!
!     SUBROUTINES AND FUNCTIONS CALLED:
!
!       APPLICATION:       RZRVEC, IZRVEC, FUNCON, 
!       LAPACK     :       DGEQPF
!       FORTRAN SUPPLIED:  ABS, MAX
!       BLAS:              DCOPY 
!

ldcj  = max(ncnln, 1)

!
!  PARTITION THE MEMEORY
!

!
!  POINTER TO ARRAY OF NONLINEAR CONSTRAINTS VALUES
!
ic    = 1
!
!  POINTER TO ARRAY WITH THE JACOBIAN OF NONLINEAR CONSTRAINTS
!
icjac = ic    + ldcj
!
!  POINTER TO MATRIX CONTAINING THE JACOBIAN  OF EQUALITY CONSTRAINTS
!
imat  = icjac + ldcj * n
!
!  POINTER TO CURRENT WORKING ARRAY
!
icurw = imat  +(nclin + ncnln + n) * n
!
!  CHECK IF MEMORY IS SUFFICIENT
!
lrs = lwrk - icurw + 1
if ( lrs < 4*n ) then
  if ( iprint >= 0 ) write( iout, 1000 ) - lrs - 4*n + 1 
  stop
endif

if ( liwrk < n ) then
  if ( iprint >= 0 ) write( iout, 1100 ) - liwrk - n + 1 
  stop
endif
!
!  SET SOME AUXILARY POINTERS
!
lmat  = imat - 1
lcjac = icjac - 1
ldmat = n + nclin + ncnln


call rzrvec( wrk(imat), ldmat*n )
!
!  IF THERE ARE NONLINEAR CONSTRAINTS, COMPUTE THEIR JACOBIAN
!
if ( ncnln > 0 ) then 
   usemerit=.true.
   call funcon(2, ncnln  , n, ldcj   , iwrk, &
               x, wrk(ic), wrk(icjac), 1     )
   usemerit=.false.
endif

!  ------------------------------------------
!  CONSTRUCT JACOBIAN OF EQUALITY CONSTRAINTS
!  ------------------------------------------

!
!  IF THE ARE ANY FIXED VARIABLES INCLUDE CORRESPONDING ROW
!  OF IDENTITY IN THE JACOBIAN
!   
m = 0 
do 20 i=1, n
  if ( ub(i)-lb(i) < pivthr ) then
    do 10 j = 1, n
      m = m + 1
      if ( j == i ) then
        wrk( lmat + (j-1)*ldmat + i )= one
      else
        wrk( lmat + (j-1)*ldmat + i )= zero
      endif
10     continue
  endif
20 continue

!
!  CHECK IF THE THERE ARE EQUALITIES AMONG LINEAR CONSTRAINTS,
!  FOR EACH EQUALITY WRITE CORRESPONDING ROW OF THE JACOBIAN
!

do 30 i=1 ,  nclin
  if ( ub(i + n) - lb(i + n) < pivthr ) then
    m = m + 1
    call dcopy(n, a( i ), lda, wrk(lmat + m), ldmat) 
  endif
30 continue

!
!  CHECK IF THE THERE ARE EQULITIES AMONG NONLINEAR CONSTRAINTS,
!  FOR EACH EQUALITY WRITE CORRESPONDING ROW OF THE JACOBIAN
!

do 40 i=1, ncnln
  if ( ub(i + n + nclin)-lb(i + n + nclin) < pivthr ) then
    m = m + 1
    call dcopy( n, wrk(lcjac + i), ldcj, wrk(lmat + m), ldmat )
  endif
40 continue
 
call izrvec( iwrk, n )


!
!  IF THE NUMBER OF EQUALITIES IS BIGGER THAN NUMBER OF VARIABLES
!  PRINT AN ERROR MESSAGE  AND QUIT (THE NUMBER OF LINEARLY INDEP.
!  CONSTRAINTS MAY STILL BE LESS THAN N, WE HAVE TO THINK HOW TO
!  TAKE CARE OF IT BETTER.
!

if ( m > n ) then
  if ( iprint >= 0 ) write( iout, 1200 ) m, n
  stop
endif
!
!  COMPUTE THE "QR" FACTORIZATION
!
call dgeqpf( m, n, wrk(imat), ldmat, iwrk, wrk(icurw), &
             wrk(icurw + n) , info )

!
!  LOOK THROUGH THE DIAGONAL ELEMENTS OF "R", FOR EACH ZERO FOUND
!  REDUCE THE NUMBER OF LIN. IND. CONSTRAINTS BY 1.
!
neqcon = m
do 50 i = m, 1, -1
  if ( abs(wrk( lmat + (i-1)*ldmat + i )) < 100*mcheps ) &
       neqcon = neqcon - 1
50 continue


return
1000 format(/ ' INTDIM: *** ERROR: LWRK IS TOO SMALL!' / &
         '             ADDITIONAL SPACE OF ',i8 ,' REQUIRED')
1100 format(/ ' INTDIM: *** ERROR: LIWRK IS TOO SMALL!' / &
         '             ADDITIONAL SPACE OF ',i8 ,' REQUIRED')
1200 format(/ ' INTDIM: *** ERROR: NUMBER OF EQUALITY CONSTRAINTS!', &
    i8,/ '             IS BIGGER THAN THE DIMENSION  ',i8 )
end  






 subroutine evalnp(val, points, j, poly, i, n, lpnt, lpoly)

!
!  *****************************************************************
!  THIS SUBROUTINE COMPUTES THE VALUE OF THE I-TH NEWTON FUNDAMENTAL
!  (STORED IN ARRAY 'POLY') AT THE J-TH POINT IN ARRAY 'POINTS'
!
!  (IF THE POLYNOMIAL IS QUADRATIC WE COMPUTE THE VALUE BY THE HORNER'S
!   RULE)
!
!  PARAMETERS
!
!  VAL    (OUTPUT) VALUE OF THE I-TH POLYNOMIAL AT THE J-TH POINT
! 
!  POINTS (INPUT)  THE ARRAY OF POINTS
!
!  J      (INPUT)  THE INDEX OF THE EVALUATION POINT
!        
!  N      (INPUT)  PROBLEM DIMENSION
!
!  POLY   (INPUT)  ARRAY OF NEWTON FUNDAMENTAL POLYNOMIALS
!          
!  I      (INPUT)  INDEX OF THE EVALUATED POLYNOMIAL
!
!  LPNT   (INPUT)  LENGTH OF THE ARRAY 'POINTS'
!  
!  LPOLY  (INPUT)  LENGTH OF THE ARRAY 'POLY'
!  *****************************************************************
!


 double precision points(lpnt), poly(lpoly), val

 integer          i, j, lpnt, lpoly, n


!
!  LOCAL VARIABLES
!

 integer          ii, ipoint, jj , k,  l,  np1
  
 double precision v

 np1    = n+1

!
!  THE POINTER TO WHERE ELEMENTS OF J-TH POINT BEGIN IN ARRAY 'POINTS'
!
 ipoint = (j-1)*n


 if (i==1) then
!
!  IF THE I-TH POLYNOMIAL IS A  CONSTANT
!
   val = poly(1)
 else if (i<=np1) then
!
!  IF THE I-TH POLYNOMIAL IS LINEAR
!
   k    = ( i - 2 )*( np1 ) + 2
   val  = poly(k)
!
!  ADD UP THE TERMS
!

   do 10 l=1,n
     val = val + poly( k + l )*points( ipoint + l )
10    continue
 else
!         
!  IF THE I-TH POLYNOMIAL IS QUADRATIC
!

   k   = ( i - np1 - 1 )*(n+2)*np1/2 + n*np1 + 2
   val = poly( k )
   jj  = k + np1
   do 30 l=1,n
     v = poly( k + l )
     do 20 ii=l, n
       v  = v + poly( jj )*points( ipoint + ii )
       jj = jj + 1
20      continue
     val=val+v*points( ipoint + l )
30    continue
 endif

 return
 end
     


!***********************************

!         NEXT SUBROUTINE          *

!***********************************


     
 subroutine evalx(val, x, poly, i, n, lpoly)


!
!  *****************************************************************
!  THIS SUBROUTINE COMPUTES THE VALUE OF THE I-TH NEWTON FUNDAMENTAL
!  (STORED IN ARRAY 'POLY') AT A GIVEN POINT X
!
!  IF THE POLYNOMIAL IS QUADRATIC WE COMPUTE THE VALUE BY THE HORNER'S
!  RULE
!
!  PARAMETERS
!
!  VAL    (OUTPUT) VALUE OF THE I-TH POLYNOMIAL AT THE J-TH POINT
! 
!  X      (INPUT)  THE EVALUATION POINT
!        
!  N      (INPUT)  PROBLEM DIMENSION
!
!  POLY   (INPUT)  ARRAY OF NEWTON FUNDAMENTAL POLYNOMIALS
!          
!  I      (INPUT)  INDEX OF THE EVALUATED POLYNOMIAL
!
!  LPNT   (INPUT)  LENGTH OF THE ARRAY 'POINTS'
!  
!  LPOLY  (INPUT)  LENGTH OF THE ARRAY 'POLY'
!  *****************************************************************
!


 double precision x(n), poly(lpoly), val

 integer          i, lpoly, n


!
!  LOCAL VARIABLES
!

 integer          ii, j , k,  l,  np1

 double precision v
 np1    = n+1

 if (i==1) then
!
!  IF THE I-TH POLYNOMIAL IS A  CONSTANT
!
   val = poly(1)
 else if (i<=np1) then
!
!  IF THE I-TH POLYNOMIAL IS LINEAR
!

   k   = ( i - 2 )*( np1 ) + 2
   val = poly(k) 
   do 10 l=1,n
     val = val + poly( k + l )*x( l )
10   continue
 else
!         
!  IF THE I-TH POLYNOMIAL IS QUADRATIC
!
   k   = ( i - np1 - 1 )*(n+2)*np1/2 + n*np1 + 2
   j   = k + np1
   val = poly( k )
   do 30 l=1,n
     v = poly( k + l )  
     do 20 ii=l,n
       v = v + poly( j )*x( ii )
       j = j + 1
20      continue
     val=val+v*x( l )
30    continue
 endif
 return
 end
     

subroutine fd( lafla , pntint, values, poly, n,  q, neqcon, &
               lptint, lpoly )


!    
!  ***********************************************************************
!  THIS SUBROUTINE COMPUTED FINITE DIFFERENCES FOR MULTIVARIATE
!  INTERPOLATION, USING 'Q' INTERPOLATION POINTS AND 'Q' NEWTON
!  FUNDAMENTAL POLYNOMIALS.
!
!  PARAMETERS:
!  
!   N      (INPUT)  DIMENTION OF THE PROBLEM
!
!   Q      (INPUT)  NUMBER OF INTERPOLATION POINTS
! 
!   PNTINT (INPUT)  LIST OF  'NIND' DATA POINTS. THE I-TH POINT OCCUPIES
!                   POSITIONS ( I - 1 ) * N + 1 TO I * N.
!   VALUES (INPUT)  VALUES OF THE OBJECTIVE FUNCTION AT THE 'NIND'
!                   DATA POINTS CONTAINED IN POINTS.
!   POLY   (INPUT)  THE ARRAY CONTAINING COEFFICIENTS OF NEWTON FUNDAMENTAL
!                   POLYNOMIALS (AS COMPUTED BY 'NBUILD')
!   NEQCON (INPUT)  NUMBER OF LINEARLY INDEP. EQUALITY CONSTRAINTS
!
!   LAFLA  (OUTPUT) THE ARRAY OF THE  FINITE DIFFERENCES
! 
!  **************************************************************************
!

 integer          n, q, neqcon, lptint, lpoly 

 double precision lafla(q), pntint(lptint), values(q), poly(lpoly)

!
!  LOCAL VARIABLES
! 
 integer          dd, i, k, l,  np1, np2, qp1
 double precision val

 np1 = n + 1
 np2 = n + 2
 dd = (np1)*(np2)/2.0
 qp1 = q + 1

!
!  INITIALISE THE LAFLA, THE FINITE DIFFERENCE OPERATORS ON F
!

 do 10 i =1,q
   lafla( i ) = values(i)
10  continue
!
!  START THE MAJOR ITERATION
!

 do 100 i=1,q
!
!  FOR ALL MULTIINDICES ALFA (WHOSE CARDINALITY IS GREATER THAN OR EQUAL
!  TO I), UPDATE THE ASSOCIATED FUNCTION L_ALFA USING
!  THE NEWTON FUNDAMENTAL POLYNOMIALS OF DEGREE EQUAL TO THE 
!  CARDINALITY OF ALFA_{I-1}
!

   if (i==2) then
!          
!  ALL FINITE DIFFERENCES CORRESPONDING TO DEGREES >= 1 ARE UPDATED.
!  THE NEWTON FUNDAMENTAL POLYNOMIALS OF DEGREE CARDINALITY OF ALFA_I = 0
!  ARE USED IN THE UPDATE. THE UPDATES ARE SKIPPED FOR DUMMY POLYNOMIALS
!  WHOSE INDEX IS ALWAYS HIGHER THAN NP1-NEQCON
!
     do 110 l = 2, min(q, np1 - neqcon) 
       call  evalnp( val, pntint, l, poly, 1, n, lptint, &
                     lpoly )
       lafla( l ) = lafla( l ) - lafla( 1 )*val
110      continue
     do 115 l = np2, q 
       call  evalnp( val, pntint, l, poly, 1, n, lptint, &
                     lpoly )
       lafla( l ) = lafla( l ) - lafla( 1 )*val
115     continue

   else
     if (i==np2) then
!
!  ALL FINITE DIFFERENCES CORRESPONDING TO DEGREES = 2 ARE UPDATED.
!  THE NEWTON FUNDAMENTAL POLYNOMIALS OF DEGREE CARDINALITY OF ALFA_I = 1
!  ARE USED IN THE UPDATE
! 
       do 120 l = np2,q
         do 130 k = 2, np1-neqcon
           call  evalnp( val, pntint, l, poly, k, n, lptint, &
                         lpoly )
           lafla( l ) = lafla( l ) - lafla( k )*val
130          continue
120        continue

     endif
   endif
100  continue


 return
 end













subroutine funcon(mode, ncnln, n, nrowj, needc, x, c, &
                  cjac, nstate)   


integer          mode, ncnln, n, nrowj, needc(*), nstate

double precision x(n), c(*), cjac(nrowj,n)

integer           i, j, k, neasy, l, ncont

double precision  half

parameter        (half=0.5d0)
include 'dfo_model_inc.inc'

neasy=nlin+nnln
!
!  IF WE ARE DEALING WITH IPOPT, THEN MODE=1,2,3 INDICATES WHETHER
!  THE VALUES OF THE CONSTRAINTS, THE JACOBIAN OF THE HESSIAN OF THE 
!  LAGRANGIAN IS EXPECTED
!
if (useipopt > 0) then
  if ( mode == 1) then
    do 5  k=1, nlin
      c(k)=0.0d0
      do 6 i=1, n
        c(k)= c(k)+amat(k, i)*x(i)
6       continue
5     continue  
    if ( nnln > 0 ) call easycon(n, x, nnln, c(nlin+1))
    if ( .not. usemerit) then
      do 10 k = 1, ncon
        l=k+neasy
        c(l)=ccon(k)
        do 30 i = 1, n
          c(l)=c(l)+lcon((k-1)*n+i)*x(i)
          do 20 j = 1, n
            c(l)=c(l)+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
20           continue
30         continue
10       continue  
    endif
  elseif (mode == 2) then    
    do 15  k=1, nlin
      do 16 i=1, n
        cjac(k,i)= amat(k, i)
16        continue
15     continue  
    if ( nnln > 0 ) call easyjac(n, x, nnln, nrowj, &
                       cjac(nlin+1,1))
    if ( .not. usemerit) then
      do 60 k = 1, ncon
        l=k+neasy
        do 50 i = 1, n
          cjac(l, i)=lcon((k-1)*n+i)
          do 40 j = 1, n
            cjac(l, i)=cjac(l, i)+qcon((k-1)*n*n+(i-1)*n+j)*x(j)
40           continue
50         continue
60       continue 
    endif 
   else if ( mode == 3 ) then
    do 65  k=1, n
      do 66 i=1, n
        cjac(k,i)= hmod(k,i)
66        continue
65      continue  
    do 70 k = nlin + 1, neasy
      call easyhess(k, n, x, nnln, nrowj, cjac, c(k))
70     continue  
    if ( .not. usemerit) then
      do 100 k = 1, ncon
        l=k+neasy
        do 90 i = 1, n
          do 80 j = 1, n
            cjac(j, i)=cjac(j, i)+c(l)*qcon((k-1)*n*n+(i-1)*n+j)
80           continue
90         continue
100       continue 
    endif           
   endif
 else
  call easycon(n, x, nnln, c(1))
  call easyjac(n, x, nnln, nrowj, cjac)
  if ( .not. usemerit ) then
    do 110 k = 1, ncon
      l=k+nnln
      c(l)=ccon(k)
      do 130 i = 1, n
      c(l)=c(l)+lcon((k-1)*n+i)*x(i)
      cjac(l, i)=lcon((k-1)*n+i)
      do 120 j = 1, n
        c(l)=c(l)+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
        cjac(l, i)=cjac(l, i)+qcon((k-1)*n*n+(i-1)*n+j)*x(j)
120       continue
130     continue
110    continue       
 endif
endif   
return
end












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
!       IF (.NOT.( F < 10.0D20 .AND. F > -10.0D20 ))
!     +            IFERR = .TRUE.
!
!       DO 10 I=1, M
!         IF (.NOT.( C(I) < 10.0D20 .AND. C(I) > -10.0D20 ))
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
 if ( method /= 4 ) then 
   do 10 i =1, m
     if (cl(i)-cnstol>c(i)) then
        val = val + pp*(cl(i)-c(i))
     elseif  (cu(i)+cnstol<c(i)) then
        val = val + pp*(c(i)-cu(i))
     endif
10    continue
 else
   do 20 i =1, m
     if (cl(i)-cnstol>c(i)) then
        val = val + pp*(cl(i)-c(i))**2
     elseif  (cu(i)+cnstol<c(i)) then
        val = val + pp*(c(i)-cu(i))**2
     endif
20    continue
 endif
 return
 end









! subroutine funmrt(mode, n, x, objf, objgrd, nstate)
    
! double precision objf, x(n), objgrd(n), val

! integer          mode, n, nstate, i, j, k

! double precision  half

! parameter        (half=0.5d0)

! if (mode/=2) then
! !  GIVE A WARNING

! endif

! objf=0
! do 10 i=1,n
!   objf=objf+gmod(i)*x(i)
!   objgrd(i)=gmod(i)
!   do 20 j=1,n
!     objf=objf+half*hmod(i,j)*x(j)*x(i)
!     objgrd(i)=objgrd(i)+hmod(i,j)*x(j)
! 20   continue
! 10 continue

! do 30 k=1,ncon
!   val=ccon(k)
!   do 40 i=1,n
!   val=val+lcon((k-1)*n+i)*x(i)
!     do 50 j=1,n
!       val=val+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
! 50     continue
! 40   continue
!   if ( val < conl(k) ) then 
!     objf = objf + penpar * ( conl(k) - val )
!     do 60 i=1,n
!       objgrd(i)=objgrd(i)-penpar*lcon((k-1)*n+i)
!       do 70 j=1,n
!         objgrd(i)=objgrd(i)-penpar*qcon((k-1)*n*n+(i-1)*n+j)*x(j)
! 70       continue
! 60     continue
!   elseif ( val > conu(k) ) then 
!     objf = objf + penpar * ( val - conu(k) )
!     do 80 i=1,n
!       objgrd(i)=objgrd(i)+penpar*lcon((k-1)*n+i)
!       do 90 j=1,n
!         objgrd(i)=objgrd(i)+penpar*qcon((k-1)*n*n+(i-1)*n+j)*x(j)
! 90       continue
! 80     continue
!   endif
! 30 continue        

! return
! end





subroutine funobj(mode, n, x, objf, objgrd, nstate)
    
double precision objf, x(n), objgrd(n)

integer          mode, n, nstate, i, j, k

double precision  half, val

parameter        (half=0.5d0)
!
!  IF IPOPT FLAG IS >0 THEN WE ARE CALLED BY IPOPT AND THEN MODE=1 MEANS
!  THAT THE OBJECTIVE FUNCTION VALUE IS EXPACTED. IF MODE = 2 THEN
!  THE GRADIENT OF THE OBJECTIVE IS REQUIRED. IF IPOPT =< 0 THEN
!  NPSOL IS CALLING, THE COMPUTE OBJF AND OBJGRAD AT ONCE.
!
include 'dfo_model_inc.inc'

if ( useipopt > 0 ) then 
  if (mode == 1 ) then      
    objf=0
    do 10 i=1,n
      objf=objf+gmod(i)*x(i)
      do 20 j=1,n
        objf=objf+half*hmod(i,j)*x(j)*x(i)
20       continue
10     continue
    if ( usemerit ) then
      do 170 k=1,ncon
      val=ccon(k)
      do 180 i=1,n
        val=val+lcon((k-1)*n+i)*x(i)
        do 190 j=1,n
          val=val+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
190         continue
180       continue
      if ( val < conl(k) ) then 
        objf = objf + penpar * ( conl(k) - val )
      elseif ( val > conu(k) ) then 
        objf = objf + penpar * ( val - conu(k) )
      endif
170       continue  
    endif
  else
    do 30 i=1,n
      objgrd(i)=gmod(i)
      do 40 j=1,n
        objgrd(i)=objgrd(i)+hmod(i,j)*x(j)
40       continue
30     continue
    if ( usemerit ) then
      do 270 k=1,ncon
        val=ccon(k)
        do 280 i=1,n
          val=val+lcon((k-1)*n+i)*x(i)
          do 290 j=1,n
            val=val+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
290           continue
280         continue
        if ( val < conl(k) ) then 
          do 200 i=1,n
            objgrd(i)=objgrd(i)-penpar*lcon((k-1)*n+i)
            do 210 j=1,n
              objgrd(i)=objgrd(i)-penpar*qcon((k-1)*n*n+(i-1)*n+j) &
              *x(j)
210             continue
200           continue
        elseif ( val > conu(k) ) then 
          do 220 i=1,n
            objgrd(i)=objgrd(i)+penpar*lcon((k-1)*n+i)
            do 230 j=1,n
              objgrd(i)=objgrd(i)+penpar*qcon((k-1)*n*n+(i-1)*n+j) &
              *x(j)
230             continue
220           continue
        endif
270       continue
    endif
  endif  
else
  objf=0
  do 50 i=1,n
    objf=objf+gmod(i)*x(i)
    objgrd(i)=gmod(i)
    do 60 j=1,n
      objf=objf+half*hmod(i,j)*x(j)*x(i)
      objgrd(i)=objgrd(i)+hmod(i,j)*x(j)
60     continue
50   continue
  if ( usemerit ) then
    do 70 k=1,ncon
      val=ccon(k)
      do 80 i=1,n
        val=val+lcon((k-1)*n+i)*x(i)
        do 90 j=1,n
          val=val+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
90         continue
80       continue
      if ( val < conl(k) ) then 
        objf = objf + penpar * ( conl(k) - val )
        do 100 i=1,n
          objgrd(i)=objgrd(i)-penpar*lcon((k-1)*n+i)
          do 110 j=1,n
            objgrd(i)=objgrd(i)-penpar*qcon((k-1)*n*n+(i-1)*n+j) &
            *x(j)
110           continue
100         continue
      elseif ( val > conu(k) ) then 
        objf = objf + penpar * ( val - conu(k) )
        do 120 i=1,n
          objgrd(i)=objgrd(i)+penpar*lcon((k-1)*n+i)
          do 130 j=1,n
            objgrd(i)=objgrd(i)+penpar*qcon((k-1)*n*n+(i-1)*n+j) &
            *x(j)
130           continue
120         continue
      endif
70     continue  
  endif
endif
return
end





!
!
!
subroutine getnp( ipoly, poly, lpoly, n, kappa, g, h)


!
!     collects  the coefficient of the polynomial IPOLY into
!     doefficient of a quadratic kappa+g'x+0.5 x'H'x


!     INPUT
!     POLY   = Newton polynomials
!     IPOLY  = The index of the polynomial we are processing

!     OUTPUT
!     KAPPA = the constant of the polynomail
!     G     = the linear coeffitients
!     H     = the n by n matrix of quadratic doeffitients
integer           ipoly, lpoly, n

double precision  poly(lpoly), kappa,  g(n), h(n,n)
!
!     Local variables
!

double precision  zero
parameter        (zero=0.0d0)
integer           i, j, kbeg, k, np1, dd

np1=n+1
dd=(np1)*(n+2)/2

kappa=zero
do 10 i=1,n
  g(i)=zero
  do 20 j=1,n
    h(i,j)=zero
20  continue
10 continue

 
if (ipoly==1) then
  kappa=poly(1)
elseif (ipoly <= np1) then
  kbeg=2+(ipoly-2)*(np1)
  kappa=poly(kbeg)
  do 30 i=1, n
   g(i)=poly(kbeg+i)
30   continue
else        
  kbeg=2+n*np1+(ipoly-np1-1)*dd
  kappa=poly(kbeg)
  do 40 i=1, n
    g(i)=poly(kbeg+i)
40   continue
  k=kbeg+n+1
  do 50 i=1, n
    h(i,i)=2*poly(k)
    k=k+1
    do 60 j=i+1,n
      h(i,j)=poly(k)
      h(j,i)=poly(k)
      k=k+1
60     continue
50   continue
endif


return
end
 subroutine mterms(kappa , b   , h, lafla , poly , &
                   pntint, nind, n, neqcon, lpoly, lptint)


!    
!  ***********************************************************************
!  THIS SUBROUTINE COMPUTES THE TERMS OF THE QUADRATIC
!  INTERPOLATING POLYNOMIAL 
!                     T         T
!      M(X)= KAPPA + G X + 0.5*X H X
!
!  USING FINITE DIFFERENCES, FOUND IN 'FG'
!
!  PARAMETERS:
!  
!   N      (INPUT)  DIMENTION OF THE PROBLEM
!
!   NIND   (INPUT)  NUMBER OF INTERPOLATION POINTS
! 
!   PNTINT (INPUT)  LIST OF  'NIND' DATA POINTS. THE I-TH POINT OCCUPIES
!                   POSITIONS ( I - 1 ) * N + 1 TO I * N.
!   POLY   (INPUT)  THE ARRAY CONTAINING COEFFICIENTS OF NEWTON FUNDAMENTAL
!                   POLYNOMIALS (AS COMPUTED BY 'NBUILD')
!   LAFLA  (INPUT)  THE ARRAY OF FINITE DIFFERENCES
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



double precision kappa, b(n), h(n,n), lafla(nind), &
                 poly(lpoly), pntint(lptint)

integer          nind,n, lpoly, lptint, neqcon

! 
!  LOCAL VARIABLES
!

integer          np1, dd, ndd, nindm1, ndmnp1, &
                 i, j, k, kk, km, jj

       


 np1    = n + 1 
 dd     = (np1)*(n+2)/2
 ndd    = dd-np1
 nindm1 = nind  - 1 
 ndmnp1 = nindm1 - n 


!
!  SET THE COEFFICIENT OF THE INTERPOLATION TO ZERO
!


 call rzrmat(h, n, n)
 call rzrvec(b, n)
!
!  COMPUTE THE CONSTANT TERM
!

!
!  INITIALIZE USING CONSTANT BLOCK
!
 kappa=lafla(1)*poly(1)

!
!  UPDATE USING THE FINITE DIFFERENCES CORRESPONDING TO LINEAR BLOCK
!

 k=min(n-neqcon,nindm1)
 do 10 j=1,k
   kappa = kappa + lafla(j+1)*poly(2 + np1*(j-1))
10  continue

 if (nind>np1) then
!
!  UPDATE USING THE FINITE DIFF. CORRESPONDING TO QUADRATIC BLOCK
!
   k=min(ndd,ndmnp1)
   do 20 j=1,k
     kappa = kappa + lafla(j+np1)*poly(2 + np1*n + dd*(j-1))
20    continue
 endif


 
 do 40 i = 1,n
!
!  COMPUTE DEGREE ONE TERMS
!
   if (nind>1) then 

!
!  UPDATE USING LINEAR BLOCK
!     
     k=min(n-neqcon,nindm1)
     do 50 j = 1,k
       b(i) = b(i) + lafla(j+1)*poly(i + 2 + np1*(j-1))
50      continue
   endif

   if (nind >np1) then
!
!  UPDATE USING QUADRATIC BLOCK
!
     k=min(ndd,ndmnp1)
     do 60 j = 1,k
       b(i) = b(i) + lafla(j+np1)*poly(i + np1*n + 2 + dd*(j-1))
60      continue
   endif


40  continue



!
!  COMPUTE QUADRATIC TERMS USING QUADRATIC BLOCK
!

 if (nind>np1) then

!
!  POINTER TO FIRST QUADRATIC  BLOCK IN  'POLY'
!

   jj = np1*np1 + 2

!
!  LOOP OVER THE NUMBER OF SECOND DEGREE POLYNOMIALS
!
   
   km=min(ndd,ndmnp1)
   do 80 j = 1,km
     k=1
     do 90 i = 1,n
!
!  UPDATE DIAGONAL ELEMENT
!
       h(i,i) = h(i,i) + 2.0d0*lafla(j+np1)* &
                         poly(jj+k-1+dd*(j-1))
       k=k+1
!
!  UPDATE OFF-DIAGONAL ELEMENTS
!
       do 91 kk = i+1,n
         h(i,kk) = h(i,kk) + lafla(j+np1)*poly(jj+k-1+dd*(j-1))
         h(kk,i) = h(i,kk)
         k=k+1
91        continue
90      continue
80    continue
 endif



 return
 end









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

 subroutine nbuild(poly  , points, value, pntint, valint, sp2in  , &    !  
                   in2sp , nind  , n    , base  , distp , delta, &
                   pivthr, pivval, neqcon)                


!
!  **********************************************************************
!  THIS FUNCTION BUILDS A SET OF NEWTON FUNDAMENTAL POLYNOMIALS (NFPs) (AS
!  COMPLETE AS POSSIBLE WITH HIGHEST DEGREE 2) FOR THE SUPPLIED Q POINTS.
!
!  THE NFPs ARE BUILT BY BLOCKS IN THE FOLLOWING WAY:
!  THE ZERO-TH (CONSTANT) BLOCK CONSIST OF ONE ELEMENT. THIS ELEMENT IS
!  A POLYNOMIAL OF DEGREE ZERO; I.E., A CONSTANT (WE CHOOSE THE CONSTANT
!  TO BE 1). WE CHOOSE THE BASE TO CORRESPOND TO THIS  FIRST NFP.
!
!  THE FIRST BLOCK (ELEMENTS OF DEGREE UP TO ONE) CONSISTS OF 'N' NFPs.
!  TO COMPUTE AN K-TH ELEMENT IN THIS BLOCK WE CHOOSE A POINT FROM
!  THE GIVEN ARRAY OF POINTS, WHICH WAS NOT CHOSEN BEFORE, WE PUT THE
!  POINT IN THE ARRAY 'PNTINT' OF "CHOSEN" POINTS. THEN WE CONSTRUCT THE
!  K-TH POLYNOMIAL SO THAT ITS VALUE AT THAT K-TH POINT IS 1 AND
!  ITS VALUE AT ALL PREVIOUSLY CHOSEN POINTS IS 0. WE ALSO
!  UPDATE THE POLYNOMIAL IN THE SAME BLOCK SO THAT THEIR VALUE IS 0
!  AT THE K-TH POINT. WE STORE THE N+1 COEFFICIENTS OF THE K-TH NFP
!  IN ARRAY POLY (AFTER THE PREVIOUSLY COMPUTED POLYNOMIALS)
!
!  FOR THE QUADRATIC BLOCK THE PROCEDURE REPEATS IN A SIMILAR WAY, EXCEPT
!  THE NUMBER OF POLYNOMIALS IN THE QUADRATIC BLOCK CAN GO UP TO
!  N(N+1)/2 AND THEIR LENGTH IS ALWAYS (N+1)(N+2)/2.
!
!  WE STOP WHEN WE EITHER RUN OUT OF POINTS, WHICH SATISFY OUR
!  REQUIREMENTS (SEE BELOW) OR COMPLETE THE BASIS;
!  I.E., BUILD (N+1)(N+2)/2 POLYNOMIALS.
!
!  WE CHOOSE THE POINT FROM THE POOL BY TWO CRITERIA - PROXIMITY
!  TO THE BASE AND "INDEPENDENCE" (GEOMETRY, WELL-POISEDNESS). 
!  WE REQUIRE POINTS TO BE SUFFICIENTLY INDEPENDENT (WELL-POISED)
!  AND TO BE WITHIN CERTAIN DISTANCE FROM THE BASE.
!  WHEN WE HAVE A FEW ELIGIBLE POINTS WE CHOSE AMONG POINT CLOSER TO THE
!  BASE A POINT WITH THE BEST 'INDEPENDENCE' PROPERTY.
!
!  PARAMETERS
!
!  POINTS (INPUT)  ARRAY (N,Q) WITH THE POOL OF POTENTIAL INTERPOLATION POINTS
!
!  PIVTHR (INPUT)  PIVOT THRESHOLD VALUE
!
!  VALUE  (INPUT)  ARRAY (Q) OF VALUES AT THE POINTS
!
!  BASE   (INPUT)  INDEX OF THE BASE POINT
!
!  DISTP  (INPUT)  ARRAY (LVALUE) OF DISTANCES OF POINTS TO THE BASE POINTS
!
!  N      (INPUT)  PROBLEM DIMENSION
!
!  DELTA  (INPUT)  TRUST REGION RADIUS
!
!  POLY   (OUTPUT) LONG ARRAY CONTAINING ALL NEWTON FUNDAMENTAL  POLYNOMIALS
!
!  NIND   (OUTPUT) NUMBER OF POINTS INCLUDED IN THE INTERPOLATION
!
!  PNTINT (OUTPUT) ARRAY (N*NIND) OF POINTS INCLUDED IN  THE  INTERPOLATION
!                  IN THE ORDER OF INCLUSION
!  VALINT (OUTPUT) ARRAY (NIND) OF VALUES OF POINTS IN PNTINT
!
!  PIVVAL (OUTPUT) ARRAY (NIND) OF VALUES OF THE PIVOTS PRODUCED IN NBUILD
!
!  SP2IN  (OUTPUT) ARRAY OF 0/1 WHICH INDICATES WHICH POINTS IN "POINTS"
!                  ARE INCLUDED IN THE INTERPOLATION
!  IN2SP  (OUTPUT) ARRAY (NIND) OF INDICES WHICH FOR EVERY POINT IN "PNTINT"
!                  INDICATES ITS POSITION IN "POINTS"
!  NEQCON (INPUT)  NUMBER OF LINEARLY INDEPENDENT EQUALITY CONSTRAINTS
!
!  ****************************************************************************
!

 double precision  poly(lpoly)   , points(lpnts) , value(lvalue), &
                   pivthr        , distp(lvalue) , delta        , &
                   pntint(lptint), valint(lvlint), pivval(lvlint)

           
 integer           nind, n, base , sp2in(lvalue) , in2sp(lvlint), &
                   neqcon        
 

!
!  COMMON VARIABLES
!   
  
 integer           lpoly, lpnts, lvalue, lptint, lvlint, lconvl
 common /rpart/    lpoly, lpnts, lvalue, lptint, lvlint, lconvl
 save /rpart/
 
 integer           iout  , iprint
 double precision  mcheps, cnstol
 common / dfocm /  iout  , iprint, mcheps, cnstol
 save / dfocm /

 integer           npmin, layer, effort
 common / opti  /  npmin, layer, effort
 save / opti  /


!
!  LOCAL VARIABLES
!

 double precision val, vmax, delu, dell

 integer          q, qp1 , dd   , ddm1 , ndd   , np1   , np2 , &
                  i, jbeg, j    , jmax , nlayer, lv    , jend, &
                  k, ki  , ll   , kk   , pbase , ipbase, &
                  l, ipol, ipend, ibeg , iend


 double precision one, zero, thous, huge
 parameter       (one=1.0d0, zero=0.0d0, thous=10.0d3, &
                  huge=1.0d20)
!
!     SUBROUTINES AND FUNCTIONS CALLED:
!
!       APPLICATION:       NEXTNP, EVALNP,  SWAPNP, IZRVEC
!       FORTRAN SUPPLIED:  ABS
!       BLAS:              DCOPY 
!


!
! COMPUTE VALUES ASSOCIATED WITH N: DD - NUMBER OF TERMS IN A COMPLETE 
! QUADRATIC POLYNOMIAL, ALSO NUMBER OF NFP IN A COMPLETE QUADRATIC BASIS.
! NDD - NUMBER OF QUADRATIC TERMS IN A QUADRATIC POLYNOMIAL, ALSO NUMBER
! OF NFP IN A COMPLETE QUADRATIC BLOCK OF THE BASIS.
!

 np1  = n + 1
 np2  = n + 2
 dd   = np1*np2/2
 ddm1 = dd-1
 ndd  = dd-np1
 qp1  = lvalue+1
 q    = lvalue


 call izrvec(sp2in, lvalue)

!
!  INITIALIZE THE POLY , THE NEWTON POLYNOMIALS WITH THE MONOMIALS
!

!
!  CONSTANT BLOCK
!
 poly(1) = 1

!
!  LINEAR BLOCK
!
 do 20 i = 1,n
   jbeg  = (i-1)*np1+2
   jend  = jbeg + n
   do 10 j = jbeg, jend
     poly(j) = zero
10    continue
   poly(jbeg+i) = one
20  continue

!
!  QUADRATIC BLOCK
!
 do 50 i = 1,ndd
   jbeg = (i-1)*dd+n*np1+2
   jend = jbeg+n
   do 30 j = jbeg, jend
     poly(j) = zero
30    continue
   jbeg = jend+1
   jend = i*dd+n*np1+1
   do 40 j = jbeg, jend
     poly(j) = zero
40    continue
   poly(i+jbeg-1) = one
50  continue


!
!  *************************
!  START THE MAJOR ITERATION
!  *************************
!

!
! PUT THE BASE POINT TO BE THE FIRST IN THE NEW INTERPOLATION SET,
! UPDATE ALL RELATED VALUES APPROPRIATELY
!

 if (base/=1) then
   pbase     = in2sp(base)
   ipbase    = (pbase-1)*n+1
   call dcopy(n, points(ipbase), 1, pntint, 1)
   valint(1)     = value(pbase)
   sp2in(pbase)= 1
   in2sp(1)  = pbase
   base      = 1
 else
   sp2in(in2sp(1))=1
 endif

 pivval(1) = one

 call izrvec(in2sp(2), lvlint-1)


 do 180 i=2,dd
!
!  IF THE WE ARE IN THE LINEAR BLOCK BUT THE NUMBER OF LINEARLY
!  INDEPENDENT POINTS HAS REACHED IT'S MAXIMUM THEN COMPLETE
!  THE THE LINEAR BLOCK TILL THE END WITH DUMMY ELEMENTS
!  (SEE COMMENTS IN SUBROUTINE COMPL)
!
   if ( i > n + 1 - neqcon .and. i <= n + 1) then
      ibeg = 2 + (n+1)*(i-2) 
      iend = 1 + (n+1)*(i-1)
      do 53 j=ibeg, iend
        poly(j)=zero
53       continue  
      valint(i)=huge
      pivval(i)=huge
      nind = n+1    
      goto 180  
   endif     
!        
!  FIND THE POINT THAT GIVES THE LARGEST PIVOT (DOING PARTIAL PIVOTING)
!  WITHIN DIFFERENT LAYERS OF PROXIMITY TO  THE BASE:
!  FIRST TRY ALL POINTS WITHIN DELTA RADIUS, IF THERE ARE NONE
!  WITH SUFFICIENTLY LARGE PIVOT VALUE, SEARCH WITHIN 2*DELTA, ETC.
!  THE LARGEST NEIGHBORHOOD OF SEARCH HAS RADIUS NLAYER*DELTA
!  

   vmax =zero
   jmax = i
!
!  IF ALL  THE POINTS ARE ALREADY USED, THEN QUIT
!
   if ( i > q ) go to 85

!
! SET THE NUMBER OF LAYERS. BY SETTING IT TO A LARGER VALUE WE ALLOW
! TO INCLUDE (IN THE INTERPOLATION) POINTS THAT ARE FURTHER AWAY FROM THE
! BASE. THUS WE CAN ACHIEVE A MORE COMPLETE MODEL BUT ON A MORE RELAXED
! REGION
!
!
!  IF WE WANT TO BUILD AT LEAST 'NPMIN' ELEMENTS OF THE BASIS (NPMIN>=2)
!  THEN UNTIL WE DO SO WE ALLOW TO LOOK FOR POINTS IN A WIDER NEIGHBORHOOD
!  DEFINE THE WIDER NEIGHBORHOOD AS THE DISTANCE TO THE FURTHEST POINT
!
   val = zero
   do 55 k = 1, q
     if ( distp(k) > val ) val = distp(k)
55    continue  
   if ( i <= npmin ) then
     nlayer = val/delta + 1
   else
     nlayer = layer   
   endif 

!
!  CYCLE FOR DIFFERENT LAYERS OF PROXIMITY
!  
   do 80 lv=1,nlayer

     if ( i <= np1 ) then
       ipend = np1
     else
       ipend = dd
     endif
!
!  CYCLE OVER ALL POLYNOMIALS IN THE SAME BLOCK FROM I-TH AND ABOVE
!  AS SOON  AS WE FIND A POLYNOMIAL WHICH GIVES US A ACCEPTABLE PIVOT
!  WE QUIT THE LOOP. (WE DO NOT DO TOTAL PIVOTING HERE, TO DO IT, WE
!  WOULD CYCLE OVER POLYNOMIAL UNTIL WE CHOOSE THE ONE WITH BEST PIVOT)
!

     do 70 ipol=i, ipend  

!
!   UPDATE THE IPOLY-TH POLYNOMIAL SO THAT IT IS ZERO AT PREVIOUS POINTS
!

       call nextnp(ipol, poly, pntint, i-1 ,n, lpoly, lptint)
!
!  SET THE BOUNDARIES OF CURRENT LAYER, DELU - OUTER, DELL - INNER BOUND
!        
       delu = lv*delta + cnstol
       dell = (lv-1)*delta + mcheps*100
       do 60 j=1,q
!
!  IF THE J-TH POINT SATISFIES DISTANCE REQUIREMENTS THEN PICK IT AS
!  A CANDIDATE
!
         if ( distp(j) <= delu .and. distp(j) > dell &
              .and. sp2in(j)==0 ) then
!
!  EVALUATE THE I-TH POLYNOMIAL AT THE PICKED POINT
!
           call evalnp( val  , points, j, poly, ipol, n, &
                        lpnts, lpoly )
!
!  IF THE VALUE OF THE PIVOT IS LARGER THAN THE LARGEST PIVOT SO FAR
!  THEN ACCEPT THIS POINT AS THE BEST CANDIDATE
!

           if ( abs(val)>abs(vmax) + 100*mcheps ) then
             vmax=val
             jmax=j
           endif

         endif

60        continue

!
!  IF WITHIN CURRENT LAYER A SUFFICIENTLY LARGE PIVOT IS FOUND
!  THEN STOP LOOKING FOR MORE POINT AND MOVE ON TO INCLUDING THE
!  CURRENT BEST.
!
       if (abs(vmax)>=pivthr) go to 90

!
!  OTHERWISE MOVE TO THE  NEXT LAYER
!

70      continue

80    continue
!
!  ALL ACCEPTABLE LAYERS ARE SEARCHED  AND STILL
!  CANNOT FIND PIVOT LARGE ENOUGH
!

85    if (abs(vmax)<pivthr  ) then
     nind=i-1

!
!  COMPLETE OR INCOMPLETE BASIS IS BUILT, EXIT
!
     return
   endif



!
!   WE GET HERE ONLY WHEN THE CURRENT PIVOT IS FOUND AND IS LARGE ENOUGH
!

90    j   = jmax
   val = vmax
   if ( iprint >= 3 ) write(iout, 8000) j, ipol, val

!
!  PLACE POLYNOMIAL WITH INDEX 'IPOLY' ON THE I-TH PLACE
!  (BY SWAPPING I-TH AND IPOLY-TH POLYNOMIALS)
!  THIS IS DONE BECAUSE OF SOME SPECIFICS OF THE CODE
!  (BY EVENTUALLY MODIFYING THE CODE WE CAN AVOID SUCH SWAPPING )
!
   if ( ipol /= i )  call swapnp( n, i, ipol, poly, lpoly )   
!
! ASSOCIATE INTERPOLATION POINT J WITH THE I-TH POLYNOMIAL
! WRITE THE POINT INTO Y AND THE CORRESPONDING VALUE INTO VALINT
!

   call dcopy(n, points((j-1)*n+1), 1, pntint((i-1)*n+1), 1)
   valint(i) = value(j)
   sp2in(j)  = 1
   in2sp(i)  = j
   pivval(i) = val/delu

!
!  NORMALIZE THE I-TH POLYNOMIAL
!

   if (i==1) then
     poly(i) = one 
     ll = 1

   else if (i<=np1) then
     k  = (i-2)*(np1)+2
     kk = k + n 
     do 100 j = k,kk
       poly(j) = poly(j)/val
100      continue

   else 
     k  = (i-np2)*dd + n*np1 + 2
     kk = k + dd - 1
     do 110 j = k,kk
       poly(j) = poly(j)/val
110      continue         
   endif    

!
!  UPDATE THE I-TH POLYNOMIAL AGAIN SO THAT IT HAS VALUE 0 AT
!  ALL PREVIOUS POINT. IN EXACT ARITHMETIC THIS IS REDUNDANT,
!  BUT WE ARE DOING IT FOR STABILITY, SINCE IF VAL IS RATHER
!  SMALL THEN AN ERROR IS INTRODUCED (THIS IS SIMILAR TO 
!  A MODIFIED GRAMM-SCHMIDT PROCEDURE) 
!

   call nextnp(i, poly, pntint, i-1 ,n, lpoly, lptint)

!
!  RENORMALIZE THE I-TH POLYNOMIAL
!

   call evalnp(val, pntint, i,  poly, i, n, lptint, lpoly)


   if (i==1) then
     poly(i) = one 
     ll = 1

   else if (i<=np1) then
     k  = (i-2)*(np1)+2
     kk = k + n
     do 120 j = k,kk
       poly(j) = poly(j)/val
120      continue

   else 
     k  = (i-np2)*dd + n*np1 + 2
     kk = k + dd - 1
     do 130 j = k,kk
       poly(j) = poly(j)/val
130      continue         
   endif    




!
!  FINALLY, UPDATE ALL THE POLYNOMIALS IN THE SAME BLOCK SO
!  THAT THEY HAVE VALUE 0 AT THE NEW I-TH POINT
!   


   if (i==1) then
!
!   POLYNOMIAL IS A CONSTANT, DO NOTING
!
     continue

   else if (i<=np1) then

! 
!  WE ARE IN LINEAR BLOCK, ALL PREVIOUS POLYNOMIALS IN LINEAR BLOCK ARE UPDATED
!
     ki  = (i-2)*(np1) + 2
     do 150 l = 2, i-1
       call evalnp(val, pntint, i, poly, l, n, lptint, lpoly)
       k  = (l-2)*(np1) + 2
       do 140 j = 0,n
         poly(k+j)  =  poly(k+j) - poly(ki+j)*val
140        continue
150      continue
    
   else
         

!
!  WE ARE IN QUADRATIC BLOCK, UPDATE POLYNOMIALS IN QUADRATIC BLOCK
!

     ki = (i-np2)*dd + n*np1 + 2
     do 170 l = np2,i-1
       call evalnp(val, pntint, i, poly, l, n, lptint, lpoly)
       k  =  (l-np2)*dd + n*np1 + 2
       do 160 j = 0,dd-1
         poly(k+j) = poly(k+j) - poly(ki+j)*val
160        continue
170      continue         

           


   endif

180  continue

 nind=dd
 return
8000  format(' NBUILD: Pivoting: point: ',i5, ', polynomial: ',i5, &
        ', pivot value: ', d14.7 )
 end





subroutine swapnp( n, i, j, poly, lpoly )

!
!  *********************************************************
!  THIS SUBROUTINE SWAPS I-TH NEWTON POLYNOMIAL WITH THE J-TH 
!  NEWTON POLYNOMIAL. BOTH POLYNOMIAL ARE STORED IN ARRAY
!  'POLY'. BOTH POLYNOMIALS SHOULD BE OF THE SAME DEGREE.
!  *********************************************************
!

!
!  SUBROUTINE PARAMETERS
!

integer           n, i, j, lpoly
double precision  poly(lpoly)
!
!  GLOBAL VARIABLES
!
integer           iout  , iprint
double precision  mcheps, cnstol
common / dfocm /  iout  , iprint, mcheps, cnstol
save / dfocm /
!
!  LOCAL VARIABLES
!
integer           np1, np2, dd, kend, k, ki, kj 
double precision  val

!
!  PARAMETERS RELATED TO LENGTH AND NUMBER OF  POLYNOMIALS
!

np1 = n + 1
np2 = n + 2
dd  = np1*np2/2
! 
!  SET BEGINNING AND END POINTERS
!

!
!  IF WE THE SWAP IS PERFORMED IN THE LINEAR BLOCK
!
if (i <= np1 ) then
  if ( j > np1 ) then
    if ( iprint >= 0 ) write( iout, 1000 ) 
    stop
  endif
  kend = np1
  ki   = 1+(i-2)*np1
  kj   = 1+(j-2)*np1
else
!
!  IF WE THE SWAP IS PERFORMED IN THE QUADRATIC BLOCK
!
  if ( j <= np1 ) then
    write( iout, 1000 ) 
    stop
  endif
  kend = dd
  ki   = 1+n*np1+(i-np2)*dd
  kj   = 1+n*np1+(j-np2)*dd
endif
!
!  PERFORM THE SWAP
!
if ( j /= i ) then
  do 60 k=1,kend
    val        = poly(kj+k)
    poly(kj+k) = poly(ki+k)
    poly(ki+k) = val
60   continue
endif
return

1000 format( ' SWAPNP:  ERROR! TRYING TO SWAP POLYNOMIALS',/ &
        '          OF DIFFERENT DEGREE. MUST BE A BUG!',/)
end  





 subroutine nextnp(ipoly , poly , pntint, nind, n, lpoly, &
                   lptint)



! 
!  *******************************************************************
!  THIS SUBROUTINE UPDATES THE IPOLY-TH  NEWTON POLYNOMIAL SO THAT
!  IT HAS ZERO VALUE AT ALL POINTS IN THE SAME BLOCK AND BELLOW
!
!  PARAMETERS 
!
!    IPOLY  (INPUT) INDEX OF THE POLYNOMIAL WHICH IS UPDATED
!
!    POLY   (INPUT/OUTPUT) THE ARRAY (LPOLY) OF NEWTON POLYNOMIALS
!
!    PNTINT (INPUT) THE ARRAY (LPTINT) OF CURRENT INTERPOLATION POINTS
!
!    NIND   (INPUT) THE NUMBER OF POINTS IN THE INTERPOLATION
!
!    N      (INPUT) PROBLEM DIMENSION
!  *******************************************************************
!


!
!  PARAMETERS
!

 double precision  poly(lpoly), pntint(lptint)
          

 integer           ipoly, nind ,n, lpoly, lptint

!
!  PROCESS CONTROL PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol
common / dfocm / iout  , iprint, mcheps, cnstol     
save / dfocm /

!
!  LOCAL VARIABLES 
!            
 double precision val

 integer          np1, dd, q, block, npbeg, npend, i, j

!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!    APPLICATION     :       EVALNP
!    FORTRAN SUPPLIED:       MAX
!

 np1 = n + 1 
 dd  = np1*(n+2)/2
 q   = max(ipoly, nind)
!
!  CHECK IN WHICH BLOCK WE ARE DOING THE UPDATES
!

 if ( ipoly <= np1 ) then
   block=1
   if ( nind > np1 ) then
     if ( iprint > 0 ) write( iout, 1000 )
     stop
   endif
 else
   block=2
 endif

 if (block==1) then

!
!  IF A POLYNOMIAL IN THE LINEAR BLOCK IS UPDATED
!
   npbeg =(ipoly-2)*np1 + 2
   npend = npbeg + n
   call evalnp(val, pntint, 1, poly, ipoly, n, lptint, lpoly)
   poly(npbeg)=poly(npbeg)-val*poly(1)

   do 20 j=2,nind
     if ( j /= ipoly ) then
       call evalnp(val, pntint, j, poly, ipoly, n, lptint, lpoly) 
       do 10 i=npbeg,npend
         poly(i)=poly(i)-val*poly(i-npbeg+2+(j-2)*np1)
10        continue 
     endif 
20    continue
 
!
!  IF A POLYNOMIAL IN THE QUADRATIC BLOCK IS UPDATED
!

 else
   npbeg = (ipoly-n-2)*dd + n*np1 + 2
   call evalnp(val, pntint, 1, poly, ipoly, n, lptint, lpoly)
   poly(npbeg) = poly(npbeg) - val*poly(1)

!
!  "ORTHOGONALIZE" THE IPOLY-TH POLYNOMIAL WITH THE LINEAR BLOCK
!

   do 40 j=2, np1
      call evalnp(val, pntint, j, poly, ipoly, n, lptint, &
                  lpoly)
      npend = npbeg + n
      do 30 i = npbeg, npend
        poly(i) = poly(i) - val*poly(i-npbeg+2+(j-2)*np1)
30       continue  
40    continue

!
!  "ORTHOGONALIZE" THE IPOLY-TH POLYNOMIAL WITH THE QUADRATIC BLOCK
!
   do 60 j = n+2, nind
     if ( j /= ipoly ) then
       call evalnp(val, pntint, j, poly, ipoly, n, &
                   lptint, lpoly)
       npend = npbeg + dd - 1
       do 50 i = npbeg,npend
         poly(i) = poly(i) - val*poly(i-npbeg+2+n*np1+(j-n-2)*dd)
50        continue  
     endif
60    continue
 endif

 return
1000  format(' NEXTNP: ERROR! UPDATING POLYNOMIAL IN LINEAR BLOCK',/ &
        '         WHEN QUADRATIC BLOCK IS PRESENT',/)  
 end



 subroutine ptexch(poly , xnew  , pntint, ipoly , nind, n, &
                   lpoly, lptint, pivthr, pivval, vmax, fail )




!
!  **********************************************************************
!  THIS SUBROUTINE INCLUDES A POINT XNEW IN THE INTERPOLATION SET BY 
!  UPDATING THE SET OF NEWTON FUNDAMENTAL POLYNOMIALS ACCORDINGLY.
!  FIRST A POINT WHICH IS  REPLACED BY XNEW IS DETERMINED, THEN
!  SUBROUTINE "PTREPL" IS CALLED, WHICH PERFORMS THE UPDATES.
!  IPOLY IS THE INDEX OF THE POINT THAT IS REPLACED BY 'XGNEW'.
!  THE SET OF NEWTON POLYNOMIALS IS UPDATED ACCORDINGLY.
!
!  IF THERE ARE QUADRATIC POLYNOMIALS IN THE INTERPOLATION BASIS THEN
!  WE EVALUATE EACH POLYNOMIAL IN THE QUADRATIC BLOCK AT 'XNEW',
!  OTHERWISE WE EVALUATE POLYNOMIALS OF THE LINEAR BLOCK AT 'XNEW'.
!  WE CHOOSE THE POLYNOMIAL  THAT HAS THE LARGEST VALUE AT XNEW.
!  IF THIS VALUE IS LARGER THEN THE PIVOT THRESHOLD THEN WE SET
!  TO THE INDEX OF THAT POLYNOMIAL.
!  OTHERWISE  WE DECLARE THAT WE FAILED TO INCLUDE XGNEW IN THE 
!  INTERPOLATION SET AND RETURN.
!
!  PARAMETERS
!    POLY   (INPUT/OUTPUT) THE SET OF NEWTON POLYNOMIALS
!
!    XNEW   (INPUT) POINT THAT SHOULD BE ADDED
!
!    PNTINT (INPUT) SET OF INTERPOLATION POINTS BEFORE XNEW IS ADDED 
!
!    NIND   (INPUT) NUMBER OF INTERPOLATION POINTS
!
!    N      (INPUT) PROBLEM DIMENSION
!
!    PIVTHR (INPUT) PIVOT THRESHOLD VALUE
!
!    PIVVAL (INPUT) ARRAY OF PIVOT VALUES 
!
!    FAIL   (OUTPUT) INDICATES IF THE POINT INTERCHANGE FAILED
!
!    IPOLY  (OUTPUT) THE INDEX OF THE POINT THAT WAS REPLACED, IF FAIL=.FALSE.
!
!    VMAX   (OUTPUT) THE VALUE OF THE PIVOT, IF FAIL=.FALSE.
!  *******************************************************************
!

!
!  PARAMETER VARIABLES
!
double precision  poly(lpoly)   , xnew(n), vmax, &
                  pntint(lptint), pivthr , pivval(nind)
     

            
integer           nind, n, lpoly, lptint, ipoly
 
logical           fail

!
!  COMMON VARIABLES
!
integer          iout, iprint

double precision mcheps, cnstol
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /

!
!  LOCAL VARIABLES 
!
integer          dd, np1, np2, block, i, base

double precision val

!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!    APPLICATION     :       EVALX, PTREPL
!    FORTRAN SUPPLIED:       ABS
!
 fail =.false.     
 base = ipoly
 np1  = n + 1 
 np2  = n + 2      
 dd   = np1*np2/2


!
!  CHECK IF INTERPOLATION SET CONTAINS ELEMENTS FROM QUADRATIC BLOCK OR NOT
!

 if ( nind <= np1 ) then
   block=1
 else
   block=2
 endif


!
!  TO IDENTIFY IPOLY WE  EVALUATE ALL NEWTON POLYNOMIALS IN THE LAST
!  BLOCK AT THIS POINT AND CHOOSE THE ONE THAT GIVES THE LARGE PIVOT.
!  IF THIS PIVOT IS LARGER THEN THE PIVOT THRESHOLD THEN WE SET IPOLY TO
!  THE INDEX OF THIS POLYNOMIAL. 
!

 vmax  = 0.0d0
 ipoly = 0
 if (block==2) then

!
!  IF WE START WITH THE QUADRATIC BLOCK
!

   do 10 i=np2, nind
     if ( i /= base ) then
       call evalx(val, xnew, poly, i, n,  lpoly)
       if ( abs(val) > abs(vmax) ) then
         vmax=val
         ipoly=i
       endif
     endif
10    continue
   if ( ipoly == 0 .or. abs(vmax*pivval(ipoly)) < pivthr ) &
     then
     fail=.true.
     return
   endif
 else

!
!  IF WE START FROM THE LINEAR BLOCK
!
   do 20 i=2, nind
     call evalx(val, xnew, poly, i, n, lpoly)
     if ( abs(val) > abs(vmax) ) then
       vmax  = val
       ipoly = i
     endif
20    continue
   if ( ipoly == 0 .or. abs(vmax*pivval(ipoly)) < pivthr ) &
     then
     fail=.true.
     return
   endif   
 endif
!
!  UPDATE THE NEWTON POLYNOMIALS, SO THAT POLYNOMIAL WITH THE
!  INDEX IPOLY CORRESPONDS TO XGNEW IN THE INTERPOLATION SET 
!

 call ptrepl( xnew , ipoly, vmax, pntint, poly, nind, n, &
              lpoly, lptint)



 return
 end










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
  if (points(i)<lb(i)) then
    points(i) = lb(i)
    bdvltd    = .true.
  elseif (points(i)>ub(i)) then
    points(i) = ub(i)
    bdvltd    = .true.
  endif
10 continue


!
!  IF PROJECTED, PRINT WARNING MESSAGE 
!
if (bdvltd) then
  if ( iprint >= 0 ) write(iout, 1000)
  call dcopy( n, points, 1, x, 1 ) 
endif

!
!  CHECK FEASIBILITY OF FIRST POINT WRT LINEAR CONSTRAINTS
!

do 20 i=1, nclin
  val=ddot(n, a(i), lda, points, 1)
  if (val > ub(n+i) .or. val < lb(n+i)) infeas=.true.
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
if ( ncnln > 0 .or. infeas) then
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



  if ( inform == 1 ) then
    if ( iprint > 0 ) write(iout,2000)
    return
  elseif (inform == 2) then
!
!  IF FEASIBLE SOLUTION WAS NOT FOUND, TRY TO INCREASE 
!  TRUST REGION RADIUS 
!
    if (del < delmax) then
      del=del*2
      goto 50
    else
!
!  IF THE TRUST REGION RADIUS IS AT IT'S MAXIMUM VALUE
!  THEN ASSUME THE PROBLEM IS INFEASIBLE AND QUIT
!
      if( iprint>0 )   write(iout,2010)
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
  if ( iprint >= 0 ) then
    if (val > 100*n*mcheps) write(iout, 1010)
  endif
endif

if (val > 100*n*mcheps .or. bdvltd) then

!
!  COMPUTE FUNCTION VALUE FOR THE FIRST POINT IF IT
!  HAD TO BE PROJECTED. IF SUCH COMPUTATION
!  FAILS, THEN QUIT THE PROGRAM
!
  nf = nf + 1
  if ( scale /= 0 ) call unscl( n, points, scal )
  call fun( n, m, points, obfval, conval, iferr )
  if ( scale /= 0 ) call scl( n, points, scal )
  if ( iferr ) then
    if( iprint>0 ) write(iout,2020)
    inform = -1
    return
  endif
  if ( nf >= maxnf ) then
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
    if ( points(inp+i-1)<lb(i)-cnstol .or. &
         points(inp+i-1)>ub(i)+cnstol ) then 
      if ( iprint >= 0 ) write(iout, 1020) ix
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
    if ( val > ub(n+i)+cnstol .or. &
         val < lb(n+i)-cnstol ) then
      if ( iprint >= 0 ) write(iout, 1030) ix
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
  if ( ncnln > 0 ) then
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



    if ( inform == 1 ) then
      if ( iprint > 0 ) write(iout,2000)
      return
    elseif (inform == 2) then
!
!  IF NO FEASIBLE SOLUTION WAS FOUND THEN SKIP THIS POINTS AND
!  MOVE TO THE NEXT ONE
!
      if ( iprint > 0 ) write( iout, 1040 ) ix
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
    if (val > 100*n*mcheps) then
      if ( iprint > 0 ) write(iout, 1040) ix
      go to 120
    endif
  endif
!
!  IF THE POINT PASSED ALL FEASIBILITY TESTS, THEN ACCEPT IT AS
!  A SAMPLE POINT AND RECORD ITS FUNCTION VALUES AND CONSTRAINT VALUES
!
  np           = np  + 1
  inp          = inp + n
  if ( np /= ix ) then
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

if ( np0 == 1 ) then
  if ( iprint >= 2 ) write( iout, 8000 ) 
  call dcopy( n, points, 1, x, 1 )
!
!  GENERATE N RANDOM NUMBERS BETWEEN 0 AND 1 AND RECORD THEM TO
!  POINTS STARTING FROM N+1-ST ENTREE
!       

  call ranlux(points(n+1), n)
  do 130 j = 1, n
    distb=zero
    distb=min(delta, ub(j)-x(j))
    if (distb>mcheps) then
      points( n + j ) = distb * points(n+j) + x( j )
    endif 
    if (distb<=mcheps) then
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
    if (val > ub(n+i) .or. val < lb(n+i)) infeas = .true.
140   continue
  val = zero
!  -------------------------------------------------------
!  FIND THE SECOND POINT FOR INTERPOLATION
!  -------------------------------------------------------

  if ( ncnln > 0 .or. infeas) then          
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

    if (inform == 1) then
      if( iprint>0 )   write(iout,2000)
      return
    elseif (inform == 2) then
      if( iprint>0 )   write(iout,2010)
      return
    endif


!  -------------------------------------------------------
!  IF FIRST AND SECOND POINTS COINCIDE, FIND A DIFFERENT SECOND POINT
!  -------------------------------------------------------

    val = zero
    do 170 i = 1, n 
      val = val + abs(x(i)-points(i))
170     continue

    if (val < 10*n*pivthr) then
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

      if ( inform == 1 ) then
        if( iprint>3 )   write(iout,2000)
        return
      elseif ( inform == 2 ) then
        if( iprint>3 )   write(iout,2010)
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


200   if ( scale /= 0 ) call unscl( n, points( n+1 ), scal )
  call fun( n, m, points( n+1 ),  obfval( 2 ), conval( m+1 ), &
            iferr )
  if ( scale /= 0 ) call scl( n, points( n+1 ), scal )
  nf = nf + 1
  if ( nf >= maxnf ) then
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
    if ( dist(2) < pivthr ) then  
      inform = -1
      return 
    endif
!
!  IF THE SECOND POINTS BECOMES NON-FEASIBLE, QUIT
!
    if ( ncnln> 0 ) then 
!            CALL FUNCON(1  , NCNLN , N  , NCNLN , IWRK,
!     *                  POINTS(N+1), WRK, WRK(NCNLN+1), 1  )
       call easycon(n, points(n+1), ncnln, wrk)
      do 220 j=1, ncnln 
        if ( wrk(j) < lb(n+nclin+j) - cnstol .or. &
             wrk(j) > ub(n+nclin+j) + cnstol     ) then
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
  if ( dist(i) > layer*delta ) bdvltd = .true.
250 continue
if ( bdvltd ) then
  if ( iprint >= 0 ) write(iout,1050)
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
  if ( values(i) < values(base) + 1.0d2*mcheps ) base=i
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



 subroutine ptrepl( xnew , ipoly, vnew, pntint, poly, nind, n, &
                    lpoly, lptint)

!
!  **********************************************************************
!  THIS SUBROUTINE REPLACES A  POINT WITH A GIVEN INDEX IN THE INTERPOLATION 
!  SET BY ANOTHER GIVEN POINT AND UPDATES NEWTON FUNDAMENTAL POLYNOMIALS
!  ACCORDINGLY. IT CAN ALSO BE USED TO ADD A POINT TO AN INCOMPLETE
!  INTERPOLATION SET, THIS CAN BE DONE BY SETTING INDEX 
!  IPOLY = NIND+1  = NUMBER OF INTERPOLATION POINTS + 1.
!  IF THE INTERPOLATION IS QUADRATIC ( NIND > N+1 ) THEN WE ARE ALLOWED
!  TO REPLACE POINTS ONLY IN THE QUADRATIC BLOCK, THUS IPOLY > N + 1.
!
!  PARAMETERS
! 
!    POLY   (INPUT/OUTPUT) THE SET OF NEWTON FUNDAMENTAL POLYNOMIALS
!
!    XNEW   (INPUT) POINT THAT IS INTRODUCED INTO THE SET
!
!    IPOLY  (INPUT) THE INDEX OF THE POINT THAT HAS TO BE REPLACED
!
!    VNEW   (INPUT) VALUE OF IPOLY-TH POLYNOMIAL AT XGNEW
!
!    PNTINT (INPUT) SET OF INTERPOLATION POINTS
!
!    NIND   (INPUT) NUMBER OF INTERPOLATION POINTS
!
!    N      (INPUT) PROBLEM DIMENSION
!  ********************************************************************
!


!
!  PARAMETERS
!


 double precision poly(lpoly)  , xnew(n), vnew, &
                  pntint(lptint)                 
                 
            
 integer          nind, n, lpoly, lptint, ipoly

!
!  COMMON VARIABLES
!

 integer          iout  , iprint
 double precision mcheps, cnstol
 common / dfocm / iout  , iprint, mcheps, cnstol
 save / dfocm /

!
!  LOCAL VARIABLES 
!

 double precision val

 integer          np1, dd, npbeg, npend, i, j, k, block

!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!    APPLICATION     :       NEXTNP, EVALX
!    FORTRAN SUPPLIED:       ABS
!


 np1 = n+1
 dd  = np1*(n+2)/2

!
!  CHECK TO WHICH BLOCK WE ARE ADDING THE POINT, CLEARLY NIND<DD
!
 if ( ipoly <= np1 ) then
   block=1
   if ( nind > np1 ) then
     if ( iprint > 0 ) write(iout,1000) ipoly
     stop 
   endif
 else
   block=2
 endif
!
!  IF WE ONLY DEAL WITH THE LINEAR BLOCK (NO QUADRATIC POLYNOMIALS YET)
!
 if (block==1) then

! 
!  "NORMALIZE" THE IPOLY-TH  POLYNOMIAL
!
   npbeg = (ipoly-2)*np1 + 2
   npend =  npbeg + n 
   do 10 i=npbeg,npend
     poly(i)=poly(i)/vnew
10    continue 

!
!  FOR NUMERICAL STABILITY  "RE-ORTHOGONALIZE" THE NEXT POLYNOMIAL
!  WITH RESPECT TO THE PREVIOUS POLYNOMIALS
! 

   if ( abs(vnew) < 1.0d2 ) then
     call nextnp(ipoly, poly, pntint, nind, n, lpoly, lptint)
     call evalx(val, xnew, poly, ipoly, n, lpoly)  
     do 20 i = npbeg, npend
       poly(i) = poly(i)/val
20      continue  
   endif 

!
!  UPDATE THE POLYNOMIALS IN THE LINEAR BLOCK SO THAT THEIR VALUE IS 
!  ZERO AT  THE POINT XNEW 
!
   do 40 j=2, nind
     if ( j /= ipoly ) then 
       call evalx( val, xnew, poly, j, n, lpoly )
       do 30 i = npbeg, npend
         k = i - npbeg + (j-2)*np1 + 2
         poly(k) = poly(k) - val*poly(i)
30        continue 
     endif 
40    continue

!
!  IF WE ARE  UPDATING POLYNOMIALS IS IN THE QUADRATIC BLOCK
!
 else


!
!  "NORMALIZE" THE IPOLY-TH  POLYNOMIAL
!

   npbeg = (ipoly - n - 2)*dd + n*np1 + 2
   npend =  npbeg + dd - 1
   do 50 i = npbeg, npend
     poly(i)=poly(i)/vnew
50    continue 

!
!  FOR NUMERICAL STABILITY  "RE-ORTHOGONALIZE" THE NEXT POLYNOMIAL
!  WITH RESPECT TO THE PREVIOUS POLYNOMIALS
! 

   if ( abs(vnew) < 1.0d2 ) then
     call nextnp(ipoly, poly, pntint, nind, n, lpoly, lptint)
     call evalx(val, xnew, poly, ipoly, n, lpoly)  
     do 60 i = npbeg, npend
       poly(i) = poly(i)/val
60     continue  
   endif

!
!  UPDATE THE POLYNOMIALS IN THE QUADRATIC BLOCK SO THAT THEIR VALUE IS 
!  ZERO AT  THE POINT XNEW  
!
   do 80 j=n+2,nind
     if ( j /= ipoly ) then
       call evalx(val, xnew, poly, j, n,  lpoly)
       do 70 i = npbeg, npend
         k = i - npbeg + (j-n-2)*dd + n*(n+1) + 2
         poly(k) = poly(k)-val*poly(i)
70        continue 
     endif 
80    continue
 endif

 return


1000  format('PTREPL  **** TRYING TO REPLACE A POINT', /, &
       '             IN THE LINEAR BLOCK, IPOLY=', i6 ) 

 end









subroutine ranlux(rvec,lenv)
!         Subtract-and-borrow random number generator proposed by
!         Marsaglia and Zaman, implemented by F. James with the name
!         RCARRY in 1991, and later improved by Martin Luescher
!         in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993
!
!       references:
!  M. Luscher, Computer Physics Communications  79 (1994) 100
!  F. James, Computer Physics Communications 79 (1994) 111
!
!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:
!
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!  Calling sequences for RANLUX:                                  ++
!!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!!                   32-bit random floating point numbers between  ++
!!!!                   zero (not included) and one (also not incl.). ++
!!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!!               which is integer between zero and MAXLEV, or if   ++
!!!!               LUX > 24, it sets p=LUX directly.  K1 and K2   ++
!!!!               should be set to zero unless restarting at a break++
!!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!!               which can be used to restart the RANLUX generator ++
!!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!!               specify how many numbers were generated since the ++
!!!!               initialization with LUX and INT.  The restarting  ++
!!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!!   A more efficient but less convenient way of restarting is by: ++
!!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!!                 32-bit integer seeds, to be used for restarting ++
!!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer,intent(in) :: lenv
double precision rvec(lenv)

double precision :: seeds(24)
integer :: iseeds(24)
integer :: isdext(25)
integer,parameter :: maxlev=4
integer,parameter :: lxdflt=3
integer :: ndskip(0:maxlev)
integer :: next(24)
double precision,parameter :: twop12=4096.
integer,parameter :: igiga=1000000000
integer,parameter :: jsdflt=314159265
integer,parameter :: itwo24=2**24
integer,parameter :: icons=2147483563

double precision :: carry , twom12 , twom24 , uni
INTEGER :: i , i24 , ilx , in24 , inner , Inout ,    &
           Ins , inseed , iouter , isd , isk ,     &
           ivec , izip , izip2 , j24 , &
           jseed , k , K1 , K2 , kount , Lout , lp , &
           Lux , mkount , nskip

save notyet, i24, j24, carry, seeds, twom24, twom12, luxlev
save nskip, ndskip, in24, next, kount, mkount, inseed
integer luxlev
logical notyet
data notyet, luxlev, in24, kount, mkount /.true., lxdflt, 0,0,0/
data i24,j24,carry/24,10,0./
!                               default
!  Luxury Level   0     1     2   *3*    4
data ndskip/0,   24,   73,  199,  365 /
!orresponds to p=24    48    97   223   389
!     time factor 1     2     3     6    10   on slow workstation
!                 1    1.5    2     3     5   on fast mainframe
!
!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential
if (notyet) then
   notyet = .false.
   jseed = jsdflt
   inseed = jseed
!         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
   luxlev = lxdflt
   nskip = ndskip(luxlev)
   lp = nskip + 24
   in24 = 0
   kount = 0
   mkount = 0
!        WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
!    +        LUXLEV,'      p =',LP
      twom24 = 1.
   do 25 i= 1, 24
      twom24 = twom24 * 0.5
   k = jseed/53668
   jseed = 40014*(jseed-k*53668) -k*12211
   if (jseed < 0)  jseed = jseed+icons
   iseeds(i) = mod(jseed,itwo24)
25    continue
   twom12 = twom24 * 4096.
   do 50 i= 1,24
   seeds(i) = real(iseeds(i))*twom24
   next(i) = i-1
50    continue
   next(1) = 24
   i24 = 24
   j24 = 10
   carry = 0.
   if (seeds(24) == 0.) carry = twom24
endif
!
!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989
!
do 100 ivec= 1, lenv
uni = seeds(j24) - seeds(i24) - carry
if (uni < 0.)  then
   uni = uni + 1.0
   carry = twom24
else
   carry = 0.
endif
seeds(i24) = uni
i24 = next(i24)
j24 = next(j24)
rvec(ivec) = uni
!  small numbers (with less than 12 "significant" bits) are "padded".
if (uni < twom12)  then
   rvec(ivec) = rvec(ivec) + twom24*seeds(j24)
!        and zero is forbidden in case someone takes a logarithm
   if (rvec(ivec) == 0.)  rvec(ivec) = twom24*twom24
endif
!        Skipping to luxury.  As proposed by Martin Luscher.
in24 = in24 + 1
if (in24 == 24)  then
   in24 = 0
   kount = kount + nskip
   do 90 isk= 1, nskip
   uni = seeds(j24) - seeds(i24) - carry
   if (uni < 0.)  then
      uni = uni + 1.0
      carry = twom24
   else
      carry = 0.
   endif
   seeds(i24) = uni
   i24 = next(i24)
   j24 = next(j24)
90    continue
endif
100 continue
kount = kount + lenv
if (kount >= igiga)  then
   mkount = mkount + 1
   kount = kount - igiga
endif
return
!
!           Entry to input and float integer seeds from previous run
entry rluxin(isdext)
!     The following IF block added by Phillip Helbig, based on conversation
!     with Fred James; an equivalent correction has been published by James.
if (notyet) then
   write(6,'(A)')  ' PROPER RESULTS ONLY WITH INITIALISATION FROM &
25 INTEGERS OBTAINED WITH RLUXUT'
   notyet = .false.
endif
   twom24 = 1.
   do 195 i= 1, 24
   next(i) = i-1
195    twom24 = twom24 * 0.5
   next(1) = 24
   twom12 = twom24 * 4096.
write(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
write(6,'(5X,5I12)') isdext
do 200 i= 1, 24
seeds(i) = real(isdext(i))*twom24
200 continue
carry = 0.
if (isdext(25) < 0)  carry = twom24
isd = iabs(isdext(25))
i24 = mod(isd,100)
isd = isd/100
j24 = mod(isd,100)
isd = isd/100
in24 = mod(isd,100)
isd = isd/100
luxlev = isd
  if (luxlev <= maxlev) then
    nskip = ndskip(luxlev)
!         WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
!    +                         LUXLEV
  else  if (luxlev >= 24) then
    nskip = luxlev - 24
!         WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
  else
    nskip = ndskip(maxlev)
!         WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
    luxlev = maxlev
  endif
inseed = -1
return
!
!                    Entry to ouput seeds as integers
entry rluxut(isdext)
do 300 i= 1, 24
   isdext(i) = int(seeds(i)*twop12*twop12)
300 continue
isdext(25) = i24 + 100*j24 + 10000*in24 + 1000000*luxlev
if (carry > 0.)  isdext(25) = -isdext(25)
return
!
!                    Entry to output the "convenient" restart point
entry rluxat(lout,inout,k1,k2)
lout = luxlev
inout = inseed
k1 = kount
k2 = mkount
return
!
!                    Entry to initialize from one or three integers
entry rluxgo(lux,ins,k1,k2)
   if (lux < 0) then
      luxlev = lxdflt
   else if (lux <= maxlev) then
      luxlev = lux
   else if (lux < 24 .or. lux > 2000) then
      luxlev = maxlev
      write (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',lux
   else
      luxlev = lux
      do 310 ilx= 0, maxlev
        if (lux == ndskip(ilx)+24)  luxlev = ilx
310       continue
   endif
if (luxlev <= maxlev)  then
   nskip = ndskip(luxlev)
!        WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
!    +        LUXLEV,'     P=', NSKIP+24
else
    nskip = luxlev - 24
!         WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
endif
in24 = 0
if (ins < 0)  write (6,'(A)') &
   ' Illegal initialization by RLUXGO, negative input seed'
if (ins > 0)  then
  jseed = ins
!       WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
!    +      JSEED, K1,K2
else
  jseed = jsdflt
!       WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
endif
inseed = jseed
notyet = .false.
twom24 = 1.
   do 325 i= 1, 24
     twom24 = twom24 * 0.5
   k = jseed/53668
   jseed = 40014*(jseed-k*53668) -k*12211
   if (jseed < 0)  jseed = jseed+icons
   iseeds(i) = mod(jseed,itwo24)
325    continue
twom12 = twom24 * 4096.
   do 350 i= 1,24
   seeds(i) = real(iseeds(i))*twom24
   next(i) = i-1
350    continue
next(1) = 24
i24 = 24
j24 = 10
carry = 0.
if (seeds(24) == 0.) carry = twom24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury > 0).
kount = k1
mkount = k2
if (k1+k2 /= 0)  then
  do 500 iouter= 1, k2+1
    inner = igiga
    if (iouter == k2+1)  inner = k1
    do 450 isk= 1, inner
      uni = seeds(j24) - seeds(i24) - carry
      if (uni < 0.)  then
         uni = uni + 1.0
         carry = twom24
      else
         carry = 0.
      endif
      seeds(i24) = uni
      i24 = next(i24)
      j24 = next(j24)
450     continue
500   continue
!         Get the right value of IN24 by direct calculation
  in24 = mod(kount, nskip+24)
  if (mkount > 0)  then
     izip = mod(igiga, nskip+24)
     izip2 = mkount*izip + in24
     in24 = mod(izip2, nskip+24)
  endif
!       Now IN24 had better be between zero and 23 inclusive
  if (in24 > 23) then
     write (6,'(A/A,3I11,A,I5)') &
    '  Error in RESTARTING with RLUXGO:','  The values', ins, &
     k1, k2, ' cannot occur at luxury level', luxlev
     in24 = 0
  endif
endif
return
end
subroutine ptnew(xgnew , ipoly , vnew  , pntint, lptint, n     , &
                 poly  , lpoly , nind  , base  , delta , lb    , &
                 ub    , a     , lda   , nclin , ncnln , pivval, &
                 pivthr, xchthr, wrk   , lwrk  , iwrk  , liwrk , &
                 fail  , info  , minmax, x0)                   

!
!  *****************************************************************************       
!  THIS SUBROUTINE TRIES TO FIND A NEW POINT THAT CAN BE INCLUDED IN 
!  INTERPOLATION TO IMPROVE GEOMETRY. IF THE INTERPOLATION SET IS NOT
!  COMPLETE, THEN WE LOOK FOR A 'GOOD' POINT TO ADD TO THE SET; IF
!  THE INTERPOLATION SET IS COMPLETE, THEN WE CHOOSE A 'BAD' POINT THAT
!  WE WOULD LIKE TO REPLACE AND TRY TO FIND A 'GOOD' POINT TO REPLACE IT.
!  IF WE FIND SUCH A POINT THEN FAIL=FASLE AND THE NEW POINT IS XGNEW.
!  IF WE DO NOT FIND IT THEN FAIL=TRUE.
!
!  PARAMETERS 
!    
!  POLY   (INPUT)  THE ARRAY OF NEWTON FUNDAMENTAL POLYNOMIALS
!
!  PNTINT (INPUT)  THE ARRAY OF INTERPOLATION POINTS
!
!  NIND   (INPUT)  CARDINALITY OF THE INTERPOLATION SET
!
!  LB     (INPUT)  LOWER BOUNDS OF THE PROBLEM
!
!  LU     (INPUT)  UPPER BOUNDS OF THE PROBLEM
!
!  DELTA  (INPUT)  TRUST REGION RADIUS
!
!  A      (INPUT)  MATRIX OF LINEAR CONSTRAINTS OF THE PROBLEM
!
!  PIVVAL (INPUT)  ARRAY OF PIVOT VALUES ASSOCIATED WITH THE
!                  INTERPOLATION POINTS                
!  PIVTHR (INPUT)  THRESHOLD FOR ACCEPTABLE PIVOT VALUE
!                
!  XCHTHR (INPUT)  THRESHOLD FOR THE MINIMUM ACCEPTABLE IMPROVEMENT 
!                  IN THE PIVOT VALUE DUE TO REPLACING A POINT
!  X0     (INPUT)  CURRENT SHIFT OF THE POINTS FROM THE ORIGINAL POSITION
!
!  XGNEW  (OUTPUT) THE NEW POINT FOUND IN ORDER TO IMPROVE GEOMETRY
!
!  VNEW   (OUTPUT) THE PIVOT VALUE CORRESPONDING TO XGNEW, IN CASE IT
!                  CAN BE ADDED               
!  IPOLY  (OUTPUT) THE INDEX WHICH SHOULD BE ASSIGNED TO XGNEW, WHEN
!                  IT IS INCLUDED IN THE INTERPOLATION SET.
!                  (IF XGNEW IS ADDED, THEN IPOLY=NIND+1, OTHERWISE
!                   IPOLY IS THE INDEX OF THE POINT WHICH IS TO BE
!                   REPLACED BY XGNEW 
!  FAIL   (OUTPUT) LOGICAL VARIABLE INDICATING IF THE SEARCH SUCCEEDED
!                  FAIL=TRUE INDICATES THAT WE COULD NEITHER FIND ANY
!                  POINT TO BE ADDED TO THE SET, NOR WE COULD FIND A
!                  POINT WHICH WOULD REPLACE ANOTHER POINT IN THE SET
!                  WITH SIGNIFICANT IMPROVEMENT IN GEOMETRY.
!                  (NOTICE: WE DID NOT FIND SUCH POINTS, DOES NOT MEAN
!                           THAT THEY DO NOT EXISTS)
!  MINMAX (INPUT)  INDICATOR THAT ENABLES USER TO MINIMIZE THE NEXT
!                  POLYNOMIAL, IF ON THE PREVIOUS CALL IT WAS MAXIMIZED
!                  0 IF WE SHOULD MAXIMIZE OF ABSOLUTE VALUE OF PIVOT
!                 -1 IF WE SHOULD MINIMIZE THE REAL VALUE OF PIVOT
!                  1 IF WE SHOULD MAXIMIZE THE REAL VALUE OF PIVOT
!         (OUTPUT) 1 IF THE INPUT VALUE WAS 0 AND ANSWER ACHIEVED
!                    BY MAXIMIZATION
!                 -1 IF THE INPUT VALUE WAS 0 AND ANSWER ACHIEVED
!                    BY MINIMIZATION
!                  0 IF THE INPUT VALUE WAS 1 OR -1
!  ********************************************************************
!






integer          nind  , n, lpoly   , lptint, lwrk, ipoly, lda  , &
                 info  , iwrk(liwrk), liwrk , base, nclin, ncnln, &
                 minmax

double precision poly(lpoly   )   , xgnew(n) , pntint(lptint), &
                 pivval(nind+1)   , wrk(lwrk), delta   , &
                 lb(n+nclin+ncnln), a(lda*n) , vnew    , &
                 ub(n+nclin+ncnln), pivthr   , xchthr  , x0(n)
 
logical          fail
!  
!  COMMON VARIABLES
!
include 'dfo_model_inc.inc'

!
!  PRINTOUT PARAMETERS
!

double precision mcheps, cnstol
integer          iout  , iprint
common / dfocm / iout  , iprint, mcheps, cnstol
save / dfocm /
!
!  EXTERNAL SUBROUTINES
!


external :: mintr


!
!  APPLICATIONS: MINTR, GETNP, NEXTNP
!



!  
!  LOCAL VARIABLES
!

double precision  vmax, kappa, mval  , pivmin, kmod
integer           lenw, icurw, ixbase, minpiv, dd, &
                  np1 , i    , j     , ig    , ih

!
!  CONSTANTS  
!
double precision  zero, one, half
parameter        (zero = 0.0d0, one = 1.0d0, half = 0.5d0)



!
!  PARTITION THE REAL SPACE
!

ig     = 1
ih     = ig+n
ixbase = ih+n*n
icurw  = ixbase+n
lenw   = lwrk-icurw+1

!
!  CHECK IF REAL SPACE IS SUFFICIENT
!
if ( lenw < 1 ) then
  if ( iprint >= 0 ) write (iout, 1000) - lenw + 1
  stop
endif


np1 = n + 1
dd  =(np1)*(n+2)/2

!
!  NO POINT WAS FOUND YET, THEREFORE FAIL IS TRUE
!
fail =.true.
vmax = zero
!
!  IF THE INTERPOLATION SET IS INCOMPLETE THEN WE LOOK FOR
!  A POINT TO BE ADDED. WE DO IT BY MAXIMIZING THE NEXT PIVOT.
!  THE NEXT PIVOT IS THE VALUE AT A CERTAIN POINT OF THE 'NEXT' 
!  (NIND+1ST) POLYNOMIAL, UPDATED SO, THAT IS IS ZERO AT ALL  
!  POINT IN THE SET. 

if ( nind < dd ) then

!
!  UPDATE THE 'NEXT POLYNOMIAL, SO THAT IT IS ZERO AT ALL POINTS
!  OF THE INTERPOLATION SET
!
 
  call nextnp(nind+1, poly, pntint, nind, n, lpoly, lptint)
!
!  PUT THE COEFFICIENTS OF THE POLYNOMIAL IN THE FORM OF QUADRATIC FORM 
!
    
  call getnp(nind+1, poly, lpoly, n, kappa, wrk(ig),wrk(ih))

  if ( minmax /= 1 ) then
!
!  SET PARAMETERS FOR TRUST-REGION MINIMIZATION
!
    call dcopy(n, pntint((base-1)*n+1), 1, wrk(ixbase), 1)
    kmod = kappa
    do 20 i=1,n
      gmod(i)= wrk(ig+i-1)
      kmod   = kmod - gmod(i)*x0(i)
      do 10 j=1,n
        hmod(i,j) = wrk(ih+(i-1)*n+j-1)
        gmod(i)   = gmod(i) - hmod(i,j)*x0(j)
        kmod      = kmod + half*x0(i)*hmod(i,j)*x0(j)
10       continue
20     continue

!
!  MINIMIZE THE 'NEXT' POLYNOMIAL OVER THE TRUST-REGION
!       
    call unshft( n, x0, wrk(ixbase) )
    call mintr( n           , wrk(ixbase), mval , delta, lb   , &
                ub          , a          , lda  , nclin, ncnln, &
                wrk( icurw ), lenw       , iwrk , liwrk, info , 1)
    call shift( n, x0, wrk(ixbase) )

!
!  STORE THE VALUE AND THE POINT
!
    if ( info == 0 ) then
      vnew = mval + kmod
      vmax = abs(vnew)
      call dcopy(n, wrk(ixbase), 1, xgnew, 1)
    endif
  endif
!
!  SET PARAMETERS FOR MAXIMIZATION OVER THE TRUST-REGION (MINIMIZATION
!  WITH NEGATIVE SIGN)
!

  if ( minmax /= -1 ) then
    kmod = -kappa
    do 40 i=1,n
      gmod(i)= -wrk(ig+i-1)
      kmod   = kmod - gmod(i)*x0(i)
      do 30 j=1,n
        hmod(i,j) = -wrk(ih+(i-1)*n+j-1)
        gmod(i)   = gmod(i) - hmod(i,j)*x0(j)
        kmod      = kmod + half*x0(i)*hmod(i,j)*x0(j)
30       continue
40     continue



    call dcopy(n, pntint((base-1)*n+1), 1, wrk(ixbase), 1)

!
!  MAXIMIZE THE 'NEXT' POLYNOMIAL OVER THE TRUST-REGION
! 
    call unshft( n, x0, wrk(ixbase) )
    call mintr( n   , wrk(ixbase), mval , delta,lb          ,ub  , &    ! 
                a   , lda        , nclin, ncnln,wrk( icurw ),lenw, &
                iwrk, liwrk      , info , 1 )
    call shift( n, x0, wrk(ixbase) )

  endif
!
!  CHOOSE THE LARGER ABSOLUTE VALUE BETWEEN THE MAXIMUM AND THE MINIMUM
!  AND PICK APPROPRIATE POINT
! 



  if ( minmax == 0 ) then
    if ( abs(mval+kmod) > vmax .and. info == 0 ) then
      vnew   = - mval - kmod
      minmax = 1
      vmax   = abs(vnew)
      call dcopy(n, wrk(ixbase), 1, xgnew, 1)
    else
      minmax = -1 
    endif
  else
    minmax = 0
    if (  abs(mval+kmod) > vmax .and. info == 0 ) then
      vnew   = - mval - kmod
      vmax   = abs(vnew)
      call dcopy(n, wrk(ixbase), 1, xgnew, 1)
    endif
  endif
!
!  IF THE PIVOT VALUE IS ACCEPTABLE, THEN WE ARE DONE
!

  if (vmax>pivthr) then
    fail  = .false.
    ipoly =  nind+1
    info  =  0
    if ( iprint >= 3 ) write(iout,8000) ipoly, vmax
  endif
endif

!
!  IF WE DID NOT MANAGE TO FIND A POINT TO ADD (BECAUSE THE
!  INTERPOLATION SET WAS FULL OR PIVOT VALUE TOO SMALL), THEN
!  WE TRY TO FIND A POINT TO INCLUDE BY REPLACING SOME OTHER POINT
!

if (fail) then
!
!  FIRST CHOOSE A POINT WHICH WE WANT TO REPLACE. IT WILL BE THE POINT
!  WITH THE SMALLEST ASSOCIATED PIVOT VALUE
!
  pivmin=one - cnstol*delta
  minpiv=0 
  do 45 i=1,nind
    if (abs(pivval(i))<abs(pivmin)) then
      pivmin=abs(pivval(i))
      minpiv=i
    endif
45   continue

!
!  IF THE THERE SMALLEST PIVOT IS REASONABLY SMALL AND THE CHOSEN POINT
!  IS NOT THE BASE, THEN WE SET IPOLY EQUAL TO MINPIV - THE INDEX OF
!  THE CANDIDATE FOR REPLACEMENT.
!



  if ( minpiv>0   .and. minpiv/=base .and. &
     ( minpiv>np1 .or.  nind<=np1)       ) then

    ipoly=minpiv
      

!
!  MAXIMIZE THE ABSOLUTE VALUE OF THE NEWTON POLYNOMIAL WITH  INDEX IPOLY
!
    call getnp(ipoly, poly, lpoly, n, kappa, wrk(ig), wrk(ih))

    call dcopy(n, pntint((base-1)*n+1), 1, wrk(ixbase), 1)


    kmod = kappa
    do 60 i=1,n
      gmod(i)= wrk(ig+i-1)
      kmod   = kmod - gmod(i)*x0(i)
      do 50 j=1,n
        hmod(i,j) = wrk(ih+(i-1)*n+j-1)
        gmod(i)   = gmod(i) - hmod(i,j)*x0(j)
        kmod      = kmod + half*x0(i)*hmod(i,j)*x0(j)
50       continue
60     continue

 
!
!  DO TRUST-REGION MINIMIZATION
!
    call unshft( n, x0, wrk(ixbase) )
    call mintr( n   , wrk(ixbase), mval, delta, lb   , &
                ub  , a          , lda , nclin, ncnln, &
                wrk( icurw )     , lenw, iwrk , liwrk, info, 1 )
    call shift( n, x0, wrk(ixbase) )

    if ( info == 0 ) then
      vnew = mval + kmod  
      vmax = abs(vnew)
      call dcopy( n, wrk(ixbase), 1, xgnew, 1 )
    endif

    kmod = -kappa
    do 80 i=1,n
      gmod(i)= -wrk(ig+i-1)
      kmod   = kmod - gmod(i)*x0(i)
      do 70 j=1,n
        hmod(i,j) = -wrk(ih+(i-1)*n+j-1)
        gmod(i)   = gmod(i) - hmod(i,j)*x0(j)
        kmod      = kmod + half*x0(i)*hmod(i,j)*x0(j)
70       continue
80     continue

     

    call dcopy(n, pntint((base-1)*n+1), 1, wrk(ixbase), 1)

!
!  DO TRUST-REGION MAXIMIZATION
! 
    call unshft( n, x0, wrk(ixbase) )
    call mintr( n   , wrk(ixbase), mval, delta, lb   , &
                ub  , a          , lda , nclin, ncnln, &
                wrk( icurw )     , lenw, iwrk , liwrk, info, 1 )
    call shift( n, x0, wrk(ixbase) )


!
!  CHOOSE THE BETTER POINT BETWEEN MAXIMIZER AND MINIMIZER
!
    if ( abs( mval+kmod ) > vmax .and. info == 0) then
      vnew = - mval - kmod  
      vmax = abs(vnew)
      call dcopy( n, wrk(ixbase), 1, xgnew, 1)
    endif


!
!  CHECK IF THE NEW PIVOT GIVES AT LEAST  'XCHTHR' TIMES IMPROVEMENT
!  OVER THE OLD PIVOT VALUE. IF IT DOES, WE ACCEPT THE POINT BY SETTING
!  FAIL=FALSE
!
    
    if ( vmax>xchthr ) then
      fail =.false.
      info = 0
      if ( iprint >= 3 ) write(iout,8020) ipoly, pivmin, vmax 
    else 
      if ( iprint >= 3 ) write(iout,8030) ipoly, pivmin, vmax
    endif
  else
    if ( iprint >= 3 ) write(iout,8010) 
  endif
endif
return

1000 format(' PTNEW:  *** ERROR: LWRK TOO SMALL!' / &
       '            IT SHOULD BE AT LEAST ',i12 ,/)
8000 format(' PTNEW: A new point is found by maximizing the pivot,',/ &
       '       polynomial: ', i4,' pivot value: ', d14.7,/ ) 
8010 format(' PTNEW: No point which can be  replaced was found',/)
8020 format(' PTNEW: A new point replaced the ',i4,'-th point',/, &
       '       old pivot: ', d14.7, ' new pivot: ', d14.7,/)
8030 format(' PTNEW: A new point DID NOT replace the ',i4,'-th point',/ &
       '       old pivot: ', d14.7, ' new pivot: ', d14.7 ,/)
end
  







 subroutine rzrmat(a,n,m)
 double precision a(n,m)
 integer          m,n,i,j

 do 10 i=1,n
   do 20 j=1,m
     a(i,j)=0.0d0
20    continue
10  continue
 return
 end


 subroutine rzrvec(b,n)
 double precision b(n)
 integer          n,i

 do 10 i=1,n        
     b(i)=0.0d0
10  continue
 return
 end



 subroutine izrmat(a,n,m)
  
 integer        a(n,m),m,n,i,j

 do 10 i=1,n
   do 20 j=1,m
     a(i,j)=0
20    continue
10  continue
 return
 end

 subroutine izrvec(b,n)

 integer           b(n),n,i

 do 10 i=1,n        
     b(i)=0
10  continue
 return
 end

end module dfo_module