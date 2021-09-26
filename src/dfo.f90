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
            wrk(icurw), lenw   , iwrk   , liwrk , inform )


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
  call bestpt( n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
             wrk(ic), nq )
  goto 48
endif
  call dfoslv( n  , m  , nclin, ncnln     , x     , fx   , &
               c  , nq , base , lp        , lv    , ld   , &
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
6200 format( a10, d12.6, / )
6300 format( a10, 3( d12.6), /) 
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


