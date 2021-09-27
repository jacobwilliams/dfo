
subroutine mintr( n   , x0   , mval  , delta, lwrbnd, uprbnd, &
                  a   , lda  , nclin , ncnln, wrk   , lwrk  , &
                  iwrk, liwrk, inform, method )
!
!  *******************************************************************
!  THIS SUBROUTINE FINDS THE MINIMUM OF THE QUADRATIC MODEL WITHIN THE
!  GIVEN REGION. THE REGION IS DEFINED BY THE INTERSECTION OF THE
!  TRUST REGION OF RADIUS DELTA AND THE ANALYTICALLY DEFINED FEASIBLE 
!  SET OF THE ORIGINAL PROBLEM.
!
!                         T       T
!                MIN [GMOD X+0.5*X HMOD X]
!         
!                  X0-DELTA <=  X   <= XO+DELTA
!        S.T.
!                               / X  \
!                        LB <= ( AX   ) <= UB
!                               \C(X)/
!  PARAMETERS:
!
!   N      (INPUT)  DIMENTION OF THE PROBLEM
!
!   X0     (INPUT)  ARRAY OF LENGTH N CONTAINING THE CENTER OF THE TRUST
!                   REGION
!          (OUTPUT) CONTAINS THE OPTIMAL POINT FOR THE MODEL
!
!   MVAL   (OUTPUT) THE VALUE OF THE MODEL AT THE OPTIMAL POINT
!
!   DELTA  (INPUT)  TRUST REGION RADIUS
!
!   LWRBND (INPUT)  ARRAY OF LENGHT N+NCLIN+NCNLN OF LOWER BOUNDS 
!
!   UPRBND (INPUT)     ''       ''         ''        UPPER   ''
!
!   NCLIN  (INPUT)  NUMBER OF LINEAR ANALYTIC CONSTRAINTS
!
!   A      (INPUT)  (LDA X N) MATRIX OF LINEAR ANALYTIC CONSTRAINTS
!  
!   NCNLN  (INPUT)  NUMBER OF NOLINEAR INEQUALITIES (DIFFICULT AND EASY)
!
!   WRK             REAL SPACE WORKING ARRAY
!
!   IWRK            INTEGER SPACE WORKING ARRAY
!
!   INFORM (OUTPUT) INFORMATION ON EXIT
!              0    SUCCESSFUL MINIMIZATION
!              1    THE DERIVATIVES OF THE CONSTRAINT OR SOME PARAMETER
!                   SET BY THE USER IS INCORRECT
!              2    MINIMIZATION FAILED FOR SOME REASON
!
!  
!   METHOD (INPUT)  METHOD FOR HANDLING CONSTRAINTS
!              1    MINIMIZE MODEL OF OBJECTIVE S.T. MODELS OF CONSTRAINTS
!              2    MINIMIZE MERIT FUNCTION OF THE MODELS OF CON AND OBJ
!             3,4   MINIMIZE MODEL OF A MERIT FUNCTION (EXACT OR QUAD)
!
!  **********************************************************************
!


integer           n ,  nclin , ncnln, liwrk, lwrk, iwrk(liwrk), &
                  lda, inform, method


double precision  x0(n), mval, delta, lwrbnd(n+nclin+ncnln), &
                  uprbnd(n+nclin+ncnln), wrk(lwrk), a(lda*n) 

!
!  COMMON VARIABLES
!

!
!  PRINTOUT PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol

common / dfocm /  iout, iprint, mcheps, cnstol
save / dfocm /
!
!  EXTERNAL SUBROUTINES
!

external          funobj, funcon

double precision ddot

external         ddot

!
!  LOCAL VARIABLES
!

double precision  val, tol, half, small, zero, one

parameter         (half = 0.5d0, small  = 1.0d-5, zero=0.0d0, &
                   one  = 1.0d0)
integer           i   , nlnd, lnneed, ibl   , ibu   , lbl  , &
                  lbu , ic  , icjac , iclmda, icuriw, icurw, &
                  ir  , inf , igrad , iistat, leniw , lenw , &
                  iter, ldcjac, ncont
intrinsic         max , min
include 'dfo_model_inc.inc'
!
!  SET THE COMMON PARAMETER 'IPOPT' (DEFINED IN DFO_MODEL_INC) TO 0
!  SINCE WE ARE USING NPSOL
!
useipopt=0
if (method == 2) then 
   usemerit=.true.
else
   usemerit=.false.
endif
ncont=ncon
if (ncnln <= nnln .and. .not.usemerit) ncon=0
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!    APPLICATION     :       FUNCON, FUNOBJ
!    NPSOL           :       NPSOL, NPOPTN 
!    BLAS            :       DCOPY, DDOT
!    FORTRAN SUPPLIED:       MIN, MAX
!


! 
!  PARTITION THE REAL SPACE
!

nlnd   = max(1,ncnln)
ibl    = 1
ibu    = ibl+n+nclin+ncnln
ic     = ibu+n+nclin+ncnln
icjac  = ic+nlnd
iclmda = icjac+nlnd*n
igrad  = iclmda+n+nclin+ncnln
ir     = igrad+n
icurw  = ir+n*n
lenw   = lwrk-icurw+1

!
!  CHECK IF THE REAL SPACE IS SUFFICIENT
!

lnneed = 2*n*n+n*nclin+2*n*ncnln+20*n+11*nclin+21*ncnln

if (lenw<lnneed) then
  if (iprint >= 0) write(iout,9000) lnneed
  stop
endif

! 
!  PARTITION THE INTEGER SPACE
!
iistat =1
icuriw  =iistat+n+nclin+ncnln
leniw=liwrk-icuriw+1

!
!  CHECK IF THE INTEGER SPACE IS SUFFICIENT
!
if (leniw<3*n+nclin+2*ncnln) then
  if (iprint >= 0) write (iout,9001) 3*n+nclin+2*ncnln
  stop
endif

!
!  SET THE JACOBIAN DIMENSION FOR NPSOL
!

ldcjac=nlnd
!
!  COMBINE TRUST REGION BOUNDS WITH THE SIMPLE BOUNDS OF THE PROBLEM
!
lbl=ibl-1
lbu=ibu-1
do 10 i=1,n
  wrk(lbl+i)=max(x0(i)-delta, lwrbnd(i))
  wrk(lbu+i)=min(x0(i)+delta, uprbnd(i))
10 continue

if (nclin+ncnln > 0) then
  call dcopy(nclin+ncnln, lwrbnd(n+1), 1, wrk(ibl+n), 1)
  call dcopy(nclin+ncnln, uprbnd(n+1), 1, wrk(ibu+n), 1)
endif
!
!  SET CONSTRAINT TOLERANCE
!

tol=cnstol
!
!  IF THERE ARE NO NONLINEAR CONSTRAINTS SET THE CONSTRAINT 
!  TOLERANCE, USED BY NPSOL, TO A DEFAULT VALUE
!
call npoptn( 'NOLIST' )
if ( ncnln == 0 ) tol = 1.0d-8

if ( tol >= 1.0d0 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D0' )
else if ( tol >= 1.0d-1 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-1' )
else if ( tol >= 1.0d-2 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-2' )
else if ( tol >= 1.0d-3 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-3' )
else if ( tol >= 1.0d-4 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-4' )
else if ( tol >= 1.0d-5 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-5' )
else if ( tol >= 1.0d-6 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-6' )
else if ( tol >= 1.0d-7 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-7' )
else if ( tol >= 1.0d-8 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-8' )
else if ( tol >= 1.0d-9 ) then 
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-9' )
else     
  call npoptn( 'FEASIBILITY TOLERANCE = 1.0D-10' )
endif
call npoptn( 'PRINT LEVEL = 0' )
   call npsol( n          , nclin       , ncnln     , lda       , &
               ldcjac     , n           , a         , wrk(ibl)  , &
               wrk(ibu)   , funcon      , funobj    , inf       , &
               iter       , iwrk(iistat), wrk(ic)   , wrk(icjac), &
               wrk(iclmda), mval        , wrk(igrad), wrk(ir)   , &
               x0         , iwrk(icuriw), leniw     , wrk(icurw), &
               lenw       )        


 if (inf==0) then
   if( iprint>=3 )    write(iout,8000) 
   inform=0
 else if (inf==1) then
   if( iprint>=3 )    write(iout,8001) 
   inform=0
 else if (inf==2) then
   if( iprint>=3 )    write(iout,8002) 
   inform=2
 else if (inf==3) then
   if( iprint>=3 )    write(iout,8003) 
   inform=3
 else if (inf==4) then
   if( iprint>=3 )    write(iout,8004) 
   inform=0
 else if (inf==6) then
   if( iprint>=3 )    write(iout,8006) 
   inform=2
 else if (inf==7) then
   if( iprint>=3 )   write(iout,8007) 
   inform=1
 else if (inf==9) then
   if( iprint>=3 )   write(iout,8009) 
   inform=1
 endif

 if ( inf == 6 .or. inf == 4  .or. inf == 1 ) then
   inform = 0
   do 40 i=1, n
     if ( x0(i) < wrk(lbl+i) - cnstol .or. &
          x0(i) > wrk(lbu+i) + cnstol     ) &
         inform = 2
40    continue
   if ( nclin > 0 .and. inform == 0 ) then
     do 70 i=1, nclin 
       val = ddot(n, a(i), lda, x0, 1 )
       if ( val < lwrbnd(n+i) - cnstol .or. &
            val > uprbnd(n+i) + cnstol     ) &
            inform = 2
70      continue
   endif
   if ( ncnln > 0 .and. inform == 0 ) then
     call funcon(1 , ncnln  , n         , ldcjac, iwrk(icuriw), &
                 x0, wrk(ic), wrk(icjac), 1     )
     do 60 i=1, ncnln 
       if ( wrk(ic+i-1) < lwrbnd(n+nclin+i) - cnstol .or. &
            wrk(ic+i-1) > uprbnd(n+nclin+i) + cnstol     ) &
            inform = 3
60      continue
   endif
 endif
 if ( ncon > 0 .and. .not. usemerit ) then
   call funcon(1 , ncnln  , n         , ldcjac, iwrk(icuriw), &
               x0, wrk(ic), wrk(icjac), 1     )
   call dcopy(ncon, wrk(ic+nnln), 1, wrk, 1)
 endif
 ncon=ncont

return


9000  format( ' DFO: MINTR: *** ERROR: LWRK TOO SMALL!' / &
        '  DFO:           IT SHOULD BE AT LEAST ',i12 )
9001  format( ' DFO: MINTR: *** ERROR: LIWRK TOO SMALL!' / &
        '  DFO:           IT SHOULD BE AT LEAST ', i12 )
8000  format( ' DFO: MINTR: SUCCESSFUL MINIMIZATION' )
8001  format( ' DFO: MINTR: NO FURTHER IMPROVEMENT CAN BE OBTAINED' )
8002  format( ' DFO: MINTR: NO FEASIBLE POINT FOUND FOR LINEAR' / &
        '  DFO:       CONSTRAINTS AND BOUNDS' )
8003  format( ' DFO: MINTR: NO FEASIBLE POINT FOUND FOR ' / &
        '  DFO:        NONLINEAR CONSTRAINTS' )
8004  format( ' DFO: MINTR: MAXIMUM NUMBER OF ITERATIONS REACHED' )
8006  format( ' DFO: MINTR: NPSOL QUIT BEFORE X SATISFIED K-T CONDIT.')
8007  format( ' DFO: MINTR: DERIVATIVES OF CONSTRAINTS ARE INCORRECT' )
8009  format( ' DFO: MINTR: AN INPUT PARAMETER IS INVALID' )

 end









