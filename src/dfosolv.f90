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

double precision mvalue, dnrmnf 
external         mvalue, dnrmnf
intrinsic        min   , max   , abs
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
if ( lrs .lt.  max( rsneed, mdblrs ) ) &
   then
   if ( iprint .ge. 0 ) write( iout, 2000 ) &
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
if ( lis .lt. max(mntris + npslis, mdblis) ) then
   if ( iprint .ge. 0 ) write( iout, 2020 ) &
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


if ( iprint .ge. 2 ) write( iout, 8000 ) 
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
if (nind.lt.2) then
  if ( iprint .ge. 0 ) write(iout, 2060)
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

   if ( iprint .ge. 2 ) write( iout, 8010 ) iter
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
     if (( i .le. n+1-neqcon .or. i .gt. n+1 ) .and. &
     ( wrk( lv + iwrk( lin2sp + i )) .le. fbase - ten * mcheps)) &
          then
       fbase = wrk(lv + iwrk( lin2sp + i ))
       base  = i
     endif
10    continue
   if ( iprint .ge. 2 ) write( iout, 8020 ) base, nind

!
!  IF THE BASE POINT HAD CHANGED, WE NEED TO RECOMPUTE THE DISTANCES
!  FROM BASE TO ALL SAMPLE POINTS IN "POINTS" 
!
   if ( obase .ne. base ) then 
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
   if ( odelta .ne. delta ) then
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
   linear=( nind .gt. n - neqcon ) 
!
!  IF THERE ARE EQUALITY CONSTRAINTS, THEN UPON REACHING THE MAXIMUM
!  POSSIBLE NUMBER OF POINTS FOR LINEAR INTERPOLATION WE COMPLETE
!  THE LINEAR INTERPOLATION TO N+1 ELEMENTS WITH DUMMY POLYNOMIALS
!  AND POINTS. THIS IS DONE BECAUSE THE PROGRAM IS WRITTEN UNDER
!  ASSUMPTION THAT QUADRATIC BLOCK IS NOT CONSTRAUCTED UNTILL THE
!  LINEAR BLOCK IS FULL. BY COMPLETETING THE LINEAR BLOCK BY DUMMY
!  ELEMENTS WE MAKE SURE THAT WE CAN MOVE ON TO QUADRATIC ELEMENTS.  
!
   
   if ( neqcon .gt. 0 .and. nind .eq. n + 1 - neqcon ) then
     call complt( wrk(inp), wrk(ivi), wrk(ipiv), n, nind, lpoly )
   endif

!  --------------------------------------
!  PRINT OUT SUMMARY OF CURRENT ITERATION
!  --------------------------------------
!
!  NINDR - NUMBER OF ELEMENTS IN THE INTERPOLATION (WITHOUT THE DUMMY)  
!
   
   if ( iprint .ge. 1 ) then
     nindr = nind
     if ( nind .gt. n + 1 - neqcon ) nindr = nind - neqcon
!           IF (ITER.EQ.( ITER/10 )*10 +1 ) WRITE( IOUT, 1000 )
     write( iout, 1000 )
     write( iout, 1010) iter, nf, fbase, rho, delta, impr, nindr, &
                        neqcon
   endif
!
!  STOPPING TEST: IF THE TRUST REGION RADIUS IS SMALL ENOUGH, WE STOP
!                 IF THE BEST VALUE HAS REACHED THE LOWER BOUND, WE STOP
!
   if ( iprint .ge. 2 ) write( iout, 8030 ) delta 
   if ( ( howstop.eq.2 .and. stopcnt.ge.2) &
         .or. (delta .lt. delmin) ) then 
     exit = 0
     call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
                 wrk(ic), nq)
     it   = iter
     return               
   endif
!
!  IN DEBUGGING MODE, CHECK THE ACCURACY OF NEWTON FUNDAMENTAL POLYNOMIALS
!
 if ( iprint .ge. 3 ) then
   do 30 i = 2, nind
     if( i .le.  n + 1 - neqcon) then
       ii = min( n + 1 - neqcon, nind )
     else 
       ii = nind
     endif
     do 20 j = 1, ii
      if (j .le. n + 1 - neqcon .or. j .gt. n + 1 ) then
        call  evalnp( val, wrk(ipi), j, wrk(inp), i, n, &
                      lptint , lpoly )
        if ( (i.eq.j .and. abs(val-one) .gt. 1d-10) .or. &
             (i.ne.j .and. abs(val)     .gt. 1d-10) ) then
          write(iout,1020) i, j, val
        endif
      endif
20      continue
30    continue
 endif

 if ( method .le. 2 ) then
!  ---------------------------------------------------------------------
!  BUILD THE QUADRATIC MODEL OF THE OBJECTIVE FUNCTION
!  ---------------------------------------------------------------------

   if ( iprint .ge. 2 ) write( iout, 8040 ) 
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

   if ( iprint .ge. 2 ) write( iout, 8040 ) 
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
 if ( iprint .ge. 3 ) then
   do 41 i=1, m
     do 40 j = 1, nind
       if ( j .le. n + 1 - neqcon .or. j .gt. n + 1 ) then
         ii = (iwrk(lin2sp + j)-1)*m
         wrk(lci + j) = wrk(lc + ii + i) 
         val = wrk(lci + j) - (wrk(lcc+i) + &
               mvalue( n, wrk( ipi+(j-1)*n), &
               wrk(icl+(i-1)*n), wrk( icq+(i-1)*n*n ), &
               wrk( icurw ), lrs))
         if ( abs(val) .gt. 1.0d-6 ) write(iout,1030)  j, val 
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

 
 if ( method .le. 2 ) then
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
 if ( iprint .ge. 3 ) then
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

 if ( iprint .ge. 2 ) write( iout, 8065 ) 
 if ( method .le. 2 ) then
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

 if ( method .eq. 2 .and. inform.eq. 3 ) then 
   usemer=.true.
   call mintr (n         , wrk( ik ), mval        , del  , lb, &
               ub        , a        , lda         , nclin, ncnln , &
               wrk(icurw), lrs      , iwrk(icuriw), lis  , inform, &
               2         )
 else
   usemer=.false.
   if ( inform .eq. 3 ) inform = 0
 endif
!
!  IF NO FEASIBLE POINT FOUND, THEN CONSIDER REDUCTION TO BE ZERO
!

 if (inform .eq. 2 ) then
   val = fbase
   goto 81
 endif 

!
!  IF THE MINIMIZATION FAILED, THEN RECORD THE BEST POINT AND EXIT
!
 if (inform.ne.0) then
   call bestpt(n, m, x, fx, c, wrk(ip), wrk(iv), wrk(iob), &
               wrk(ic), nq)
   it    = iter
   exit=-7
   return
 endif
!
!  COMPUTE THE PREDICTED VALUE OF THE MERIT FUNCTION
!
 if ( method .eq. 1 .or. (method .eq. 2 .and. .not. usemer) ) then
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

 if ( prered .lt. max( noisef,delmin*1.0d-2) .or. &
      prered .lt. delmin*abs(fbase) .and. snorm .lt. delmin .or. &
      prered .lt. cnstol*abs(fbase)*1.0d-3 ) then
   if ( iprint .ge. 2 ) write( iout, 8067 ) prered
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
   if ( nf .gt. maxnf ) then
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
   if ( scale .ne. 0 ) call unscl( n, wrk(ik), scal )
   call fun( n, m, wrk( ik ), wrk(iob+nq), wrk(ic + nq*m), iferr )
   if ( scale .ne. 0 ) call scl( n, wrk(ik), scal )
!
!  IF FUNCTION COMPUTATION FAILS, THEN REDUCE THE RADIUS FOR 
!  MINIMIZATION TO HALF OF THE DISTANCE BETWEEN THE NEW POINT 
!  AND THE BASE. REPEAT THE MINIMIZATION
!
   if (iferr) then
     del  = snorm/2
     if ( del .ge. delmin ) then 
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
   if ( iprint .ge. 2 ) write( iout, 8080 ) rho

!  ------------------------------------------------------------
!  COMPUTE RELATIVE ACHIEVED REDUCTION AND CHECK IF IT IS SMALL
!  ------------------------------------------------------------
   val = ( fbase - fnq ) /(one+abs(fbase))
   if ( ( val .gt. 0 ).and. ( val .lt. whenstop).and. &
        ( rho .gt. rhomin )) then
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

   if ( rho .ge. rhojmp  .and. rho .le. two - rhojmp &
            .and. abs(snorm-delta).lt.cnstol ) then
     if ( iprint .ge. 2 ) write( iout, 8085 )
     theta = zero
     do 90 i = 1, nq - 1
!
!  FOR ALL SAMPLE POINTS (EXCEPT LAST) WHICH ARE NOT IN INTERPOLATION
!  AND ARE FURTHER THAT DELTA AWAY  FROM THE BASE COMPUTE THE AGREEMENT
!  BETWEEN THE MODEL AND FUNCTION
!
       if ( iwrk( lsp2in+ i ) .eq. 0   .and. &
            wrk(  ld + i ) .ge. delta ) then
         ibeg = lp + (i-1) * n 
!
!  COMPUTE THE PREDICTED VALUE OF THE MERIT FUNCTION
!
         mval= mvalue( n, wrk( ibeg+1 ), wrk( ig ), &
                       wrk( ih ), wrk( icurw), lrs )
         if ( method .le. 2  ) then
           call fmerit(m, val, mval+kappa, wrk(ici+(i-1)*m), &
                       ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1), &
                       pp, method )
         else
           val = mval + kappa
         endif
         val=val-fbase
         if ( abs(val).gt.mcheps ) then
           rhow = ( wrk( lv + i ) - fbase ) / val
!
!  IF FOR A GIVEN POINT THE AGREEMENT IF GOOD, THEN WE TRY TO JUMP 
!  AT LEAST AS FAR AS THIS POINT FROM THE BASE
!
           if ( rhow .ge. one - rhomin .and. &
                rhow .le. one + rhomin      ) then
             theta = max( theta, wrk( ld + i ) )
           endif 
         endif                  
       endif
90      continue

!
!  MAKE SURE THAT THE SIZE OF THE JUMP DOES NOT EXCEED THE LIMIT
!          
     theta = min( theta, maxjmp * delta )
     if ( iprint .ge. 2 ) write( iout, 8090 ) theta
       
     del = theta

!
!  IF THE POSSIBLE SIZE OF THE JUMP BIG ENOUGH (SO IT COVERS AREA LARGER
!  THAN THE AREA THAT WOULD BE COVERED BY THE NEXT ITERATION) THEN
!  COMPUTE THE JUMP
!
100      if ( del .gt. snorm + ratio * delta) then

!
!  SET INITIAL POINT TO THE BASE POINT AND SHIFT IT TO ORIGINAL POSITION
!
       call dcopy(n, wrk(ibase), 1, wrk( ik), 1)        
       call unshft(n, x, wrk(ik) )
!
!  CALL MINIMIZATION
!
       if ( method .le. 2 ) then 
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
       if ( inform .eq. 3 ) inform=0
       if ( inform .eq. 2 ) goto 115
!
!  COMPUTE THE PREDICTED VALUE OF THE MERIT FUNTION
!
       if ( method .le. 2 ) then
         call fmerit(m, val, mval+kmod, wrk(icurw), &
          ub(n+nclin+ncnln+1), lb(n+nclin+ncnln+1), pp, &
          method )
       else
         val=mval+kmod
       endif
!
!  IF MINIMIZATION HAD FAILED, THEN RECORD THE BEST POINT AND EXIT
!
       if (inform.ne.0) then
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
       if ( iprint .ge. 2 ) write( iout, 8092 ) stnorm
!
!  IF THE ACTUAL JUMP IS LARGE ENOUGH THEN SAMPLE THE FUNCTION VALUE
!  AT THE NEW POINT
!
       if ( stnorm .gt.  snorm + ratio * delta ) then
         call dcopy(n, wrk(ik), 1, wrk( inq + n ), 1)
         nf = nf + 1
!
!  IF THE NUMBER OF FUNCTION CALLS EXCEEDS MAXIMUM, THEN RECORD THE
!  BEST POINT AND EXIT
!
         if ( nf .ge. maxnf ) then
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
         if ( scale .ne. 0 ) call unscl( n, wrk(ik), scal )
         call fun( n, m, wrk( ik ), wrk(iob+nq), wrk(ic + nq*m), &
                   iferr )
         if ( scale .ne. 0 ) call scl( n, wrk(ik), scal )
!
!  IF THE FUNCTION COMPUTATION FAILS THEN REDUCE THE SIZE OF THE JUMP
!  AND TRY TO COMPUTE A NEW JUMP
!
         if ( iferr ) then 
           del=stnorm/2
           if ( del .ge. delta ) then 
             goto 100
           else
             goto 115
           endif
         endif

!
!  IF THE FUNCTION VALUE IS COMPUTED THEN RECORD IT AND DO APPROPRIATE
!  BOOKKEEPING AND UPDATES
!
         if ( iprint .ge. 2 ) write( iout, 8095 )
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
            
         if ( ft .lt.  fnq ) then
         
!
!  COMPUTE 'RHO' - ACHIEVED REDUCTION VS. PREDICTED REDUCTION
!
           rho   = ( ft - fbase ) / (val - fbase)
           next  = nq
           snorm = stnorm
           impr  = 20
           if ( iprint .ge. 2 ) write( iout, 8097 ) rho
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

   if ( rho .gt. rhomin ) then
     if ( iprint .ge. 2 ) write( iout, 8100 ) 
!
!  IF THE MODEL IS INCOMPLETE, THEN WE ADD THE NEXT POINT TO THE
!  MODEL ( IF EFFORT >= 4, THEN IT WILL BE DONE LATER )
!
     if ( nind .lt. dd .and. effort .lt. 4 ) then
!
!  DETERMINE THE INDEX OF LAST POLYNOMIAL IN THE CURRENT BLOCK
!
       if ( nind .le. n ) then
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
          
         if (abs(val) .gt. pivthr ) go to 130

120        continue
       go to 140
!
!  PLACE THE IPOLY-TH POLYNOMIAL IN THE PLACE OF NIND+1-ST
!  POLYNOMIAL, BY SWAPPING THEM IN ARRAY POLY
!
130        if ( ipoly .ne. nind+1 ) then
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
       if ( effort .eq. 4 ) then
         base              = nind 
         lbase = lpi + ( base - 1 ) * n
         ibase = lbase + 1
         obase = base
         call getdis( n, nq, 0 , wrk( ip ) , iwrk(lin2sp+base), &
                      wrk( id ), wrk(icurw), lrs )
       endif
       if ( iprint .ge. 2 ) write( iout, 8110 ) val
     endif


140      oldnq = nq          

!
!  IF THE MODEL IS FULL OR THE PIVOT VALUE FOR ADDING THE POINT
!  IS TOO SMALL, THEN TRY TO INCLUDE THE NEW POINT BY REPLACING
!  ANOTHER POINT BY IT. IF SUCCEED, UPDATE THE NEWTON POLYNOMIALS
!  ( IF EFFORT >= 4, THEN IT WILL BE DONE LATER )
!

     if ( ifadd.eq.0 ) then
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
         if ( effort .eq. 4 ) then
           base              = ipoly 
           lbase = lpi + ( base - 1 ) * n
           ibase = lbase + 1
           obase = base
           call getdis( n, nq, 0 , wrk( ip ) , iwrk(lin2sp+base), &
                        wrk( id ), wrk(icurw), lrs )
         endif
         if ( iprint .ge. 2 ) write( iout, 8120 ) ipoly,  val
       endif
     endif
!
!  IF THE MODEL REDUCTION IS NOT SUFFICIENTLY GOOD, WE STILL
!  TRY TO ADD THE NEW POINT (SINCE IT CONTAINS NEW INFORMATION)
!

   else
     if ( iprint .ge. 2 ) write( iout, 8105 )
!
!  IF THE MODEL IS INCOMPLETE TRY ADDING THE NEW POINT
!  IF THE EFFORT >= 3 TRY TO DO IT LATER
!

     if ( nind.lt.dd .and. effort .lt. 3 ) then
!
!  DETERMINE THE INDEX OF LAST POLYNOMIAL IN THE CURRENT BLOCK
!
       if ( nind .le. n ) then
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
         if (abs(val) .gt. addthr ) go to 160

150        continue
!
!  IF THERE IS NO POLYNOMIAL WHICH GIVES ACCEPTABLE PIVOT, THEN
!  MOVE ON TO MODEL IMPROVEMENT

       go to 170
!
!  PLACE THE IPOLY-TH POLYNOMIAL IN THE PLACE OF NIND+1-ST
!  POLYNOMIAL, BY SWAPPING THEM IN ARRAY POLY
!
160        if (ipoly .ne. nind+1) then
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
       if ( iprint .ge. 2 ) write( iout, 8110 ) val
     endif
 
!
!  END OF ATTEMPT OF INCLUDING A NEW POINT
!
   endif
 endif
170  continue

 if ( iprint .ge. 3 ) then
   write( iout, 8170 )
   do 175 i = 1, nind
     write( iout, 8055 ) wrk(ipiv + i - 1)
175    continue
 endif 
!
!  IF THE IMPROVEMENT OF THE MODEL IS REQUIRED, SET IMPR TO 0
! 
 if ( (ifadd .eq. 0) .or. (ifxg .eq. 1) ) impr = 0
!
!  IF THE NEW POINT HASN'T BEEN ADDED, BUT IT GIVES REDUCTION, THEN
!  SET A VALUE OF IMPR, SO THAT THIS POINT IS SET TO BASE IN 'IMPMOD' 
! 
 if ( ifadd.eq.0 .and. rho.gt. rhomin ) impr=-next
!
!  TRY TO IMPROVE THE MODEL (BY FINDING ANOTHER POINT OR  RECOMPUTING
!  THE BASIS OF NEWTON POLYNOMIALS AND THE INTERPOLATION SET IF:
!      1. NOTHING WAS ADDED TO THE POINT
!      2. IF A NEW 'GEOMETRY' POINT IS REQUIRED 
!      3. THE MODEL IS TOO OLD, THAT IS THE CURRENT BASE IS LAYER*DELTA
!         AWAY FROM THE FIRST INTERPOLATION POINT
!      4. IF THE EFFORT LEVEL REQUIRES INTERPOLATION TO BE RECOMPUTED      
!
 if (  impr .le. 0                            .or. &
      (wrk(ld+iwrk(iin2sp)) .gt. layer*delta) .or. &
      (effort .eq. 3 .and. rho .le. rhomin)   .or. &
       effort .eq. 4 ) then

!
!  CHECK IF THE NUMBER IF FUNCTION CALLS REACHED MAXIMUM
!  IF IT DOES RECORD THE BEST POINT AND EXIT
!
   if ( nf .ge. maxnf ) then
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

   if ( neqcon .gt. 0 .and. nind .eq. n + 1 - neqcon ) then
     call complt( wrk(inp), wrk(ivi), wrk(ipiv), n, nind, lpoly )
   endif              
   if ( iprint .ge. 2 ) write( iout, 8140 )
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
   if ( nf .ge. maxnf ) then
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

   if (nind.lt.2) then
     if ( iprint .ge. 0 ) write(iout, 2060)
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

   if ( oldnq .ne. nq ) then
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
 if ((( impr .ge. 4  .or. rho .eq. - thous) .and. linear) &
       .or. impr .eq. 0 ) okgeom=.true.


!  --------------------------------------------------------------
!
!  UPDATING TRUST REGION RADIUS
!
!  --------------------------------------------------------------


!
!  IF THE REDUCTION IF GOOD, INCREASE DELTA ACCORDING TO 
!  THE STEP TAKEN BY TRUST REGION  MINIMIZATION  
!
 if ( rho .ge. rhoinc .or. impr .eq. 20 ) then
   if ( iprint .ge. 2 ) write( iout, 8130 )
   delta = max( delta, min( delmax, ratio*snorm ))

!
!  IF THE REDUCTION WAS BAD, AND THE MODEL WAS ALREADY GOOD ENOUGH
!  THEN REDUCE DELTA
!
 elseif ( rho .lt. rhomin .and. okgeom ) then
   if ( iprint .ge. 2 ) write( iout, 8150 )
   delthr =  ratio  * delmin
   delta = max( delta / ratio,  delthr )
 endif

!
!  IF DELTA IS JUST ABOVE THE MINIMUM THRESHOLD, MAKE SURE THAT IT
!  DOES NOT GET DECREASED UNTIL WE ARE SURE THAT WE CANNOT GET
!  FURTHER IMPROVEMENT
!

 if ( delta .le. ratio * delmin + 10*mcheps .and. &
      rho  .lt. rhomin  ) then
   if ( iprint .ge. 2 ) write( iout, 8160 )
   nstop = nstop + 1
   if (  nstop .ge. 5 .and. okgeom ) then
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
              i4,1x, 'ERROR IS:', d12.6)
1030 format('DFO: INEXACT MODEL AT POINT:', 1x, i4,  1x,'ERROR IS:', &
             d12.6)
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
  if ( values(i) .lt. val - mcheps ) then
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

intrinsic        abs, max


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
if ( lrs .lt. 4*n ) then
  if ( iprint .ge. 0 ) write( iout, 1000 ) - lrs - 4*n + 1 
  stop
endif

if ( liwrk .lt. n ) then
  if ( iprint .ge. 0 ) write( iout, 1100 ) - liwrk - n + 1 
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
if ( ncnln .gt. 0 ) then 
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
  if ( ub(i)-lb(i) .lt. pivthr ) then
    do 10 j = 1, n
      m = m + 1
      if ( j .eq. i ) then
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
  if ( ub(i + n) - lb(i + n) .lt. pivthr ) then
    m = m + 1
    call dcopy(n, a( i ), lda, wrk(lmat + m), ldmat) 
  endif
30 continue

!
!  CHECK IF THE THERE ARE EQULITIES AMONG NONLINEAR CONSTRAINTS,
!  FOR EACH EQUALITY WRITE CORRESPONDING ROW OF THE JACOBIAN
!

do 40 i=1, ncnln
  if ( ub(i + n + nclin)-lb(i + n + nclin) .lt. pivthr ) then
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

if ( m .gt. n ) then
  if ( iprint .ge. 0 ) write( iout, 1200 ) m, n
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
  if ( abs(wrk( lmat + (i-1)*ldmat + i )) .lt. 100*mcheps ) &
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






