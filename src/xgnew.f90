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


external          mintr, funobj, funcon


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
if ( lenw .lt. 1 ) then
  if ( iprint .ge. 0 ) write (iout, 1000) - lenw + 1
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

if ( nind .lt. dd ) then

!
!  UPDATE THE 'NEXT POLYNOMIAL, SO THAT IT IS ZERO AT ALL POINTS
!  OF THE INTERPOLATION SET
!
 
  call nextnp(nind+1, poly, pntint, nind, n, lpoly, lptint)
!
!  PUT THE COEFFICIENTS OF THE POLYNOMIAL IN THE FORM OF QUADRATIC FORM 
!
    
  call getnp(nind+1, poly, lpoly, n, kappa, wrk(ig),wrk(ih))

  if ( minmax .ne. 1 ) then
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
    if ( info .eq. 0 ) then
      vnew = mval + kmod
      vmax = abs(vnew)
      call dcopy(n, wrk(ixbase), 1, xgnew, 1)
    endif
  endif
!
!  SET PARAMETERS FOR MAXIMIZATION OVER THE TRUST-REGION (MINIMIZATION
!  WITH NEGATIVE SIGN)
!

  if ( minmax .ne. -1 ) then
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



  if ( minmax .eq. 0 ) then
    if ( abs(mval+kmod) .gt. vmax .and. info .eq. 0 ) then
      vnew   = - mval - kmod
      minmax = 1
      vmax   = abs(vnew)
      call dcopy(n, wrk(ixbase), 1, xgnew, 1)
    else
      minmax = -1 
    endif
  else
    minmax = 0
    if (  abs(mval+kmod) .gt. vmax .and. info .eq. 0 ) then
      vnew   = - mval - kmod
      vmax   = abs(vnew)
      call dcopy(n, wrk(ixbase), 1, xgnew, 1)
    endif
  endif
!
!  IF THE PIVOT VALUE IS ACCEPTABLE, THEN WE ARE DONE
!

  if (vmax.gt.pivthr) then
    fail  = .false.
    ipoly =  nind+1
    info  =  0
    if ( iprint .ge. 3 ) write(iout,8000) ipoly, vmax
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
    if (abs(pivval(i)).lt.abs(pivmin)) then
      pivmin=abs(pivval(i))
      minpiv=i
    endif
45   continue

!
!  IF THE THERE SMALLEST PIVOT IS REASONABLY SMALL AND THE CHOSEN POINT
!  IS NOT THE BASE, THEN WE SET IPOLY EQUAL TO MINPIV - THE INDEX OF
!  THE CANDIDATE FOR REPLACEMENT.
!



  if ( minpiv.gt.0   .and. minpiv.ne.base .and. &
     ( minpiv.gt.np1 .or.  nind.le.np1)       ) then

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

    if ( info .eq. 0 ) then
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
    if ( abs( mval+kmod ) .gt. vmax .and. info .eq. 0) then
      vnew = - mval - kmod  
      vmax = abs(vnew)
      call dcopy( n, wrk(ixbase), 1, xgnew, 1)
    endif


!
!  CHECK IF THE NEW PIVOT GIVES AT LEAST  'XCHTHR' TIMES IMPROVEMENT
!  OVER THE OLD PIVOT VALUE. IF IT DOES, WE ACCEPT THE POINT BY SETTING
!  FAIL=FALSE
!
    
    if ( vmax.gt.xchthr ) then
      fail =.false.
      info = 0
      if ( iprint .ge. 3 ) write(iout,8020) ipoly, pivmin, vmax 
    else 
      if ( iprint .ge. 3 ) write(iout,8030) ipoly, pivmin, vmax
    endif
  else
    if ( iprint .ge. 3 ) write(iout,8010) 
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
  







