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
 intrinsic        abs
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





