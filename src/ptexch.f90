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

intrinsic        abs

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










