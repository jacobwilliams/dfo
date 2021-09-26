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

 intrinsic        abs
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









