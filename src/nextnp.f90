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

 intrinsic        max
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



