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

 if (nind.gt.np1) then
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
   if (nind.gt.1) then 

!
!  UPDATE USING LINEAR BLOCK
!     
     k=min(n-neqcon,nindm1)
     do 50 j = 1,k
       b(i) = b(i) + lafla(j+1)*poly(i + 2 + np1*(j-1))
50      continue
   endif

   if (nind .gt.np1) then
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

 if (nind.gt.np1) then

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









