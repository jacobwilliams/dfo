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


 if (i.eq.1) then
!
!  IF THE I-TH POLYNOMIAL IS A  CONSTANT
!
   val = poly(1)
 else if (i.le.np1) then
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

 if (i.eq.1) then
!
!  IF THE I-TH POLYNOMIAL IS A  CONSTANT
!
   val = poly(1)
 else if (i.le.np1) then
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
     

