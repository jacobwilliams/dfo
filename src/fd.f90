subroutine fd( lafla , pntint, values, poly, n,  q, neqcon, &
               lptint, lpoly )


!    
!  ***********************************************************************
!  THIS SUBROUTINE COMPUTED FINITE DIFFERENCES FOR MULTIVARIATE
!  INTERPOLATION, USING 'Q' INTERPOLATION POINTS AND 'Q' NEWTON
!  FUNDAMENTAL POLYNOMIALS.
!
!  PARAMETERS:
!  
!   N      (INPUT)  DIMENTION OF THE PROBLEM
!
!   Q      (INPUT)  NUMBER OF INTERPOLATION POINTS
! 
!   PNTINT (INPUT)  LIST OF  'NIND' DATA POINTS. THE I-TH POINT OCCUPIES
!                   POSITIONS ( I - 1 ) * N + 1 TO I * N.
!   VALUES (INPUT)  VALUES OF THE OBJECTIVE FUNCTION AT THE 'NIND'
!                   DATA POINTS CONTAINED IN POINTS.
!   POLY   (INPUT)  THE ARRAY CONTAINING COEFFICIENTS OF NEWTON FUNDAMENTAL
!                   POLYNOMIALS (AS COMPUTED BY 'NBUILD')
!   NEQCON (INPUT)  NUMBER OF LINEARLY INDEP. EQUALITY CONSTRAINTS
!
!   LAFLA  (OUTPUT) THE ARRAY OF THE  FINITE DIFFERENCES
! 
!  **************************************************************************
!

 integer          n, q, neqcon, lptint, lpoly 

 double precision lafla(q), pntint(lptint), values(q), poly(lpoly)

 intrinsic        min
!
!  LOCAL VARIABLES
! 
 integer          dd, i, k, l,  np1, np2, qp1
 double precision val

 np1 = n + 1
 np2 = n + 2
 dd = (np1)*(np2)/2.0
 qp1 = q + 1

!
!  INITIALISE THE LAFLA, THE FINITE DIFFERENCE OPERATORS ON F
!

 do 10 i =1,q
   lafla( i ) = values(i)
10  continue
!
!  START THE MAJOR ITERATION
!

 do 100 i=1,q
!
!  FOR ALL MULTIINDICES ALFA (WHOSE CARDINALITY IS GREATER THAN OR EQUAL
!  TO I), UPDATE THE ASSOCIATED FUNCTION L_ALFA USING
!  THE NEWTON FUNDAMENTAL POLYNOMIALS OF DEGREE EQUAL TO THE 
!  CARDINALITY OF ALFA_{I-1}
!

   if (i==2) then
!          
!  ALL FINITE DIFFERENCES CORRESPONDING TO DEGREES >= 1 ARE UPDATED.
!  THE NEWTON FUNDAMENTAL POLYNOMIALS OF DEGREE CARDINALITY OF ALFA_I = 0
!  ARE USED IN THE UPDATE. THE UPDATES ARE SKIPPED FOR DUMMY POLYNOMIALS
!  WHOSE INDEX IS ALWAYS HIGHER THAN NP1-NEQCON
!
     do 110 l = 2, min(q, np1 - neqcon) 
       call  evalnp( val, pntint, l, poly, 1, n, lptint, &
                     lpoly )
       lafla( l ) = lafla( l ) - lafla( 1 )*val
110      continue
     do 115 l = np2, q 
       call  evalnp( val, pntint, l, poly, 1, n, lptint, &
                     lpoly )
       lafla( l ) = lafla( l ) - lafla( 1 )*val
115     continue

   else
     if (i==np2) then
!
!  ALL FINITE DIFFERENCES CORRESPONDING TO DEGREES = 2 ARE UPDATED.
!  THE NEWTON FUNDAMENTAL POLYNOMIALS OF DEGREE CARDINALITY OF ALFA_I = 1
!  ARE USED IN THE UPDATE
! 
       do 120 l = np2,q
         do 130 k = 2, np1-neqcon
           call  evalnp( val, pntint, l, poly, k, n, lptint, &
                         lpoly )
           lafla( l ) = lafla( l ) - lafla( k )*val
130          continue
120        continue

     endif
   endif
100  continue


 return
 end













