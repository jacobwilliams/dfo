!     ( Last modified on Tue May 28 14:07:02 CDT 2002 )
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
implicit none

integer           n ,  nclin , ncnln, liwrk, lwrk, iwrk(liwrk), &
                  lda, inform, method


double precision  x0(n), mval, delta, lwrbnd(n+nclin+ncnln), &
                  uprbnd(n+nclin+ncnln), wrk(lwrk), a(lda,n)

double precision ddot

external         ddot


! HERE WE ASSUME THAT IPOPT HAS BEEN COMPILED WITH DYNAMIC MEMORY ALLOCATION
integer lrw, liw
parameter( lrw=0, liw=0 )
integer iw
double precision rw

! USE THE FOLLOWING, IF IPOPT HAS NOT BEEN COMPILED TO USE DYNAMIC MEMORY
! ALLOCATION
!      INTEGER LRW, LIW
!      PARAMETER (LRW = 100000, LIW = 100000)
!      INTEGER IW( LIW )
!      DOUBLE PRECISION RW( LRW )

!
!     Parameters for IPOPT
!
include 'dfo_model_inc.inc'
!
!  PRINTOUT PARAMETERS
!
integer          iout  , iprint
double precision mcheps, cnstol

common / dfocm /  iout, iprint, mcheps, cnstol
save / dfocm /

integer dfo_nmax, dfo_mmax
parameter( dfo_mmax = nconmx + nnlnmx + nlinmx )
parameter( dfo_nmax = nvarmx + dfo_mmax )
double precision dfo_inf
parameter( dfo_inf = 1.d20 )

integer nipopt, mipopt, nlb, nub
double precision x( dfo_nmax )
integer ilb( dfo_nmax )
integer iub( dfo_nmax )
double precision bnds_l( dfo_nmax )
double precision bnds_u( dfo_nmax )
double precision v_l( dfo_nmax )
double precision v_u( dfo_nmax )
double precision lam( dfo_mmax )
double precision c( dfo_mmax )
integer idat(3)
double precision dat(1)
!
!     Algorithmic Parameters (INITPARAMS)
!
integer nargs
double precision args( 50 )
character*20 cargs( 50 )

integer iter
integer ierr

integer i, idummy, ncont, nnlnt, nlint, m
double precision dummy, val

external eval_f, eval_g, eval_c, eval_a, eval_h
external ev_hlv_dummy, ev_hov_dummy, ev_hcv_dummy

useipopt=1
if (method == 2) then 
   usemerit=.true.
else
   usemerit=.false.
endif
ncont=ncon
if ( ncnln <= nnln .and. .not.usemerit) ncon=0
do i = 1, n
  call dcopy(nclin, a(1,i), 1, amat(1,i), 1)
enddo
!
!     Consider bounds for original variables
!
mipopt = nclin + ncnln
nipopt = n + mipopt
if( nipopt>dfo_nmax ) then
   write(*,*) 'DFO_NMAX too small. Must be at least ',nipopt
   inform = 1
   goto 9999
endif
if( mipopt>dfo_mmax ) then
   write(*,*) 'DFO_MMAX too small. Must be at least ',mipopt
   inform = 1
   goto 9999
endif
do i = 1, n
   ilb(i)    = i
   bnds_l(i) = max(x0(i)-delta, lwrbnd(i))
   iub(i)    = i
   bnds_u(i) = min(x0(i)+delta, uprbnd(i))
enddo
nlb = n
nub = n
!
!     We create slacks for all constraints.  For the equality constraints,
!     IPOPT will then treat them as fixed variables (and this way it also
!     knows about the right hand side.)
!
do i = n+1, n+mipopt
   if( lwrbnd(i)>-dfo_inf ) then
      nlb         = nlb + 1
      ilb(nlb)    = i
      bnds_l(nlb) = lwrbnd(i)
   endif
   if( uprbnd(i)< dfo_inf ) then
      nub         = nub + 1
      iub(nub)    = i
      bnds_u(nub) = uprbnd(i)
   endif
enddo
!
!     Compute initial point for IPOPT
!
call dcopy(n, x0, 1, x, 1)
call funcon(1, ncnln, n, ncnln, 0, x0, x(n+1), dummy, idummy)
!
!     Initialize parameters
!
args(1) = -1.d0
cargs(1) = 'IPRINT'
args(2) = 1.d3
cargs(2) = 'DMU0'
args(3) = 1.d-1
cargs(3) = 'DBNDPUSH'
args(4) = 1.d-1
cargs(4) = 'DBNDFRAC'
args(5) = 1.d-2
cargs(5) = 'DRESTOREDFACT'
nargs = 5
!
!     The following data is communicated for the EVAL_* routines
!
idat(1) = n
idat(2) = nclin
idat(3) = ncnln

call ipopt(nipopt, x, mipopt, nlb, ilb, bnds_l, nub, iub, &
     bnds_u, v_l, v_u, lam, c, lrw, rw, liw, iw, iter, ierr, &
     eval_f, eval_c, eval_g, eval_a, eval_h, ev_hlv_dummy, &
     ev_hov_dummy, ev_hcv_dummy, dat, idat, nargs, args, cargs)
call dcopy(n, x, 1, x0, 1)
 if (ierr==0) then
   if( iprint>=3 )    write(iout,8000) 
   inform=0
 else if (ierr==1) then
   if( iprint>=3 )    write(iout,8004) 
   inform=3
 else if (ierr>=97 .or. ( ierr >= 4 .and. ierr <= 8) &
          .or. ierr == 15 ) then
   if( iprint>=3 )   write(iout,8009) 
   inform=1
 else
   if( iprint>=3 )    write(iout,8001) 
   inform=3        
 endif

 if ( inform == 3 ) then
   inform = 0
   do 40 i=1, n
     if ( x0(i) < bnds_l(i) - cnstol .or. &
          x0(i) > bnds_u(i) + cnstol     ) &
         inform = 2
40    continue
   if ( nclin > 0 .and. inform == 0 ) then
     do 70 i=1, nclin 
       val = ddot(n, a(i, 1), lda, x0, 1 )
       if ( val < lwrbnd(n+i) - cnstol .or. &
            val > uprbnd(n+i) + cnstol     ) &
            inform = 2
70      continue
   endif
   if ( ncnln > 0 .and. inform == 0 ) then
     call funcon(1, ncnln+nclin, n, ncnln+nclin, iwrk, x, &
                 wrk, wrk(ncnln+nclin+1), 1)
     do 60 i=1, nnln 
       if ( wrk(nlin+i) < lwrbnd(n+nlin+i) - cnstol .or. &
            wrk(nlin+i) > uprbnd(n+nlin+i) + cnstol     ) &
            inform = 2
60      continue
     if ( .not. usemerit ) then 
       do 80 i=1, ncon 
         if ( wrk(nlin+nnln+i) < lwrbnd(n+nlin+nnln+i) - cnstol &    ! 
              .or. wrk(nlin+nnln+i) > uprbnd(n+nlin+nnln+i) &
              + cnstol     ) &
              inform = 3
80        continue
     endif
   endif
 endif
call funobj(1, n, x, mval, wrk, 1)
if ( ncon > 0 .and. .not. usemerit) then
   if ( inform == 0 ) then
     call dcopy(ncon, x(n+nlin+nnln+1), 1, wrk, 1)
   else
     nlint=nlin
     nnlnt=nnln
     nnln=0
     nlin=0
     call funcon(1, ncon, n, ncon, iwrk, x, wrk, wrk(ncon+1), 1)
     nnln=nnlnt
     nlin=nlint
   endif
endif
9999 continue
ncon=ncont
return
8000  format( ' DFO: MINTR: SUCCESSFUL MINIMIZATION' )
8001  format( ' DFO: MINTR: NO FURTHER IMPROVEMENT CAN BE OBTAINED' )
8004  format( ' DFO: MINTR: MAXIMUM NUMBER OF ITER. REACHED IN IPOPT' )
8009  format( ' DFO: MINTR: AN INPUT PARAMETER TO IPOPT IS INVALID' )
end



! Copyright (C) 2002, Carnegie Mellon University and others.
! All Rights Reserved.
! This code is published under the Common Public License.
!*******************************************************************************
!
subroutine eval_a(task, n, x, nz, a, acon, avar, dat, idat)
!
!*******************************************************************************
!
!    $Id: mintr_ipopt.f 20 2004-04-20 18:41:53Z andreasw $
!
!-------------------------------------------------------------------------------
!                                 Title
!-------------------------------------------------------------------------------
!
!T    Compute Jacobian of constraints for DFO problem
!
!-------------------------------------------------------------------------------
!                          Programm description
!-------------------------------------------------------------------------------
!
!B    
!
!-------------------------------------------------------------------------------
!                             Author, date
!-------------------------------------------------------------------------------
!
!A    Andreas Waechter      03/27/03
!
!-------------------------------------------------------------------------------
!                             Documentation
!-------------------------------------------------------------------------------
!
!D
!
!-------------------------------------------------------------------------------
!                             Parameter list    
!-------------------------------------------------------------------------------
!
!    Name     I/O   Type   Meaning
!
!P   TASK      I    INT     =0: Obtain NZ
!P                         <>0: Compute Jacobian
!P   N         I    INT    number of variables in problem statement
!P                            (including slacks for inequality constraints)
!P   X         I    DP     point where A is to be evaluated
!P   NZ       I/O   INT    TASK = 0: O: number of nonzero elements
!P                         otherwise: number of nonzero elements
!P                                     (size of A, AVAR, ACON)
!P   A         O    DP     (only TASK<>0) values in Jacobian
!P   ACON      O    INT    (only TASK<>0) row indices
!P   AVAR      O    INT    (only TASK<>0) column indices
!P   DAT       I    DP     ignored
!P   IDAT      I    INT    IDAT(1) = N
!P                         IDAT(2) = NCLIN
!P                         IDAT(3) = NCNLN
!
!-------------------------------------------------------------------------------
!                             local variables
!-------------------------------------------------------------------------------
!
!L
!
!-------------------------------------------------------------------------------
!                             used subroutines
!-------------------------------------------------------------------------------
!
!S    FUNCON
!
!*******************************************************************************
!
!                              Declarations
!
!*******************************************************************************
!
implicit none
!
!*******************************************************************************
!
!                              Include files
!
!*******************************************************************************
!
!
!
!-------------------------------------------------------------------------------
!                             Parameter list
!-------------------------------------------------------------------------------
!
integer task
integer n
double precision x(n)
integer nz
double precision a(nz)
integer acon(nz)
integer avar(nz)
double precision dat(*)
integer idat(3)
!
!-------------------------------------------------------------------------------
!                            Local varibales
!-------------------------------------------------------------------------------
!
double precision dummy
integer idummy, i, j
!
!*******************************************************************************
!
!                           Executable Statements
!
!*******************************************************************************
!
if( task==0 ) then
!
!     We assume dense Jacobian for original variables, and have the slacks
!     for all constraints
!
   nz = (idat(1)+1)*(idat(2)+idat(3))

else
!
!     Call FUNCON to get Jacobian (for nonslack variables)
!
   call funcon(2, idat(3), idat(1), idat(2)+idat(3), 0, x, dummy, &
        a, idummy)
   nz = 0
   do i = 1, idat(1)
      do j = 1, idat(2)+idat(3)
         nz = nz + 1
         acon(nz) = j
         avar(nz) = i
      enddo
   enddo
!
!     Now add entries for the slacks
!
   do i = 1, idat(2)+idat(3)
      acon(nz+i) = i
      avar(nz+i) = idat(1) + i
   enddo
   call dcopy(idat(2)+idat(3), -1.d0, 0, a(nz+1), 1)
   nz = nz + idat(2)+idat(3)

endif

return
end
! Copyright (C) 2002, Carnegie Mellon University and others.
! All Rights Reserved.
! This code is published under the Common Public License.
!*******************************************************************************
!
subroutine eval_c(n, x, m, c, dat, idat)
!
!*******************************************************************************
!
!    $Id: mintr_ipopt.f 20 2004-04-20 18:41:53Z andreasw $
!
!-------------------------------------------------------------------------------
!                                 Title
!-------------------------------------------------------------------------------
!
!T    Compute values of constraints for DFO problem
!
!-------------------------------------------------------------------------------
!                          Programm description
!-------------------------------------------------------------------------------
!
!B    
!
!-------------------------------------------------------------------------------
!                             Author, date
!-------------------------------------------------------------------------------
!
!A    Andreas Waechter      03/27/03
!
!-------------------------------------------------------------------------------
!                             Documentation
!-------------------------------------------------------------------------------
!
!D
!
!-------------------------------------------------------------------------------
!                             Parameter list    
!-------------------------------------------------------------------------------
!
!    Name     I/O   Type   Meaning
!
!P   N         I    INT    number of variables in problem statement
!P                            (including slacks for inequality constraints)
!P   X         I    DP     point where G is to be evaluated
!P   M         I    INT    number of constraints
!P   C         O    DP     values of constraints
!P   DAT       I    DP     ignored
!P   IDAT      I    INT    IDAT(1) = N
!P                         IDAT(2) = NCLIN
!P                         IDAT(3) = NCNLN
!
!-------------------------------------------------------------------------------
!                             local variables
!-------------------------------------------------------------------------------
!
!L
!
!-------------------------------------------------------------------------------
!                             used subroutines
!-------------------------------------------------------------------------------
!
!S    FUNOBJ
!
!*******************************************************************************
!
!                              Declarations
!
!*******************************************************************************
!
implicit none
!
!*******************************************************************************
!
!                              Include files
!
!*******************************************************************************
!
!
!
!-------------------------------------------------------------------------------
!                             Parameter list
!-------------------------------------------------------------------------------
!                        
integer n
double precision x(n)
integer m
double precision c(m)
double precision dat(*)
integer idat(3)
!
!-------------------------------------------------------------------------------
!                            Local varibales
!-------------------------------------------------------------------------------
!
double precision dummy
integer idummy
!
!*******************************************************************************
!
!                           Executable Statements
!
!*******************************************************************************
!

!
!     First get constraint values
!
call funcon(1, idat(3), idat(1), idat(2)+idat(3), 0, x, c, &
     dummy, idummy)
!
!     Now substract values of slacks
!
call daxpy(idat(2)+idat(3), -1.d0, x(idat(1)+1), 1, c, 1)

return
end
! Copyright (C) 2002, Carnegie Mellon University and others.
! All Rights Reserved.
! This code is published under the Common Public License.
!*******************************************************************************
!
subroutine eval_f(n, x, f, dat, idat)
!
!*******************************************************************************
!
!    $Id: mintr_ipopt.f 20 2004-04-20 18:41:53Z andreasw $
!
!-------------------------------------------------------------------------------
!                                 Title
!-------------------------------------------------------------------------------
!
!T    Compute objective function value for DFO problem
!
!-------------------------------------------------------------------------------
!                          Programm description
!-------------------------------------------------------------------------------
!
!B    
!
!-------------------------------------------------------------------------------
!                             Author, date
!-------------------------------------------------------------------------------
!
!A    Andreas Waechter      03/27/03
!
!-------------------------------------------------------------------------------
!                             Documentation
!-------------------------------------------------------------------------------
!
!D
!
!-------------------------------------------------------------------------------
!                             Parameter list    
!-------------------------------------------------------------------------------
!
!    Name     I/O   Type   Meaning
!
!P   N         I    INT    number of variables in problem statement
!P                            (including slacks for inequality constraints)
!P   X         I    DP     point where F is to be evaluated
!P   F         O    DP     objective function value
!P   DAT       I    DP     ignored
!P   IDAT      I    INT    IDAT(1) = N
!P                         IDAT(2) = NCLIN
!P                         IDAT(3) = NCNLN
!
!-------------------------------------------------------------------------------
!                             local variables
!-------------------------------------------------------------------------------
!
!L
!
!-------------------------------------------------------------------------------
!                             used subroutines
!-------------------------------------------------------------------------------
!
!S    FUNOBJ
!
!*******************************************************************************
!
!                              Declarations
!
!*******************************************************************************
!
implicit none
!
!*******************************************************************************
!
!                              Include files
!
!*******************************************************************************
!
!
!
!-------------------------------------------------------------------------------
!                             Parameter list
!-------------------------------------------------------------------------------
!                        
integer n
double precision x(n)
double precision f
double precision dat(*)
integer idat(3)
!
!-------------------------------------------------------------------------------
!                            Local varibales
!-------------------------------------------------------------------------------
!
double precision dummy
integer idummy
!
!*******************************************************************************
!
!                           Executable Statements
!
!*******************************************************************************
!

call funobj(1, idat(1), x, f, dummy, idummy)

9999 continue
return
end
! Copyright (C) 2002, Carnegie Mellon University and others.
! All Rights Reserved.
! This code is published under the Common Public License.
!*******************************************************************************
!
subroutine eval_g(n, x, g, dat, idat)
!
!*******************************************************************************
!
!    $Id: mintr_ipopt.f 20 2004-04-20 18:41:53Z andreasw $
!
!-------------------------------------------------------------------------------
!                                 Title
!-------------------------------------------------------------------------------
!
!T    Compute gradient of objective function for DFO problem
!
!-------------------------------------------------------------------------------
!                          Programm description
!-------------------------------------------------------------------------------
!
!B    
!
!-------------------------------------------------------------------------------
!                             Author, date
!-------------------------------------------------------------------------------
!
!A    Andreas Waechter      03/27/03
!
!-------------------------------------------------------------------------------
!                             Documentation
!-------------------------------------------------------------------------------
!
!D
!
!-------------------------------------------------------------------------------
!                             Parameter list    
!-------------------------------------------------------------------------------
!
!    Name     I/O   Type   Meaning
!
!P   N         I    INT    number of variables in problem statement
!P                            (including slacks for inequality constraints)
!P   X         I    DP     point where G is to be evaluated
!P   G         O    DP     gradient of objective function
!P   DAT       I    DP     ignored
!P   IDAT      I    INT    IDAT(1) = N
!P                         IDAT(2) = NCLIN
!P                         IDAT(3) = NCNLN
!
!-------------------------------------------------------------------------------
!                             local variables
!-------------------------------------------------------------------------------
!
!L
!
!-------------------------------------------------------------------------------
!                             used subroutines
!-------------------------------------------------------------------------------
!
!S    FUNOBJ
!
!*******************************************************************************
!
!                              Declarations
!
!*******************************************************************************
!
implicit none
!
!*******************************************************************************
!
!                              Include files
!
!*******************************************************************************
!
!
!
!-------------------------------------------------------------------------------
!                             Parameter list
!-------------------------------------------------------------------------------
!                        
integer n
double precision x(n)
double precision g(*)
double precision dat(*)
integer idat(3)
!
!-------------------------------------------------------------------------------
!                            Local varibales
!-------------------------------------------------------------------------------
!
double precision f
integer idummy
!
!*******************************************************************************
!
!                           Executable Statements
!
!*******************************************************************************
!
call funobj(2, idat(1), x, f, g, idummy)
!
!     Add entries for slack variables
!
call dcopy(idat(2)+idat(3), 0.d0, 0, g(idat(1)+1), 1)

9999 continue
return
end
! Copyright (C) 2002, Carnegie Mellon University and others.
! All Rights Reserved.
! This code is published under the Common Public License.
!*******************************************************************************
!

subroutine eval_h(task, n, x, m, lam, nnzh, hess, irnh, icnh, &
     dat, idat)
!
!*******************************************************************************
!
!    $Id: mintr_ipopt.f 20 2004-04-20 18:41:53Z andreasw $
!
!-------------------------------------------------------------------------------
!                                 Title
!-------------------------------------------------------------------------------
!
!T    Compute Hessian of Lagrangian for DFO problem
!
!-------------------------------------------------------------------------------
!                          Programm description
!-------------------------------------------------------------------------------
!
!B    
!
!-------------------------------------------------------------------------------
!                             Author, date
!-------------------------------------------------------------------------------
!
!A    Andreas Waechter      03/27/03
!
!-------------------------------------------------------------------------------
!                             Documentation
!-------------------------------------------------------------------------------
!
!D
!
!-------------------------------------------------------------------------------
!                             Parameter list    
!-------------------------------------------------------------------------------
!
!    Name     I/O   Type   Meaning
!
!P   TASK      I    INT     =0: Obtain NNZH
!P                         <>0: Compute Hessian
!P   N         I    INT    number of variables
!P   X         I    DP     point at which constraints are to be evaluated
!P   M         I    INT    number of equality constraints
!P   LAM       I    DP     vector of Lagrangian multipliers
!P   NNZH     I/O   INT    size of HESS, IRNH, ICNH
!P   HESS      O    DP     nonzero elements of Hessian
!P   IRNH      O    INT    row indices of nonzero elements
!P   ICNH      O    INT    column indices of nonzero elements
!P                         (for each i=1,..,NNZH the nonzero element HESS(i)
!P                         is the element of Hessian in row IRNH(i) and column
!P                         ICNH(i) as well as in row ICNH(i) and column IRNH(i).
!P                         For non-diagonal elements provide only one of them.)
!P   DAT       I    DP     ignored
!P   IDAT      I    INT    IDAT(1) = N
!P                         IDAT(2) = NCLIN
!P                         IDAT(3) = NCNLN
!
!-------------------------------------------------------------------------------
!                             local variables
!-------------------------------------------------------------------------------
!
!L
!
!-------------------------------------------------------------------------------
!                             used subroutines
!-------------------------------------------------------------------------------
!
!S    FUNCON
!
!*******************************************************************************
!
!                              Declarations
!
!*******************************************************************************
!
implicit none
!
!*******************************************************************************
!
!                              Include files
!
!*******************************************************************************
!
include 'dfo_model_inc.inc'
!
!
!-------------------------------------------------------------------------------
!                             Parameter list
!-------------------------------------------------------------------------------
!
integer task
integer n
integer m 
integer nnzh
double precision lam(m)
double precision x(n)
double precision hess(nnzh)
integer irnh(nnzh)
integer icnh(nnzh)
double precision dat(*)
integer idat(3)
!
!-------------------------------------------------------------------------------
!                            Local varibales
!-------------------------------------------------------------------------------
!
integer idummy, i, j
double precision h(nvarmx, nvarmx)
!
!*******************************************************************************
!
!                           Executable Statements
!
!*******************************************************************************
!
if( task==0 ) then
!
!     Get number of nonzeros in Hessian of the Lagrangian.  Assume dense
!     Hessian (for original variables and lower triangular part only)
!
   nnzh = (idat(1)*(idat(1)+1))/2
else
   if( idat(1)> nvarmx) then
      write(*,*) &
           'In eval_h: NVARMX too small.  Must be at least', &
           idat(1)
      stop
   endif
!
!     Call FUNCON to get values for Hessian
!
   call funcon(3, idat(3), idat(1), nvarmx, 0, x, lam, &
        h, idummy)
!
!     Now copy out the lower triangular part
!
   nnzh = 0
   do i = 1, idat(1)
      do j = 1, i
         nnzh = nnzh + 1
         hess(nnzh) = h(i,j)
         irnh(nnzh) = i
         icnh(nnzh) = j
      enddo
   enddo
         
endif

9999 continue
return
end

subroutine ev_hlv_dummy()
end
subroutine ev_hov_dummy()
end
subroutine ev_hcv_dummy()
end
