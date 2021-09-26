
program dfotest

!  K. Scheinberg
!  2003.
!
integer          nmax   , maxfun, lwrk  , liwrk , ddmax , nclmax, &
                 ncnmax , lmax  , rsmax0, rsmax1, rsmax2, rsmax3, &
                 rsmax4 , ismax0, ismax1, ismax2, nfpmax, nxmax , &
                 mmax
parameter      ( nmax   =  30 ,  nclmax = 10, ncnmax=10, &
                 maxfun = 4000,  nxmax  = 10, mmax  = 15)
parameter      ( lmax   = nmax+nclmax+ncnmax+mmax )

parameter      ( ddmax  = (nmax+1)*(nmax+2)/2)

parameter      ( nfpmax = ddmax*(nmax+1)*nmax/2+(nmax+1)*nmax+1 )

parameter      ( rsmax0 = maxfun*(nmax+2) + ddmax*(nmax+2) &
                          + 2*nmax + nmax**2  + nfpmax)
parameter      ( rsmax1 = 3*(lmax) + (nmax+1)*ncnmax &
                          + nmax*nmax + 2*nmax + 1 )
parameter      ( rsmax2 = 2*nmax*nmax + nmax*nclmax +2*nmax*ncnmax &
                          + 20*nmax + 11*nclmax + 21*ncnmax )
parameter      ( rsmax3 = ddmax*ddmax + (ddmax-1)*nmax + &
                          ddmax-1 )

parameter      ( rsmax4 = 3*nmax+nmax**2 )

parameter      ( ismax0 = maxfun + ddmax )
parameter      ( ismax1 = 3*nmax + nclmax+2*ncnmax )
parameter      ( ismax2 = lmax )
parameter      ( lwrk   = rsmax0 + rsmax1 +  rsmax2 +  rsmax4+ &
                          rsmax3  )
parameter      ( liwrk  = ismax0 + ismax1 +  ismax2 + ddmax )


integer          n    , it , nf, exit, input, &
                 nclin        , ncnln, lda, nx, ldx , iout , &
                 majit        , itest, ln, maxit, maxnf, scale, &
                 iprint, stpcrtr          

      
double precision x( nmax*nxmax ), fx(nxmax)        , &
                 lb(lmax)       , ub(lmax)   , alin(nclmax*nmax), &
                 x0( nmax*nxmax ), v(mmax)   , conx(nxmax*mmax), &
                 pp, delmin, delta, stpthr, cnstol
logical          linear(mmax), equatn(mmax), efirst, lfirst, &
                 nvfrst   

parameter      ( iout = 6 )

character * 256    vname( nmax ), pname, cname(mmax+nclmax+ncnmax)

character * 256   prbdat
parameter      ( input = 55  )

integer ntest
common  /test/    ntest
save /test/
character*80 value

integer          m
double precision c(mmax), ldcj

!
!  SET THE TEST PROBLEM NUMBER
!     
call getarg(1, value)
read(value,*) ntest
!
!  SET UP THE DATA STRUCTURES AND SOME DATA FOR DFO
!
call setup ( iout  , n, x, nmax, nclmax, ncnmax, nclin, &
             ncnln, lb , ub, &
              alin, lda, m, mmax, &
              pname, vname,  cname )

ldcj=max(1, ncnln)
ln= n+nclin+ncnln


!
!  SET THE NUMBER OF POINTS PROVIDED AND THEIR DIMENSION
!      
call dcopy(n, x, 1, x0, 1)
nx = 1
ldx= n
maxit=1000
maxnf=1000
stpcrtr=1
delmin=0.001d0
stpthr=0.001d0
cnstol=0.00001d0
delta=1d0
pp=1000d0
scale=0
iprint=1
  do 5555 itest=1,1
    call dcopy(n, x0, 1, x, 1)
     call dfo(n    , nx  , x    , ldx  , fx, conx, .false., m, &
              c    , nclin, ncnln, lb, ub  , &
              alin , lda , vname, pname, cname, it, nf  , exit, &
              maxit,  maxnf, stpcrtr, delmin, stpthr, &
              cnstol, delta, pp, scale,  iout , iprint)
     
5555   continue  
stop
end

!  ***************************************************
!
!  SUBROUTINE COMPUTING INITIAL DATA
!
!  ***************************************************


subroutine setup( iout  , n    , x     , nmax, nclmax, ncnmax   , &
                  nclin , ncnln, bl    , bu  , alin  , lda, m, &
                  mmax, pname, xnames,  gnames )
!
!  TEST PROBLEMs FOR DFO:
!  SETUP NUMBER OF VARIABLES, A STARTING POINT AND BOUNDS FOR 
!  PROBLEMS FROM HOCK AND SCHITTKOWSKI
!
!  NTEST = 1 --  HOCK AND SCHITTKOWSKI # 16
!  
!  NTEST = 2 --  HOCK AND SCHITTKOWSKI # 113
!
!  NTEST = 3 --  HOCK AND SCHITTKOWSKI # 112
!
!  NTEST = 4 --  ROSENBROCK BANANA FUNCTION
!
!  NTEST = 5 --  HOCK AND SCHITTKOWSKI # 111
!
!  NTEST = 6 --  HOCK AND SCHITTKOWSKI # 26
!
!  NTEST = 7 --  HOCK AND SCHITTKOWSKI # 50
!
!  NTEST = 8 --  HOCK AND SCHITTKOWSKI # 112
!
integer            iout  , n     , nmax, nclin, ncnln, &
                   nclmax, ncnmax, mmax, m, lda, j

character * 256     xnames( nmax ), pname, &
                   gnames( mmax )

double precision   x( nmax ), bl( nmax+nclmax+ncnmax+mmax ), &
                   bu(nmax+nclmax+ncnmax+mmax), &
                   alin(nclmax, nmax)
!
!  COMMON VARIABLE FOR THE TEST NUMBER
!
integer           ntest
common  /test/    ntest
save /test/
!
!  LOCAL VARIABLES
!
double precision  optval
integer           i
!
!  PROBLEM DIMENSION
!
if ( ntest == 1 ) then
  n=2
  optval=0.25
elseif ( ntest == 2 ) then
  n=10
  optval=24.30621
elseif ( ntest  == 3 ) then
  n=10
  optval=-47.76065
elseif ( ntest  == 4 ) then
  n=2
  optval=0.0
elseif ( ntest  == 5 ) then
  n=10
  optval=-47.75997
elseif ( ntest  == 6 ) then
  n=3
  optval=0.0
elseif ( ntest  == 7 ) then
  n=5
  optval=0.0
elseif ( ntest  == 8 ) then
  n=10
  optval=-47.76065   
!      ELSEIF ( NTEST == 9 ) THEN
!        N=2
!        OPTVAL=0.0
!      ELSEIF ( NTEST == 10 ) THEN
!        N=9
!        OPTVAL=0.0
!      ELSEIF ( NTEST  == 11 ) THEN
!        N=10
!        OPTVAL=-47.75997
!      ELSEIF ( NTEST == 12 ) THEN
!        N=10
!        OPTVAL=24.30621
else
   if ( iout > 0 ) then
      write( iout, 2100 )
   end if
   stop        
endif

!
!  CHECK IF THE MAXIMUM ALLOWED DIMENSION IS BIG ENOUGH
!
if ( n > nmax ) then
   if ( iout > 0 ) then
      write( iout, 2000 ) 'X   ', 'NMAX  ', n - nmax
   end if
   stop
end if



!
!  NUMBER OF CONSTRAINT (LINEAR AND NONLINEAR) 
!

   ncnln=0
if ( ntest == 1 ) then
   m=2
   nclin=0
elseif ( ntest == 2 ) then
   m=5
   nclin=3
elseif ( ntest == 3 ) then
   m=0
   nclin=3
elseif ( ntest == 4 ) then
   nclin=0
   m=0
elseif ( ntest == 5 ) then
   nclin=0
   m=3
elseif ( ntest == 6 ) then
   nclin=0
   m=1
elseif ( ntest == 7 ) then
   nclin=3
   m=0
elseif ( ntest == 8 ) then
   nclin=0
   m=3
elseif ( ntest == 9 ) then
   nclin=0
   m=0
elseif ( ntest == 10 ) then
   nclin=0
   m=0
elseif ( ntest == 11 ) then
   nclin=0
   m=1
   ncnln=2
elseif ( ntest == 12 ) then
   m=3
   nclin=3
   ncnln=2
endif
lda = max(1, nclmax)    
!
!  CHECK THAT THE NUMBER OF CONSTRAINTS DO NOT EXCEED MAXIMUM
!
if ( nclin > nclmax ) then
   if ( iout > 0 ) then
      write( iout, 2000 ) 'ALIN  ', 'NCLMAX  ', nclin - nclmax
   end if
   stop
end if

if ( ncnln > ncnmax ) then
   if ( iout > 0 ) then
      write( iout, 2000 ) 'WRK  ', 'NCLMAX  ', ncnln - ncnmax
   end if
   stop
end if


!
!  SET THE STARTING POINT
!
if ( ntest == 1 ) then
  x(1)=-2.0d0
  x(2)=1.0d0
elseif ( ntest == 2) then
  x(1)=2.0d0
  x(2)=3.0d0
  x(3)=5.0d0
  x(4)=5.0d0
  x(5)=1.0d0
  x(6)=2.0d0
  x(7)=7.0d0
  x(8)=3.0d0
  x(9)=6.0d0
  x(10)=10.0d0
elseif ( ntest == 3) then
  do 30 i = 1, n
    x(i)=0.1d0
30   continue  
elseif ( ntest == 4) then
  x(1)=-1.2d0
  x(2)=1.0d0
elseif ( ntest == 5) then
  do 50 i = 1, n
    x(i)=-2.3d0
50   continue        
elseif ( ntest == 6) then
  x(1)=-2.6d0
  x(2)=2d0
  x(3)=2d0
elseif ( ntest == 7) then
  x(1)=35d0
  x(2)=-31d0
  x(3)=11d0
  x(4)=5d0
  x(5)=-5d0
elseif ( ntest == 8) then
  do 80 i = 1, n
    x(i)=0.1d0
80   continue     
elseif ( ntest == 9 ) then
  x(1)=0.5
  x(2)=0.5
elseif ( ntest == 10 ) then
  x(1)=0.5
  x(2)=0.5
  x(3)=0.5
  x(4)=-0.5
  x(5)=-0.5
  x(6)=-0.5
  x(7)=0.5
  x(8)=-0.5
  x(9)=0.5
elseif ( ntest == 11) then
  do 90 i = 1, n
    x(i)=-2.3d0
90   continue  
elseif ( ntest == 12) then
  x(1)=2.0d0
  x(2)=3.0d0
  x(3)=5.0d0
  x(4)=5.0d0
  x(5)=1.0d0
  x(6)=2.0d0
  x(7)=7.0d0
  x(8)=3.0d0
  x(9)=6.0d0
  x(10)=10.0d0
endif

!
!  SET THE BOUNDS
!
if ( ntest == 1 ) then 
  bl(1)=-0.5d0
  bl(2)=-1.0d20
  bl(3)= 0.0d0
  bl(4)= 0.0d0
  bu(1)=0.5d0
  bu(2)=1.0d0
  bu(3)=1.0d20
  bu(4)=1.0d20
elseif ( ntest == 2 ) then
  do 200 i=1,n
    bl(i)=-1.0d20 
    bu(i)=1.0d20 
200   continue
    bl(n+1)=-105.0d0
    bl(n+2)=0.0d0
    bl(n+3)=-12.0d0
    bu(n+1)=1.0d20
    bu(n+2)=1.0d20
    bu(n+3)=1.0d20
  do 210 i = 1, m
    bl(n+nclin+i)=0.0d0 
    bu(n+nclin+i)=1.0d20 
210   continue
elseif ( ntest == 3 ) then
  do 300 i=1,n
    bl(i)=1.0d-6
    bu(i)=1.0d20 
300   continue
  bl(11)=2.0d0 
  bu(11)=2.0d0 
  bl(12)=1.0d0 
  bu(12)=1.0d0 
  bl(13)=1.0d0 
  bu(13)=1.0d0 
elseif ( ntest == 4 ) then
  bl(1)=-1.0d20
  bl(2)=-1.0d20
  bu(1)=1.0d20
  bu(2)=1.0d20
elseif ( ntest == 5 ) then
  do 500 i=1,n
    bl(i)=-1.0d2
    bu(i)=1.0d2
500   continue
  bl(11)=2.0d0 
  bu(11)=2.0d0 
  bl(12)=1.0d0 
  bu(12)=1.0d0 
  bl(13)=1.0d0 
  bu(13)=1.0d0 
elseif ( ntest == 6 ) then
  do 600 i=1,n
    bl(i)=-1.0d20
    bu(i)=1.0d20
600   continue
    bl(4)=0.0d0
    bu(4)=0.0d0 
elseif ( ntest == 7 ) then
  do 700 i=1,n
    bl(i)=-1.0d20
    bu(i)=1.0d20
700   continue
  do 710 i=1,n
    bl(i+n)=6.0d0
    bu(i+n)=6.0d0
710   continue
elseif ( ntest == 8 ) then
  do 800 i=1,n
    bl(i)=1.0d-6
    bu(i)=1.0d20 
800   continue
  bl(11)=2.0d0 
  bu(11)=2.0d0 
  bl(12)=1.0d0 
  bu(12)=1.0d0 
  bl(13)=1.0d0 
  bu(13)=1.0d0 
elseif ( ntest == 9 ) then
  bl(1)=-1.0d20
  bu(1)=1.0d20
  bl(2)=-1.0d20
  bu(2)=1.0d20
elseif ( ntest == 10 ) then
  bl(1)=-1.0d0
  bu(1)=1.0d0
  bl(2)=-1.0d0
  bu(2)=1.0d0
  bl(3)=-1.0d0
  bu(3)=1.0d0
  bl(4)=-1.0d0
  bu(4)=1.0d0
  bl(5)=-1.0d0
  bu(5)=1.0d0
  bl(6)=-1.0d0
  bu(6)=1.0d0
  bl(7)=-1.0d0
  bu(7)=1.0d0
  bl(8)=-1.0d0
  bu(8)=1.0d0
  bl(9)=-1.0d0
  bu(9)=1.0d0
elseif ( ntest == 11 ) then
  do 900 i=1,n
    bl(i)=-1.0d2
    bu(i)=1.0d2
900   continue
  bl(11)=2.0d0 
  bu(11)=2.0d0 
  bl(12)=1.0d0 
  bu(12)=1.0d0 
  bl(13)=1.0d0 
  bu(13)=1.0d0 
elseif ( ntest == 12 ) then
  do 1200 i=1,n
    bl(i)=-1.0d20 
    bu(i)=1.0d20 
1200   continue
    bl(n+1)=-105.0d0
    bl(n+2)=0.0d0
    bl(n+3)=-12.0d0
    bu(n+1)=1.0d20
    bu(n+2)=1.0d20
    bu(n+3)=1.0d20
  do 1210 i = 1, ncnln+m
    bl(n+nclin+i)=0.0d0 
    bu(n+nclin+i)=1.0d20 
1210   continue
endif
!
!  MATRIX OF LINEAR CONSTRAINTS
!
if ( ntest == 1 ) then
  continue
elseif ( ntest == 2 ) then
  do 230 i = 1, n
    do 240 j = 1, nclin
      alin(j, i)=0.0d0
240     continue
230   continue  
  alin(1,1)=-4.0d0
  alin(1,2)=-5.0d0
  alin(1,7)=3.0d0
  alin(1,8)=-9.0d0
  alin(2,1)=-10.0d0
  alin(2,2)=8.0d0
  alin(2,7)=17.0d0
  alin(2,8)=-2.0d0
  alin(3,1)=8.0d0
  alin(3,2)=-2.0d0
  alin(3,9)=-5.0d0
  alin(3,10)=2.0d0
elseif ( ntest == 3 ) then
  do 330 i = 1, n
    do 340 j = 1, nclin
      alin(j, i)=0.0d0
340     continue
330   continue  
  alin(1,1)=1.0d0
  alin(1,2)=2.0d0
  alin(1,3)=2.0d0
  alin(1,6)=1.0d0
  alin(1,10)=1.0d0
  alin(2,4)=1.0d0
  alin(2,5)=2.0d0
  alin(2,6)=1.0d0
  alin(2,7)=1.0d0
  alin(3,3)=1.0d0
  alin(3,7)=1.0d0
  alin(3,8)=1.0d0
  alin(3,9)=2.0d0
  alin(3,10)=1.0d0
elseif ( ntest == 4 ) then
  continue
elseif ( ntest == 5 ) then
  continue
elseif ( ntest == 6 ) then
  continue
elseif ( ntest == 7 ) then
  do 730 i = 1, n
    do 740 j = 1, nclin
      alin(j, i)=0.0d0
740     continue
730   continue  
  alin(1,1)=1.0d0
  alin(1,2)=2.0d0
  alin(1,3)=3.0d0
  alin(2,2)=1.0d0
  alin(2,3)=2.0d0
  alin(2,4)=3.0d0
  alin(3,3)=1.0d0
  alin(3,4)=2.0d0
  alin(3,5)=3.0d0
elseif ( ntest == 8 ) then
  continue
elseif ( ntest == 9 ) then
  continue
elseif ( ntest == 10 ) then
  continue
elseif ( ntest == 11 ) then
  continue
elseif ( ntest == 12 ) then
  do 830 i = 1, n
    do 840 j = 1, nclin
      alin(j, i)=0.0d0
840     continue
830   continue  
  alin(1,1)=-4.0d0
  alin(1,2)=-5.0d0
  alin(1,7)=3.0d0
  alin(1,8)=-9.0d0
  alin(2,1)=-10.0d0
  alin(2,2)=8.0d0
  alin(2,7)=17.0d0
  alin(2,8)=-2.0d0
  alin(3,1)=8.0d0
  alin(3,2)=-2.0d0
  alin(3,9)=-5.0d0
  alin(3,10)=2.0d0
endif

!
!  OBTAIN NAMES
!
call names( n, m, pname, xnames, gnames )

!
!  PRINT THE EXPECTED OPTIMAL VALUE
! 
write(iout, 3000) pname, optval

return

2000 format( /, ' ** Program USETUP: array length ', a6, ' too small.', &
        /, ' -- Miminimization abandoned.', &
        /, ' -- Increase the parameter ', a6, ' by at least ', i8, &
           ' and restart.'  )
2100 format( /, ' ** Program DFOTES: invalid test problem number;' , &
        /, '1 to 8 are valid.')
3000 format( ' STARTING TEST PROBLEM ', a10, &
        /, ' OPTIMAL VALUE SHOULD BE ABOUT ', d12.4 )
end


!  ***************************************************
!
!  SUBROUTINE FOR NAMES
!
!  ***************************************************


subroutine names( n, m, pname, vnames, gnames )
integer n, m
character * 256 pname, vnames( n ), gnames( * )
!
!  COMMON VARIABLE FOR THE TEST NUMBER
!
integer           ntest
common  /test/    ntest
save /test/
!
!  SETTING NAMES FOR VARIABLES AND PROBLEM
!
if ( ntest == 1 ) then
  pname='HS16'
  vnames(1)='X1'
  vnames(2)='X2'
  gnames(1)='CON1'
  gnames(2)='CON2'
elseif ( ntest == 2 ) then
  pname='HS113'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  vnames(4)='X4'
  vnames(5)='X5'
  vnames(6)='X6'
  vnames(7)='X7'
  vnames(8)='X8'
  vnames(9)='X9'
  vnames(10)='X10'
  gnames(1)='CON1'
  gnames(2)='CON2'
  gnames(3)='CON3'
  gnames(4)='CON4'
  gnames(5)='CON5'
elseif ( ntest == 3 ) then
  pname='HS112'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  vnames(4)='X4'
  vnames(5)='X5'
  vnames(6)='X6'
  vnames(7)='X7'
  vnames(8)='X8'
  vnames(9)='X9'
  vnames(10)='X10'
  gnames(1)='CON1'
  gnames(2)='CON2'
  gnames(3)='CON3'
elseif ( ntest == 4 ) then
  pname='ROSENBR'
  vnames(1)='X1'
  vnames(2)='X2'
elseif ( ntest == 5 ) then
  pname='HS111'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  vnames(4)='X4'
  vnames(5)='X5'
  vnames(6)='X6'
  vnames(7)='X7'
  vnames(8)='X8'
  vnames(9)='X9'
  vnames(10)='X10'
  gnames(1)='CON1'
  gnames(2)='CON2'
  gnames(3)='CON3'
elseif ( ntest == 6 ) then
  pname='HS26'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  gnames(1)='CON1'
elseif ( ntest == 7 ) then
  pname='HS50'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  vnames(4)='X4'
  vnames(5)='X5'
  gnames(1)='CON1'
  gnames(2)='CON2'
  gnames(3)='CON3'
elseif ( ntest == 8 ) then
  pname='HS112'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  vnames(4)='X4'
  vnames(5)='X5'
  vnames(6)='X6'
  vnames(7)='X7'
  vnames(8)='X8'
  vnames(9)='X9'
  vnames(10)='X10'
  gnames(1)='CON1'
  gnames(2)='CON2'
  gnames(3)='CON3'    
elseif ( ntest == 9 ) then
  pname='RANDOM-NOISE'
  vnames(1)='X1'
  vnames(2)='X2'
elseif ( ntest == 10 ) then
  pname='UGLY-PROBLEM'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  vnames(4)='X1'
  vnames(5)='X2'
  vnames(6)='X3'
  vnames(7)='X1'
  vnames(8)='X2'
  vnames(9)='X3'
elseif ( ntest == 11 ) then
  pname='HS111'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  vnames(4)='X4'
  vnames(5)='X5'
  vnames(6)='X6'
  vnames(7)='X7'
  vnames(8)='X8'
  vnames(9)='X9'
  vnames(10)='X10'
  gnames(1)='CON1'
  gnames(2)='CON2'
  gnames(3)='CON3'
elseif ( ntest == 12 ) then
  pname='HS113'
  vnames(1)='X1'
  vnames(2)='X2'
  vnames(3)='X3'
  vnames(4)='X4'
  vnames(5)='X5'
  vnames(6)='X6'
  vnames(7)='X7'
  vnames(8)='X8'
  vnames(9)='X9'
  vnames(10)='X10'
  gnames(1)='CON1'
  gnames(2)='CON2'
  gnames(3)='CON3'
  gnames(4)='CON4'
  gnames(5)='CON5'
endif     
return
end


!  ***************************************************
!
!  SUBROUTINE FOR CONSTRAINTS
!
!  ***************************************************
!
!
!      SUBROUTINE FUNEASYCON(MODE, NCNLN, N, NROWJ, NEEDC, X, C,
!     *                  CJAC, NSTATE)   
!C
!C  CONSTRAINTS OF THE 7 TEST PROBLEMS 
!C
!      INTEGER          MODE, NCNLN, N, NROWJ, NEEDC(NCNLN), NSTATE, I, J
!      
!      DOUBLE PRECISION X(N), C(NROWJ), CJAC(NROWJ,N)
!C
!C  COMMON VARIABLE FOR THE TEST NUMBER
!C
!      INTEGER           NTEST
!      COMMON  /TEST/    NTEST
!      IF ( NTEST == 1 ) THEN 
!        DO 10 I=1, NCNLN
!          IF( NEEDC(I)/=0 ) THEN
!            IF ( MODE == 0 .OR. MODE == 2 )  THEN    
!              IF(I==1) THEN
!                C(1)=X(1)+X(2)*X(2)
!              ENDIF
!              IF(I==2) THEN
!                 C(2)=X(1)*X(1)+X(2)
!              ENDIF
!            ENDIF
!            IF ( MODE == 1 .OR. MODE == 2 ) THEN
!              IF(I==1) THEN
!                CJAC(1,1)=1
!                CJAC(1,2)=2*X(2)
!              ENDIF
!              IF(I==2) THEN
!                CJAC(2,1)=2*X(1)
!                CJAC(2,2)=1      
!              ENDIF
!            ENDIF
!          ENDIF      
! 10     CONTINUE
!      ELSEIF ( NTEST == 2 ) THEN
!       DO 20 I=1, NCNLN
!          IF( NEEDC(I)/=0 ) THEN
!            IF ( MODE == 0 .OR. MODE == 2 )  THEN    
!              IF(I==1) THEN
!                C(1)=-3*(X(1)-2)**2-4*(X(2)-3)**2-2*X(3)**2+
!     +                7*X(4)+120
!              ENDIF
!              IF(I==2) THEN
!                C(2)=-5*X(1)**2-8*X(2)
!     +                - (X(3)-6)**2+2*X(4)+ 40
!              ENDIF
!              IF(I==3) THEN
!                C(3)=-0.5*(X(1)-8)**2-2*(X(2)-4)**2-3*X(5)**2+
!     +                X(6) + 30
!              ENDIF
!              IF(I==4) THEN
!                C(4)=-X(1)**2-2*(X(2)-2)**2+2*X(1)*X(2)-
!     +                14*X(5)+6*X(6)
!              ENDIF
!              IF(I==5) THEN
!                C(5)=3*X(1)-6*X(2)-12*(X(9)-8)**2+7*X(10)
!              ENDIF
!            ENDIF
!            IF ( MODE == 1 .OR. MODE == 2 ) THEN
!              DO 25 J = 1, N
!                CJAC(I, J) = 0.0D0
! 25           CONTINUE  
!              IF(I==1) THEN
!                CJAC(1,1)=-6*X(1)+12D0
!                CJAC(1,2)=-8*X(2)+24D0    
!                CJAC(1,3)=-4*X(3)
!                CJAC(1,4)=7D0    
!              ENDIF
!              IF(I==2) THEN
!                CJAC(2,1)=-10D0*X(1)
!                CJAC(2,2)=-8D0  
!                CJAC(2,3)=-2D0*X(3)+12
!                CJAC(2,4)=2D0     
!              ENDIF
!              IF(I==3) THEN
!                CJAC(3,1)=-X(1)+8D0
!                CJAC(3,2)=-4*X(2)+16D0   
!                CJAC(3,5)=-6*X(5)
!                CJAC(3,6)=1D0   
!              ENDIF
!              IF(I==4) THEN
!                CJAC(4,1)=-2*X(1)+2D0*X(2)
!                CJAC(4,2)=-4*X(2)+8D0+2*X(1)   
!                CJAC(4,5)=-14D0
!                CJAC(4,6)=6D0  
!              ENDIF
!              IF(I==5) THEN
!                CJAC(5,1)=3D0
!                CJAC(5,2)=-6D0  
!                CJAC(5,9)=-24D0*X(9)+192
!                CJAC(5,10)=7D0   
!              ENDIF
!            ENDIF
!          ENDIF      
! 20     CONTINUE        
!      ELSEIF ( NTEST == 3 ) THEN
!        CONTINUE
!      ELSEIF ( NTEST == 4 ) THEN
!        CONTINUE
!      ELSEIF ( NTEST == 5 ) THEN
!        CONTINUE
!      ELSEIF ( NTEST == 6 ) THEN
!        IF( NEEDC(1)/=0 ) THEN
!          IF ( MODE == 0 .OR. MODE == 2 )  THEN    
!                C(1)=(1D0+X(2)**2)*X(1)+X(3)**4-3
!          ENDIF
!          IF ( MODE == 1 .OR. MODE == 2 ) THEN
!            CJAC(1,1)=1D0+X(2)**2
!            CJAC(1,2)=2D0*X(1)*X(2)
!            CJAC(1,3)=4*X(3)**3
!          ENDIF
!        ENDIF
!      ELSEIF ( NTEST == 7 ) THEN 
!        CONTINUE
!      ELSEIF ( NTEST == 8 ) THEN
!        DO 80 I=1, NCNLN
!          IF( NEEDC(I)/=0 ) THEN
!            IF ( MODE == 0 .OR. MODE == 2 )  THEN    
!              IF(I==1) THEN
!                C(1)=X(1)+2*X(2)+2*X(3)+X(6)**2+X(10)
!              ENDIF
!              IF(I==2) THEN
!                C(2)=X(4)+2*X(5)+X(6)**2+X(7)
!              ENDIF
!              IF(I==3) THEN
!                C(3)=X(3)+X(7)+X(8)+2*X(9)+X(10)
!              ENDIF
!            ENDIF
!            IF ( MODE == 1 .OR. MODE == 2 ) THEN
!              DO 85 J = 1, N
!                CJAC(I, J) = 0.0D0
! 85           CONTINUE  
!              IF(I==1) THEN
!                CJAC(1,1)=1.0d0
!                CJAC(1,2)=2.0d0    
!                CJAC(1,3)=2.0d0
!                CJAC(1,6)=2.0d0*X(6)
!                CJAC(1,10)=1.0d0    
!              ENDIF
!              IF(I==2) THEN
!                CJAC(2,4)=1.0d0
!                CJAC(2,5)=2.0d0  
!                CJAC(2,6)=2.0d0*X(6)
!                CJAC(2,7)=1.0d0   
!              ENDIF
!              IF(I==3) THEN
!                CJAC(3,3)=1.0d0
!                CJAC(3,7)=1.0d0 
!                CJAC(3,8)=1.0d0
!                CJAC(3,9)=2.0d0
!                CJAC(3,10)=1.0d0 
!              ENDIF
!            ENDIF
!          ENDIF      
! 80     CONTINUE 
!      ELSEIF ( NTEST == 9 ) THEN
!        CONTINUE
!      ELSEIF ( NTEST == 10 ) THEN
!        CONTINUE
!      ELSEIF ( NTEST == 11 ) THEN
!        DO 50 I=1, NCNLN
!          IF( NEEDC(I)/=0 ) THEN
!            IF ( MODE == 0 .OR. MODE == 2 )  THEN    
!              IF(I==1) THEN
!                C(1)=EXP(X(1))+2*EXP(X(2))+2*EXP(X(3))+EXP(X(6))+
!     +               EXP(X(10))
!              ENDIF
!              IF(I==2) THEN
!                C(2)=EXP(X(4))+2*EXP(X(5))+EXP(X(6))+EXP(X(7))
!              ENDIF
!              IF(I==3) THEN
!                C(3)=EXP(X(3))+EXP(X(7))+EXP(X(8))+2*EXP(X(9))+
!     +               EXP(X(10))
!              ENDIF
!            ENDIF
!            IF ( MODE == 1 .OR. MODE == 2 ) THEN
!              DO 55 J = 1, N
!                CJAC(I, J) = 0.0D0
! 55           CONTINUE  
!              IF(I==1) THEN
!                CJAC(1,1)=EXP(X(1))
!                CJAC(1,2)=2*EXP(X(2))    
!                CJAC(1,3)=2*EXP(X(3))
!                CJAC(1,6)=EXP(X(6))
!                CJAC(1,10)=EXP(X(10))    
!              ENDIF
!              IF(I==2) THEN
!                CJAC(2,4)=EXP(X(4))
!                CJAC(2,5)=2*EXP(X(5))  
!                CJAC(2,6)=EXP(X(6))
!                CJAC(2,7)=EXP(X(7))   
!             ENDIF
!              IF(I==3) THEN
!                CJAC(3,3)=EXP(X(3))
!                CJAC(3,7)=EXP(X(7))  
!                CJAC(3,8)=EXP(X(8))
!                CJAC(3,9)=2*EXP(X(9))
!                CJAC(3,10)=EXP(X(10)) 
!              ENDIF
!            ENDIF
!          ENDIF      
! 50     CONTINUE
!      ENDIF
!      RETURN
!      END



!  *********************************************************
!
!  ROUTINE COMPUTING FUNCTION VALUE
!
!  *********************************************************
subroutine fun(n, m, x, val, c, iferr)

!
!  OBJECTIVE FUNCTION FOR 7  TEST PROBLEM FOR DFO
!
!     
integer n, m, i
double precision x(n), val, c(m), sum
logical          iferr
!
!  COMMON VARIABLE FOR THE TEST NUMBER
!
integer           ntest
common  /test/    ntest
save /test/
!
!  LOCAL VARIABLES
!
double precision chs112(10)
real             gennor
external         gennor

intrinsic        log, exp, sqrt

iferr=.false.
if ( ntest == 1 ) then 
  val = 100.0d0*(x(2)-x(1)**2)**2+(1.0d0-x(1))**2
  c(1)=x(1)+x(2)*x(2)
  c(2)=x(1)*x(1)+x(2)
elseif ( ntest == 2 ) then
  val = x(1)**2 + x(2)**2 + x(1)*x(2) - 14d0*x(1) - 16d0*x(2)+ &
        (x(3)-10)**2 + 4d0*(x(4)-5d0)**2+(x(5)-3d0)**2+ &
        2.0d0*(x(6)-1d0)**2+ &
        5.0d0*x(7)**2+7.0d0*(x(8)-11d0)**2+2d0*(x(9)-10d0)**2+ &
        (x(10)-7)**2+45d0
          c(1)=-3*(x(1)-2)**2-4*(x(2)-3)**2-2*x(3)**2+ &
                7*x(4)+120
          c(2)=-5*x(1)**2-8*x(2) &
                - (x(3)-6)**2+2*x(4)+ 40
          c(3)=-0.5*(x(1)-8)**2-2*(x(2)-4)**2-3*x(5)**2+ &
                x(6) + 30
          c(4)=-x(1)**2-2*(x(2)-2)**2+2*x(1)*x(2)- &
                14*x(5)+6*x(6)
          c(5)=3*x(1)-6*x(2)-12*(x(9)-8)**2+7*x(10)
elseif ( ntest == 3 ) then
  chs112(1)=-6.089d0
  chs112(2)=-17.164d0
  chs112(3)=-34.054d0
  chs112(4)=-5.914d0
  chs112(5)=-24.721d0
  chs112(6)=-14.986d0
  chs112(7)=-24.100d0
  chs112(8)=-10.708d0
  chs112(9)=-26.662d0
  chs112(10)=-22.179d0
  sum = 0.0d0
  val = 0.0d0
  do 30 i=1, n
    sum = sum + x(i)
30   continue
  do 31 i=1, n  
    val = val + x(i)*(chs112(i)+log(x(i)/sum))
31   continue  
elseif ( ntest == 4 ) then
  val = 100.0d0*(x(2)-x(1)**2)**2+(1.0d0-x(1))**2        
elseif ( ntest == 5 ) then
  chs112(1)=-6.089d0
  chs112(2)=-17.164d0
  chs112(3)=-34.054d0
  chs112(4)=-5.914d0
  chs112(5)=-24.721d0
  chs112(6)=-14.986d0
  chs112(7)=-24.100d0
  chs112(8)=-10.708d0
  chs112(9)=-26.662d0
  chs112(10)=-22.179d0
  sum = 0.0d0
  val = 0.0d0
  do 50 i=1, n
    sum = sum + exp(x(i))
50   continue
  do 51 i=1, n  
    val = val + exp(x(i))*(chs112(i)+x(i)-log(sum))
51   continue  
          c(1)=exp(x(1))+2*exp(x(2))+2*exp(x(3))+exp(x(6))+ &
               exp(x(10))
          c(2)=exp(x(4))+2*exp(x(5))+exp(x(6))+exp(x(7))
          c(3)=exp(x(3))+exp(x(7))+exp(x(8))+2*exp(x(9))+ &
               exp(x(10))
elseif ( ntest == 6 ) then
  val = (x(1)-x(2))**2+(x(2)-x(3))**4
  c(1)=(1d0+x(2)**2)*x(1)+x(3)**4-3
elseif ( ntest == 7 ) then
  val = (x(1)-x(2))**2+(x(2)-x(3))**2+(x(3)-x(4))**4+ &
        (x(4)-x(5))**2
elseif ( ntest == 8 ) then
  chs112(1)=-6.089d0
  chs112(2)=-17.164d0
  chs112(3)=-34.054d0
  chs112(4)=-5.914d0
  chs112(5)=-24.721d0
  chs112(6)=-14.986d0
  chs112(7)=-24.100d0
  chs112(8)=-10.708d0
  chs112(9)=-26.662d0
  chs112(10)=-22.179d0
  sum = 0.0d0
  val = 0.0d0
  do 80 i=1, n
    sum = sum + x(i)
80   continue
  do 81 i=1, n  
    val = val + x(i)*(chs112(i)+log(x(i)/sum))
81   continue  
          c(1)=x(1)+2*x(2)+2*x(3)+x(6)**2+x(10)
          c(2)=x(4)+2*x(5)+x(6)**2+x(7)
          c(3)=x(3)+x(7)+x(8)+2*x(9)+x(10)
elseif ( ntest == 9 ) then
!        val=gennor(0.0, 0.0001)
  val=x(1)**2+x(2)**2
elseif ( ntest == 10 ) then
  val=sqrt((x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2)
  sum=exp((-12.0d0)*log(val))-2.0d0*exp((-6.0d0)*log(val))
  val=sqrt((x(1)-x(7))**2+(x(2)-x(8))**2+(x(3)-x(9))**2)
  sum=sum+exp((-12.0d0)*log(val))-2.0d0*exp((-6.0d0)*log(val))
  val=sqrt((x(7)-x(4))**2+(x(8)-x(5))**2+(x(9)-x(6))**2)
  sum=sum+exp((-12.0d0)*log(val))-2.0d0*exp((-6.0d0)*log(val))
  val=sum
elseif ( ntest == 11 ) then
  chs112(1)=-6.089d0
  chs112(2)=-17.164d0
  chs112(3)=-34.054d0
  chs112(4)=-5.914d0
  chs112(5)=-24.721d0
  chs112(6)=-14.986d0
  chs112(7)=-24.100d0
  chs112(8)=-10.708d0
  chs112(9)=-26.662d0
  chs112(10)=-22.179d0
  sum = 0.0d0
  val = 0.0d0
  do 90 i=1, n
    sum = sum + exp(x(i))
90   continue
  do 91 i=1, n  
    val = val + exp(x(i))*(chs112(i)+x(i)-log(sum))
91   continue  
!                C(1)=EXP(X(1))+2*EXP(X(2))+2*EXP(X(3))+EXP(X(6))+
!     +               EXP(X(10))
!                C(2)=EXP(X(4))+2*EXP(X(5))+EXP(X(6))+EXP(X(7))
          c(1)=exp(x(3))+exp(x(7))+exp(x(8))+2*exp(x(9))+ &
               exp(x(10))
elseif ( ntest == 12 ) then
  val = x(1)**2 + x(2)**2 + x(1)*x(2) - 14d0*x(1) - 16d0*x(2)+ &
        (x(3)-10)**2 + 4d0*(x(4)-5d0)**2+(x(5)-3d0)**2+ &
        2.0d0*(x(6)-1d0)**2+ &
        5.0d0*x(7)**2+7.0d0*(x(8)-11d0)**2+2d0*(x(9)-10d0)**2+ &
        (x(10)-7)**2+45d0
!                C(1)=-3*(X(1)-2)**2-4*(X(2)-3)**2-2*X(3)**2+
!     +                7*X(4)+120
!                C(2)=-5*X(1)**2-8*X(2)
!     +                - (X(3)-6)**2+2*X(4)+ 40
          c(1)=-0.5*(x(1)-8)**2-2*(x(2)-4)**2-3*x(5)**2+ &
                x(6) + 30
          c(2)=-x(1)**2-2*(x(2)-2)**2+2*x(1)*x(2)- &
                14*x(5)+6*x(6)
          c(3)=3*x(1)-6*x(2)-12*(x(9)-8)**2+7*x(10)
endif     
return
end


