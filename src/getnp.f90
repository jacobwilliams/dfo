!
!
!
subroutine getnp( ipoly, poly, lpoly, n, kappa, g, h)


!
!     collects  the coefficient of the polynomial IPOLY into
!     doefficient of a quadratic kappa+g'x+0.5 x'H'x


!     INPUT
!     POLY   = Newton polynomials
!     IPOLY  = The index of the polynomial we are processing

!     OUTPUT
!     KAPPA = the constant of the polynomail
!     G     = the linear coeffitients
!     H     = the n by n matrix of quadratic doeffitients
integer           ipoly, lpoly, n

double precision  poly(lpoly), kappa,  g(n), h(n,n)
!
!     Local variables
!

double precision  zero
parameter        (zero=0.0d0)
integer           i, j, kbeg, k, np1, dd

np1=n+1
dd=(np1)*(n+2)/2

kappa=zero
do 10 i=1,n
  g(i)=zero
  do 20 j=1,n
    h(i,j)=zero
20  continue
10 continue

 
if (ipoly.eq.1) then
  kappa=poly(1)
elseif (ipoly .le. np1) then
  kbeg=2+(ipoly-2)*(np1)
  kappa=poly(kbeg)
  do 30 i=1, n
   g(i)=poly(kbeg+i)
30   continue
else        
  kbeg=2+n*np1+(ipoly-np1-1)*dd
  kappa=poly(kbeg)
  do 40 i=1, n
    g(i)=poly(kbeg+i)
40   continue
  k=kbeg+n+1
  do 50 i=1, n
    h(i,i)=2*poly(k)
    k=k+1
    do 60 j=i+1,n
      h(i,j)=poly(k)
      h(j,i)=poly(k)
      k=k+1
60     continue
50   continue
endif


return
end
