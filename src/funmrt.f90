subroutine funmrt(mode, n, x, objf, objgrd, nstate)
    
double precision objf, x(n), objgrd(n), val

integer          mode, n, nstate, i, j, k

double precision  half

parameter        (half=0.5d0)

if (mode/=2) then
!  GIVE A WARNING

endif

objf=0
do 10 i=1,n
  objf=objf+gmod(i)*x(i)
  objgrd(i)=gmod(i)
  do 20 j=1,n
    objf=objf+half*hmod(i,j)*x(j)*x(i)
    objgrd(i)=objgrd(i)+hmod(i,j)*x(j)
20   continue
10 continue

do 30 k=1,ncon
  val=ccon(k)
  do 40 i=1,n
  val=val+lcon((k-1)*n+i)*x(i)
    do 50 j=1,n
      val=val+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
50     continue
40   continue
  if ( val < conl(k) ) then 
    objf = objf + penpar * ( conl(k) - val )
    do 60 i=1,n
      objgrd(i)=objgrd(i)-penpar*lcon((k-1)*n+i)
      do 70 j=1,n
        objgrd(i)=objgrd(i)-penpar*qcon((k-1)*n*n+(i-1)*n+j)*x(j)
70       continue
60     continue
  elseif ( val > conu(k) ) then 
    objf = objf + penpar * ( val - conu(k) )
    do 80 i=1,n
      objgrd(i)=objgrd(i)+penpar*lcon((k-1)*n+i)
      do 90 j=1,n
        objgrd(i)=objgrd(i)+penpar*qcon((k-1)*n*n+(i-1)*n+j)*x(j)
90       continue
80     continue
  endif
30 continue        

return
end





