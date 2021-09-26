subroutine funobj(mode, n, x, objf, objgrd, nstate)
    
double precision objf, x(n), objgrd(n)

integer          mode, n, nstate, i, j, k

double precision  half, val

parameter        (half=0.5d0)
!
!  IF IPOPT FLAG IS >0 THEN WE ARE CALLED BY IPOPT AND THEN MODE=1 MEANS
!  THAT THE OBJECTIVE FUNCTION VALUE IS EXPACTED. IF MODE = 2 THEN
!  THE GRADIENT OF THE OBJECTIVE IS REQUIRED. IF IPOPT =< 0 THEN
!  NPSOL IS CALLING, THE COMPUTE OBJF AND OBJGRAD AT ONCE.
!
include 'dfo_model_inc.inc'

if ( useipopt > 0 ) then 
  if (mode == 1 ) then      
    objf=0
    do 10 i=1,n
      objf=objf+gmod(i)*x(i)
      do 20 j=1,n
        objf=objf+half*hmod(i,j)*x(j)*x(i)
20       continue
10     continue
    if ( usemerit ) then
      do 170 k=1,ncon
      val=ccon(k)
      do 180 i=1,n
        val=val+lcon((k-1)*n+i)*x(i)
        do 190 j=1,n
          val=val+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
190         continue
180       continue
      if ( val < conl(k) ) then 
        objf = objf + penpar * ( conl(k) - val )
      elseif ( val > conu(k) ) then 
        objf = objf + penpar * ( val - conu(k) )
      endif
170       continue  
    endif
  else
    do 30 i=1,n
      objgrd(i)=gmod(i)
      do 40 j=1,n
        objgrd(i)=objgrd(i)+hmod(i,j)*x(j)
40       continue
30     continue
    if ( usemerit ) then
      do 270 k=1,ncon
        val=ccon(k)
        do 280 i=1,n
          val=val+lcon((k-1)*n+i)*x(i)
          do 290 j=1,n
            val=val+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
290           continue
280         continue
        if ( val < conl(k) ) then 
          do 200 i=1,n
            objgrd(i)=objgrd(i)-penpar*lcon((k-1)*n+i)
            do 210 j=1,n
              objgrd(i)=objgrd(i)-penpar*qcon((k-1)*n*n+(i-1)*n+j) &
              *x(j)
210             continue
200           continue
        elseif ( val > conu(k) ) then 
          do 220 i=1,n
            objgrd(i)=objgrd(i)+penpar*lcon((k-1)*n+i)
            do 230 j=1,n
              objgrd(i)=objgrd(i)+penpar*qcon((k-1)*n*n+(i-1)*n+j) &
              *x(j)
230             continue
220           continue
        endif
270       continue
    endif
  endif  
else
  objf=0
  do 50 i=1,n
    objf=objf+gmod(i)*x(i)
    objgrd(i)=gmod(i)
    do 60 j=1,n
      objf=objf+half*hmod(i,j)*x(j)*x(i)
      objgrd(i)=objgrd(i)+hmod(i,j)*x(j)
60     continue
50   continue
  if ( usemerit ) then
    do 70 k=1,ncon
      val=ccon(k)
      do 80 i=1,n
        val=val+lcon((k-1)*n+i)*x(i)
        do 90 j=1,n
          val=val+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
90         continue
80       continue
      if ( val < conl(k) ) then 
        objf = objf + penpar * ( conl(k) - val )
        do 100 i=1,n
          objgrd(i)=objgrd(i)-penpar*lcon((k-1)*n+i)
          do 110 j=1,n
            objgrd(i)=objgrd(i)-penpar*qcon((k-1)*n*n+(i-1)*n+j) &
            *x(j)
110           continue
100         continue
      elseif ( val > conu(k) ) then 
        objf = objf + penpar * ( val - conu(k) )
        do 120 i=1,n
          objgrd(i)=objgrd(i)+penpar*lcon((k-1)*n+i)
          do 130 j=1,n
            objgrd(i)=objgrd(i)+penpar*qcon((k-1)*n*n+(i-1)*n+j) &
            *x(j)
130           continue
120         continue
      endif
70     continue  
  endif
endif
return
end





