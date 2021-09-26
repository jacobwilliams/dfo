subroutine funcon(mode, ncnln, n, nrowj, needc, x, c, &
                  cjac, nstate)   


integer          mode, ncnln, n, nrowj, needc(*), nstate

double precision x(n), c(*), cjac(nrowj,n)

integer           i, j, k, neasy, l, ncont

double precision  half

parameter        (half=0.5d0)
include 'dfo_model_inc.inc'

neasy=nlin+nnln
!
!  IF WE ARE DEALING WITH IPOPT, THEN MODE=1,2,3 INDICATES WHETHER
!  THE VALUES OF THE CONSTRAINTS, THE JACOBIAN OF THE HESSIAN OF THE 
!  LAGRANGIAN IS EXPECTED
!
if (useipopt .gt. 0) then
  if ( mode .eq. 1) then
    do 5  k=1, nlin
      c(k)=0.0d0
      do 6 i=1, n
        c(k)= c(k)+amat(k, i)*x(i)
6       continue
5     continue  
    if ( nnln .gt. 0 ) call easycon(n, x, nnln, c(nlin+1))
    if ( .not. usemerit) then
      do 10 k = 1, ncon
        l=k+neasy
        c(l)=ccon(k)
        do 30 i = 1, n
          c(l)=c(l)+lcon((k-1)*n+i)*x(i)
          do 20 j = 1, n
            c(l)=c(l)+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
20           continue
30         continue
10       continue  
    endif
  elseif (mode .eq. 2) then    
    do 15  k=1, nlin
      do 16 i=1, n
        cjac(k,i)= amat(k, i)
16        continue
15     continue  
    if ( nnln .gt. 0 ) call easyjac(n, x, nnln, nrowj, &
                       cjac(nlin+1,1))
    if ( .not. usemerit) then
      do 60 k = 1, ncon
        l=k+neasy
        do 50 i = 1, n
          cjac(l, i)=lcon((k-1)*n+i)
          do 40 j = 1, n
            cjac(l, i)=cjac(l, i)+qcon((k-1)*n*n+(i-1)*n+j)*x(j)
40           continue
50         continue
60       continue 
    endif 
   else if ( mode .eq. 3 ) then
    do 65  k=1, n
      do 66 i=1, n
        cjac(k,i)= hmod(k,i)
66        continue
65      continue  
    do 70 k = nlin + 1, neasy
      call easyhess(k, n, x, nnln, nrowj, cjac, c(k))
70     continue  
    if ( .not. usemerit) then
      do 100 k = 1, ncon
        l=k+neasy
        do 90 i = 1, n
          do 80 j = 1, n
            cjac(j, i)=cjac(j, i)+c(l)*qcon((k-1)*n*n+(i-1)*n+j)
80           continue
90         continue
100       continue 
    endif           
   endif
 else
  call easycon(n, x, nnln, c(1))
  call easyjac(n, x, nnln, nrowj, cjac)
  if ( .not. usemerit ) then
    do 110 k = 1, ncon
      l=k+nnln
      c(l)=ccon(k)
      do 130 i = 1, n
      c(l)=c(l)+lcon((k-1)*n+i)*x(i)
      cjac(l, i)=lcon((k-1)*n+i)
      do 120 j = 1, n
        c(l)=c(l)+half*qcon((k-1)*n*n+(i-1)*n+j)*x(i)*x(j)
        cjac(l, i)=cjac(l, i)+qcon((k-1)*n*n+(i-1)*n+j)*x(j)
120       continue
130     continue
110    continue       
 endif
endif   
return
end









