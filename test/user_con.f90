subroutine easycon(n, x, ncnln, c)   


integer           ncnln, n

double precision x(n), c(*)

integer           i, j
!
!  COMMON VARIABLE FOR THE TEST NUMBER
!
integer           ntest
common  /test/    ntest
save /test/
!
!  CAN BE A DUMMY ROUTINE. IF CONSTRAINTS ARE PRESENT SHOULD COMPUTE
!  THE VALUES OF THE CONSTRAINTS AT POINT X AND STORE IT IN C
!
if ( ntest == 11 ) then
  if ( ncnln > 0 ) then
    do 50 i=1, ncnln
        if(i==1) then
          c(1)=exp(x(1))+2*exp(x(2))+2*exp(x(3))+exp(x(6))+ &
               exp(x(10))
        endif
        if(i==2) then
          c(2)=exp(x(4))+2*exp(x(5))+exp(x(6))+exp(x(7))
        endif
50     continue
  endif
elseif (ntest == 12 ) then
  if ( ncnln > 0 ) then
          c(1)=-3*(x(1)-2)**2-4*(x(2)-3)**2-2*x(3)**2+ &
                7*x(4)+120
          c(2)=-5*x(1)**2-8*x(2) &
                - (x(3)-6)**2+2*x(4)+ 40
  endif         
endif 
return
end

subroutine easyjac(n, x, ncnln, nrowc, cjac)   
integer           ntest
common  /test/    ntest 
save /test/

integer           ncnln, n, nrowc

double precision x(n), cjac(nrowc,n)

integer           i, j
!
!  CAN BE A DUMMY ROUTINE. IF CONSTRAINTS ARE PRESENT SHOULD COMPUTE
!  THE JACOBIAN OF THE CONSTRAINTS AT POINT X AND STORE IT IN CJAC
!  SO THAT THE FIST ROW IS THE GRADIENT OF THE FIRST CONSTRAINT
!
if ( ntest == 11 ) then
  do 50 i=1, ncnln
        do 55 j = 1, n
          cjac(i, j) = 0.0d0
55         continue  
        if(i==1) then
          cjac(1,1)=exp(x(1))
          cjac(1,2)=2*exp(x(2))    
          cjac(1,3)=2*exp(x(3))
          cjac(1,6)=exp(x(6))
          cjac(1,10)=exp(x(10))    
        endif
        if(i==2) then
          cjac(2,4)=exp(x(4))
          cjac(2,5)=2*exp(x(5))  
          cjac(2,6)=exp(x(6))
          cjac(2,7)=exp(x(7))   
        endif
        if(i==3) then
          cjac(3,3)=exp(x(3))
          cjac(3,7)=exp(x(7))  
          cjac(3,8)=exp(x(8))
          cjac(3,9)=2*exp(x(9))
          cjac(3,10)=exp(x(10)) 
        endif
50   continue
elseif (ntest == 12 ) then
  if ( ncnln > 0 ) then
   do 65 j = 1, n
     cjac(1, j) = 0.0d0
65    continue  
   do 66 j = 1, n
     cjac(2, j) = 0.0d0
66    continue  
    cjac(1,1)=-6d0*(x(1)-2d0)
    cjac(1,2)=-8d0*(x(2)-3d0)
    cjac(1,3)=-4d0*x(3)
    cjac(1,4)=7d0
    cjac(2,1)=-10d0*x(1)
    cjac(2,2)=-8d0
    cjac(2,3)=-2d0*(x(3)-6d0)
    cjac(2,4)=2d0
  endif         
endif 
return
end

subroutine easyhess(k, n, x, ncnln, nrowc, chess, lambda)   
integer           ntest
common  /test/    ntest
save /test/
integer           ncnln, n, nrowc, k

double precision x(n), chess(nrowc,n), lambda

integer           i, j
!
!  CAN BE A DUMMY ROUTINE. IF CONSTRAINTS ARE PRESENT SHOULD COMPUTE
!  THE HESSIAN OF THE K-TH CONSTRAINTS AT POINT X AND ADD IT TO
!  EXISTING CHESS WITH COEFFICIENT LAMBDA 
!      CHESS=CHESS+LAMBDA*HESS(C_K(X))

if ( ntest == 11 ) then
  if ( k== 1) then 
    chess(1,1)=chess(1,1)+exp(x(1))*lambda
    chess(2,2)=chess(2,2)+2*exp(x(2))*lambda
    chess(3,3)=chess(3,3)+exp(x(3))*(2*lambda)
    chess(6,6)=chess(6,6)+exp(x(6))*(lambda)
    chess(10,10)=chess(10,10)+exp(x(10))*(lambda)  
  else if ( k== 2) then   
    chess(4,4)=chess(4,4)+exp(x(4))*lambda
    chess(5,5)=chess(5,5)+2*exp(x(5))*lambda
    chess(6,6)=chess(6,6)+exp(x(6))*lambda
    chess(7,7)=chess(7,7)+exp(x(7))*(lambda)
  endif
elseif (ntest == 12 ) then
  if ( ncnln > 0 ) then
   if ( k == 1 ) then
    chess(1,1)=chess(1,1)-lambda*6d0
    chess(2,2)=chess(2,2)-lambda*8d0
    chess(3,3)=chess(3,3)-lambda*4d0
   elseif ( k== 2) then 
    chess(1,1)=chess(1,1)-lambda*10d0
    chess(3,3)=chess(3,3)-lambda*2d0
  endif 
 endif           
endif  
return
end




