 subroutine rzrmat(a,n,m)
 double precision a(n,m)
 integer          m,n,i,j

 do 10 i=1,n
   do 20 j=1,m
     a(i,j)=0.0d0
20    continue
10  continue
 return
 end


 subroutine rzrvec(b,n)
 double precision b(n)
 integer          n,i

 do 10 i=1,n        
     b(i)=0.0d0
10  continue
 return
 end



 subroutine izrmat(a,n,m)
  
 integer        a(n,m),m,n,i,j

 do 10 i=1,n
   do 20 j=1,m
     a(i,j)=0
20    continue
10  continue
 return
 end

 subroutine izrvec(b,n)

 integer           b(n),n,i

 do 10 i=1,n        
     b(i)=0
10  continue
 return
 end
