      subroutine orthreg(nx,n,p,x,in, y,res)
      implicit double precision (a-h,o-z)
      integer n,nx,p, in(p)
      double precision x(nx,p),y(n),res(n)
      do 23000 i=1,n 
      res(i)=y(i)
23000 continue
      do 23002 j=1,p 
      if(.not.(in(j).eq.1))goto 23004
      temp1=0
      temp2=0
      do 23006 i=1,n 
      temp1=temp1+res(i)*x(i,j)
      temp2=temp2+x(i,j)*x(i,j)
23006 continue
      beta=temp1/temp2
      do 23008 i=1,n 
      res(i)=res(i)-beta*x(i,j)
23008 continue
23004 continue
23002 continue
      return
      end
