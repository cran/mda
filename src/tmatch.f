      subroutine tpack(n,p,match,x,xbar)
      integer n,p,match(n)
      double precision x(n),xbar(n)
      do 23000 i=1,p
      xbar(i)=0e0
23000 continue
      do 23002 i=1,n
      xbar(match(i))=xbar(match(i))+x(i)
23002 continue
      return
      end
      subroutine suff(n,p,match,x,y,w,xbar,ybar,wbar,work)
      integer n,p,match(n)
      double precision x(n),xbar(n),y(n),ybar(n),w(n),wbar(n),work(n)
      call tpack(n,p,match,w,wbar)
      do 23004 i=1,n
      xbar(match(i))=x(i)
23004 continue
      do 23006 i=1,n
      work(i)=y(i)*w(i)
23006 continue
      call tpack(n,p,match,work,ybar)
      do 23008 i=1,p
      if(.not.(wbar(i).gt.0e0))goto 23010
      ybar(i)=ybar(i)/wbar(i)
      goto 23011
23010 continue
      ybar(i)=0e0
23011 continue
23008 continue
      return
      end
      subroutine untpack(n,p,match,xbar,x)
      integer n,p, match(n)
      double precision x(n),xbar(p+1)
      if(.not.(p.lt.n))goto 23012
      xbar(p+1)=0e0
23012 continue
      do 23014 i = 1,n
      x(i)=xbar(match(i))
23014 continue
      return
      end
