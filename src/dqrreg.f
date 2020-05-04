      subroutine qrreg(nx,n,px,p,nclass,x,xsc,in,y,qpivot,qrank,beta,
     &     res,rss,cvar,var,varsc,scr1,work)
      implicit double precision (a-h,o-z)
      integer nx,n,p,px, qpivot(p),qrank,nclass,in(p)
      double precision x(nx,p), xsc(n,p), y(n,nclass),res(nx,nclass),
     &     beta(px,nclass),work(*),scr1(p),var(px,p),varsc(px,p)
      logical cvar
      ii=0
      do 23000 j=1,p 
      if(.not.(in(j).eq.1))goto 23002
      ii=ii+1
      do 23004 i=1,n 
      xsc(i,ii)=x(i,j)
23004 continue
23002 continue
23000 continue
      nt=ii
      ijob=101
      info=1
      temp3=1d-2
      do 23006 i=1,p 
      qpivot(i)=i
23006 continue
      call dqrdc2(xsc,n,n,nt,temp3,qrank,scr1,qpivot,work)
      rss=0.0
      do 23008 k=1,nclass
      call dqrsl(xsc,n,n,qrank,scr1,y(1,k),work(1),work(1),beta(1,k),
     &     work(1),res(1,k),ijob,info)
      do 23010 i=1,n 
      res(i,k)=y(i,k)-res(i,k)
      rss=rss+res(i,k)*res(i,k)
23010 continue
23008 continue
      if(.not.(cvar))goto 23012
c$$$ Naras fix
c$$$      call calcvar(nx,n,px,xsc,qrank,qpivot,var,varsc,work)
      call calcvar(nx,px,xsc,qrank,var,varsc)
23012 continue
      return
      end
