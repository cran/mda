c$$$ Naras fix
c$$$      subroutine calcvar(nx,n,px,qr,qrank,qpivot,cov,tmpcov,work)
      subroutine calcvar(nx,px,qr,qrank,cov,tmpcov)
      implicit double precision (a-h,o-z)
c$$$ Naras fix
c$$$      integer n,px,qrank,qpivot(px)
c$$$      double precision qr(nx,px),cov(px,px), tmpcov(px,px),work(1)
      integer px,qrank
      double precision qr(nx,px),cov(px,px), tmpcov(px,px)
      double precision dsum
      integer i,j,km
      do 23000 i=1,qrank
      do 23002 j=1,qrank
      tmpcov(i,j)=0d0
      cov(i,j)=qr(i,j)
23002 continue
      tmpcov(i,i)=1e0
23000 continue
      info=0
c R version has different args
c      call dbksl(cov,px,qrank,tmpcov,px,info)
      do 20 j = 1, qrank
         call dtrsl(cov, px, qrank, tmpcov(1,j), 01, info)
 20   continue
      do 23004 i=1,qrank
      do 23006 j=i,qrank
      dsum=0e0
      km=max(i,j)
      k=km
23008 if(.not.(k.le.qrank))goto 23010
      dsum=dsum+tmpcov(i,k)*tmpcov(j,k)
      k=k+1
      goto 23008
23010 continue
      tmpcov(i,j)=dsum
      tmpcov(j,i)=dsum
23006 continue
23004 continue
      do 23011 i=1,qrank
      do 23013 j=1,qrank
      cov(i,j)=tmpcov(i,j)
23013 continue
23011 continue
      return
      end
