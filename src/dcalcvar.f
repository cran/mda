      subroutine calcvar(nx,n,px,qr,qrank,qpivot,cov,tmpcov,work)
      implicit double precision (a-h,o-z)
      integer n,px,qrank,qpivot(px)
      double precision qr(nx,px),cov(px,px), tmpcov(px,px),work(1)
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
 20      call dtrsl(cov, px, qrank, tmpcov(1,j), 01, info)
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
