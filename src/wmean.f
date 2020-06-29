      double precision function wmean(n,y,w)
      integer n
      double precision y(n),w(n),wtot,wsum
      wtot=0e0
      wsum=0e0
      do 23028 i=1,n
      wsum=wsum+y(i)*w(i)
      wtot=wtot+w(i)
23028 continue
      if(.not.(wtot .gt. 0e0))goto 23030
      wmean=wsum/wtot
      goto 23031
23030 continue
      wmean=0e0
23031 continue
      return
      end
