      subroutine sknotl(x,n,knot,k)
      implicit double precision(a-h,o-z)
      integer n,k,ndk,j
      double precision x(n),knot(n+6),a1,a2,a3,a4
      a1 = log(50e0)/log(2e0) 
      a2 = log(100e0)/log(2e0)
      a3 = log(140e0)/log(2e0) 
      a4 = log(200e0)/log(2e0)
      if(.not.(n.lt.50))goto 23000
      ndk = n 
      goto 23001
23000 continue
      if(.not.(n.ge.50 .and. n.lt.200))goto 23002
c$$$ Naras fix
c$$$      ndk = 2.**(a1+(a2-a1)*(n-50.)/150.) 
      ndk = int(2.**(a1+(a2-a1)*(n-50.)/150.))
      goto 23003
23002 continue
      if(.not.(n.ge.200 .and. n.lt.800))goto 23004
c$$$ Naras fix
c$$$      ndk = 2.**(a2+(a3-a2)*(n-200.)/600.) 
      ndk = int(2.**(a2+(a3-a2)*(n-200.)/600.))
      goto 23005
23004 continue
      if(.not.(n.ge.800 .and. n.lt.3200))goto 23006
c$$$ Naras fix
c$$$      ndk = 2.**(a3+(a4-a3)*(n-800.)/2400.) 
      ndk = int(2.**(a3+(a4-a3)*(n-800.)/2400.))
      goto 23007
23006 continue
      if(.not.(n.ge.3200))goto 23008
c$$$ Naras fix
c$$$      ndk = 200. + (n-3200)**.2 
      ndk = int(200. + (n-3200)**.2)
23008 continue
23007 continue
23005 continue
23003 continue
23001 continue
      k = ndk + 6
      do 23010 j=1,3 
      knot(j) = x(1) 
23010 continue
      do 23012 j=1,ndk 
      knot(j+3) = x( 1 + (j-1)*(n-1)/(ndk-1) ) 
23012 continue
      do 23014 j=1,3 
      knot(ndk+3+j) = x(n) 
23014 continue
      return
      end
