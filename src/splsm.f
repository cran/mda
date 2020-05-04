c$$$ Naras fix
c$$$      subroutine splsm(x,y,w,n,match,nef,spar,dof,smo,s0,cov,ifcov,work,
c$$$     &      lenw)
c$$$      implicit double precision(a-h,o-z)
c$$$      integer n,match(n),nef,lenw
c$$$      double precision x(n),y(n),w(n),spar,dof,smo(n),s0,cov(n),work(
c$$$     &      lenw)
c$$$      logical ifcov
c$$$      call splsm1(x, y, w, n, match, nef, spar, dof, smo, s0, cov, 
c$$$     &      ifcov, work(1), work(nef+2), work(2*nef+3), work(3*nef+4), 
c$$$     &      work(3*nef+n+10), lenw)
      subroutine splsm(x,y,w,n,match,nef,spar,dof,smo,s0,cov,work,
     &      lenw)
      implicit double precision(a-h,o-z)
      integer n,match(n),nef,lenw
      double precision x(n),y(n),w(n),spar,dof,smo(n),s0,cov(n),work(
     &      lenw)
      call splsm1(x, y, w, n, match, nef, spar, dof, smo, s0, cov, 
     &      work(1), work(nef+2), work(2*nef+3), work(3*nef+4), 
     &      work(3*nef+n+10), lenw)
      return
      end
c$$$ Naras fix
c$$$      subroutine splsm1(x,y,w,n,match,nef,spar,dof,smo,s0,lev,ifcov,xin,
c$$$     &      yin,win,knot,work,lenw)
c$$$      implicit double precision(a-h,o-z)
c$$$      integer n,match(n),nef,lenw
c$$$      double precision x(n),y(n),w(n),spar,dof,smo(n),s0,lev(n),work(
c$$$     &      lenw)
c$$$      logical ifcov
      subroutine splsm1(x,y,w,n,match,nef,spar,dof,smo,s0,lev,xin,
     &      yin,win,knot,work,lenw)
      implicit double precision(a-h,o-z)
      integer n,match(n),nef,lenw
      double precision x(n),y(n),w(n),spar,dof,smo(n),s0,lev(n),work(
     &      lenw)
      double precision xin(nef+1),yin(nef+1),win(nef+1),knot(nef+4)
      integer nk,ldnk,ld4,k
      double precision xmin,xrange
      call suff(n,nef,match,x,y,w,xin,yin,win,work(1))
      xmin=xin(1)
      xrange=xin(nef)-xin(1)
      do 23000 i=1,nef 
      xin(i)=(xin(i)-xmin)/xrange
23000 continue
      call sknotl(xin,nef,knot,k)
      nk=k-4
      ld4=4
      ldnk=1
c$$$ Naras fix
c$$$      call splsm2(x, y, w, n, match, nef, spar, dof, smo, s0, lev, 
c$$$     &    ifcov, xin, yin,  win, knot, work(1),  work(nk+1),  
c$$$     &    work(nk+nef+2), work(nk+2*nef+3),  work(2*nk+2*nef+3), 
c$$$     &    work(3*nk+2*nef+3), work(4*nk+2*nef+3), work(5* nk+2*nef+3), 
c$$$     &    work(6*nk+2*nef+3), work(7*nk+2*nef+3), work(8*nk+2*nef+ 3), 
c$$$     &    work(9*nk+2*nef+3), work(10*nk+2*nef+3), 
c$$$     &    work((10+ld4)*nk+2*nef+ 3), work((10+2*ld4)*nk+2*nef+3), 
c$$$     &    ld4, ldnk, nk)
      call splsm2(y, w, n, match, nef, spar, dof, smo, s0, lev, 
     &    xin, yin,  win, knot, work(1),  work(nk+1),  
     &    work(nk+nef+2), work(nk+2*nef+3),  work(2*nk+2*nef+3), 
     &    work(3*nk+2*nef+3), work(4*nk+2*nef+3), work(5* nk+2*nef+3), 
     &    work(6*nk+2*nef+3), work(7*nk+2*nef+3), work(8*nk+2*nef+ 3), 
     &    work(9*nk+2*nef+3), work(10*nk+2*nef+3), 
     &    work((10+ld4)*nk+2*nef+ 3), work((10+2*ld4)*nk+2*nef+3), 
     &    ld4, ldnk, nk)
      return
      end
c$$$ Naras fix
c$$$      subroutine splsm2(x, y, w, n, match, nef, spar, dof, smo, s0, 
c$$$     &    lev, ifcov, xin, yin, win, knot, coef, sout, levout, xwy, 
c$$$     &    hs0, hs1, hs2, hs3, sg0, sg1, sg2, sg3, abd, p1ip, p2ip, 
c$$$     &    ld4, ldnk, nk)
c$$$      implicit double precision(a-h,o-z)
c$$$      integer n,match(n),nef
c$$$      double precision x(n),y(n),w(n),spar,dof,smo(n),s0,lev(n)
c$$$      integer nk,ldnk,ld4
c$$$      logical ifcov
      subroutine splsm2(y, w, n, match, nef, spar, dof, smo, s0, 
     &    lev, xin, yin, win, knot, coef, sout, levout, xwy, 
     &    hs0, hs1, hs2, hs3, sg0, sg1, sg2, sg3, abd, p1ip, p2ip, 
     &    ld4, ldnk, nk)
      implicit double precision(a-h,o-z)
      integer n,match(n),nef
      double precision y(n),w(n),spar,dof,smo(n),s0,lev(n)
      integer nk,ldnk,ld4
      double precision xin(nef+1),yin(nef+1),win(nef+1),knot(nk+4)
      double precision coef(nk), sout(nef+1), levout(nef+1), xwy(nk), 
     &   hs0(nk), hs1(nk), hs2(nk), hs3(nk), sg0(nk), sg1(nk), 
     &   sg2(nk), sg3(nk), abd(ld4, nk), p1ip(ld4, nk), p2ip(ldnk, 1)
      integer ispar,icrit,isetup,ier
      double precision lspar,uspar,tol,penalt,sumwin,dofoff,crit,xbar,
     &      dsum
      double precision wmean
      crit=0e0
      if(.not.(spar.eq.0e0))goto 23002
      ispar=0
      dofoff=0e0
      if(.not.(dof.eq.0e0))goto 23004
      icrit=2 
      goto 23005
23004 continue
      dofoff=dof+1e0
      icrit=3
23005 continue
      goto 23003
23002 continue
      ispar=1
      dofoff=0e0
      icrit=3
23003 continue
      isetup=0
      ier=1
      penalt=1e0
      lspar=1e-10
      uspar=1.5
      tol=1e-3
      call sbart(penalt, dofoff, xin, yin, win, nef, knot, nk, coef, 
     &   sout, levout, crit, icrit, spar, ispar, lspar, uspar, tol, 
     &   isetup, xwy, hs0, hs1, hs2, hs3, sg0, sg1, sg2, sg3, abd, 
     &   p1ip, p2ip, ld4, ldnk, ier)
      nit=0
23006 if(.not.((ier .ne. 0) .and. (nit .lt.7)))goto 23007
      crit=0e0
      spar=0
      if(.not.(spar.eq.0e0))goto 23008
      ispar=0
      dofoff=0e0
      if(.not.(dof.eq.0e0))goto 23010
      icrit=2 
      goto 23011
23010 continue
      dofoff=dof+1e0
      icrit=3
23011 continue
      goto 23009
23008 continue
      ispar=1
      dofoff=0e0
      icrit=3
23009 continue
      isetup=0
      ier=1
      penalt=1e0
      lspar=1e-10
      tol=1e-3
      nit=nit+1
      uspar=.9+(uspar-.9)*.5
      call sbart(penalt, dofoff, xin, yin, win, nef, knot, nk, coef, 
     &     sout, levout, crit, icrit, spar, ispar, lspar, uspar, 
     &     tol, isetup, xwy, hs0, hs1, hs2, hs3, sg0, sg1, sg2, 
     &     sg3, abd, p1ip, p2ip, ld4, ldnk, ier)
      goto 23006
23007 continue
      dof=0e0
      sumwin=0e0
      do 23012 i=1,nef 
      win(i)=win(i)*win(i)
23012 continue
      sbar=wmean(nef,sout,win)
      xbar=wmean(nef,xin,win)
      do 23014 i=1,nef
      sumwin=sumwin+win(i)
23014 continue
      s0=wmean(n,y,w)
      do 23016 i=1,nef 
      lev(i)=(xin(i)-xbar)*sout(i)
23016 continue
      xsbar=wmean(nef,lev,win)
      do 23018 i=1,nef 
      lev(i)=(xin(i)-xbar)**2
23018 continue
      dsum=wmean(nef,lev,win)
      do 23020 i=1,nef 
      if(.not.(win(i).gt.0e0))goto 23022
      lev(i)=levout(i)/win(i)-1e0/sumwin -lev(i)/(sumwin*dsum)
      goto 23023
23022 continue
      lev(i)=0e0
23023 continue
23020 continue
      do 23024 i=1,nef 
      dof=dof+lev(i)*win(i)
23024 continue
      dof=dof+1e0
      do 23026 i=1,nef
      sout(i)=sout(i)-sbar -(xin(i)-xbar)*xsbar/dsum
23026 continue
      call untpack(n,nef,match,sout,smo)
      return
      end
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
