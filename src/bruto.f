      subroutine bruto(x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,dfmax,
     *cost,lambda,df,coef,type,xrange,gcvsel,gcvbak,dfit,maxit,nit,eta,
     *resid,thresh,work,iwork,trace)
      implicit double precision(a-h,o-z)
      integer n,q,p,nkmax,nk(q),match(n,q),nef(q),type(q),maxit(2),nit(
     *2),iwork(n)
      double precision x(n,q),y(n,p),w(n),knot(nkmax+4,q),wp(p),dfmax(q)
     *,cost,lambda(q),df(q),coef(nkmax*p,q),xrange(2,q),gcvsel(q,*),
     *gcvbak(q,*),dfit(q,*),eta(n,p),resid(n,p),
     *thresh,work(*)
      logical trace, select
      do 23000 j=1,p
      do 23002 i=1,n
      resid(i,j)=y(i,j)-eta(i,j)
23002 continue
23000 continue
      select=.true.
      call bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,dfmax,
     *cost,lambda,df,coef,type,xrange,gcvsel,dfit,maxit(1),nit(1),eta,
     *resid,thresh*10d0,work,iwork,trace)
      select=.false.
      call bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,dfmax,
     *cost,lambda,df,coef,type,xrange,gcvbak,dfit,maxit(2),nit(2),eta,
     *resid,thresh,work,iwork,trace)
      do 23004 j=1,p
      do 23006 i=1,n
      eta(i,j)=y(i,j)-resid(i,j)
23006 continue
23004 continue
      return
      end
      subroutine bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,
     *dfmax,cost,lambda,df,coef,type,xrange,gcv,dfit,maxit,nit,s,resid,
     *thresh,work,iwork,trace)
      implicit double precision(a-h,o-z)
      integer n,q,p,nkmax,nk(q),match(n,q),nef(q),type(q),maxit,nit,ier,
     *ntype,iwork(n)
      double precision x(n,q),y(n,p),w(n),knot(nkmax+4,q),wp(p),dfmax(q)
     *,cost,lambda(q),df(q),coef(nkmax*p,q),xrange(2,q),gcv(q,maxit),
     *dfit(q,maxit),s(n,p),resid(n,p),thresh,work(*)
      double precision dfoff, gcv0, ndf,gcv1,gcvrat,ndfoff,wmean,sbar,
     *rss,tol
      logical center, select, trace
      center=.true.
      tol=1d-3
      rss=0d0
      do 23008 j=1,p
      sbar=wmean(n,resid(1,j),w)
      do 23010 i=1,n
      resid(i,j)=resid(i,j) - sbar
      work(i)=resid(i,j)**2
23010 continue
      rss=rss+wp(j)*wmean(n,work,w)
23008 continue
      dfoff=0
      do 23012 k=1,q 
      dfoff=dfoff+df(k)
23012 continue
      gcv1=rss/((1-(1+dfoff*cost)/n)**2)
      gcvrat=1d0
      nit=0
23014 if(.not.(nit.lt.maxit.and.gcvrat .gt.thresh ))goto 23015
      gcv0=gcv1
      nit=nit+1
      do 23016 k=1,q
      gcv(k,nit)=gcv1
      if(.not.(.not.select.and.type(k).eq.1))goto 23018
      goto 23016
23018 continue
      if(.not.(type(k).gt.1))goto 23020
      call psspl2(x(1,k),n,p,knot(1,k),nk(k),xrange(1,k),coef(1,k),coef(
     *1,k),s,0,type(k))
      do 23022 j=1,p
      sbar=wmean(n,s(1,j),w)
      do 23024 i=1,n
      resid(i,j)=resid(i,j) + s(i,j)-sbar
23024 continue
23022 continue
23020 continue
      ndfoff=dfoff-df(k)
      if(.not.(select))goto 23026
      ntype=0
      goto 23027
23026 continue
      ntype=type(k)
23027 continue
      call sspl2(x(1,k),resid,w,n,p,knot(1,k),nk(k),wp,match(1,k),nef(k)
     *,ndfoff,dfmax(k),cost,lambda(k),ndf,gcv1,coef(1,k),s,ntype,center,
     *xrange(1,k),work,iwork,tol,ier)
      gcv(k,nit)=gcv1
      if(.not.(select))goto 23028
      dfit(k,nit)=ndf
      df(k)=ndf
      dfoff=ndfoff+ndf
      type(k)=ntype
23028 continue
      if(.not.(type(k).gt.1))goto 23030
      do 23032 j=1,p
      do 23034 i=1,n
      resid(i,j)=resid(i,j) - s(i,j)
23034 continue
23032 continue
23030 continue
23016 continue
      gcvrat=dabs(gcv1-gcv0)/gcv0
      if(.not.(trace))goto 23036
      call intpr("outer iteration",15,nit,1)
      call dblepr("gcv  ",5,gcv1,1)
      call dblepr("ratio",5,gcvrat,1)
23036 continue
      goto 23014
23015 continue
      return
      end
      subroutine pbruto(x,n,q,ybar,p,knot,nkmax,nk,coef,type,xrange,eta,
     *work)
      implicit double precision(a-h,o-z)
      integer n,q,p,nkmax,nk(q),type(q)
      double precision x(n,q),ybar(p),knot(nkmax+4,q),coef(nkmax*p,q),
     *xrange(2,q),eta(n,p),work(n,p)
      do 23038 j=1,p
      do 23040 i =1,n 
      eta(i,j)=ybar(j)
23040 continue
23038 continue
      do 23042 k=1,q
      if(.not.(type(k).eq.1))goto 23044
      goto 23042
23044 continue
      call psspl2(x(1,k),n,p,knot(1,k),nk(k),xrange(1,k),coef(1,k),coef(
     *1,k),work,0,type(k))
      do 23046 j=1,p
      do 23048 i=1,n
      eta(i,j)=eta(i,j) + work(i,j)
23048 continue
23046 continue
23042 continue
      return
      end
