      subroutine bruto(x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,dfmax,
     *     cost,lambda,df,coef,type,xrange,gcvsel,gcvbak,dfit,maxit,nit,
     *     eta,resid,thresh,work,iwork,trace)
      implicit double precision(a-h,o-z)
      integer n,q,p,nkmax,nk(q),match(n,q),nef(q),type(q),maxit(2),
     *     nit(2),iwork(n)
      double precision x(n,q),y(n,p),w(n),knot(nkmax+4,q),wp(p),
     *     dfmax(q),cost,lambda(q),df(q),coef(nkmax*p,q),xrange(2,q),
     *     gcvsel(q,*),gcvbak(q,*),dfit(q,*),eta(n,p),resid(n,p),
     *     thresh,work(*)
      logical trace, select

      do 01 j=1,p
         do 02 i=1,n
            resid(i,j)=y(i,j)-eta(i,j)
 02      continue
 01   continue

      select=.true.
c$$$ Naras fix
c$$$      call bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,dfmax,
c$$$     *cost,lambda,df,coef,type,xrange,gcvsel,dfit,maxit(1),nit(1),eta,
c$$$     *resid,thresh*10d0,work,iwork,trace)
c$$$      select=.false.
c$$$      call bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,dfmax,
c$$$     *cost,lambda,df,coef,type,xrange,gcvbak,dfit,maxit(2),nit(2),eta,
c$$$     *resid,thresh,work,iwork,trace)

      call bakssp(select,x,n,q,p,w,knot,nkmax,nk,wp,match,nef,dfmax,
     *cost,lambda,df,coef,type,xrange,gcvsel,dfit,maxit(1),nit(1),eta,
     *resid,thresh*10d0,work,iwork,trace)
      select=.false.
      call bakssp(select,x,n,q,p,w,knot,nkmax,nk,wp,match,nef,dfmax,
     *cost,lambda,df,coef,type,xrange,gcvbak,dfit,maxit(2),nit(2),eta,
     *resid,thresh,work,iwork,trace)
      do 04 j=1,p
         do 06 i=1,n
            eta(i,j)=y(i,j)-resid(i,j)
 06      continue
 04   continue
      return
      end

c$$$ Naras fix      
c$$$      subroutine bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,
c$$$     *     dfmax,cost,lambda,df,coef,type,xrange,gcv,dfit,maxit,nit,s,
c$$$     *     resid,thresh,work,iwork,trace)
c$$$      implicit double precision(a-h,o-z)
c$$$      integer n,q,p,nkmax,nk(q),match(n,q),nef(q),type(q),maxit,nit,ier,
c$$$     *     ntype,iwork(n)
c$$$      double precision x(n,q),y(n,p),w(n),knot(nkmax+4,q),wp(p),dfmax(q)
c$$$     *     ,cost,lambda(q),df(q),coef(nkmax*p,q),xrange(2,q),
c$$$     *     gcv(q,maxit),dfit(q,maxit),s(n,p),resid(n,p),thresh,work(*)
c$$$
      subroutine bakssp(select,x,n,q,p,w,knot,nkmax,nk,wp,match,nef,
     *     dfmax,cost,lambda,df,coef,type,xrange,gcv,dfit,maxit,nit,s,
     *     resid,thresh,work,iwork,trace)
      implicit double precision(a-h,o-z)
      integer n,q,p,nkmax,nk(q),match(n,q),nef(q),type(q),maxit,nit,ier,
     *     ntype,iwork(n)
      double precision x(n,q),w(n),knot(nkmax+4,q),wp(p),dfmax(q)
     *     ,cost,lambda(q),df(q),coef(nkmax*p,q),xrange(2,q),
     *     gcv(q,maxit),dfit(q,maxit),s(n,p),resid(n,p),thresh,work(*)
      double precision dfoff, gcv0, ndf,gcv1,gcvrat,ndfoff,wmean,sbar,
     *     rss,tol
      logical center, select, trace
      center=.true.
      tol=1d-3
      rss=0d0
      do 08 j=1,p
         sbar=wmean(n,resid(1,j),w)
         do 10 i=1,n
            resid(i,j)=resid(i,j) - sbar
            work(i)=resid(i,j)**2
 10      continue
         rss=rss+wp(j)*wmean(n,work,w)
 08   continue
      dfoff=0
      do 12 k=1,q
         dfoff=dfoff+df(k)
 12   continue
      gcv1=rss/((1-(1+dfoff*cost)/n)**2)
      gcvrat=1d0
      nit=0
 14   if(.not.(nit.lt.maxit.and.gcvrat .gt.thresh ))goto 15
      gcv0=gcv1
      nit=nit+1
      do 16 k=1,q
         gcv(k,nit)=gcv1
         if(.not.select.and.type(k).eq.1)goto 16

         if(type(k) .gt. 1) then
            call psspl2(x(1,k),n,p,knot(1,k),nk(k),xrange(1,k),
     *           coef(1,k),coef(1,k),s,0,type(k))
            do 22 j=1,p
               sbar=wmean(n,s(1,j),w)
               do 24 i=1,n
                  resid(i,j)=resid(i,j) + s(i,j)-sbar
 24            continue
 22         continue
         endif

         ndfoff=dfoff-df(k)
         if(select) then
            ntype=0
         else
            ntype=type(k)
         endif

         call sspl2(x(1,k),resid,w,n,p,knot(1,k),nk(k),wp,match(1,k),
     *        nef(k),ndfoff,dfmax(k),cost,lambda(k),ndf,gcv1,coef(1,k),
     *        s,ntype,center,xrange(1,k),work,iwork,tol,ier)
         gcv(k,nit)=gcv1
         if(select) then
            dfit(k,nit)=ndf
            df(k)=ndf
            dfoff=ndfoff+ndf
            type(k)=ntype
         endif

         if(type(k) .gt. 1) then
            do 32 j=1,p
               do 34 i=1,n
                  resid(i,j)=resid(i,j) - s(i,j)
 34            continue
 32         continue
         endif

 16   continue

      gcvrat=dabs(gcv1-gcv0)/gcv0
      if(trace) then
         call intpr("outer iteration",15,nit,1)
         call dblepr("gcv  ",5,gcv1,1)
         call dblepr("ratio",5,gcvrat,1)
      endif
      continue
      goto 14
 15   continue
      return
      end

      subroutine pbruto(x,n,q,ybar,p,knot,nkmax,nk,coef,type,xrange,eta,
     *     work)
      implicit double precision(a-h,o-z)
      integer n,q,p,nkmax,nk(q),type(q)
      double precision x(n,q),ybar(p),knot(nkmax+4,q),coef(nkmax*p,q),
     *     xrange(2,q),eta(n,p),work(n,p)

      do 38 j=1,p
         do 40 i =1,n
            eta(i,j)=ybar(j)
 40      continue
 38   continue
      do 42 k=1,q
         if(type(k) .eq .1) goto 42
         call psspl2(x(1,k),n,p,knot(1,k),nk(k),xrange(1,k),coef(1,k),
     *        coef(1,k),work,0,type(k))
         do 46 j=1,p
            do 48 i=1,n
               eta(i,j)=eta(i,j) + work(i,j)
 48         continue
 46      continue
 42   continue
      return
      end
