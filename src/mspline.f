      subroutine sspl(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,dfoff,
     *dfmax,cost,lambda,df,cv,gcv,coef,s,lev,xwy,hs,sg,abd,p1ip,ier)
      implicit double precision(a-h,o-z)
      integer n,p,ldy,nk,method,ier
      double precision x(n),y(ldy,p),w(n),knot(nk+4),tol,wp(p),ssy(p),
     *dfoff,dfmax,cost,lambda,df,cv,gcv,coef(nk,p),s(ldy,p),lev(n),xwy(
     *nk,p),hs(nk,4),sg(nk,4),abd(4,nk),p1ip(4,nk)
      call sgram(sg(1,1),sg(1,2),sg(1,3),sg(1,4),knot,nk)
      call stxwx2(x,y,w,n,ldy,p,knot,nk,xwy,hs(1,1),hs(1,2),hs(1,3),hs(
     *1,4))
      if(.not.(method.eq.1))goto 23000
      call sslvr2(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,dfoff,cost,
     *lambda,df,cv,gcv,coef,s,lev,xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,4),
     *sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
      goto 23001
23000 continue
      call fmm(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,dfoff,cost,
     *lambda,df,cv,gcv,coef,s,lev,xwy,hs,sg,abd,p1ip,ier)
      if(.not.(method.gt.2.and.df.gt.dfmax))goto 23002
      df=dfmax
      call fmm(x,y,w,n,ldy,p,knot,nk,2,tol,wp,ssy,dfoff,cost,lambda,df,
     *cv,gcv,coef,s,lev,xwy,hs,sg,abd,p1ip,ier)
23002 continue
23001 continue
      return
      end
      subroutine fmm(xs,ys,ws,n,ldy,nvar,knot,nk,method,tol,wp,ssy,
     *dfoff,cost,lambda,df,cv,gcv,coef,s,lev,xwy,hs,sg,abd,p1ip,ier)
      integer n,ldy,nvar,nk,method,ier
      double precision xs(n),ys(ldy,nvar),ws(n),knot(nk+4),tol,wp(nvar),
     *ssy(nvar),dfoff,cost,lambda,df,cv,gcv,coef(nk,nvar),s(ldy,nvar),
     *lev(n),xwy(nk,nvar),hs(nk,4),sg(nk,4),abd(4,nk),p1ip(4,nk)
      double precision t1,t2,ratio, a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,
     *v,w, fu,fv,fw,fx,x,targdf,ax,bx
      integer i
      ax=1d-10
      bx=1.5
      t1=0. 
      t2=0.
      targdf=df
      do 23004 i=3,nk-3 
      t1 = t1 + hs(i,1) 
23004 continue
      do 23006 i=3,nk-3 
      t2 = t2 + sg(i,1) 
23006 continue
      ratio = t1/t2
      c = 0.5*(3. - dsqrt(5d0))
      eps = 1d0
10    eps = eps/2d0
      tol1 = 1d0 + eps
      if(.not.(tol1 .gt. 1d0))goto 23008
      go to 10
23008 continue
      eps = dsqrt(eps)
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0
      lambda = ratio*16.**(-2. + x*(6.))
      call sslvr2(xs,ys,ws,n,ldy,nvar,knot,nk,method,tol,wp,ssy,dfoff,
     *cost,lambda,df,cv,gcv,coef,s,lev,xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,
     *4),sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
      I23010 = (method)
      if(.not.(I23010.eq.( 2)))goto 23011
      fx=3d0+(targdf-df)**2
      goto 23010
23011 continue
      if(.not.(I23010.eq.( 3)))goto 23012
      fx=gcv
      goto 23010
23012 continue
      if(.not.(I23010.eq.( 4)))goto 23013
      fx=cv
23013 continue
23010 continue
      fv = fx
      fw = fx
20    xm = 0.5*(a + b)
      tol1 = eps*dabs(x) + tol/3d0
      tol2 = 2d0*tol1
      if(.not.(dabs(x - xm) .le. (tol2 - 0.5*(b - a))))goto 23014
      go to 90
23014 continue
      if(.not.(dabs(e) .le. tol1))goto 23016
      go to 40
23016 continue
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.00*(q - r)
      if(.not.(q .gt. 0.0))goto 23018
      p = -p
23018 continue
      q = dabs(q)
      r = e
      e = d
30    if(.not.(dabs(p) .ge. dabs(0.5*q*r)))goto 23020
      go to 40
23020 continue
      if(.not.(p .le. q*(a - x)))goto 23022
      go to 40
23022 continue
      if(.not.(p .ge. q*(b - x)))goto 23024
      go to 40
23024 continue
      d = p/q
      u = x + d
      if(.not.((u - a) .lt. tol2))goto 23026
      d = dsign(tol1, xm - x)
23026 continue
      if(.not.((b - u) .lt. tol2))goto 23028
      d = dsign(tol1, xm - x)
23028 continue
      go to 50
40    if(.not.(x .ge. xm))goto 23030
      e = a - x
23030 continue
      if(.not.(x .lt. xm))goto 23032
      e = b - x
23032 continue
      d = c*e
50    if(.not.(dabs(d) .ge. tol1))goto 23034
      u = x + d
23034 continue
      if(.not.(dabs(d) .lt. tol1))goto 23036
      u = x + dsign(tol1, d)
23036 continue
      lambda = ratio*16.**(-2. + u*(6.))
      call sslvr2(xs,ys,ws,n,ldy,nvar,knot,nk,method,tol,wp,ssy,dfoff,
     *cost,lambda,df,cv,gcv,coef,s,lev,xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,
     *4),sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
      I23038 = (method)
      if(.not.(I23038.eq.( 2)))goto 23039
      fu=3d0+(targdf-df)**2
      goto 23038
23039 continue
      if(.not.(I23038.eq.( 3)))goto 23040
      fu=gcv
      goto 23038
23040 continue
      if(.not.(I23038.eq.( 4)))goto 23041
      fu=cv
23041 continue
23038 continue
      if(.not.(fu .gt. fx))goto 23042
      go to 60
23042 continue
      if(.not.(u .ge. x))goto 23044
      a = x
23044 continue
      if(.not.(u .lt. x))goto 23046
      b = x
23046 continue
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
60    if(.not.(u .lt. x))goto 23048
      a = u
23048 continue
      if(.not.(u .ge. x))goto 23050
      b = u
23050 continue
      if(.not.(fu .le. fw))goto 23052
      go to 70
23052 continue
      if(.not.(w .eq. x))goto 23054
      go to 70
23054 continue
      if(.not.(fu .le. fv))goto 23056
      go to 80
23056 continue
      if(.not.(v .eq. x))goto 23058
      go to 80
23058 continue
      if(.not.(v .eq. w))goto 23060
      go to 80
23060 continue
      go to 20
70    v = w
      fv = fw
      w = u
      fw = fu
      go to 20
80    v = u
      fv = fu
      go to 20
90    continue
      if(.not.(method.eq.2))goto 23062
      call sslvr2(xs,ys,ws,n,ldy,nvar,knot,nk,1,tol,wp,ssy,dfoff,cost,
     *lambda,df,cv,gcv,coef,s,lev,xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,4),
     *sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
23062 continue
      return
      end
      subroutine stxwx2(x,z,w,k,ldy,pz,xknot,n,y,hs0,hs1,hs2,hs3)
      implicit double precision(a-h,o-z)
      integer k,n,pz,ldy,j,i,pp,ileft,mflag
      double precision z(ldy,pz),w(k),x(k),xknot(n+4),y(n,pz),hs0(n),
     *hs1(n),hs2(n),hs3(n),eps,vnikx(4,1),work(16)
      do 23064 i=1,n 
      hs0(i)=0d0
      hs1(i)=0d0
      hs2(i)=0d0
      hs3(i)=0d0
      do 23066 j=1,pz 
      y(i,j)=0d0 
23066 continue
23064 continue
      eps = .1d-9
      do 23068 i=1,k 
      call interv(xknot(1),(n+1),x(i),ileft,mflag)
      if(.not.(mflag.eq. 1))goto 23070
      if(.not.(x(i).le.(xknot(ileft)+eps)))goto 23072
      ileft=ileft-1
      goto 23073
23072 continue
      return
23073 continue
23070 continue
      call bsplvd (xknot,n+8,4,x(i),ileft,work,vnikx,1)
      j= ileft-4+1
      do 23074 pp=1,pz 
      y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(1,1)
23074 continue
      hs0(j)=hs0(j)+w(i)*vnikx(1,1)**2
      hs1(j)=hs1(j)+w(i)*vnikx(1,1)*vnikx(2,1)
      hs2(j)=hs2(j)+w(i)*vnikx(1,1)*vnikx(3,1)
      hs3(j)=hs3(j)+w(i)*vnikx(1,1)*vnikx(4,1)
      j= ileft-4+2
      do 23076 pp=1,pz 
      y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(2,1)
23076 continue
      hs0(j)=hs0(j)+w(i)*vnikx(2,1)**2
      hs1(j)=hs1(j)+w(i)*vnikx(2,1)*vnikx(3,1)
      hs2(j)=hs2(j)+w(i)*vnikx(2,1)*vnikx(4,1)
      j= ileft-4+3
      do 23078 pp=1,pz 
      y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(3,1)
23078 continue
      hs0(j)=hs0(j)+w(i)*vnikx(3,1)**2
      hs1(j)=hs1(j)+w(i)*vnikx(3,1)*vnikx(4,1)
      j= ileft-4+4
      do 23080 pp=1,pz 
      y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(4,1)
23080 continue
      hs0(j)=hs0(j)+w(i)*vnikx(4,1)**2
23068 continue
      return
      end
      subroutine sslvr2(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,dfoff,
     *cost,lambda,df,cv,gcv,coef,sz,lev,xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,
     *sg3,abd,p1ip,info)
      implicit double precision(a-h,o-z)
      integer n,p,ldy,nk,method,info
      double precision x(n),y(ldy,p),w(n),knot(nk+4),tol,wp(p),ssy(p),
     *dfoff,cost,lambda,df,cv,gcv,coef(nk,p),sz(ldy,p),lev(n),xwy(nk,p),
     *hs0(nk),hs1(nk),hs2(nk),hs3(nk),sg0(nk),sg1(nk),sg2(nk),sg3(nk),
     *abd(4,nk),p1ip(4,nk)
      double precision b0,b1,b2,b3,eps,vnikx(4,1),work(16),xv,bvalue,
     *rss,tssy
      integer ld4,i,icoef,ileft,ilo,j,mflag
      logical fittoo
      fittoo= (method.ne.2)
      ilo = 1 
      eps = .1d-10 
      ld4=4
      if(.not.(fittoo))goto 23082
      do 23084 i=1,nk 
      do 23086 j=1,p 
      coef(i,j) = xwy(i,j) 
23086 continue
23084 continue
23082 continue
      do 23088 i=1,nk 
      abd(4,i) = hs0(i)+lambda*sg0(i) 
23088 continue
      do 23090 i=1,(nk-1) 
      abd(3,i+1) = hs1(i)+lambda*sg1(i) 
23090 continue
      do 23092 i=1,(nk-2) 
      abd(2,i+2) = hs2(i)+lambda*sg2(i) 
23092 continue
      do 23094 i=1,(nk-3) 
      abd(1,i+3) = hs3(i)+lambda*sg3(i) 
23094 continue
      call dpbfa(abd,ld4,nk,3,info)
      if(.not.(info.ne.0))goto 23096
      return
23096 continue
      if(.not.(fittoo))goto 23098
      do 23100 j=1,p
      call dpbsl(abd,ld4,nk,3,coef(1,j)) 
23100 continue
      icoef = 1
      do 23102 i=1,n 
      xv = x(i)
      do 23104 j=1,p
      sz(i,j) = bvalue(knot,nk+8,coef(1,j),nk,4,xv,0) 
23104 continue
23102 continue
23098 continue
      call sinrp2(abd,ld4,nk,p1ip)
      do 23106 i=1,n 
      xv = x(i)
      call interv(knot(1),(nk+1),xv,ileft,mflag)
      if(.not.(mflag.eq.-1))goto 23108
      ileft = 4 
      xv = knot(4)+eps 
23108 continue
      if(.not.(mflag.eq.1))goto 23110
      ileft = nk 
      xv = knot(nk+1)-eps 
23110 continue
      j=ileft-3
      call bsplvd(knot,nk+8,4,xv,ileft,work,vnikx,1)
      b0=vnikx(1,1)
      b1=vnikx(2,1)
      b2=vnikx(3,1)
      b3=vnikx(4,1)
      lev(i) = (p1ip(4,j)*b0**2 + 2.*p1ip(3,j)*b0*b1 +2.*p1ip(2,j)*b0*
     *b2 + 2.*p1ip(1,j)*b0*b3 +p1ip(4,j+1)*b1**2 + 2.*p1ip(3,j+1)*b1*b2 
     *+2.*p1ip(2,j+1)*b1*b3 +p1ip(4,j+2)*b2**2 + 2.*p1ip(3,j+2)*b2*b3 +
     *p1ip(4,j+3)*b3**2 )*w(i)
23106 continue
      rss = 0d0 
      df = 0d0 
      sumw=0d0
      gcv=0d0
      cv=0d0
      do 23112 i=1,n 
      df = df + lev(i)
23112 continue
      if(.not.(fittoo))goto 23114
      do 23116 i=1,n 
      sumw = sumw + w(i)
      do 23118 j=1,p
      rss = rss + w(i)*wp(j)*(y(i,j)-sz(i,j))**2
      cv = cv +w(i)*wp(j)*((y(i,j)-sz(i,j))/(1-lev(i)))**2
23118 continue
23116 continue
      tssy=0d0
      do 23120 j=1,p
      tssy=tssy+wp(j)*ssy(j)
23120 continue
      gcv=((rss+tssy)/sumw)/((1d0-((dfoff+df-1)*cost+1)/sumw)**2)
      cv=(cv+tssy)/sumw
23114 continue
      return
      end
      subroutine sinrp2(abd,ld4,nk,p1ip)
      implicit double precision(a-h,o-z)
      integer ld4,nk,i,j,k
      double precision abd(ld4,nk),p1ip(ld4,nk), wjm3(3),wjm2(2),wjm1(1)
     *,c0,c1,c2,c3
      wjm3(1)=0d0
      wjm3(2)=0d0
      wjm3(1)=0d0
      wjm2(1)=0d0
      wjm2(2)=0d0
      wjm1(1)=0d0
      do 23122 i=1,nk 
      j=nk-i+1
      c0 = 1d0/abd(4,j)
      if(.not.(j.le.nk-3))goto 23124
      c1 = abd(1,j+3)*c0
      c2 = abd(2,j+2)*c0
      c3 = abd(3,j+1)*c0
      goto 23125
23124 continue
      if(.not.(j.eq.nk-2))goto 23126
      c1 = 0d0
      c2 = abd(2,j+2)*c0
      c3 = abd(3,j+1)*c0
      goto 23127
23126 continue
      if(.not.(j.eq.nk-1))goto 23128
      c1 = 0d0
      c2 = 0d0
      c3 = abd(3,j+1)*c0
      goto 23129
23128 continue
      if(.not.(j.eq.nk))goto 23130
      c1 = 0d0
      c2 = 0d0
      c3 = 0d0
23130 continue
23129 continue
23127 continue
23125 continue
      p1ip(1,j) = 0d0- (c1*wjm3(1)+c2*wjm3(2)+c3*wjm3(3))
      p1ip(2,j) = 0d0- (c1*wjm3(2)+c2*wjm2(1)+c3*wjm2(2))
      p1ip(3,j) = 0d0- (c1*wjm3(3)+c2*wjm2(2)+c3*wjm1(1))
      p1ip(4,j) = c0**2 +c1**2*wjm3(1)+2.*c1*c2*wjm3(2)+2.*c1*c3*wjm3(3)
     * +c2**2*wjm2(1)+2.*c2*c3*wjm2(2) +c3**2*wjm1(1)
      wjm3(1)=wjm2(1) 
      wjm3(2)=wjm2(2) 
      wjm3(3)=p1ip(2,j)
      wjm2(1)=wjm1(1) 
      wjm2(2)=p1ip(3,j)
      wjm1(1)=p1ip(4,j)
23122 continue
      return
      end
      subroutine suff2(n,p,ny,match,y,w,ybar,wbar,ssy,work)
      integer match(n),n,ny,p,i,j
      double precision y(n,ny),ybar(p+1,ny),w(n),wbar(p+1),ssy(ny),work(
     *n)
      double precision tsum
      call tpack(n,p,match,w,wbar)
      do 23132 j=1,ny
      do 23134 i=1,n
      work(i)=y(i,j)*w(i)
23134 continue
      call tpack(n,p,match,work,ybar(1,j))
      do 23136 i=1,p
      if(.not.(wbar(i).gt.0d0))goto 23138
      ybar(i,j)=ybar(i,j)/wbar(i)
      goto 23139
23138 continue
      ybar(i,j)=0d0
23139 continue
23136 continue
      call untpack(n,p,match,ybar(1,j),work)
      tsum=0d0
      do 23140 i=1,n
      tsum=tsum+ w(i)*(y(i,j)-work(i))**2
23140 continue
      ssy(j)=tsum
23132 continue
      return
      end
      subroutine sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,
     *dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,xrange,work,iwork,
     *ier)
      integer n,p,nk,method,ier,nef, nefp1, n2,match(*),iwork(n)
      double precision x(n),y(n,p),w(n),knot(nk+4),tol,wp(p),dfoff,
     *dfmax,cost,lambda,df,cv,gcv,coef(*),s(n,p),lev(nef),xrange(2),
     *work(*)
      logical center
      double precision xmiss,sigtol,xdiff,temp
      if(.not.(nef.eq.0))goto 23142
      xmiss=1d20
      sigtol=1d-5
      call namat(x,match,n,nef,work,iwork,xmiss,sigtol)
      xrange(1)=work(1)
      xrange(2)=work(nef)
      goto 23143
23142 continue
      do 23144 i=1,n 
      work(match(i))=x(i)
23144 continue
23143 continue
      xdiff=xrange(2)-xrange(1)
      do 23146 i=1,nef 
      work(i)=(work(i)-xrange(1))/xdiff
23146 continue
      if(.not.(nk.eq.0))goto 23148
      call sknotl(work,nef,knot,nk)
      nk=nk-4
23148 continue
      if(.not.(dfmax .gt. dble(nk)))goto 23150
      dfmax=dble(nk)
23150 continue
      if(.not.(cost.gt.0))goto 23152
      if(center) then
         icenter = 1
      else
         icenter = 0
      endif
      temp=dble(n-icenter)/cost - dfoff
      if(.not.(dfmax.gt.temp))goto 23154
      dfmax=temp
23154 continue
23152 continue
      nefp1=nef+1
      n2=nefp1*(2*p+2)+1
      call sspl1(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,nefp1,center,
     *dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,work(1),work(nefp1+1)
     *,work(nefp1*(p+1)+1),work(nefp1*(p+2)+1),work(n2),work(n2+p*nk),
     *work(n2+(p+4)*nk),work(n2+(p+8)*nk),work(n2+(p+12)*nk),work(n2+(p+
     *16)*nk),work(n2+(p+16)*nk+p),ier)
      return
      end
      subroutine sspl1(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,nefp1,
     *center,dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,xin,yin,win,
     *sout,xwy,hs,sg,abd,p1ip,ssy,work,ier)
      integer n,p,nefp1,nk,method,ier,match(n),nef
      double precision x(n),y(n,p),w(n),knot(nk+4),tol,wp(p),dfoff,
     *dfmax,cost,lambda,df,cv,gcv,coef(nk,p),s(n,p),lev(nef),xin(nefp1),
     *yin(nefp1,p),win(nefp1),sout(nefp1,p),xwy(nk,p),hs(nk,4),sg(nk,4),
     *abd(4,nk),p1ip(4,nk),ssy(p),work(n)
      logical center
      double precision sbar, wmean,tdfoff
      call suff2(n,nef,p,match,y,w,yin,win,ssy,work)
      if(.not.(center))goto 23156
      if(.not.(cost.gt.0))goto 23158
      tdfoff=dfoff-1/cost
23158 continue
23156 continue
      call sspl(xin,yin,win,nef,nefp1,p,knot,nk,method,tol,wp,ssy,
     *tdfoff,dfmax,cost,lambda,df,cv,gcv,coef,sout,lev,xwy,hs,sg,abd,
     *p1ip,ier)
      do 23160 j=1,p
      call untpack(n,nef,match,sout(1,j),s(1,j))
      if(.not.(center))goto 23162
      sbar=wmean(nef,sout(1,j),win)
      do 23164 i=1,n
      s(i,j)=s(i,j)-sbar
23164 continue
23162 continue
23160 continue
      if(.not.(center))goto 23166
      df=df-1
23166 continue
      return
      end
      subroutine namat(x,match,n,nef,work,iwork,xmiss,tol)
      implicit double precision(a-h,o-z)
      integer match(*),n,nef,iwork(n),index
      double precision x(*),xmiss,work(*),tol,xend,xstart
      do 23168 i=1,n 
      work(i)=x(i)
      iwork(i)=i
23168 continue
      call sortdi(work,n,iwork,1,n)
      xstart=x(iwork(1))
      index=n
      xend=x(iwork(n))
23170 if(.not.(xend .ge. xmiss .and. index .gt. 1))goto 23171
      index=index-1
      xend=x(iwork(index))
      goto 23170
23171 continue
      tol=tol*(xend-xstart)
      index=1
      work(1)=xstart
      i=1
23172 if(.not.(i.le.n))goto 23174
23175 if(.not.((x(iwork(i))-xstart).lt.tol))goto 23176
      match(iwork(i))=index
      i=i+1
      if(.not.(i.gt.n))goto 23177
      goto 10
23177 continue
      goto 23175
23176 continue
      xstart= x(iwork(i))
      index=index+1
      match(iwork(i))=index
      work(index)=xstart
      i=i+1
      goto 23172
23174 continue
10    if(.not.(xstart .ge. xmiss))goto 23179
      nef=index-1
      goto 23180
23179 continue
      nef=index
23180 continue
      return
      end
      subroutine simfit(x,y,w,n,p,dfoff,cost,wp,gcv,coef,s,type,center,
     *work)
      implicit double precision (a-h,o-z)
      integer n,p,type
      double precision x(n),y(n,p),w(n),cost,dfoff,wp(p),gcv,coef(2,p),
     *s(n,p),work(p)
      logical center
      double precision sx,sy,sumw, dcent,sxx,syy,sxy
      if(center) then
         icenter = 1
      else
         icenter = 0
      endif
      dcent=1-dble(icenter)
      sumw=0d0
      gcvc=0d0
      do 23181 i=1,n 
      sumw=sumw+w(i)
23181 continue
      if(.not.(type.ne.1))goto 23183
      sx=0.0 
      sxx=0.0
      gcvl=0d0
      do 23185 i=1,n 
      sx=sx+w(i)*x(i)
23185 continue
      xbar=sx/sumw
      do 23187 i=1,n 
      sxx=sxx+w(i)*(x(i)-xbar)*x(i)
23187 continue
23183 continue
      do 23189 j=1,p
      sy=0d0
      syy=0d0
      do 23191 i=1,n
      sy=sy+w(i)*y(i,j)
23191 continue
      work(j)=sy/sumw
      do 23193 i=1,n
      syy=syy+w(i)*(y(i,j)-work(j))*y(i,j)
23193 continue
      gcvc=gcvc+wp(j)*syy
      if(.not.(type.ne.1))goto 23195
      sxy=0.0
      do 23197 i=1,n
      sxy=sxy+w(i)*(x(i)-xbar)*y(i,j)
23197 continue
      coef(2,j)=sxy/sxx
      gcvl=gcvl+wp(j)*(syy -sxy*coef(2,j))
23195 continue
23189 continue
      if(.not.(type.eq.0))goto 23199
      gcvc =gcvc/ (sumw* (1 - (dfoff*cost + dcent)/sumw)**2 )
      gcvl=gcvl/(sumw* (1 - (dcent + (dfoff +1)* cost)/sumw)**2)
      if(.not.(gcvc.le.gcvl))goto 23201
      type=1
      gcv=gcvc
      goto 23202
23201 continue
      type=2
      gcv=gcvl
23202 continue
      goto 23200
23199 continue
      if(.not.(type.eq.1))goto 23203
      gcv=gcvc/(sumw* (1 - (dfoff*cost + dcent)/sumw)**2 )
      goto 23204
23203 continue
      gcv=gcvl/(sumw* (1 - (dcent + (dfoff + 1)*cost)/sumw)**2)
23204 continue
23200 continue
      if(.not.(type.eq.1))goto 23205
      do 23207 j=1,p
      coef(1,j)=work(j)*dcent
      coef(2,j)=0d0
      do 23209 i=1,n 
      s(i,j)=coef(1,j)
23209 continue
23207 continue
      goto 23206
23205 continue
      do 23211 j=1,p
      coef(1,j)=work(j)*dcent-xbar*coef(2,j)
      do 23213 i=1,n 
      s(i,j)=coef(1,j)+coef(2,j)*x(i)
23213 continue
23211 continue
23206 continue
      return
      end
      subroutine sspl2(x,y,w,n,p,knot,nk,wp,match,nef,dfoff,dfmax,cost,
     *lambda,df,gcv,coef,s,type,center,xrange,work,iwork,tol,ier)
      integer n,p,nk,nef,type,match(n),ier,method,iwork(n)
      double precision x(n),y(n,p),w(n),knot(nk+4),wp(p),dfoff,dfmax,
     *cost,lambda,df,gcv,coef(*),s(n,p),xrange(2),work(*),tol
      double precision coef1,coef2,cv
      logical center
      if(.not.(type.eq.3))goto 23215
      method=1
      call sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,dfoff,
     *dfmax,cost,lambda,df,cv,gcv,coef,s,work(1),xrange,work(n+1),iwork,
     *ier)
      return
23215 continue
      if(.not.(type.gt.0))goto 23217
      call simfit(x,y,w,n,p,dfoff,cost,wp,gcv,coef,s,type,center,work)
      if(center) then
         icenter = 1
      else
         icenter = 0
      endif
      df=dble(type) - dble(icenter)
      return
23217 continue
      method=3
      call sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,dfoff,
     *dfmax,cost,lambda,df,cv,gcv,coef,s,work(1),xrange,work(n+1),iwork,
     *ier)
      gcv1=gcv
      call simfit(x,y,w,n,p,dfoff,cost,wp,gcv,work,work(2*p+1),type,
     *center,work((n+2)*p+1))
      if(.not.(gcv.le.gcv1))goto 23219
      if(center) then
         icenter = 1
      else
         icenter = 0
      endif
      df=dble(type) - dble(icenter)
      do 23221 j=1,p
      coef1=work(1+(j-1)*2)
      coef2=work(2+(j-1)*2)
      if(.not.(type.eq.1))goto 23223
      do 23225 i=1,n 
      s(i,j)=coef1
23225 continue
      goto 23224
23223 continue
      do 23227 i=1,n 
      s(i,j) =coef1+coef2*x(i)
23227 continue
23224 continue
      coef(1+(j-1)*2)=coef1
      coef(2+(j-1)*2)=coef2
23221 continue
      goto 23220
23219 continue
      type=3
      gcv=gcv1
23220 continue
      return
      end
      subroutine psspl2(x,n,p,knot,nk,xrange,coef,coefl,s,order,type)
      implicit double precision(a-h,o-z)
      integer n,p,nk,order, type
      double precision x(n),knot(nk+4),xrange(2),coef(nk,p),coefl(2,1),
     *s(n,p)
      double precision ytemp
      I23229 = (type)
      if(.not.(I23229.eq.( 1)))goto 23230
      do 23231 j=1,p
      if(.not.(order.ge.1))goto 23233
      ytemp=0d0
      goto 23234
23233 continue
      ytemp=coefl(1,j)
23234 continue
      do 23235 i=1,n 
      s(i,j)=ytemp
23235 continue
23231 continue
      goto 23229
23230 continue
      if(.not.(I23229.eq.( 2)))goto 23237
      if(.not.(order.ge.1))goto 23238
      do 23240 j=1,p
      if(.not.(order.eq.1))goto 23242
      ytemp=coefl(2,j)
      goto 23243
23242 continue
      ytemp=0d0
23243 continue
      do 23244 i =1,n 
      s(i,j)=ytemp
23244 continue
23240 continue
      goto 23239
23238 continue
      do 23246 j=1,p
      do 23248 i=1,n 
      s(i,j)=coefl(1,j)+coefl(2,j)*x(i)
23248 continue
23246 continue
23239 continue
      goto 23229
23237 continue
      if(.not.(I23229.eq.( 3)))goto 23250
      call psspl(x,n,p,knot,nk,xrange,coef,s,order)
23250 continue
23229 continue
      return
      end
      subroutine psspl(x,n,p,knot,nk,xrange,coef,s,order)
      implicit double precision(a-h,o-z)
      integer n,p,nk,order
      double precision x(n),knot(nk+4),xrange(2),coef(nk,p),s(n,p)
      double precision xcs,xmin,xdif, endv(2),ends(2),xends(2),stemp
      double precision bvalue
      integer where
      if(.not.(order.gt.2.or.order.lt.0))goto 23251
      return
23251 continue
      xdif=xrange(2)-xrange(1)
      xmin=xrange(1)
      xends(1)=0d0
      xends(2)=1d0
      do 23253 j=1,p
      if(.not.(order.eq.0))goto 23255
      endv(1)=bvalue(knot,nk+8,coef(1,j),nk,4,0d0,0)
      endv(2)=bvalue(knot,nk+8,coef(1,j),nk,4,1d0,0)
23255 continue
      if(.not.(order.le.1))goto 23257
      ends(1)=bvalue(knot,nk+8,coef(1,j),nk,4,0d0,1)
      ends(2)=bvalue(knot,nk+8,coef(1,j),nk,4,1d0,1)
23257 continue
      do 23259 i=1,n
      xcs=(x(i)-xmin)/xdif
      where=0
      if(.not.(xcs.lt.0d0))goto 23261
      where=1
23261 continue
      if(.not.(xcs.gt.1d0))goto 23263
      where=2
23263 continue
      if(.not.(where.gt.0))goto 23265
      I23267 = (order)
      if(.not.(I23267.eq.( 0)))goto 23268
      stemp=endv(where)+(xcs-xends(where))*ends(where)
      goto 23267
23268 continue
      if(.not.(I23267.eq.( 1)))goto 23269
      stemp=ends(where)
      goto 23267
23269 continue
      if(.not.(I23267.eq.( 2)))goto 23270
      stemp=0d0
23270 continue
23267 continue
      goto 23266
23265 continue
      stemp=bvalue(knot,nk+8,coef(1,j),nk,4,xcs,order)
23266 continue
      if(.not.(order.gt.0))goto 23271
      s(i,j)=stemp/(xdif**dble(order))
      goto 23272
23271 continue
      s(i,j)=stemp
23272 continue
23259 continue
23253 continue
      return
      end
