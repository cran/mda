      subroutine sspl(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,dfoff,
     1     dfmax,cost,lambda,df,cv,gcv,coef,s,lev,xwy,hs,sg,abd,p1ip,
     2     ier)
      implicit double precision(a-h,o-z)
C A Cubic B-spline Smoothing routine.
C
C          The algorithm minimises:
C
C      (1/n) * sum w(i)* (y(i)-s(i))**2 + lambda* int ( s"(x) )**2 dx
C
C	for each of p response variables in y
C Input
C  x(n)		vector containing the ordinates of the observations
C  y(ldy,p)	matrix (n x p) of responses (ldy can be greater than n)
C  w(n)		vector containing the weights given to each data point
C  n    	number of data points
C  ldy		leading dimension of y
C  p		number of columns in y
C  knot(nk+4)	vector of knot points defining the cubic b-spline basis.
C  nk		number of b-spline coefficients to be estimated
C			nk <= n+2
C  method 	method for selecting amount of smoothing, lambda
C		1 = fixed lambda
C		2 = fixed df
C		3 = gcv
C		4 = cv
C  tol		used in Golden Search routine
C  wp(p)	weights, length p,  used to combine cv or gcv in 3 or 4 above
C  ssy(p)	offsets for weighted sum of squares for y; can be all zero,
C			else should be the variability lost due to collapsing
C			onto unique values
C  dfoff	offset df used in gcv calculations (0 is good default)
C  dfmax	maximum value for df allowed when gcv or cv are used
C  		routine simply returns the value at dfmax if it was exceeded
C  cost		cost per df (1 is good default)
C Input/Output
C  lambda	penalised likelihood smoothing parameter
C  df		trace(S)
C Output
C  cv		omnibus cv criterion
C  gcv		omnibus gcv criterion (including penalty and offset)
C  coef(nk,p)	vector of spline coefficients
C  s(ldy,p)	matrix of smoothed y-values
C  lev(n)	vector of leverages
C Working arrays/matrix
C  xwy(nk,p)	X'Wy
C  hs(nk,4)	the diagonals of the X'WX matrix
C  sg(nk,4)   	the diagonals of the Gram matrix
C  abd(4,nk)	[ X'WX+lambda*SIGMA] in diagonal form
C  p1ip(4,nk)	inner products between columns of L inverse
C  ier          error indicator
C                  ier = 0 ___  everything fine
C                  ier = 1 ___  spar too small or too big
C                               problem in cholesky decomposition
      integer n,p,ldy,nk,method,ier
      double precision x(n),y(ldy,p),w(n),knot(nk+4),tol,wp(p),ssy(p),
     *     dfoff,dfmax,cost,lambda,df,cv,gcv,coef(nk,p),s(ldy,p),
     *     lev(n),xwy(nk,p),hs(nk,4),sg(nk,4),abd(4,nk),p1ip(4,nk)

C  Compute SIGMA, X'WX, X'WY, trace, ratio, s0, s1.
C SIGMA-> sg[]
C X'WX -> hs[]
C X'WY -> xwy[]
      call sgram(sg(1,1),sg(1,2),sg(1,3),sg(1,4),knot,nk)
      call stxwx2(x,y,w,n,ldy,p,knot,nk,xwy,hs(1,1),hs(1,2),
     *     hs(1,3),hs(1,4))
C Compute estimate
      if(method.eq.1)then
C     Value of lambda supplied
         call sslvr2(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,dfoff,cost,
     1        lambda,df,cv,gcv,coef,s,lev,xwy,
     2        hs(1,1),hs(1,2),hs(1,3),hs(1,4),
     3        sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
      else
C     Use Forsythe, Malcom and Moler routine to minimise criterion
         call fmm(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,dfoff,cost,
     *        lambda,df,cv,gcv,coef,s,lev,xwy,hs,sg,abd,p1ip,ier)

         if(method.gt.2.and.df.gt.dfmax)then
            df=dfmax
            call fmm(x,y,w,n,ldy,p,knot,nk,2,tol,wp,ssy,dfoff,cost,
     *           lambda,df,cv,gcv,coef,s,lev,xwy,hs,sg,abd,p1ip,ier)
         endif
      endif
      return
      end

      subroutine fmm(xs,ys,ws,n,ldy,nvar,knot,nk,method,tol,wp,ssy,
     1     dfoff,cost,lambda,df,cv,gcv,coef,s,lev,xwy,hs,sg,abd,p1ip,
     2     ier)
      integer n,ldy,nvar,nk,method,ier
      double precision xs(n),ys(ldy,nvar),ws(n),knot(nk+4),tol,wp(nvar),
     *     ssy(nvar),dfoff,cost,lambda,df,cv,gcv,coef(nk,nvar),
     *     s(ldy,nvar),lev(n),xwy(nk,nvar),hs(nk,4),sg(nk,4),abd(4,nk),
     *     p1ip(4,nk)
      double precision t1,t2,ratio, a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,
     *     v,w, fu,fv,fw,fx,x,targdf,ax,bx
      integer i
C used to be lspar
      ax=1d-10
C used to be uspar
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
C
C  an approximation  x  to the point where f  attains a minimum  on
C  the interval  (ax,bx)  is determined.
C
C
C  input..
C
C  ax	 left endpoint of initial interval
C  bx	 right endpoint of initial interval
C  f	 function subprogram which evaluates  f(x)  for any  x
C	 in the interval  (ax,bx)
C  tol	 desired length of the interval of uncertainty of the final
C	 result ( .ge. 0.0)
C
C
C  output..
C
C  fmin  abcissa approximating the point where	f  attains a minimum
C
C
C      the method used is a combination of  golden  section  search  and
C  successive parabolic interpolation.	convergence is never much slower
C  than  that  for  a  fibonacci search.  if  f  has a continuous second
C  derivative which is positive at the minimum (which is not  at  ax  or
C  bx),  then  convergence  is	superlinear, and usually of the order of
C  about  1.324....
C      the function  f	is never evaluated at two points closer together
C  than  eps*dabs(fmin) + (tol/3), where eps is	approximately the square
C  root  of  the  relative  machine  precision.   if   f   is a unimodal
C  function and the computed values of	 f   are  always  unimodal  when
C  separated by at least  eps*dabs(x) + (tol/3), then  fmin  approximates
C  the abcissa of the global minimum of  f  on the interval  ax,bx  with
C  an error less than  3*eps*dabs(fmin) + tol.  if   f	is not unimodal,
C  then fmin may approximate a local, but perhaps non-global, minimum to
C  the same accuracy.
C      this function subprogram is a slightly modified	version  of  the
C  algol  60 procedure	localmin  given in richard brent, algorithms for
C  minimization without derivatives, prentice - hall, inc. (1973).
C
C
C      double precision  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
C      double precision  fu,fv,fw,fx,x
C
C  c is the squared inverse of the golden ratio
C
      c = 0.5*(3. - dsqrt(5d0))
C
C  eps is approximately the square root of the relative machine
C  precision.
C
      eps = 1d0
 10   eps = eps/2d0
      tol1 = 1d0 + eps
      if(tol1 .gt. 1d0)then
      go to 10
      endif
      eps = dsqrt(eps)
C
C  initialization
C
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0
      lambda = ratio*16.**(-2. + x*(6.))
      call sslvr2(xs,ys,ws,n,ldy,nvar,knot,nk,method,tol,wp,ssy,dfoff,
     1     cost,lambda,df,cv,gcv,coef,s,lev,xwy,
     2     hs(1,1),hs(1,2),hs(1,3),hs(1,4),
     3     sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
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
C
C  main loop starts here
C
20    xm = 0.5*(a + b)
      tol1 = eps*dabs(x) + tol/3d0
      tol2 = 2d0*tol1
C
C  check stopping criterion
C
      if(dabs(x - xm) .le. (tol2 - 0.5*(b - a)))then
      go to 90
C
C is golden-section necessary
C
      endif
      if(dabs(e) .le. tol1)then
         go to 40
C
C     fit parabola
C
      endif
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.00*(q - r)
      if(q .gt. 0.0)then
         p = -p
      endif
      q = dabs(q)
      r = e
      e = d
C
C is parabola acceptable
C
30    if(dabs(p) .ge. dabs(0.5*q*r))then
         go to 40
      endif
      if(p .le. q*(a - x))then
         go to 40
      endif
      if(p .ge. q*(b - x))then
         go to 40
      endif
C
C  a parabolic interpolation step
C
      d = p/q
      u = x + d
C
C  f must not be evaluated too close to ax or bx
C
      if((u - a) .lt. tol2)then
      d = dsign(tol1, xm - x)
      endif
      if((b - u) .lt. tol2)then
      d = dsign(tol1, xm - x)
      endif
      go to 50
C
C  a golden-section step
C
40    if(x .ge. xm)then
      e = a - x
      endif
      if(x .lt. xm)then
      e = b - x
      endif
      d = c*e
C
C  f must not be evaluated too close to x
C
50    if(dabs(d) .ge. tol1)then
      u = x + d
      endif
      if(dabs(d) .lt. tol1)then
      u = x + dsign(tol1, d)
      endif
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
C
C  update  a, b, v, w, and x
C
      if(fu .gt. fx)then
      go to 60
      endif
      if(u .ge. x)then
      a = x
      endif
      if(u .lt. x) then
      b = x
      endif
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
60    if(u .lt. x)then
      a = u
      endif
      if(u .ge. x)then
      b = u
      endif
      if(fu .le. fw)then
      go to 70
      endif
      if(w .eq. x)then
      go to 70
      endif
      if(fu .le. fv)then
      go to 80
      endif
      if(v .eq. x)then
      go to 80
      endif
      if(v .eq. w)then
      go to 80
      endif
      go to 20
70    v = w
      fv = fw
      w = u
      fw = fu
      go to 20
80    v = u
      fv = fu
      go to 20
C
C  end of main loop
C
90    continue
      if(method.eq.2)then
      call sslvr2(xs,ys,ws,n,ldy,nvar,knot,nk,1,tol,wp,ssy,dfoff,cost,
     1        lambda,df,cv,gcv,coef,s,lev,xwy,
     2        hs(1,1),hs(1,2),hs(1,3),hs(1,4),
     3        sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
      endif
      return
      end

      subroutine stxwx2(x,z,w,k,ldy,pz,xknot,n,y,hs0,hs1,hs2,hs3)
      implicit double precision(a-h,o-z)
      integer k,n,pz,ldy,j,i,pp,ileft,mflag
      double precision z(ldy,pz),w(k),x(k),xknot(n+4),y(n,pz),hs0(n),
     *     hs1(n),hs2(n),hs3(n),eps,vnikx(4,1),work(16)
C Initialise the output vectors
      do 23064 i=1,n
         hs0(i)=0d0
         hs1(i)=0d0
         hs2(i)=0d0
         hs3(i)=0d0
         do 23066 j=1,pz
            y(i,j)=0d0
23066    continue
23064 continue
C Compute X'WX -> hs0,hs1,hs2,hs3  and X'WZ -> y
      ileft = 1
      eps = .1d-9
      do 23068 i=1,k
         ileft = interv(xknot(1),(n+1),x(i), 0,0, ileft,mflag)
         if(mflag.eq. 1)then
            if(x(i).le.(xknot(ileft)+eps))then
               ileft=ileft-1
            else
               return
            endif
         endif
         call bsplvd (xknot,n+8,4,x(i),ileft,work,vnikx,1)
         j= ileft-4+1
         do 23074 pp=1,pz
            y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(1,1)
23074    continue
         hs0(j)=hs0(j)+w(i)*vnikx(1,1)**2
         hs1(j)=hs1(j)+w(i)*vnikx(1,1)*vnikx(2,1)
         hs2(j)=hs2(j)+w(i)*vnikx(1,1)*vnikx(3,1)
         hs3(j)=hs3(j)+w(i)*vnikx(1,1)*vnikx(4,1)
         j= ileft-4+2
         do 23076 pp=1,pz
            y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(2,1)
23076    continue
         hs0(j)=hs0(j)+w(i)*vnikx(2,1)**2
         hs1(j)=hs1(j)+w(i)*vnikx(2,1)*vnikx(3,1)
         hs2(j)=hs2(j)+w(i)*vnikx(2,1)*vnikx(4,1)
         j= ileft-4+3
         do 23078 pp=1,pz
            y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(3,1)
23078    continue
         hs0(j)=hs0(j)+w(i)*vnikx(3,1)**2
         hs1(j)=hs1(j)+w(i)*vnikx(3,1)*vnikx(4,1)
         j= ileft-4+4
         do 23080 pp=1,pz
            y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(4,1)
23080    continue
         hs0(j)=hs0(j)+w(i)*vnikx(4,1)**2
23068 continue
      return
      end

      subroutine sslvr2(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,dfoff,
     *     cost,lambda,df,cv,gcv,coef,sz,lev,xwy,hs0,hs1,hs2,hs3,
     *     sg0,sg1,sg2,sg3,abd,p1ip,info)
      implicit double precision(a-h,o-z)
      integer n,p,ldy,nk,method,info
      double precision x(n),y(ldy,p),w(n),knot(nk+4),tol,wp(p),ssy(p),
     *     dfoff,cost,lambda,df,cv,gcv,coef(nk,p),sz(ldy,p),lev(n),
     $     xwy(nk,p), hs0(nk),hs1(nk),hs2(nk),hs3(nk),
     $     sg0(nk),sg1(nk),sg2(nk),sg3(nk), abd(4,nk),p1ip(4,nk)
C local storage
      double precision b0,b1,b2,b3,eps,vnikx(4,1),work(16),xv,bvalue,
     *     rss,tssy
      integer ld4,i,icoef,ileft,j,mflag
      logical fittoo
      fittoo= (method.ne.2)
      ileft = 1
      eps = .1d-10
      ld4=4
C Purpose : Solves the smoothing problem and computes the
C           criterion functions (CV and GCV).
C The coeficients of estimated smooth
      if(fittoo)then
         do 23084 i=1,nk
            do 23086 j=1,p
               coef(i,j) = xwy(i,j)
23086       continue
23084    continue
      endif
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
      if(info.ne.0)then
         return
      endif
      if(fittoo)then
         do 23100 j=1,p
            call dpbsl(abd,ld4,nk,3,coef(1,j))
23100    continue
C     Value of smooths at the data points
         icoef = 1
         do 23102 i=1,n
            xv = x(i)
            do 23104 j=1,p
               sz(i,j) = bvalue(knot,nk+8,coef(1,j),nk,4,xv,0)
23104       continue
23102    continue
23098    continue
      endif
C     Compute the criterion functions
C     Get Leverages First
      call sinrp2(abd,ld4,nk,p1ip)
      do 23106 i=1,n
         xv = x(i)
         ileft = interv(knot(1),(nk+1),xv, 0,0, ileft,mflag)
         if(mflag.eq.-1)then
            ileft = 4
            xv = knot(4)+eps
         endif
         if(mflag.eq.1)then
            ileft = nk
            xv = knot(nk+1)-eps
         endif
         j=ileft-3
         call bsplvd(knot,nk+8,4,xv,ileft,work,vnikx,1)
         b0=vnikx(1,1)
         b1=vnikx(2,1)
         b2=vnikx(3,1)
         b3=vnikx(4,1)
         lev(i) = (p1ip(4,j)*b0**2 + 2.*p1ip(3,j)*b0*b1 +
     1          2.*p1ip(2,j)*b0*b2 + 2.*p1ip(1,j)*b0*b3 +
     2           p1ip(4,j+1)*b1**2 + 2.*p1ip(3,j+1)*b1*b2 +
     3        2.*p1ip(2,j+1)*b1*b3 +    p1ip(4,j+2)*b2**2 +
     4        2.*p1ip(3,j+2)*b2*b3 +    p1ip(4,j+3)*b3**2 )*w(i)
23106 continue
C     Evaluate Criteria
      rss = 0d0
      df = 0d0
      sumw=0d0
      gcv=0d0
      cv=0d0
      do 23112 i=1,n
         df = df + lev(i)
23112 continue

      if(fittoo)then
         do 23116 i=1,n
            sumw = sumw + w(i)
            do 23118 j=1,p
               rss = rss + w(i)*wp(j)*(y(i,j)-sz(i,j))**2
               cv = cv +w(i)*wp(j)*((y(i,j)-sz(i,j))/(1-lev(i)))**2
23118       continue
23116    continue
         tssy=0d0
         do 23120 j=1,p
            tssy=tssy+wp(j)*ssy(j)
23120    continue
         gcv=((rss+tssy)/sumw)/((1d0-((dfoff+df-1)*cost+1)/sumw)**2)
C note: the weights should sum to n (the number of original observations)
         cv=(cv+tssy)/sumw
C     lev includes the weights
C     Note that this version of cv omits ALL observations at
C     tied x values, since the data are already collapsed here
      endif
      return
      end

      subroutine sinrp2(abd,ld4,nk,p1ip)
C Purpose :  Computes Inner Products between columns of L(-1)
C	     where L = abd is a Banded Matrix with 3 subdiagonals
C		A refinement of Elden's trick is used.
C	Coded by Finbarr O'Sullivan

      implicit double precision(a-h,o-z)
      integer ld4,nk,i,j
      double precision abd(ld4,nk),p1ip(ld4,nk),
     *     wjm3(3),wjm2(2),wjm1(1),c0,c1,c2,c3
      wjm3(1)=0d0
      wjm3(2)=0d0
C     wjm3(1)=0d0   spotted by Brian Ripley 10/8/2015 (earlier actually)
      wjm3(3)=0d0
      wjm2(1)=0d0
      wjm2(2)=0d0
      wjm1(1)=0d0
      do 100 i=1,nk
         j=nk-i+1
         c0 = 1d0/abd(4,j)
         if(j.le.nk-3)then
            c1 = abd(1,j+3)*c0
            c2 = abd(2,j+2)*c0
            c3 = abd(3,j+1)*c0
         else if(j.eq.nk-2)then
            c1 = 0d0
            c2 = abd(2,j+2)*c0
            c3 = abd(3,j+1)*c0
         else if(j.eq.nk-1)then
            c1 = 0d0
            c2 = 0d0
            c3 = abd(3,j+1)*c0
         else if(j.eq.nk)then
            c1 = 0d0
            c2 = 0d0
            c3 = 0d0
         endif

         p1ip(1,j) = 0d0- (c1*wjm3(1)+c2*wjm3(2)+c3*wjm3(3))
         p1ip(2,j) = 0d0- (c1*wjm3(2)+c2*wjm2(1)+c3*wjm2(2))
         p1ip(3,j) = 0d0- (c1*wjm3(3)+c2*wjm2(2)+c3*wjm1(1))
         p1ip(4,j) = c0**2 +c1**2*wjm3(1)+2.*c1*c2*wjm3(2)+
     *        2.*c1*c3*wjm3(3)+ c2**2*wjm2(1)+2.*c2*c3*wjm2(2) +
     2        c3**2*wjm1(1)
         wjm3(1)=wjm2(1)
         wjm3(2)=wjm2(2)
         wjm3(3)=p1ip(2,j)
         wjm2(1)=wjm1(1)
         wjm2(2)=p1ip(3,j)
         wjm1(1)=p1ip(4,j)
 100  continue
      return
      end

      subroutine suff2(n,p,ny,match,y,w,ybar,wbar,ssy,work)
      integer match(n),n,ny,p,i,j
      double precision y(n,ny),ybar(p+1,ny),w(n),wbar(p+1),ssy(ny),
     *     work(n)
      double precision tsum

      call tpack(n,p,match,w,wbar)
      do 100 j=1,ny
         do 34 i=1,n
            work(i)=y(i,j)*w(i)
 34      continue
         call tpack(n,p,match,work,ybar(1,j))
         do 36 i=1,p
            if(wbar(i).gt.0d0)then
               ybar(i,j)=ybar(i,j)/wbar(i)
            else
               ybar(i,j)=0d0
            endif
 36      continue
         call untpack(n,p,match,ybar(1,j),work)
         tsum=0d0
         do 40 i=1,n
            tsum=tsum+ w(i)*(y(i,j)-work(i))**2
 40      continue
         ssy(j)=tsum
 100  continue
      return
      end

      subroutine sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,
     *     dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,xrange,
     $     work,iwork, ier)
      integer n,p,nk,method,ier,nef, nefp1, n2,match(*),iwork(n)
      double precision x(n),y(n,p),w(n),knot(nk+4),tol,wp(p),dfoff,
     *     dfmax,cost,lambda,df,cv,gcv,coef(*),s(n,p),lev(nef),
     *     xrange(2), work(*)
C     workspace must be (2*p+2)*nefp1 + (p+16)*nk + n +p

      logical center
      double precision xmiss,sigtol,xdiff,temp

      if(nef.eq.0)then
C     match has not been initialized
         xmiss=1d20
         sigtol=1d-5
         call namat(x,match,n,nef,work,iwork,xmiss,sigtol)
         xrange(1)=work(1)
         xrange(2)=work(nef)
      else
         do 44 i=1,n
            work(match(i))=x(i)
 44      continue
      endif
      xdiff=xrange(2)-xrange(1)
      do 46 i=1,nef
         work(i)=(work(i)-xrange(1))/xdiff
 46   continue
      if(nk.eq.0)then
         call sknotl(work,nef,knot,nk)
         nk=nk-4
      endif
      if(dfmax .gt. dble(nk))then
         dfmax=dble(nk)
      endif
      if(cost.gt.0)then
         if(center) then
            icenter = 1
         else
            icenter = 0
         endif
         temp=dble(n-icenter)/cost - dfoff
         if(dfmax.gt.temp)then
            dfmax=temp
         endif
      endif
      nefp1=nef+1
      n2=nefp1*(2*p+2)+1
      call sspl1(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,nefp1,center,
     *     dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,
     *     work(1),work(nefp1+1),
c          xin,    yin
     $     work(nefp1*(p+1)+1),work(nefp1*(p+2)+1),
c          win,                sout
     *     work(n2),work(n2+p*nk),work(n2+(p+4)*nk),
c          xwy,     hs,           sg
     $     work(n2+(p+8)*nk),work(n2+(p+12)*nk),work(n2+(p+16)*nk),
     $     work(n2+(p+16)*nk+p), ier)
      return
      end

C Memory management subroutine
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
      if(center)then
         tdfoff=dfoff
         if(cost.gt.0)then
            tdfoff=dfoff-1/cost
         endif
      endif
      call sspl(xin,yin,win,nef,nefp1,p,knot,nk,method,tol,wp,ssy,
     *     tdfoff,dfmax,cost,lambda,df,cv,gcv,coef,sout,lev,xwy,
     $     hs,sg,abd, p1ip,ier)

C now unpack the results
      do 60 j=1,p
         call untpack(n,nef,match,sout(1,j),s(1,j))
         if(center) then
            sbar=wmean(nef,sout(1,j),win)
            do 64 i=1,n
               s(i,j)=s(i,j)-sbar
 64         continue
         endif
 60   continue

      if(center) df=df-1
      return
      end

      subroutine namat(x,match,n,nef,work,iwork,xmiss,tol)
C returns match (order) and work(1:nef) is the sorted unique x
      implicit double precision(a-h,o-z)
      integer match(*),n,nef,iwork(n),index
      double precision x(*),xmiss,work(*),tol,xend,xstart

      do 68 i=1,n
         work(i)=x(i)
         iwork(i)=i
 68   continue
      call sortdi(work,n,iwork,1,n)
      xstart=x(iwork(1))
      index=n
      xend=x(iwork(n))
c     while( . )
 70   if(xend .ge. xmiss .and. index .gt. 1)then
         index=index-1
         xend=x(iwork(index))
         goto 70
      endif

      tol=tol*(xend-xstart)
      index=1
      work(1)=xstart
c     for(i = 1:n) -- not allowed to increment i inside do-loop!
      i=1
 100  if(i .le. n) then
c        while() :
 75      if((x(iwork(i))-xstart).lt.tol) then
            match(iwork(i))=index
            i=i+1
            if(i.gt.n) goto 10
c      this next line was missing; added by TH in Nov/2017!            
            goto 75
         endif
         xstart= x(iwork(i))
         index=index+1
         match(iwork(i))=index
         work(index)=xstart
         i=i+1
         goto 100
      endif

 10   if(xstart .ge. xmiss)then
         nef=index-1
      else
         nef=index
      endif
      return
      end

      subroutine simfit(x,y,w,n,p,dfoff,cost,wp,gcv,coef,s,type,center,
     *     work)
C
C computes constant and linear fits, and selects the best using gcv
C
      implicit double precision (a-h,o-z)
      integer n,p,type
      double precision x(n),y(n,p),w(n),cost,dfoff,wp(p),gcv,coef(2,p),
     *     s(n,p),work(p)
      logical center
      double precision sx,sy,sumw, dcent,sxx,syy,sxy
C center is F for no centering, else T
C Note: if type enters 1 or 2, no selection is made
      if(center) then
         icenter = 1
      else
         icenter = 0
      endif
      dcent=1-dble(icenter)
      sumw=0d0
      gcvc=0d0
      do 81 i=1,n
         sumw=sumw+w(i)
 81   continue

      if(type.ne.1)then
C     either 0 or 2 in which case the linear is needed as well
         sx=0.0
         sxx=0.0
         gcvl=0d0
         do 85 i=1,n
            sx=sx+w(i)*x(i)
 85      continue
         xbar=sx/sumw
         do 87 i=1,n
            sxx=sxx+w(i)*(x(i)-xbar)*x(i)
 87      continue
      endif

      do 100 j=1,p
         sy=0d0
         syy=0d0
         do 91 i=1,n
            sy=sy+w(i)*y(i,j)
 91      continue
         work(j)=sy/sumw
         do 93 i=1,n
            syy=syy+w(i)*(y(i,j)-work(j))*y(i,j)
 93      continue
         gcvc=gcvc+wp(j)*syy
         if(type.ne.1) then
C     once again, do for linear as well
            sxy=0.0
            do 97 i=1,n
               sxy=sxy+w(i)*(x(i)-xbar)*y(i,j)
 97         continue
            coef(2,j)=sxy/sxx
            gcvl=gcvl+wp(j)*(syy -sxy*coef(2,j))
         endif
 100  continue
      if(type.eq.0) then
         gcvc =gcvc/ (sumw* (1 - (dfoff*cost + dcent)/sumw)**2 )
         gcvl=gcvl/(sumw* (1 - (dcent + (dfoff +1)* cost)/sumw)**2)
         if(gcvc.le.gcvl) then
            type=1
            gcv=gcvc
         else
            type=2
            gcv=gcvl
         endif
      else if(type.eq.1) then
         gcv=gcvc/(sumw* (1 - (dfoff*cost + dcent)/sumw)**2 )
      else
         gcv=gcvl/(sumw* (1 - (dcent + (dfoff + 1)*cost)/sumw)**2)
      endif

      if(type.eq.1) then
         do 107 j=1,p
            coef(1,j)=work(j)*dcent
            coef(2,j)=0d0
            do 109 i=1,n
               s(i,j)=coef(1,j)
 109        continue
 107     continue
      else
         do 111 j=1,p
            coef(1,j)=work(j)*dcent-xbar*coef(2,j)
            do 113 i=1,n
               s(i,j)=coef(1,j)+coef(2,j)*x(i)
 113        continue
 111     continue
      endif
      return
      end

      subroutine sspl2(x,y,w,n,p,knot,nk,wp,match,nef,dfoff,dfmax,cost,
     *     lambda,df,gcv,coef,s,type,center,xrange,work,iwork,tol,ier)
C this routine selects from the linear and constant model as well
C see documentation for sspl
C workspace must be (2*p+2)*nefp1 + (p+16)*nk + 2*n
C if type>0 then no selection is performed; the fit is simply computed.
      integer n,p,nk,nef,type,match(n),ier,method,iwork(n)
      double precision x(n),y(n,p),w(n),knot(nk+4),wp(p),dfoff,dfmax,
     *     cost,lambda,df,gcv,coef(*),s(n,p),xrange(2),work(*),tol
      double precision coef1,coef2,cv
      logical center
Ccenter is F for no centering, else T
      if(type.eq.3)then
         method=1
         call sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,
     $        dfoff, dfmax,cost,lambda,df,cv,gcv,coef,s,work(1),
     $        xrange,work(n+1),iwork, ier)
         return
      endif

      if(type.gt.0)then
         call simfit(x,y,w,n,p,dfoff,cost,wp,gcv,coef,s,type,center,
     $        work)
         if(center) then
            icenter = 1
         else
            icenter = 0
         endif
         df=dble(type) - dble(icenter)
         return
      endif

      method=3
      call sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,dfoff,
     *     dfmax,cost,lambda,df,cv,gcv,coef,s,work(1),xrange,work(n+1),
     $     iwork, ier)
      gcv1=gcv
      call simfit(x,y,w,n,p,dfoff,cost,wp,gcv,work,work(2*p+1),type,
     *     center,work((n+2)*p+1))
      if(gcv.le.gcv1) then
C     the coef swapping is so as not to destroy the spline coefs if needed
         if(center) then
            icenter = 1
         else
            icenter = 0
         endif
         df=dble(type) - dble(icenter)
         do 30 j=1,p
            coef1=work(1+(j-1)*2)
            coef2=work(2+(j-1)*2)
            if(type.eq.1)then
               do 25 i=1,n
                  s(i,j)=coef1
 25            continue
            else
               do 27 i=1,n
                  s(i,j) =coef1+coef2*x(i)
 27            continue
            endif
            coef(1+(j-1)*2)=coef1
            coef(2+(j-1)*2)=coef2
 30      continue
      else
         type=3
         gcv=gcv1
      endif
      return
      end

      subroutine psspl2(x,n,p,knot,nk,xrange,coef,coefl,s,order,type)
      implicit double precision(a-h,o-z)
      integer n,p,nk, order, type
      double precision x(n),knot(nk+4),xrange(2),coef(nk,p),coefl(2,1),
     *     s(n,p)
c Var
      double precision ytemp

C switch(type){ case 1, 2, 3 } :
      if(type .eq. 1) then
         do 30 j=1,p
            if(order.ge.1)then
               ytemp=0d0
            else
               ytemp=coefl(1,j)
            endif
            do 35 i=1,n
               s(i,j)=ytemp
 35         continue
 30      continue

      else if(type .eq. 2) then

         if(order.ge.1)then
            do 40 j=1,p
               if(order.eq.1)then
                  ytemp=coefl(2,j)
               else
                  ytemp=0d0
               endif
               do 44 i =1,n
                  s(i,j)=ytemp
 44            continue
 40         continue
         else
            do 46 j=1,p
               do 48 i=1,n
                  s(i,j)=coefl(1,j)+coefl(2,j)*x(i)
 48            continue
 46         continue
         endif

      else if(type .eq. 3) then

         call psspl(x,n,p,knot,nk,xrange,coef,s,order)

      endif
      return
      end

      subroutine psspl(x,n,p,knot,nk,xrange,coef,s,order)
C     make predictions from a fitted smoothing spline, linear or constant
      implicit double precision(a-h,o-z)
      integer n,p,nk,order
      double precision x(n),knot(nk+4),xrange(2),coef(nk,p),s(n,p)
      double precision xcs,xmin,xdif, endv(2),ends(2),xends(2),stemp
      double precision bvalue
      integer where

      if(order.gt.2.or.order.lt.0)  return

      xdif=xrange(2)-xrange(1)
      xmin=xrange(1)
      xends(1)=0d0
      xends(2)=1d0
      do 23253 j=1,p
         if(order.eq.0)then
            endv(1)=bvalue(knot,nk+8,coef(1,j),nk,4,0d0,0)
            endv(2)=bvalue(knot,nk+8,coef(1,j),nk,4,1d0,0)
         endif
         if(order.le.1)then
            ends(1)=bvalue(knot,nk+8,coef(1,j),nk,4,0d0,1)
            ends(2)=bvalue(knot,nk+8,coef(1,j),nk,4,1d0,1)
         endif
         do 23259 i=1,n
            xcs=(x(i)-xmin)/xdif
            if(xcs.lt.0d0)then
               where=1
            else if(xcs.gt.1d0)then
               where=2
            else
               where=0
            endif

            if(where.gt.0)then
C     beyond extreme knots
c 	switch(order){ case 0: / 1 / 2  }  :
               if(order .eq. 0) then
                  stemp=endv(where)+(xcs-xends(where))*ends(where)
               else if(order .eq. 1) then
                  stemp=ends(where)
               else if(order .eq. 2) then
                  stemp=0d0
               endif
            else
               stemp=bvalue(knot,nk+8,coef(1,j),nk,4,xcs,order)
            endif

            if(order.gt.0)then
               s(i,j)=stemp/(xdif**dble(order))
            else
               s(i,j)=stemp
            endif

23259    continue
23253 continue
      return
      end
