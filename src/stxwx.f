C Output from Public domain Ratfor, version 1.0
      subroutine stxwx(x,z,w,k,xknot,n,y,hs0,hs1,hs2,hs3)

      implicit none
      integer k,n
      DOUBLE precision x(k),z(k),w(k), xknot(n+4),y(n),
     &	   hs0(n),hs1(n),hs2(n),hs3(n)
C local
      DOUBLE precision eps,vnikx(4,1),work(16)
      integer lenxk,ilo, i,j, ileft,mflag

      lenxk=n+4
C     Initialise the output vectors
      do 1 i=1,n
	 y(i)=0d0
	 hs0(i)=0d0
	 hs1(i)=0d0
	 hs2(i)=0d0
	 hs3(i)=0d0
 1    continue

C     Compute X'WX -> hs0,hs1,hs2,hs3  and X'WZ -> y
      ilo=1
      eps= .1d-9

      do 100 i=1,k
C	 call intrv (xknot(1),(n+1),x(i),ilo,ileft,mflag)
	 call interv(xknot(1),(n+1),x(i),ileft,mflag)
C        if(mflag==-1) {write(6,'("Error in hess ",i2)')mflag;stop}
C        if(mflag==-1) {return}
	 if(mflag.eq. 1)then
	    if(x(i).le.(xknot(ileft)+eps))then
	       ileft=ileft-1
	    else
	       return
	    endif
C        else{write(6,'("Error in hess ",i2)')mflag;stop}}
	 endif
C	 call bspvd (xknot,4,1,x(i),ileft,4,vnikx,work)
	 call bsplvd (xknot,lenxk,4,x(i),ileft,work,vnikx,1)

C     check 2
C     print 2999;2999 format(" got through check2 stxwx ")
	 j= ileft-4+1
	 y(j) = y(j)+w(i)**2*z(i)*vnikx(1,1)
	 hs0(j)=hs0(j)+w(i)**2*vnikx(1,1)**2
	 hs1(j)=hs1(j)+w(i)**2*vnikx(1,1)*vnikx(2,1)
	 hs2(j)=hs2(j)+w(i)**2*vnikx(1,1)*vnikx(3,1)
	 hs3(j)=hs3(j)+w(i)**2*vnikx(1,1)*vnikx(4,1)
	 j= ileft-4+2
	 y(j) = y(j)+w(i)**2*z(i)*vnikx(2,1)
	 hs0(j)=hs0(j)+w(i)**2*vnikx(2,1)**2
	 hs1(j)=hs1(j)+w(i)**2*vnikx(2,1)*vnikx(3,1)
	 hs2(j)=hs2(j)+w(i)**2*vnikx(2,1)*vnikx(4,1)
	 j= ileft-4+3
	 y(j) = y(j)+w(i)**2*z(i)*vnikx(3,1)
	 hs0(j)=hs0(j)+w(i)**2*vnikx(3,1)**2
	 hs1(j)=hs1(j)+w(i)**2*vnikx(3,1)*vnikx(4,1)
	 j= ileft-4+4
	 y(j) = y(j)+w(i)**2*z(i)*vnikx(4,1)
	 hs0(j)=hs0(j)+w(i)**2*vnikx(4,1)**2
 100  continue

      return
      end
