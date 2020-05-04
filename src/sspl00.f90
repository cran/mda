subroutine sspl00(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,dfoff&
  &,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,xrange,work,ier)
! Naras fix
! integer         n,p,nk,method,ier,nef, nefp1, n2,match(1),center
integer         n,p,nk,method,ier,nef,match(1),center
double precision x(n),y(n,p),w(n),knot(nk+4),tol,wp(p),dfoff,dfmax,cost&
  &,lambda,df,cv,gcv,coef(1),s(n,p),lev(nef),xrange(2),work(1)
!workspace must be (2*p+2)*nefp1 + (p+16)*nk + n +p
logical center2
if (center.eq.0) then
  center2=.false.
else
  center2=.true.
end if 
call sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center2,dfoff,dfmax&
  &,cost,lambda,df,cv,gcv,coef,s,lev,xrange,work,ier)
return
end  
