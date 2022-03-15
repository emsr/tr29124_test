c
c
      subroutine mccheb(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,
     *finl,finr,endl,endr,xfer,wfer,a,b,fnu,alpha,beta,ncap,kount,
     *ierr,be,x,w,xm,wm,s,s0,s1,s2)
c
c This is a multiple-component discretized modified Chebyshev
c algorithm, basically a modified Chebyshev algorithm in which the
c modified moments are discretized in the same manner as the inner
c product in the discretization procedure  mcdis. The input and
c output parameters are as in  mcdis. In addition, the arrays  a,b
c must be filled with the recursion coefficients  a(k-1),b(k-1),
c k=1,2,...,2*n-1, defining the modified moments. The arrays
c be,x,w,xm,wm,s,s0,s1,s2  are used for working space. The routine
c calls upon the subroutine  cheb. The routine exits immediately with
c ierr=-1  if  n  is not in range.
c
      dimension xp(*),yp(*),endl(mc),endr(mc),xfer(ncapm),
     *wfer(ncapm),a(*),b(*),fnu(*),alpha(n),beta(n),be(n),x(ncapm),
     *w(ncapm),xm(*),wm(*),s(n),s0(*),s1(*),s2(*)
      logical finl,finr
c
c The arrays  xp,yp  are assumed to have dimension  mp  if mp > 0, 
c the arrays  a,b  dimension 2*n-1, the arrays  fnu,s0,s1,s2  dimension 
c 2*n, and the arrays  xm,wm  dimension  mc*ncapm+mp.
c
      nd=2*n
      if(idelta.le.0) idelta=1
      if(n.lt.1) then
        ierr=-1
        return
      end if
c
c Initialization
c
      incap=1
      kount=-1
      ierr=0
      do 10 k=1,n
        beta(k)=0.
   10 continue
      ncap=(nd-1)/idelta
   20 do 30 k=1,n
        be(k)=beta(k)
   30 continue
      kount=kount+1
      if(kount.gt.1) incap=2**(kount/5)*n
      ncap=ncap+incap
      if(ncap.gt.ncapm) then
        ierr=ncapm
        return
      end if
c
c Discretization of the modified moments
c
      mtncap=mc*ncap
      do 50 i=1,mc
        im1tn=(i-1)*ncap
        if(iq.eq.1) then
          call quad(ncap,x,w,i,ierr)
        else
          call qgp(ncap,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)
        end if
        if(ierr.ne.0) then
          ierr=i
          return
        end if
        do 40 k=1,ncap
          xm(im1tn+k)=x(k)
          wm(im1tn+k)=w(k)
   40   continue
   50 continue
      if(mp.ne.0) then
        do 60 k=1,mp
          xm(mtncap+k)=xp(k)
          wm(mtncap+k)=yp(k)
   60   continue
      end if
      mtnpmp=mtncap+mp
      do 90 k=1,nd
        km1=k-1
        sum=0.
        do 80 i=1,mtnpmp
          p1=0.
          p=1.
          if(k.gt.1) then
            do 70 l=1,km1
              pm1=p1
              p1=p
              p=(xm(i)-a(l))*p1-b(l)*pm1
   70       continue
          end if
          sum=sum+wm(i)*p
   80   continue
        fnu(k)=sum
   90 continue
c
c Computation of the desired recursion coefficients
c
      call cheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c In the following statement, the absolute value of the beta's is
c used to guard against failure in cases where the routine is applied
c to variable-sign weight functions and hence the positivity of the
c beta's is not guaranteed.
c
      do 100 k=1,n
        if(abs(beta(k)-be(k)).gt.eps*abs(beta(k))) goto 20
  100 continue
      return
      end

