c
c
      subroutine qgp(n,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)
c
c This is a general-purpose discretization routine that can be used
c as an alternative to the routine  quad  in the multiple-component 
c discretization procedure  mcdis. It takes no account of the special
c nature of the weight function involved and hence may result in slow 
c convergence of the discretization procedure. This routine, therefore, 
c should be used only as a last resort, when no better, more natural 
c discretization can be found.
c
c It is assumed that there are  mc.ge.1  disjoint component intervals.
c The discretization is effected by the Fejer quadrature rule,
c suitably transformed to the respective interval. An interval that
c extends to minus infinity has to be indexed by  1; one that extends
c to plus infinity has to be indexed by  mc.
c
c The output variable  ierr  is given the value  0. Additional input 
c parameters and working space used by this routine are as follows:
c
c          mc      - the number of component intervals; type integer
c          finl    - a logical variable to be set .true. if the
c                    extreme left interval is finite, and .false.
c                    otherwise
c          finr    - a logical variable to be set .true. if the
c                    extreme right interval is finite, and .false.
c                    otherwise
c          endl    - an array of dimension  mc  containing the left
c                    endpoints of the component intervals; if the
c                    first of these extends to minus infinity,  endl(1)
c                    can be set to an arbitrary value
c          endr    - an array of dimension  mc  containing the right
c                    endpoints of the component intervals; if the
c                    last of these extends to plus infinity,  endr(mc)
c                    can be set to an arbitrary value
c          xfer,wfer-working arrays holding the Fejer nodes and
c                    weights, respectively, for the interval [-1,1].
c
c The user has to supply the routine 
c
c                     function wf(x,i),
c
c which evaluates the weight function at the point  x  on the i-th 
c component interval. The routine also uses the subroutines  fejer,
c symtr  and  tr, which are appended.
c
      dimension x(n),w(n),endl(mc),endr(mc),xfer(*),wfer(*)
      logical finl,finr
c
c The arrays  xfer,wfer  are dimensioned in the routine  mcdis.
c
      ierr=0
      if(i.eq.1) call fejer(n,xfer,wfer)
      if(i.gt.1 .and. i.lt.mc) goto 60
      if(mc.eq.1) then
        if(finl.and.finr) goto 60
        if(finl) goto 20
        if(finr) goto 40
        do 10 k=1,n
          call symtr(xfer(k),phi,phi1)
          x(k)=phi
          w(k)=wfer(k)*wf(phi,i)*phi1
   10   continue
        return
      else
        if((i.eq.1.and.finl).or.(i.eq.mc.and.finr)) goto 60
        if(i.eq.1) goto 40
      end if
   20 do 30 k=1,n
        call tr(xfer(k),phi,phi1)
        x(k)=endl(mc)+phi
        w(k)=wfer(k)*wf(x(k),mc)*phi1
   30 continue
      return
   40 do 50 k=1,n
        call tr(-xfer(k),phi,phi1)
        x(k)=endr(1)-phi
        w(k)=wfer(k)*wf(x(k),1)*phi1
   50 continue
      return
   60 do 70 k=1,n
        x(k)=.5*((endr(i)-endl(i))*xfer(k)+endr(i)+endl(i))
        w(k)=.5*(endr(i)-endl(i))*wfer(k)*wf(x(k),i)
   70 continue
      return
      end

      subroutine symtr(t,phi,phi1)
c
c This implements a particular transformation  x=phi(t)  mapping
c the t-interval [-1,1] to the x-interval [-oo,oo].
c
c        input:   t
c        output:  phi=phi(t)
c                 phi1=derivative of phi(t)
c
      t2=t*t
      phi=t/(1.-t2)
      phi1=(t2+1.)/(t2-1.)**2
      return
      end

      subroutine tr(t,phi,phi1)
c
c This implements a particular transformation  x=phi(t)  mapping
c the t-interval [-1,1] to the x-interval [0,oo].
c
c         input:   t
c         output:  phi=phi(t)
c                  phi1=derivative of phi(t)
c
      phi=(1.+t)/(1.-t)
      phi1=2./(t-1.)**2
      return
      end

      subroutine fejer(n,x,w)
c
c This routine generates the n-point Fejer quadrature rule.
c
c         input:   n   - the number of quadrature nodes
c         output:  x,w - arrays of dimension  n  holding the quadrature
c                        nodes and weights, respectively; the nodes
c                        are ordered increasingly
c
      dimension x(n),w(n)
      pi=4.*atan(1.)
      nh=n/2
      np1h=(n+1)/2
      fn=real(n)
      do 10 k=1,nh
        x(n+1-k)=cos(.5*real(2*k-1)*pi/fn)
        x(k)=-x(n+1-k)
   10 continue
      if(2*nh.ne.n) x(np1h)=0.
      do 30 k=1,np1h
        c1=1.
        c0=2.*x(k)*x(k)-1.
        t=2.*c0
        sum=c0/3.
        do 20 m=2,nh
          c2=c1
          c1=c0
          c0=t*c1-c2
          sum=sum+c0/real(4*m*m-1)
   20   continue
        w(k)=2.*(1.-2.*sum)/fn
        w(n+1-k)=w(k)
   30 continue
      return
      end

