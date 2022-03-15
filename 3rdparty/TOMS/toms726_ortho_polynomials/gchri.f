c
c
      subroutine gchri(n,iopt,nu0,numax,eps,a,b,x,y,alpha,beta,nu,
     *  ierr,ierrc,fnu,rho,rold,s,s0,s1,s2)
c
c This routine implements the generalized Christoffel theorem, using
c the method of modified moments (cf. Section 4 of W. Gautschi,
c Minimal solutions of three-term recurrence relations and orthogonal
c polynomials'', Math. Comp. 36, 1981, 547-554). Given the recursion
c coefficients  a(k), b(k), k=0,1,...n, for the (monic) orthogonal
c polynomials with respect to some measure  dlambda(t), it generates
c the recursion coefficients  alpha(k), beta(k), k=0,1,2,...,n-1 for
c the measure
c
c         dlambda(t)/(t-x)        if iopt=1
c         dlambda(t)/{(t-x)**2+y**2} if iopt=2
c
c   Input:  n   - - the number of recurrence coefficients desired;
c                   type integer
c           iopt  - an integer selecting the desired weight distribution
c           nu0   - an integer estimating the starting backward 
c                   recurrence index; in the absence of any better
c                   choice, take  nu0 = 3*n
c           numax - an integer controlling termination of backward
c                   recursion in case of nonconvergence; a conservative 
c                   choice is  numax = 500
c           eps - - a relative error tolerance; type real
c           a,b - - arrays of dimension numax to be supplied with the
c                   recursion coefficients a(k)=alpha(k-1),b(k)=beta(k),
c                   k=1,2,...,numax, for the measure  dlambda
c           x,y - - real parameters defining the linear and quadratic
c                   divisors of  dlambda
c
c   Output: alpha,beta - arrays of dimension  n  containing the desired
c                   recursion coefficients  alpha(k-1), beta(k-1), k=1,
c                   2,...,n
c           nu  - - the backward recurrence index yielding convergence; 
c                   in case of nonconvergence,  nu  will have the value 
c                   numax
c           ierr  - an error flag, where
c                   ierr=0     on normal return
c                   ierr=1     if  iopt  is neither 1 nor 2
c                   ierr=nu0   if  nu0 > numax
c                   ierr=numax if the backward recurrence algorithm does
c                              not converge
c                   ierr=-1    if  n  is not in range
c           ierrc - an error flag inherited from the routine  cheb
c
c The arrays  fnu,s,s0,s1,s2  are working space. The routine calls
c upon the routines  knum  and  cheb.
c
      complex rho,rold,z
      dimension a(numax),b(numax),alpha(n),beta(n),fnu(*),rho(*),
     *rold(*),s(n),s0(*),s1(*),s2(*)
c
c The arrays  fnu,rho,rold,s0,s1,s2  are assumed to have dimension  2*n.
c
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      nd=2*n
      ndm1=nd-1
c
c Linear divisor
c
      if(iopt.eq.1) then
c
c Generate the modified moments of  dlambda.
c
        z=cmplx(x,0.)
        call knum(ndm1,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
        do 10 k=1,nd
          fnu(k)=-real(rho(k))
   10   continue
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm.
c
        call cheb(n,a,b,fnu,alpha,beta,s,ierrc,s0,s1,s2)
        return
c
c Quadratic divisor
c
      else if(iopt.eq.2) then
c
c Generate the modified moments of  dlambda.
c
        y=abs(y)
        z=cmplx(x,y)
        call knum(ndm1,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
        do 20 k=1,nd
          fnu(k)=-aimag(rho(k))/y
   20   continue
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm.
c
        call cheb(n,a,b,fnu,alpha,beta,s,ierrc,s0,s1,s2)
        return
      else
        ierr=1
        return
      end if
      end

