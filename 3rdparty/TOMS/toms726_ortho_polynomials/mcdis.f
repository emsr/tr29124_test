c
c
      subroutine mcdis(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,irout,
     *finl,finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,
     *be,x,w,xm,wm,p0,p1,p2)
c
c This is a multiple-component discretization procedure as described in
c Section 4.3 of the companion paper. It generates to a relative
c accuracy of  eps  the recursion coefficients  alpha(k), beta(k),
c k=0,1,...,n-1, for the polynomials orthogonal with respect to a
c weight distribution consisting of the sum of  mc  continuous
c components and a discrete component with  mp  points. The continuous
c part of the spectrum is made up of  mc  weight functions, each
c supported on its own interval. These intervals may or may not be
c disjoint. The discretization of the inner product on the i-th
c interval is furnished either by a user-supplied subroutine  quad,
c or by the general-purpose subroutine  qgp  provided in this package,
c depending on whether  iq  is equal, or not equal, to  1, respectively.
c The user-supplied routine must have the form  quad(n,x,w,i,ierr)  and
c is assumed to supply the abscissas  x(k)  and weights  w(k), k=1,2,
c ...,n, to be used in approximating the i-th inner product
c
c               integral of p(x)*q(x)*wf(x,i)dx
c
c by the
c
c       sum over k from 1 to n of w(k)*p(x(k))*q(x(k)),
c
c                                        i=1,2,...,mc.
c
c The desired recurrence coefficients are then approximated by the
c recursion coefficients of the discrete orthogonal polynomials 
c belonging to the discretized inner product, which in turn are 
c computed by either the Stieltjes procedure or the Lanczos algorithm 
c according as  irout  is equal to, or not equal to  1, respectively.
c Two error flags  ierr,ie  are provided which signal the occurrence 
c of an error condition in the quadrature process, or in the routine 
c sti  or  lancz  (whichever is used), respectively. The point spectrum 
c is given through its abscissas  xp  and jumps  yp.
c
c If the quadrature routine  quad  has polynomial degree of exactness
c at least  id(n)  for each i, and if  id(n)/n = idelta + O(1/n)  as 
c n  goes to infinity, then the procedure is designed to converge after
c one iteration, provided  idelta  is set with the appropriate
c integer. Normally,  idelta=1 (for interpolatory rules) or  idelta=2
c (for Gaussian rules). The default value is  idelta=1.
c
c    Input:  n    - - the number of recursion coefficients desired;
c                     type integer
c            ncapm  - a discretization parameter indicating an upper
c                     limit of the fineness of the discretization;
c                     ncapm=500  will usually be satisfactory; type
c                     integer
c            mc  - -  the number of disjoint intervals in the
c                     continuous part of the spectrum; type integer
c            mp  - -  the number of points in the discrete part of
c                     the spectrum; type integer. If there is no
c                     point spectrum, set  mp=0.
c            xp  - -  an array of dimension  mp  containing the
c                     abscissas of the point spectrum
c            yp  - -  an array of dimension  mp  containing the jumps
c                     of the point spectrum
c            quad  -  a subroutine determining the discretization of
c                     the inner product on each component interval,
c                     or a dummy routine if  iq  is not equal to  1
c                     (see below)
c            eps  - - the desired relative accuracy of the nonzero
c                     recursion coefficients; type real
c            iq   - - an integer selecting a user-supplied quadrature
c                     routine  quad  if  iq=1  or the ORTHPOL routine 
c                     qgp  otherwise
c            idelta - a nonzero integer, typically  1  or  2, inducing
c                     fast convergence in the case of special quadrature
c                     routines
c            irout  - an integer selecting the routine for generating
c                     the recursion coefficients from the discrete
c                     inner product. Specifically,  irout=1  selects the
c                     routine  sti, whereas any other value selects the
c                     routine  lancz
c            
c The logical variables  finl,finr, the arrays  endl,endr  of 
c dimension  mc, and the arrays  xfer,wfer  of dimension  ncapm  are 
c input variables to the subroutine  qgp  and are used (and hence need
c to be properly dimensioned) only if  iq  is not equal to  1.
c
c    Output:  alpha,beta - arrays of dimension n, holding as k-th
c                     element  alpha(k-1), beta(k-1), k=1,2,...,n,
c                     respectively
c             ncap  - an integer indicating the fineness of the
c                     discretization that yields convergence within
c                     the eps-tolerance
c             kount - the number of iterations used
c             ierr  - an error flag, equal to  0  on normal return,
c                     equal to  -1  if  n  is not in the proper range,
c                     equal to  i  if there is an error condition in
c                     the discretization of the i-th interval,
c                     and equal to  ncapm  if the discretized 
c                     Stieltjes procedure does not converge within the
c                     discretization resolution specified by  ncapm
c             ie - -  an error flag inherited from the routine  sti
c                     or  lancz  (whichever is used)
c
c The array  be  of dimension  n, the arrays  x,w  of dimension  ncapm,
c and the arrays  xm,wm,p0,p1,p2  of dimension mc*ncapm+mp  are used 
c for working space. The routine calls upon the subroutine  sti  or
c lancz, depending on the choice of  irout.
c
      dimension xp(*),yp(*),endl(mc),endr(mc),xfer(ncapm),wfer(ncapm),
     *alpha(n),beta(n),be(n),x(ncapm),w(ncapm),xm(*),wm(*),p0(*),p1(*),
     *p2(*)
      logical finl,finr
c
c The arrays  xp,yp  are assumed to have dimension  mp  if mp > 0, and
c the arrays  xm,wm,p0,p1,p2  dimension  mc*ncapm+mp.
c
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
      ncap=(2*n-1)/idelta
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
c Discretization of the inner product
c
      mtncap=mc*ncap
      do 50 i=1,mc
        im1tn=(i-1)*ncap
        if(iq.eq.1) then
          call quad(ncap,x,w,i,ierr)
        else
          call qgp(ncap,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,
     *      wfer)
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
c
c Computation of the desired recursion coefficients
c
      if(irout.eq.1) then
        call sti(n,mtncap+mp,xm,wm,alpha,beta,ie,p0,p1,p2)
      else
        call lancz(n,mtncap+mp,xm,wm,alpha,beta,ie,p0,p1)
      end if 
c
c In the following statement, the absolute value of the beta's is
c used to guard against failure in cases where the routine is applied
c to variable-sign weight functions and hence the positivity of the
c beta's is not guaranteed.
c
      do 70 k=1,n
        if(abs(beta(k)-be(k)).gt.eps*abs(beta(k))) goto 20
   70 continue
      return
      end

