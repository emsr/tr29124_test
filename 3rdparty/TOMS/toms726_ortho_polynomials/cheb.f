c
c
      subroutine cheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c Given a set of polynomials  p(0),p(1),...,p(2*n-1)  satisfying
c
c        p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
c                        k=0,1,...,2*n-2,
c
c        p(-1)(x)=0,  p(0)(x)=1,
c
c and associated modified moments
c
c           fnu(k)=integral of p(k)(x)*dlambda(x),
c                        k=0,1,...,2*n-1,
c
c this subroutine uses the modified Chebyshev algorithm (see, e.g.,
c Section 2.4 of W. Gautschi,On generating orthogonal polynomials'',
c SIAM J. Sci. Statist. Comput. 3, 1982, 289-317) to generate the
c recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the
c polynomials  pi(k)  orthogonal with respect to the integration
c measure  dlambda(x), i.e.,
c
c        pi(k+1)(x)=(x-alpha(k))*pi(k)(x)-beta(k)*pi(k-1)(x),
c                          k=0,1,...,n-1,
c
c        pi(-1)(x)=0,  pi(0)(x)=1.
c
c     Input:    n - - the number of recursion coefficients desired
c               a,b-- arrays of dimension 2*n-1 to be filled with the
c                     values of  a(k-1),b(k-1), k=1,2,...,2*n-1
c               fnu-- array of dimension  2*n  to be filled with the 
c                     values of the modified moments  fnu(k-1), k=1,2,
c                     ...,2*n
c     Output:   alpha,beta-- arrays containing, respectively, the
c                     recursion coefficients  alpha(k-1),beta(k-1),
c                     k=1,2,...,n, where  beta(0)  is the total mass.
c               s - - array containing the normalization factors
c                     s(k)=integral [pi(k)(x)]**2 dlambda(x), k=0,1,
c                     2,...,n-1.
c               ierr- an error flag, equal to  0  on normal return, 
c                     equal to  1  if  abs(fnu(0))  is less than the 
c                     machine zero, equal to  2  if  n  is out of range,
c                     equal to  -k  if  s(k), k=0,1,2,...,n-1, is about 
c                     to underflow, and equal to  +k  if it is about to 
c                     overflow.
c
c The arrays  s0,s1,s2  are needed for working space.
c
c On machines with limited exponent range, the occurrence of underflow
c [overflow] in the computation of the  alpha's  and  beta's  can often 
c be avoided by multiplying all modified moments by a sufficiently large
c [small] scaling factor and dividing the new  beta(0)  by the same 
c scaling factor.
c
c The routine uses the function subroutine  r1mach.
c
      dimension a(*),b(*),fnu(*),alpha(n),beta(n),s(n),s0(*),s1(*),
     *s2(*)
c
c The arrays  a,b  are assumed to have dimension  2*n-1, the arrays
c fnu,s0,s1,s2  dimension  2*n.
c
      nd=2*n
      tiny=10.*r1mach(1)
      huge=.1*r1mach(2)
      ierr=0
      if(abs(fnu(1)).lt.tiny) then
        ierr=1
        return
      end if
      if(n.lt.1) then
        ierr=2
        return
      end if
c
c Initialization
c
      alpha(1)=a(1)+fnu(2)/fnu(1)
      beta(1)=fnu(1)
      if(n.eq.1) return
      s(1)=fnu(1)
      do 10 l=1,nd
        s0(l)=0.
        s1(l)=fnu(l)
   10 continue
c
c Continuation
c
      do 40 k=2,n
        lk=nd-k+1
        do 20 l=k,lk
c
c The quantities  s2(l)  for l > k are auxiliary quantities which may
c be zero or may become so small as to underflow, without however
c causing any harm.
c
          s2(l)=s1(l+1)-(alpha(k-1)-a(l))*s1(l)-beta(k-1)*s0(l)
     *      +b(l)*s1(l-1)
          if(l.eq.k) s(k)=s2(k)
c
c Check impending underflow or overflow
c
          if(abs(s(k)).lt.tiny) then
            ierr=-(k-1)
            return
          else if(abs(s(k)).gt.huge) then
            ierr=k-1
            return
          end if
   20   continue
c
c Compute the alpha- and beta-coefficient
c
        alpha(k)=a(k)+(s2(k+1)/s2(k))-(s1(k)/s1(k-1))
        beta(k)=s2(k)/s1(k-1)
        do 30 l=k,lk
          s0(l)=s1(l)
          s1(l)=s2(l)
   30   continue
   40 continue
      return
      end

