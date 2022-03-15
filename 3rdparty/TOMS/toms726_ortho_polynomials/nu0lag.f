c
c
      function nu0lag(n,z,al,eps)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Laguerre measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      complex z
      pi=4.*atan(1.)
      x=real(z)
      y=aimag(z)
      phi=.5*pi
      if(y.lt.0.) phi=1.5*pi
      if(x.eq.0.) goto 10
      phi=atan(y/x)
      if(y.gt.0. .and. x.gt.0.) goto 10
      phi=phi+pi
      if(x.lt.0.) goto 10
      phi=phi+pi
   10 nu0lag=(sqrt(real(n+1)+.5*(al+1.))+alog(1./eps)/(4.*(x*x+
     *  y*y)**.25*cos(.5*(phi-pi))))**2-.5*(al+1.)
      return
      end

