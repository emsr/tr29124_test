c
c
      function nu0jac(n,z,eps)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Jacobi measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      complex z
      pi=4.*atan(1.)
      x=real(z)
      y=abs(aimag(z))
      if(x.lt.1.) then
        if(x.lt.-1.) angle=.5*(2.*pi+atan(y/(x-1.))+atan(y/(x+1.)))
        if(x.eq.-1.) angle=.5*(1.5*pi-atan(.5*y))
        if(x.gt.-1.) angle=.5*(pi+atan(y/(x-1.))+atan(y/(x+1.)))
      else
        if(x.eq.1.) angle=.5*(.5*pi+atan(.5*y))
        if(x.gt.1.) angle=.5*(atan(y/(x-1.))+atan(y/(x+1.)))
      end if
      x2=x*x
      y2=y*y
      r=((x2-y2-1.)**2+4.*x2*y2)**.25
      r=sqrt((x+r*cos(angle))**2+(y+r*sin(angle))**2)
      nu0jac=real(n+1)+.5*alog(1./eps)/alog(r)
      return
      end 

