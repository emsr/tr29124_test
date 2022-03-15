c
c
      function nu0her(n,z,eps)
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Hermite measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      complex z
      nu0her=2.*(sqrt(.5*real(n+1))+.25*alog(1./eps)/
     *  abs(aimag(z)))**2
      return
      end

