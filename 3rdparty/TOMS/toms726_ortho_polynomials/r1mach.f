c
c Check this for modern constants:
c https://people.sc.fsu.edu/~jburkardt/f_src/machine/machine.html
c
      real function r1mach(i)
c
c  Single-precision machine constants
c
c  r1mach(1) = b**(emin-1), the smallest positive magnitude.
c
c  r1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c  r1mach(3) = b**(-t), the smallest relative spacing.
c
c  r1mach(4) = b**(1-t), the largest relative spacing.
c
c  r1mach(5) = log10(b)
c
      real rmach(5)
      real x

c     standard machine constants
      x = 1.0
      rmach(1) = tiny(x)
      rmach(2) = huge(x)
      rmach(3) = epsilon(x) / radix(x)
      rmach(4) = epsilon(x)
      rmach(5) = log(x * radix(x))

      if (i .lt. 1  .or.  i .gt. 5) goto 999
      r1mach = rmach(i)
      return

  999 write(*,1999) i
 1999 format(' r1mach - i out of bounds',i10)
      stop
      end

