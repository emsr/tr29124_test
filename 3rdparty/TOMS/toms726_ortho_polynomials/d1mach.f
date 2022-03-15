c
c Check this for modern constants:
c https://people.sc.fsu.edu/~jburkardt/f_src/machine/machine.html
c
      double precision function d1mach(i)
c
c  Double-precision machine constants
c
c  d1mach( 1) = b**(emin-1), the smallest positive magnitude.
c
c  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c  d1mach( 3) = b**(-t), the smallest relative spacing.
c
c  d1mach( 4) = b**(1-t), the largest relative spacing.
c
c  d1mach( 5) = log10(b)
c
c  To alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  removing the c from column 1.
c  On rare machines a static statement may need to be added.
c  (But probably more systems prohibit it than require it.)
c
c  For IEEE-arithmetic machines (binary standard), one of the second
c  two sets of constants below should be appropriate.
c
c  Where possible, decimal, octal or hexadecimal constants are used
c  to specify the constants exactly.  Sometimes this requires using
c  equivalent integer arrays.  If your compiler uses half-word
c  integers by default (sometimes called integer*2), you may need to
c  change integer to integer*4 or otherwise instruct your compiler
c  to use full-word integers in the next 5 declarations.
c
      double precision dmach(5)
      double precision x

c     standard machine constants
      x = 1.0d0
      dmach(1) = tiny(x)
      dmach(2) = huge(x)
      dmach(3) = epsilon(x) / radix(x)
      dmach(4) = epsilon(x)
      dmach(5) = log(x * radix(x))

      if (i .lt. 1  .or.  i .gt. 5) goto 999
      d1mach = dmach(i)
      return

  999 write(*,1999) i
 1999 format(' d1mach - i out of bounds',i10)
      stop
      end

