      program main

c***********************************************************************
c
cc TOMS468_PRB tests TOMS468.
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS468_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 468, automatic'
      write ( *, '(a)' ) '  numerical integration.'
      write ( *, '(a)' ) ' '

      call test01
      call test02
      call test03

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS468_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test01

c***********************************************************************
c
cc TEST01 tests QUAD.
c
      implicit none

      real a
      real b
      real epsil
      real exact
      real f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13
      external f01
      external f02
      external f03
      external f04
      external f05
      external f06
      external f07
      external f08
      external f09
      external f10
      external f11
      external f12
      external f13
      integer icheck
      integer k
      integer npts
      real result(8)

      epsil = 0.001

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test QUAD, for simple quadrature.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Error tolerance EPSIL = ', epsil
      write ( *, '(a)' ) ' '

      write ( *, '(a,a)' ) 
     &  '      A         B   ICHECK K     NFUNC     ',
     &  'RESULT(K)        EXACT'
      write ( *, '(a)' ) ' '

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f01 )

      exact = 2.0E+00 / 3.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = -1.0E+00
      b =  1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f02 )

      exact = 0.4794282267E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = -1.0E+00
      b =  1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f03 )

      exact = 1.582232964E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f04 )

      exact = 2.0E+00 / 5.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f05 )

      exact = 0.8669729873E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f06 )

      exact = 1.154700669E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f07 )

      exact = 0.7775046341E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.1E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f08 )

      exact = 0.009098645256E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a =  0.0E+00
      b = 10.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f09 )

      exact = 0.4993638029E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 3.1415927E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f10 )

      exact = 0.8386763234E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f11 )

      exact = -1.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f12 )

      exact = -0.6346651825E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f13 )
c
c  The reference lists an exact value of 0.0013492485650E+00 but this is
c  apparently a typo.
c
      exact = 0.013492485650E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &  a, b, icheck, k, npts, result(k), exact

      return
      end
      subroutine test02

c***********************************************************************
c
cc TEST02 tests QSUB.
c
      implicit none

      real a
      real b
      real epsil
      real exact
      real f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13
      external f01
      external f02
      external f03
      external f04
      external f05
      external f06
      external f07
      external f08
      external f09
      external f10
      external f11
      external f12
      external f13
      integer icheck
      integer npts
      real qsub
      real relerr
      real result

      epsil = 0.001

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test QSUB, for quadrature with subdivision.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Error tolerance EPSIL = ', epsil
      write ( *, '(a)' ) ' '

      write ( *, '(a,a)' ) 
     &  '      A         B   ICHECK   NFUNC     ',
     &  '   RESULT        EXACT          RELERR'
      write ( *, '(a)' ) ' '

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f01 )

      exact = 2.0E+00 / 3.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = -1.0E+00
      b =  1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f02 )

      exact = 0.4794282267E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = -1.0E+00
      b =  1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f03 )

      exact = 1.582232964E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f04 )

      exact = 2.0E+00 / 5.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f05 )

      exact = 0.8669729873E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f06 )

      exact = 1.154700669E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f07 )

      exact = 0.7775046341E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.1E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f08 )

      exact = 0.009098645256E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a =  0.0E+00
      b = 10.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f09 )

      exact = 0.4993638029E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 3.1415927E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f10 )

      exact = 0.8386763234E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f11 )

      exact = -1.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f12 )

      exact = -0.6346651825E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f13 )
c
c  The reference lists an exact value of 0.0013492485650E+00 but this is
c  apparently a typo.
c
      exact = 0.013492485650E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      return
      end
      subroutine test03

c***********************************************************************
c
cc TEST03 tests QSUBA.
c
      implicit none

      real a
      real b
      real epsil
      real exact
      real f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13
      external f01
      external f02
      external f03
      external f04
      external f05
      external f06
      external f07
      external f08
      external f09
      external f10
      external f11
      external f12
      external f13
      integer icheck
      integer npts
      real qsuba
      real relerr
      real result

      epsil = 0.001

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test QSUBA,'
      write ( *, '(a)' ) '  for adaptive quadrature with subdivision.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Error tolerance EPSIL = ', epsil
      write ( *, '(a)' ) ' '

      write ( *, '(a,a)' ) 
     &  '      A         B   ICHECK   NFUNC     ',
     &  '   RESULT        EXACT          RELERR'
      write ( *, '(a)' ) ' '

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f01 )

      exact = 2.0E+00 / 3.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = -1.0E+00
      b =  1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f02 )

      exact = 0.4794282267E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = -1.0E+00
      b =  1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f03 )

      exact = 1.582232964E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f04 )

      exact = 2.0E+00 / 5.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f05 )

      exact = 0.8669729873E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f06 )

      exact = 1.154700669E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f07 )

      exact = 0.7775046341E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.1E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f08 )

      exact = 0.009098645256E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a =  0.0E+00
      b = 10.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f09 )

      exact = 0.4993638029E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 3.1415927E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f10 )

      exact = 0.8386763234E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f11 )

      exact = -1.0E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f12 )

      exact = -0.6346651825E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f13 )
c
c  The reference lists an exact value of 0.0013492485650E+00 but this is
c  apparently a typo.
c
      exact = 0.013492485650E+00

      write ( *, 
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' ) 
     &  a, b, icheck, npts, result, exact, relerr

      return
      end
      function f01 ( x )

c***********************************************************************
c
cc F01 evaluates test function 1.
c
      implicit none

      real f01
      real x

      f01 = sqrt ( x )

      return
      end
      function f02 ( x )

c***********************************************************************
c
cc F02 evaluates test function 2.
c
      implicit none

      real f02
      real x

      f02 = 0.92E+00 * cosh ( x ) - cos ( x )

      return
      end
      function f03 ( x )

c***********************************************************************
c
cc F03 evaluates test function 3.
c
      implicit none

      real f03
      real x

      f03 = 1.0 / ( x**4 + x**2 + 0.9 )

      return
      end
      function f04 ( x )

c***********************************************************************
c
cc F04 evaluates test function 4.
c
      implicit none

      real f04
      real x

      f04 = sqrt ( x**3 )

      return
      end
      function f05 ( x )

c***********************************************************************
c
cc F05 evaluates test function 5.
c
      implicit none

      real f05
      real x

      f05 = 1.0E+00 / ( 1.0E+00 + x**4 )

      return
      end
      function f06 ( x )

c***********************************************************************
c
cc F06 evaluates test function 6.
c
      implicit none

      real f06
      real x

      f06 = 1.0E+00 / ( 1.0E+00 + 0.5E+00 * sin ( 31.4159E+00 * x ) )

      return
      end
      function f07 ( x )

c***********************************************************************
c
cc F07 evaluates test function 7.
c
      implicit none

      real f07
      real x

      f07 = x / ( exp ( x ) - 1.0E+00 )

      return
      end
      function f08 ( x )

c***********************************************************************
c
cc F08 evaluates test function 8.
c
      implicit none

      real f08
      real x

      f08 = sin ( 314.159E+00 * x ) / ( 3.14159E+00 * x )

      return
      end
      function f09 ( x )

c***********************************************************************
c
cc F09 evaluates test function 9.
c
      implicit none

      real f09
      real x

      f09 = 50.0E+00 / ( 2500.0E+00 * x**2 + 1.0E+00 ) / 3.14159E+00

      return
      end
      function f10 ( x )

c***********************************************************************
c
cc F10 evaluates test function 10.
c
      implicit none

      real arg
      real f10
      real x

      arg =         cos (           x ) 
     &  + 3.0E+00 * sin (           x ) 
     &  + 2.0E+00 * cos ( 2.0E+00 * x )
     &  + 3.0E+00 * cos ( 3.0E+00 * x ) 
     &  + 3.0E+00 * sin ( 2.0E+00 * x )

      f10 = cos ( arg )

      return
      end
      function f11 ( x )

c***********************************************************************
c
cc F11 evaluates test function 11.
c
      implicit none

      real f11
      real x

      f11 = log ( x )

      return
      end
      function f12 ( x )

c***********************************************************************
c
cc F12 evaluates test function 12.
c
      implicit none

      real f12
      real pi
      real x

      pi = 3.141592653589793E+00

      f12 = 4.0E+00 * pi**2 * x * sin ( 20.0E+00 * pi * x ) 
     &  * cos ( 2.0E+00 * pi * x )

      return
      end
      function f13 ( x )

c***********************************************************************
c
cc F13 evaluates test function 13.
c
      implicit none

      real f13
      real x

      f13 = 1.0E+00 / ( 1.0E+00 + ( 230.0E+00 * x - 30.0E+00 )**2 )

      return
      end
