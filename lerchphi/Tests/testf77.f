c     ---------------------------
c     Test program for lerchphi() 
c     ---------------------------
c
c     This program is copyright by
c
c     Sergej V. Aksenov (http://www.geocities.com/saksenov) and 
c     Ulrich D. Jentschura (jentschura@physik.tu-dresden.de), 2002.
c
c     Version 1.00 (May 1, 2002)
c
c     -------------------------------------------------------------------
c     BRIEF DESCRIPTION:
c     performs a simple calculation with lerchphi (double precision)
c     and demonstrates how lerchphi.c, compiled with the -DADD_UNDERSCORE
c     option, can be linked with a Fortran program.
c     -------------------------------------------------------------------
      program lerchphitest
      real*8 z,s,v,res,acc
      integer it
      z = -0.99999d0
      s = 2.0d0
      v = 1.0d0
      acc = 1.0d-14
      write(*,*) 'lerchphi() test'
      write (*,10) z
      write (*,11) s
      write (*,12) v
      write (*,13) acc
10    format ("  z=",e22.16)
11    format ("  s=",e22.16)
12    format ("  v=",e22.16)
13    format ("acc=",e22.16)
      call lerchphi(z,s,v,acc,res,it)
      write (*,14) res,it
14    format ('Phi=',e22.16,2x,'Iterations=',i3)
      stop
      end
