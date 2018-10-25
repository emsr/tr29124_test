c
c gfortran -std=f95 -g -o run_coulfg RUN_COULFG.FOR COULFG.FOR
c ./run_coulfg > run_coulfg.txt
c
      program runcoulfg

      implicit none
      integer mode, kfn, ifail
      integer ih, ir
      double precision eta, lam, rho
      double precision fc(100), gc(100), fcp(100), gcp(100)

      integer nfp, npq, iexp, m1
      double precision paccq
      common /steed/ paccq, nfp, npq, iexp, m1

      mode = 1
      kfn = 0
      lam = 0
      do ih = -2, 10, 2
        eta = ih * 1.0
        write(6, *) ''
        write(6, *) ''
        write(6, *) 'eta = ', eta
        do ir = 1, 200
          rho = ir * 0.1
          fc(1) = 0
          gc(1) = 0
          fcp(1) = 0
          gcp(1) = 0
          call coulfg(rho, eta, lam, lam + 1, fc, gc, fcp, gcp, 
     *                mode, kfn, ifail)
          write(6, '(G16.6,G16.6,G16.6,G16.6,G16.6,I4,I4)')
     *          rho, fc(1), gc(1), fcp(1), gcp(1), ifail, m1
        end do
      end do

      end program
