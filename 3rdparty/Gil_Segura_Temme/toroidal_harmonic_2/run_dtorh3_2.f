c gfortran -o run run_dtorh3_2.f dtorh3_2.f  Rout.f adrout.f
c ./run > run.txt

      program run_dtorh3_2

        implicit none
        integer ix, m, n;
        integer mmax, nmax
        integer mnew, nnew
        integer mdim, ndim
        double precision x, pl, ql
        dimension x(3)
        dimension pl(0:11,0:11), ql(0:11,0:11)

        x(1) = 1.01d0
        x(2) = 2.0d0
        x(3) = 5.0d0

        do ix=1,3
          mmax = 10
          nmax = 10
          write(6,*) 'x = ', x(ix)
          call dtorh3(x(ix),11,11,mmax,nmax,pl,ql,mnew,nnew)
          do m=0,mmax
            do n=0,nmax
              write(6, '(f24.6)', advance='no') pl(m,n)
            end do
            write(6,*)
          end do
          write(6,*)
        end do
      end program
