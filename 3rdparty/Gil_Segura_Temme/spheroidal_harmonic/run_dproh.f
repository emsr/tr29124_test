c gfortran -o run_dproh run_dproh.f dproh.f fact.f 
c ./run_dproh > dproh.txt

      program run_dproh

      implicit none
      integer ix, nmax, m, mode, nuevo, n;
      double precision, dimension(5) :: x
      double precision, dimension(0:11) :: rl
      double precision, dimension(0:11) :: tl

      x(1) = 1.1d0
      x(2) = 1.5d0
      x(3) = 2.0d0
      x(4) = 5.0d0
      x(5) = 10.0d0
      nmax = 10
      m = 3
      mode = 0

      do ix=1,5
        write(6,*) 'x = ', x(ix)
        call dproh(x(ix), m, nmax, mode, rl, tl, nuevo)

        write(6, '(a4)', advance='no') 'Rl:'
        do n=0, nmax
          write(6, '(f24.6)', advance='no') rl(n)
        end do
        write(6,*)

        write(6, '(a4)', advance='no') 'Tl:'
        do n=0, nmax
          write(6, '(f24.6)', advance='no') tl(n)
        end do
        write(6,*)
      end do

      end program