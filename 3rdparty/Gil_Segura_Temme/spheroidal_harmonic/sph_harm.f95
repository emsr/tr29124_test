
module sph_harm

contains

  subroutine obl_sph_harm(x, m, nmax, mode, rl, tl, nuevo) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), value :: x
    integer(c_int), value :: m, nmax, mode, nuevo
    real(c_double) :: rl(0:nmax+1), tl(0:nmax+1)
    call doblh(x, m, nmax, mode, rl, tl, nuevo)
  end subroutine

  subroutine pro_sph_harm(x, m, nmax, mode, pl, ql, nuevo) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), value :: x
    integer(c_int), value :: m, nmax, mode, nuevo
    real(c_double) :: pl(0:nmax+1), ql(0:nmax+1)
    call dproh(x, m, nmax, mode, pl, ql, nuevo)
  end subroutine

end module
