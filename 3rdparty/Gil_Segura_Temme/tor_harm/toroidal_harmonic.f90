
module toroidal_harmonic

contains

  subroutine tor_harmonic(z, m, nmax, pl, ql, newn) bind(c)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value :: m, nmax
    real(c_double), intent(in), value :: z
    real(c_double), intent(out) :: pl(0:nmax+1), ql(0:nmax+1)
    integer(c_int), intent(inout) :: newn
    call dtorh1(z, m, nmax, pl, ql, newn)
  end subroutine

end module
