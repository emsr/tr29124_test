
module incomplete_gamma

contains

  subroutine inc_gamma(a, x, p, q, ierr) bind(c)
    use iso_c_binding
    use IncgamFI
    implicit none
    real(c_double), intent(in), value :: a
    real(c_double), intent(in), value :: x
    real(c_double), intent(out) :: p
    real(c_double), intent(out) :: q 
    integer(c_int), intent(out) :: ierr
    call incgam(a, x, p, q, ierr)
  end subroutine

  subroutine inv_inc_gamma(a, p, q, x, ierr) bind(c)
    use iso_c_binding
    use IncgamFI
    implicit none
    real(c_double), intent(in), value :: a
    real(c_double), intent(in), value :: p
    real(c_double), intent(in), value :: q 
    real(c_double), intent(out) :: x
    integer(c_int), intent(out) :: ierr
    call invincgam(a, p, q, x, ierr)
  end subroutine

end module
