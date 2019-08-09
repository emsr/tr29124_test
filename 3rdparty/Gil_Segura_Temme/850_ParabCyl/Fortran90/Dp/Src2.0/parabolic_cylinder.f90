
module parabolic_cylinder

contains

  subroutine parab_cyl(a, x, mode, uaxx, vaxx, ierr) bind(c)
    use iso_c_binding
    use Parabolic
    implicit none
    integer(c_int), intent(in), value :: mode
    real(c_double), intent(in), value :: a
    real(c_double), intent(in), value :: x
    real(c_double), intent(out) :: uaxx(2)
    real(c_double), intent(out) :: vaxx(2) 
    integer(c_int), intent(out) :: ierr
    call parab(a, x, mode, uaxx, vaxx, ierr)
  end subroutine

end module
