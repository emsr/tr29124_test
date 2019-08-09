
module scorer

contains

  subroutine scorer_gi(ifacg, x, y, reg, img, regp, imgp, ierrog) bind(c)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value :: ifacg
    real(c_double), intent(in), value :: x
    real(c_double), intent(out) :: y, reg, img, regp, imgp
    integer(c_int), intent(out) :: ierrog
    call giz(ifacg, x, y, reg, img, regp, imgp, ierrog)
  end subroutine

  subroutine scorer_hi(ifach, x, y, reh, imh, rehp, imhp, ierroh) bind(c)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value :: ifach
    real(c_double), intent(in), value :: x
    real(c_double), intent(out) :: y, reh, imh, rehp, imhp
    integer(c_int), intent(out) :: ierroh
    call hiz(ifach, x, y, reh, imh, rehp, imhp, ierroh)
  end subroutine

end module
