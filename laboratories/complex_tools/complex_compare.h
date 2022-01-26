#ifndef COMPLEX_COMPARE_H
#define COMPLEX_COMPARE_H 1

#include <complex>

/**
 * A sane operator< for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<(const std::complex<_Tp>& y, const std::complex<_Up>& z)
  {
    if (y.real() < z.real())
      return true;
    else if (z.real() < y.real())
      return false;
    else
      return y.imag() < z.imag();
  }

/**
 * A sane operator< for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<(const std::complex<_Tp>& y, _Up z)
  {
    if (y.real() < z)
      return true;
    else if (z < y.real())
      return false;
    else
      return y.imag() < _Up{};
  }

/**
 * A sane operator< for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<(_Tp y, const std::complex<_Up>& z)
  {
    if (y < z.real())
      return true;
    else if (z.real() < y)
      return false;
    else
      return _Tp{} < z.imag();
  }

/**
 * A sane operator> for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>(const std::complex<_Tp>& y, const std::complex<_Up>& z)
  { return operator<(z, y); }

/**
 * A sane operator> for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>(const std::complex<_Tp>& y, _Up z)
  { return operator<(z, y); }

/**
 * A sane operator> for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>(_Tp y, const std::complex<_Up>& z)
  { return operator<(z, y); }

/**
 * A sane operator>= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>=(const std::complex<_Tp>& y, const std::complex<_Up>& z)
  { return !operator<(y, z); }

/**
 * A sane operator>= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>=(const std::complex<_Tp>& y, _Up z)
  { return !operator<(y, z); }

/**
 * A sane operator>= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>=(_Tp y, const std::complex<_Up>& z)
  { return !operator<(y, z); }

/**
 * A sane operator<= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<=(const std::complex<_Tp>& y, const std::complex<_Up>& z)
  { return !operator>(y, z); }

/**
 * A sane operator<= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<=(const std::complex<_Tp>& y, _Up z)
  { return !operator>(y, z); }

/**
 * A sane operator<= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<=(_Tp y, const std::complex<_Up>& z)
  { return !operator>(y, z); }

#endif // COMPLEX_COMPARE_H
