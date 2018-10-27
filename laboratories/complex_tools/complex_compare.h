#ifndef COMPLEX_COMPARE_H
#define COMPLEX_COMPARE_H 1

#include <complex>

/**
 * A sane operator< for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<(const std::complex<_Tp>& __y, const std::complex<_Up>& __z)
  {
    if (__y.real() < __z.real())
      return true;
    else if (__z.real() < __y.real())
      return false;
    else
      return __y.imag() < __z.imag();
  }

/**
 * A sane operator< for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<(const std::complex<_Tp>& __y, _Up __z)
  {
    if (__y.real() < __z)
      return true;
    else if (__z < __y.real())
      return false;
    else
      return __y.imag() < _Up{};
  }

/**
 * A sane operator< for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<(_Tp __y, const std::complex<_Up>& __z)
  {
    if (__y < __z.real())
      return true;
    else if (__z.real() < __y)
      return false;
    else
      return _Tp{} < __z.imag();
  }

/**
 * A sane operator> for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>(const std::complex<_Tp>& __y, const std::complex<_Up>& __z)
  { return operator<(__z, __y); }

/**
 * A sane operator> for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>(const std::complex<_Tp>& __y, _Up __z)
  { return operator<(__z, __y); }

/**
 * A sane operator> for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>(_Tp __y, const std::complex<_Up>& __z)
  { return operator<(__z, __y); }

/**
 * A sane operator>= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>=(const std::complex<_Tp>& __y, const std::complex<_Up>& __z)
  { return !operator<(__y, __z); }

/**
 * A sane operator>= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>=(const std::complex<_Tp>& __y, _Up __z)
  { return !operator<(__y, __z); }

/**
 * A sane operator>= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator>=(_Tp __y, const std::complex<_Up>& __z)
  { return !operator<(__y, __z); }

/**
 * A sane operator<= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<=(const std::complex<_Tp>& __y, const std::complex<_Up>& __z)
  { return !operator>(__y, __z); }

/**
 * A sane operator<= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<=(const std::complex<_Tp>& __y, _Up __z)
  { return !operator>(__y, __z); }

/**
 * A sane operator<= for complex numbers.
 */
template<typename _Tp, typename _Up>
  constexpr bool
  operator<=(_Tp __y, const std::complex<_Up>& __z)
  { return !operator>(__y, __z); }

#endif // COMPLEX_COMPARE_H
