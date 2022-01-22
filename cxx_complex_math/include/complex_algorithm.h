#ifndef COMPLEX_ALGORITHM_H
#define COMPLEX_ALGORITHM_H 1

#include <numeric> // For midpoint.
#include <cmath> // For lerp.

#if __cplusplus >= 201703L
#if __cpp_lib_interpolate >=  201902L
namespace emsr //std
{
  // Overload midpoint.
  template<typename _Tp>
    constexpr std::complex<_Tp>
    midpoint(const std::complex<_Tp>& u, const std::complex<_Tp>& v)
    {
      return std::complex<_Tp>{std::midpoint(u.real(), v.real()),
			       std::midpoint(u.imag(), v.imag())};
    }

  // Overload lerp.
  template<typename _Tp>
    constexpr std::complex<_Tp>
    lerp(const std::complex<_Tp>& u, const std::complex<_Tp>& v, _Tp t)
    {
      return std::complex<_Tp>{std::lerp(u.real(), v.real(), t),
			       std::lerp(u.imag(), v.imag(), t)};
    }
}
#endif // __cpp_lib_interpolate
#endif // C++17

#endif // COMPLEX_ALGORITHM_H
