/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -I. -o test_sincos test_sincos.cpp wrap_boost.cpp gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_sincos > test_sincos.txt

g++ -std=c++14 -o test_sincos test_sincos.cpp
./test_sincos > test_sincos.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>


namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Tp>
    struct __sincos_t
    {
      _Tp sin_value;
      _Tp cos_value;
    };

  template<typename _Tp>
    __sincos_t<_Tp>
    __sincos(_Tp __x)
    { return __sincos_t<_Tp>{std::sin(__x), std::cos(__x)}; }

  template<>
    __sincos_t<float>
    __sincos(float __x)
    {
      float __sin, __cos;
      __builtin_sincosf(__x, &__sin, &__cos);
      return __sincos_t<float>{__sin, __cos};
    }

  template<>
    __sincos_t<double>
    __sincos(double __x)
    {
      double __sin, __cos;
      __builtin_sincos(__x, &__sin, &__cos);
      return __sincos_t<double>{__sin, __cos};
    }

  template<>
    __sincos_t<long double>
    __sincos(long double __x)
    {
      long double __sin, __cos;
      __builtin_sincosl(__x, &__sin, &__cos);
      return __sincos_t<long double>{__sin, __cos};
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

int
main()
{
  auto a = std::__detail::__sincos(1.5f);
  auto b = std::__detail::__sincos(1.5);
  auto c = std::__detail::__sincos(1.5l);
}
