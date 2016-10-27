/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_sincos test_sincos.cpp wrap_boost.cpp -lquadmath
./test_sincos > test_sincos.txt

$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -I. -o test_sincos test_sincos.cpp
./test_sincos > test_sincos.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <bits/float128_io.h>

//namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
//{

  /**
   * A struct to store a cosine and a sine value.
   */
  template<typename _Tp>
    struct __sincos_t
    {
      _Tp sin_value;
      _Tp cos_value;
    };

//} // namespace __gnu_cxx

//namespace std _GLIBCXX_VISIBILITY(default)
//{
//// Implementation-space details.
//namespace __detail
//{
//_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Default implementation of sincos.
   */
  template<typename _Tp>
    inline __gnu_cxx::__sincos_t<_Tp>
    __sincos(_Tp __x)
    { return __gnu_cxx::__sincos_t<_Tp>{std::sin(__x), std::cos(__x)}; }

  template<>
    inline __gnu_cxx::__sincos_t<float>
    __sincos(float __x)
    {
      float __sin, __cos;
      __builtin_sincosf(__x, &__sin, &__cos);
      return __gnu_cxx::__sincos_t<float>{__sin, __cos};
    }

  template<>
    inline __gnu_cxx::__sincos_t<double>
    __sincos(double __x)
    {
      double __sin, __cos;
      __builtin_sincos(__x, &__sin, &__cos);
      return __gnu_cxx::__sincos_t<double>{__sin, __cos};
    }

  template<>
    inline __gnu_cxx::__sincos_t<long double>
    __sincos(long double __x)
    {
      long double __sin, __cos;
      __builtin_sincosl(__x, &__sin, &__cos);
      return __gnu_cxx::__sincos_t<long double>{__sin, __cos};
    }

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  template<>
    inline __gnu_cxx::__sincos_t<__float128>
    __sincos(__float128 __x)
    {
      __float128 __sin, __cos;
      ::sincosq(__x, &__sin, &__cos);
      return __gnu_cxx::__sincos_t<__float128>{__sin, __cos};
    }
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

  /**
   * Reperiodized sincos.
   */
  template<typename _Tp>
    __gnu_cxx::__sincos_t<_Tp>
    __sincos_pi(_Tp __x)
    {
      const auto _S_pi = __gnu_cxx::__const_pi<_Tp>(__x);
      const auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>(__x);
      if (std::isnan(__x))
	return __gnu_cxx::__sincos_t<_Tp>{_S_NaN, _S_NaN};
      else if (__x < _Tp{0})
	{
	  __gnu_cxx::__sincos_t<_Tp> __tempsc = __sincos_pi(-__x);
	  return __gnu_cxx::__sincos_t<_Tp>{-__tempsc.sin_value,
					     __tempsc.cos_value};
	}
      else if (__x < _Tp{0.5L})
	return __sincos(_S_pi * __x);
      else if (__x < _Tp{1})
	{
	  __gnu_cxx::__sincos_t<_Tp>
	    __tempsc = __sincos(_S_pi * (_Tp{1} - __x));
	  return __gnu_cxx::__sincos_t<_Tp>{__tempsc.sin_value,
					   -__tempsc.cos_value};
	}
      else
	{
	  auto __nu = std::floor(__x);
	  auto __arg = __x - __nu;
	  auto __sign = (int(__nu) & 1) == 1 ? _Tp{-1} : _Tp{+1};

	  auto __sinval = (__arg < _Tp{0.5L})
			? std::sin(_S_pi * __arg)
			: std::sin(_S_pi * (_Tp{1} - __arg));
	  auto __cosval = std::cos(_S_pi * __arg);
	  return __gnu_cxx::__sincos_t<_Tp>{__sign * __sinval,
					    __sign * __cosval};
	}
    }

  /**
   * Reperiodized complex constructor.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __polar_pi(_Tp __rho, _Tp __phi_pi)
    {
      __gnu_cxx::__sincos_t<_Tp> __sc = __sincos_pi(__phi_pi);
      return std::complex<_Tp>(__rho * __sc.cos_value, __rho * __sc.sin_value);
    }

//_GLIBCXX_END_NAMESPACE_VERSION
//} // namespace __detail
//} // namespace std


template<typename _Tp>
  void
  test_sincos(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto pi = __gnu_cxx::__const_pi<_Tp>(proto);

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sincos.sin"
	      << std::setw(width) << "sincos_pi.sin"
	      << std::setw(width) << "delta sin"
	      << std::setw(width) << "delta sin"
	      << std::setw(width) << "sincos.cos"
	      << std::setw(width) << "sincos_pi.cos"
	      << std::setw(width) << "delta cos"
	      << std::setw(width) << "delta cos"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    for (int i = -40; i <= +40; ++i)
      {
	auto x = _Tp(0.1Q * i);
	auto sc = std::__detail::__sincos(pi * x);
	auto scpi = std::__detail::__sincos_pi(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << sc.sin_value
		  << std::setw(width) << scpi.sin_value
		  << std::setw(width) << sc.sin_value - scpi.sin_value
		  << std::setw(width) << scpi.sin_value - std::sin(pi * x)
		  << std::setw(width) << sc.cos_value
		  << std::setw(width) << scpi.cos_value
		  << std::setw(width) << sc.cos_value - scpi.cos_value
		  << std::setw(width) << scpi.cos_value - std::cos(pi * x)
		  << '\n';
      }
  }

int
main()
{
  constexpr auto pif = __gnu_cxx__const_pi<float>();
  constexpr auto pi = __gnu_cxx::__const_pi<double>();
  constexpr auto pil = __gnu_cxx::__const_pi<long double>();
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  constexpr auto piq = __gnu_cxx::__const_pi<__float128>();
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

  auto a1 = /*std::__detail::*/__sincos(pif * 1.5f);
  auto a2 = /*std::__detail::*/__sincos_pi(1.5f);

  auto b1 = /*std::__detail::*/__sincos(pi * 1.5);
  auto b2 = /*std::__detail::*/__sincos_pi(1.5);

  auto c1 = /*std::__detail::*/__sincos(pil * 1.5l);
  auto c2 = /*std::__detail::*/__sincos_pi(1.5l);

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  auto d1 = /*std::__detail::*/__sincos(piq * 1.5q);
  auto d2 = /*std::__detail::*/__sincos_pi(1.5q);
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

  test_sincos<double>();
}
