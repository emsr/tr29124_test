/**
 *
 */

#include <wrap_boost.h>
#include <wrap_gsl.h>
#include <bits/numeric_limits.h>
#include <cmath>
#include <ext/math_constants.h>

  /**
   * Try complex (and negative real) values. Nope.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __arith_geom_mean(std::complex<_Tp> __x, std::complex<_Tp> __y)
    {
      if (__x == std::complex<_Tp>{} || __y == std::complex<_Tp>{})
	return std::complex<_Tp>{};
      else
        {
	  const auto _S_eps = __gnu_cxx::__epsilon(std::real(__x)) / _Tp{2};
	  while (true)
  	    {
	      __y = std::sqrt(__y * std::exchange(__x, (__x + __y) / _Tp{2}));
	      if (std::abs(__x - __y) < _S_eps * std::abs(__x + __y))
		break;
	    }
	  return (__x + __y) / _Tp{2};
	}
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __arith_geom_mean(_Tp __x, _Tp __y)
    {
      if (__x == _Tp{0} || __y == _Tp{0})
	return _Tp{0};
      else
        {
	  const auto _S_eps = __gnu_cxx::__epsilon(__x) / _Tp{2};
	  while (true)
  	    {
	      __y = std::sqrt(__y * std::exchange(__x, (__x + __y) / _Tp{2}));
	      if (std::abs(__x - __y) < _S_eps * (__x + __y))
		break;
	    }
	  return (__x + __y) / _Tp{2};
	}
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __log_agm(_Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::math::__pi_v<_Tp>;
      constexpr auto _S_log_10 = __gnu_cxx::math::__ln_10_v<_Tp>;
      constexpr auto _S_log_2 = __gnu_cxx::math::__ln_2_v<_Tp>;
      const auto __p = std::numeric_limits<_Tp>::digits;
      const auto __b = std::numeric_limits<_Tp>::radix;
      const auto __n = std::ilogb(__x);
      const auto __m = __p / 2 - __n;
      const auto __s = __x * std::ldexp(_Tp{1}, __m);
      const auto __logb = (__b == 2 ? _S_log_2 : _S_log_10);
      return _S_pi / __arith_geom_mean(_Tp{1}, _Tp{4} / __s) /  _Tp{2}
	   - __m * __logb;
    }

/**
 * 
 */
template<typename _Tp>
  void
  test_arith_geom_mean(_Tp __proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(__proto));
    auto w = 8 + std::cout.precision();

    for (int i = 0; i < 100; ++i)
      {
	auto x = i * _Tp{1};
	for (int j = 0; j < 100; ++j)
	  {
	    auto y = j * _Tp{1};
	    auto agm = __arith_geom_mean(x, y);
	    std::cout << ' ' << x
		      << ' ' << y
		      << ' ' << std::setw(w) << agm << '\n';
	  }
      }
  }

template<typename _Tp>
  void
  test_arith_geom_mean_cmplx(_Tp __proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(__proto));
    auto w = 4 + 2 * (6 + std::cout.precision());

    const auto del = _Tp{0.0625};
    for (int ire = -10; ire <= +10; ++ire)
      for (int iim = -10; iim <= +10; ++iim)
	{
	  std::cout << '\n';
	  const auto x = std::complex<_Tp>(del * ire, del * iim);
	  for (int jre = -10; jre <= +10; ++jre)
	    for (int jim = -10; jim <= +10; ++jim)
	      {
		const auto y = std::complex<_Tp>(del * jre, del * jim);
		auto agm = __arith_geom_mean(x, y);
		std::cout << ' ' << x
			  << ' ' << y
			  << ' ' << std::setw(w) << agm << '\n';
	      }
	}
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_log(_Tp __proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(__proto));
    auto w = 8 + std::cout.precision();

    for (int i = 1; i < 100; ++i)
      {
	auto x = i * _Tp{1};
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << __log_agm(x)
		  << ' ' << std::setw(w) << std::log(x) << '\n';
      }
  }

int
main()
{
  test_arith_geom_mean(1.0);
  //test_arith_geom_mean_cmplx(1.0); Nope.
  test_log(1.0);
}

