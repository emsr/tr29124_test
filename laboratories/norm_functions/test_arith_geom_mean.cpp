/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

#include <wrap_boost.h>
#include <wrap_gsl.h>
#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

  /**
   * Try complex (and negative real) values. Nope.
   */
  template<typename _Tp>
    std::complex<_Tp>
    arith_geom_mean(std::complex<_Tp> x, std::complex<_Tp> y)
    {
      if (x == std::complex<_Tp>{} || y == std::complex<_Tp>{})
	return std::complex<_Tp>{};
      else
        {
	  const auto _S_eps = emsr::epsilon(std::real(x)) / _Tp{2};
	  while (true)
  	    {
	      y = std::sqrt(y * std::exchange(x, (x + y) / _Tp{2}));
	      if (std::abs(x - y) < _S_eps * std::abs(x + y))
		break;
	    }
	  return (x + y) / _Tp{2};
	}
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    arith_geom_mean(_Tp x, _Tp y)
    {
      if (x == _Tp{0} || y == _Tp{0})
	return _Tp{0};
      else
        {
	  const auto _S_eps = emsr::epsilon(x) / _Tp{2};
	  while (true)
  	    {
	      y = std::sqrt(y * std::exchange(x, (x + y) / _Tp{2}));
	      if (std::abs(x - y) < _S_eps * (x + y))
		break;
	    }
	  return (x + y) / _Tp{2};
	}
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    log_agm(_Tp x)
    {
      constexpr auto _S_pi = emsr::pi_v<_Tp>;
      constexpr auto _S_log_10 = emsr::ln10_v<_Tp>;
      constexpr auto _S_log_2 = emsr::ln2_v<_Tp>;
      const auto p = std::numeric_limits<_Tp>::digits;
      const auto b = std::numeric_limits<_Tp>::radix;
      const auto n = std::ilogb(x);
      const auto m = p / 2 - n;
      const auto s = x * std::ldexp(_Tp{1}, m);
      const auto logb = (b == 2 ? _S_log_2 : _S_log_10);
      return _S_pi / arith_geom_mean(_Tp{1}, _Tp{4} / s) /  _Tp{2}
	   - m * logb;
    }

/**
 * 
 */
template<typename _Tp>
  void
  test_arith_geom_mean(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = 8 + std::cout.precision();

    for (int i = 0; i < 100; ++i)
      {
	auto x = i * _Tp{1};
	for (int j = 0; j < 100; ++j)
	  {
	    auto y = j * _Tp{1};
	    auto agm = arith_geom_mean(x, y);
	    std::cout << ' ' << x
		      << ' ' << y
		      << ' ' << std::setw(w) << agm << '\n';
	  }
      }
  }

template<typename _Tp>
  void
  test_arith_geom_mean_cmplx(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
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
		auto agm = arith_geom_mean(x, y);
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
  test_log(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = 8 + std::cout.precision();

    for (int i = 1; i < 100; ++i)
      {
	auto x = i * _Tp{1};
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << log_agm(x)
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

