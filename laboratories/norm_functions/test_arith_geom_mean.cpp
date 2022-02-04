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
  template<typename Tp>
    std::complex<Tp>
    arith_geom_mean(std::complex<Tp> x, std::complex<Tp> y)
    {
      if (x == std::complex<Tp>{} || y == std::complex<Tp>{})
	return std::complex<Tp>{};
      else
        {
	  const auto _S_eps = emsr::epsilon(std::real(x)) / Tp{2};
	  while (true)
  	    {
	      y = std::sqrt(y * std::exchange(x, (x + y) / Tp{2}));
	      if (std::abs(x - y) < _S_eps * std::abs(x + y))
		break;
	    }
	  return (x + y) / Tp{2};
	}
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    arith_geom_mean(Tp x, Tp y)
    {
      if (x == Tp{0} || y == Tp{0})
	return Tp{0};
      else
        {
	  const auto _S_eps = emsr::epsilon(x) / Tp{2};
	  while (true)
  	    {
	      y = std::sqrt(y * std::exchange(x, (x + y) / Tp{2}));
	      if (std::abs(x - y) < _S_eps * (x + y))
		break;
	    }
	  return (x + y) / Tp{2};
	}
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    log_agm(Tp x)
    {
      constexpr auto _S_pi = emsr::pi_v<Tp>;
      constexpr auto _S_log_10 = emsr::ln10_v<Tp>;
      constexpr auto _S_log_2 = emsr::ln2_v<Tp>;
      const auto p = std::numeric_limits<Tp>::digits;
      const auto b = std::numeric_limits<Tp>::radix;
      const auto n = std::ilogb(x);
      const auto m = p / 2 - n;
      const auto s = x * std::ldexp(Tp{1}, m);
      const auto logb = (b == 2 ? _S_log_2 : _S_log_10);
      return _S_pi / arith_geom_mean(Tp{1}, Tp{4} / s) /  Tp{2}
	   - m * logb;
    }

/**
 * 
 */
template<typename Tp>
  void
  test_arith_geom_mean(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = 8 + std::cout.precision();

    for (int i = 0; i < 100; ++i)
      {
	auto x = i * Tp{1};
	for (int j = 0; j < 100; ++j)
	  {
	    auto y = j * Tp{1};
	    auto agm = arith_geom_mean(x, y);
	    std::cout << ' ' << x
		      << ' ' << y
		      << ' ' << std::setw(w) << agm << '\n';
	  }
      }
  }

template<typename Tp>
  void
  test_arith_geom_mean_cmplx(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = 4 + 2 * (6 + std::cout.precision());

    const auto del = Tp{0.0625};
    for (int ire = -10; ire <= +10; ++ire)
      for (int iim = -10; iim <= +10; ++iim)
	{
	  std::cout << '\n';
	  const auto x = std::complex<Tp>(del * ire, del * iim);
	  for (int jre = -10; jre <= +10; ++jre)
	    for (int jim = -10; jim <= +10; ++jim)
	      {
		const auto y = std::complex<Tp>(del * jre, del * jim);
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
template<typename Tp>
  void
  test_log(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = 8 + std::cout.precision();

    for (int i = 1; i < 100; ++i)
      {
	auto x = i * Tp{1};
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

