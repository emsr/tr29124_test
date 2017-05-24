/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_arith_geom_mean test_arith_geom_mean.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_arith_geom_mean > test_arith_geom_mean.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_arith_geom_mean test_arith_geom_mean.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
PATH=wrappers/debug:$PATH ./test_arith_geom_mean > test_arith_geom_mean.txt
*/

#include <wrap_boost.h>
#include <wrap_gsl.h>
#include <bits/numeric_limits.h>

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

int
main()
{
  test_arith_geom_mean(1.0);
  //test_arith_geom_mean_cmplx(1.0); Nope.
}

