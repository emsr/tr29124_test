/**
 *
 */

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

#include <emsr/float128_io.h>
#include <emsr/sf_polylog.h>

  /**
   * Return the Bose-Einstein probability distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename _Tp>
    _Tp
    bose_einstein_pdf(_Tp mu, _Tp beta, _Tp x)
    {
      _Tp{1} / (std::exp(beta * (x - mu)) - _Tp{1});
    }

  /**
   * Return the Bose-Einstein cumulative distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename _Tp>
    _Tp
    bose_einstein_p(_Tp mu, _Tp beta, _Tp x)
    {
      ;
    }

  /**
   * Return the Bose-Einstein cumulative distribution function (continuous).
   * @f[
   *    \int_{0}^{\infty}dx \frac{x^{p-1}}{e^x-1} = \Gamma(p) \zeta(p)
   * @f]
   */
  template<typename _Tp>
    _Tp
    bose_einstein_integral(_Tp mu, _Tp beta, _Tp x)
    {
      ;
    }


template<typename _Tp>
  void
  run_bose_einstein()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<_Tp> svec{_Tp{1}/_Tp{2}, _Tp{1}, _Tp{3}/_Tp{2},
			  _Tp{2}, _Tp{3}, _Tp{4}, _Tp{5}};

    for (auto s : svec)
      {
	std::cout << "\n\n\n s = " << std::setw(width) << s << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int i = -250; i <= 250; ++i)
	  {
	    auto x = del * i;
	    auto G = emsr::detail::bose_einstein(s, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << G << '\n';

	  }
      }
  }


int
main()
{
  run_bose_einstein<float>();

  run_bose_einstein<double>();

  run_bose_einstein<long double>();

// This works but takes too long.
#if 0 && !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //run_bose_einstein<__float128>();
#endif
}
