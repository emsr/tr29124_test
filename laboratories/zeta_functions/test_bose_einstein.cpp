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
  template<typename Tp>
    Tp
    bose_einstein_pdf(Tp mu, Tp beta, Tp x)
    {
      Tp{1} / (std::exp(beta * (x - mu)) - Tp{1});
    }

  /**
   * Return the Bose-Einstein cumulative distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename Tp>
    Tp
    bose_einstein_p(Tp mu, Tp beta, Tp x)
    {
      ;
    }

  /**
   * Return the Bose-Einstein cumulative distribution function (continuous).
   * @f[
   *    \int_{0}^{\infty}dx \frac{x^{p-1}}{e^x-1} = \Gamma(p) \zeta(p)
   * @f]
   */
  template<typename Tp>
    Tp
    bose_einstein_integral(Tp mu, Tp beta, Tp x)
    {
      ;
    }


template<typename Tp>
  void
  run_bose_einstein()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<Tp> svec{Tp{1}/Tp{2}, Tp{1}, Tp{3}/Tp{2},
			  Tp{2}, Tp{3}, Tp{4}, Tp{5}};

    for (auto s : svec)
      {
	std::cout << "\n\n\n s = " << std::setw(width) << s << '\n';
	const auto del = Tp{1} / Tp{10};
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
#ifdef EMSR_HAVE_FLOAT128
  //run_bose_einstein<__float128>();
#endif
}
