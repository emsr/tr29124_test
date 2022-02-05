/**
 *
 */

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <wrap_gsl.h>

#include <emsr/float128_io.h>
#include <emsr/sf_polylog.h>

  /**
   * Return the Fermi-Dirac probability distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename Tp>
    Tp
    fermi_dirac_pdf(Tp mu, Tp beta, Tp x)
    {
      Tp{1} / (std::exp(beta * (x - mu)) + Tp{1});
    }

  /**
   * return the Fermi-Dirac cumulative distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename Tp>
    Tp
    fermi_dirac_p(Tp mu, Tp beta, Tp x)
    {
      ;
    }


template<typename Tp>
  void
  run_fermi_dirac()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<Tp> svec{Tp{-1}/Tp{2}, Tp{0}, Tp{1}/Tp{2}, Tp{1},
			  Tp{3}/Tp{2}, Tp{2}, Tp{3}, Tp{4}, Tp{5}};

    for (auto s : svec)
      {
	std::cout << "\n\n\n s = " << std::setw(width) << s << '\n';
	const auto del = Tp{1} / Tp{10};
	for (int i = -250; i <= 250; ++i)
	  {
	    auto x = del * i;
	    auto F = emsr::detail::fermi_dirac(s, x);
	    auto F_GSL = gsl::fermi_dirac(s, x);
	    auto err = (F - F_GSL) / std::abs(F_GSL);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << F
		      << ' ' << std::setw(width) << F_GSL
		      << ' ' << std::setw(width) << err << '\n';

	  }
      }
  }


int
main()
{
  run_fermi_dirac<float>();

  run_fermi_dirac<double>();

  run_fermi_dirac<long double>();

// This works but takes too long.
#ifdef EMSR_HAVE_FLOAT128
  //run_fermi_dirac<__float128>();
#endif
}
