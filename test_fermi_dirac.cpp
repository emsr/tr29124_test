/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_fermi_dirac test_fermi_dirac.cpp -lquadmath -Lwrappers -lwrap_gsl
./test_fermi_dirac > test_fermi_dirac.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_fermi_dirac test_fermi_dirac.cpp -lquadmath -Lwrappers -lwrap_gsl
./test_fermi_dirac > test_fermi_dirac.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include "wrap_gsl.h"
#include <bits/float128_io.h>


  /**
   * Return the Fermi-Dirac probability distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename _Tp>
    _Tp
    __fermi_dirac_pdf(_Tp __mu, _Tp __beta, _Tp __x)
    {
      _Tp{1} / (std::exp(__beta * (__x - __mu)) + _Tp{1});
    }

  /**
   * return the Fermi-Dirac cumulative distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename _Tp>
    _Tp
    __fermi_dirac_cdf(_Tp __mu, _Tp __beta, _Tp __x)
    {
      ;
    }


template<typename _Tp>
  void
  run_fermi_dirac()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<_Tp> svec{_Tp{-1}/_Tp{2}, _Tp{0}, _Tp{1}/_Tp{2}, _Tp{1},
			  _Tp{3}/_Tp{2}, _Tp{2}, _Tp{3}, _Tp{4}, _Tp{5}};

    for (auto s : svec)
      {
	std::cout << "\n\n\n s = " << std::setw(width) << s << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int i = -250; i <= 250; ++i)
	  {
	    auto x = del * i;
	    auto F = std::__detail::__fermi_dirac(s, x);
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
  std::cout << "\nfloat\n=====\n\n";
  run_fermi_dirac<float>();

  std::cout << "\ndouble\n======\n";
  run_fermi_dirac<double>();

  std::cout << "\nlong double\n===========\n";
  run_fermi_dirac<long double>();

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\n__float128\n==========\n";
  run_fermi_dirac<__float128>();
#endif
}
