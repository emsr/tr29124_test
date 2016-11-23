/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_fermi_dirac test_fermi_dirac.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
./test_fermi_dirac > test_fermi_dirac.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_fermi_dirac test_fermi_dirac.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
./test_fermi_dirac > test_fermi_dirac.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include "wrap_gsl.h"
#include <bits/float128_io.h>


template<typename _Tp>
  void
  run_fermi_dirac()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<_Tp> svec{_Tp{-0.5Q}, _Tp{0}, _Tp{0.5Q}, _Tp{1}, _Tp{1.5Q},
			  _Tp{2}, _Tp{3}, _Tp{4}, _Tp{5}};

    for (auto s : svec)
      {
	std::cout << " s = " << std::setw(width) << s << '\n';
	std::cout << '\n';
	for (int i = -250; i <= 250; ++i)
	  {
	    auto x = _Tp{0.1Q} * i;
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
