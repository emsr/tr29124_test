/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_bose_einstein test_bose_einstein.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
./test_bose_einstein > test_bose_einstein.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_bose_einstein test_bose_einstein.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
./test_bose_einstein > test_bose_einstein.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include <complex>
#include <bits/float128_io.h>


template<typename _Sp, typename _Tp>
  std::enable_if<std::is_floating_point_v<_Sp>, void>
  run_bose_einstein()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<_Sp> svec{_Sp{0.5L}, _Sp{1}, _Sp{1.5L},
			  _Sp{2}, _Sp{3}, _Sp{4}, _Sp{5}};

    for (auto s : svec)
      {
	std::cout << " s = " << std::setw(width) << s << '\n';
	std::cout << '\n';
	for (int i = -250; i <= 250; ++i)
	  {
	    auto x = _Tp{0.1L} * i;
	    auto G = std::__detail::__bose_einstein(s, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << G << '\n';

	  }
      }
  }


int
main()
{
  std::cout << "\nfloat\n=====\n\n";
  run_bose_einstein<float, float>();

  std::cout << "\ndouble\n======\n";
  run_bose_einstein<double, double>();

  std::cout << "\nlong double\n===========\n";
  run_bose_einstein<long double, long double>();

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\n__float128\n==========\n";
  run_bose_einstein<__float128, __float128>();
#endif
}
