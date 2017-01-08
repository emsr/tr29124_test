/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bose_einstein test_bose_einstein.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
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


  /**
   * Return the Bose-Einstein probability distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename _Tp>
    _Tp
    __bose_einstein_pdf(_Tp __mu, _Tp __beta, _Tp __x)
    {
      _Tp{1} / (std::exp(__beta * (__x - __mu)) - _Tp{1});
    }

  /**
   * Return the Bose-Einstein cumulative distribution function (continuous).
   *
   * @see Kieth Oldham, Jan Myland, and Jerome Spanier,
   * An Atlas of Functions, 2nd Edition p. 266
   */
  template<typename _Tp>
    _Tp
    __bose_einstein_cdf(_Tp __mu, _Tp __beta, _Tp __x)
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
    __bose_einstein_integral(_Tp __mu, _Tp __beta, _Tp __x)
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
  run_bose_einstein<float>();

  std::cout << "\ndouble\n======\n";
  run_bose_einstein<double>();

  std::cout << "\nlong double\n===========\n";
  run_bose_einstein<long double>();

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\n__float128\n==========\n";
  run_bose_einstein<__float128>();
#endif
}
