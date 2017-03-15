/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_legendre_ellint test_legendre_ellint.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_legendre_ellint > test_legendre_ellint.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_legendre_ellint test_legendre_ellint.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
PATH=wrappers/debug:$PATH ./test_legendre_ellint > test_legendre_ellint.txt
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include "wrap_gsl.h"

template<typename _Tp>
  void
  test_comp_ellint_1(_Tp __proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(__proto));
    auto w = 8 + std::cout.precision();

    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i < 100; ++i)
      {
	_Tp k = i * (0.01L);
	auto K_gcc = std::comp_ellint_1(k);
	auto K_gsl = gsl::comp_ellint_1(k);
	std::cout << std::setw(5) << k
                  << ' ' << std::setw(w) << K_gcc
                  << ' ' << std::setw(w) << K_gsl
                  << ' ' << std::setw(w) << K_gcc - K_gsl
                  << '\n';
      }
  }

int
main()
{
  test_comp_ellint_1(1.0);
}
