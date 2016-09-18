/*
$HOME/bin_tr29124/bin/g++ -std=c++17 -Wall -Wextra -o test_little_airy test_little_airy.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas 2> err.txt
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_little_airy > test_little_airy.txt

g++ -std=c++14 -DNO_CBRT -DNO_LOGBQ -Wall -Wextra -o test_little_airy test_little_airy.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas 2> err.txt
./test_little_airy > test_little_airy.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "airy.tcc"

#include "wrap_gsl.h"

int
main()
{
  std::complex<double> t;
  std::complex<double> Ai, Aip, Bi, Bip;
  std::complex<double> w1, w1p, w2, w2p;

  std::cout.precision(16);
  std::cout.flags(std::ios::showpoint);

  for (auto i = -200; i <= +200; ++i)
    {
      t = i * 0.05;
      auto Ai_gsl = gsl::airy_ai(std::real(t));
      auto Bi_gsl = gsl::airy_bi(std::real(t));
      airy(t, Ai, Aip, Bi, Bip, w1, w1p, w2, w2p);
      std::cout << ' ' << std::setw(8) << std::real(t)
                << ' ' << std::setw(20) << std::real(Ai)
                << ' ' << std::setw(20) << std::real(Bi)
                << ' ' << std::setw(20) << Ai_gsl
                << ' ' << std::setw(20) << Bi_gsl
                << ' ' << std::setw(20) << std::real(Ai - Ai_gsl)
                << ' ' << std::setw(20) << std::real(Bi - Bi_gsl)
                << '\n';
    }
}
