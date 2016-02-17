// $HOME/bin_specfun/bin/g++ -std=c++1z -Wall -Wextra -o test_little_airy test_little_airy.cpp gsl_wrap.cpp -lgsl -lgslcblas 2> err.txt

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_little_airy > test_little_airy.txt

// g++ -std=c++14 -DNO_CBRT -DNO_LOGBQ -Wall -Wextra -o test_little_airy test_little_airy.cpp gsl_wrap.cpp -lgsl -lgslcblas 2> err.txt

// ./test_little_airy > test_little_airy.txt

#include <iostream>
#include <iomanip>
#include <cmath>

#include "airy.tcc"

#include "gsl_wrap.h"

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
