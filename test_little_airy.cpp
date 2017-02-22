/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_little_airy test_little_airy.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_little_airy > test_little_airy.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_little_airy test_little_airy.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
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
