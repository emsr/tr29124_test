/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_sincos_tables build_sincos_tables.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./build_sincos_tables > build_sincos_tables.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_sincos_tables build_sincos_tables.cpp -lquadmath -lmpfr
*/

#include <mpreal.h>
#include "bits/numeric_limits_mpreal.h"

int
main()
{
  mpfr::mpreal p(0, 128);
  std::cout.precision(__gnu_cxx::__max_digits10(p));
  auto w = 8 + std::cout.precision();
  mpfr::mpreal pi4 = mpfr::const_pi(128) / mpfr::mpreal(4 * 46, 128);

  std::cout << "\n// sin table\n";
  for (int i = 0; i <= 46; ++i)
    {
      auto x = i * pi4;
      std::cout << "  " << std::setw(w) << mpfr::sin(x) << '\n';
    }

  std::cout << "\n// cos table\n";
  for (int i = 0; i <= 46; ++i)
    {
      auto x = i * pi4;
      std::cout << "  " << std::setw(w) << mpfr::cos(x) << '\n';
    }
}
