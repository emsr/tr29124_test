/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_sincos_tables build_sincos_tables.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./build_sincos_tables > build_sincos_tables.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_sincos_tables build_sincos_tables.cpp -lquadmath -lmpfr
./build_sincos_tables
*/

#include <mpreal.h>
#include "bits/numeric_limits_mpreal.h"

int
main()
{
  std::size_t prec = 128;
  mpfr::mpreal p(0, prec);
  std::cout.precision(__gnu_cxx::__max_digits10(p));
  auto w = 8 + std::cout.precision();
  mpfr::mpreal pi4 = mpfr::const_pi(prec) / mpfr::mpreal(4 * 46, prec);

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
