/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_etam1_table build_etam1_table.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./build_etam1_table > build_etam1_table.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_etam1_table build_etam1_table.cpp -lquadmath -lmpfr
*/

#include <mpreal.h>
#include "bits/numeric_limits_mpreal.h"

int
main()
{
  int prec = 128;
  mpfr::mpreal p(0, prec);
  std::cout.precision(__gnu_cxx::__max_digits10(p));
  auto w = 8 + std::cout.precision();

  for (int i = 0; i <= 120; ++i)
    {
      auto x = mpfr::mpreal(i, prec);
      auto pow2 = mpfr::pow(mpfr::mpreal(2, prec), 1 - x);
      std::cout << "  " << std::setw(w) << (1 - pow2) * mpfr::zeta(x) - 1 << '\n';
    }
}
