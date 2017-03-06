/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_etam1_table build_etam1_table.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./build_etam1_table > build_etam1_table.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_etam1_table build_etam1_table.cpp -lquadmath -lmpfr
*/

#include <mpreal.h>
#include "bits/numeric_limits_mpreal.h"

int
main()
{
  mpfr::mpreal p(0, 128);
  std::cout.precision(__gnu_cxx::__max_digits10(p));
  auto w = 8 + std::cout.precision();

  for (int i = 0; i <= 120; ++i)
    {
      auto x = mpfr::mpreal(i, 128);
      std::cout << "  " << std::setw(w) << (1 - mpfr::pow(2, 1 - x)) * mpfr::zeta(x) - 1 << '\n';
    }
}
