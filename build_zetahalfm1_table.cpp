/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_zetahalfm1_table build_zetahalfm1_table.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./build_zetahalfm1_table > build_zetahalfm1_table.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_zetahalfm1_table build_zetahalfm1_table.cpp -lquadmath -lmpfr
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
  mpfr::mpreal half = mpfr::mpreal(1, prec) / mpfr::mpreal(2, prec);

  for (int i = 0; i <= 120; ++i)
    {
      auto x = mpfr::mpreal(i, prec) + half;
      std::cout << "  " << std::setw(w) << mpfr::zeta(x) - 1 << '\n';
    }
}
