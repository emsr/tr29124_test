/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_log_table build_log_table.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./build_log_table > build_log_table.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_log_table build_log_table.cpp -lquadmath -lmpfr
*/

#include <mpreal.h>
#include "bits/numeric_limits_mpreal.h"

int
main()
{
  mpfr::mpreal p(0, 128);
  std::cout.precision(__gnu_cxx::__max_digits10(p));
  auto w = 8 + std::cout.precision();

  for (int i = 0; i <= 100; ++i)
    {
      auto x = mpfr::mpreal(i, 128);
      std::cout << "  " << std::setw(w) << mpfr::log(x) << '\n';
    }
}
