/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_log_table build_log_table.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./build_log_table > build_log_table.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_log_table build_log_table.cpp -lquadmath -lmpfr
./build_log_table
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

  for (int i = 0; i <= 100; ++i)
    {
      auto x = mpfr::mpreal(i, prec);
      std::cout << "  " << std::setw(w) << mpfr::log(x) << '\n';
    }
}
