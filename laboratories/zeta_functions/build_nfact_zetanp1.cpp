/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_nfact_zetanp1 build_nfact_zetanp1.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./build_nfact_zetanp1 > build_nfact_zetanp1.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_nfact_zetanp1 build_nfact_zetanp1.cpp -lquadmath -lmpfr
*/

#include <mpreal.h>
#include "bits/numeric_limits_mpreal.h"

// List of n! zeta(n+1)
int
main()
{
  std::size_t prec = 128;
  mpfr::mpreal p(0, prec);
  std::cout.precision(__gnu_cxx::__max_digits10(p));
  auto w = 8 + std::cout.precision();

  mpfr::mpreal fact(1, prec);
  for (int i = 0; i <= 120; ++i)
    {
      auto x = mpfr::mpreal(i + 1, prec);
      std::cout << "  " << std::setw(w) << fact * mpfr::zeta(x) << '\n';
      fact *= (i + 1);
    }
}
