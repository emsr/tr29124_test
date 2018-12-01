/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_zeta_trig_tables build_zeta_trig_tables.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./build_zeta_trig_tables > build_zeta_trig_tables.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../mpreal -o build_zeta_trig_tables build_zeta_trig_tables.cpp -lquadmath -lmpfr
PATH=$HOME/bin/lib64:$PATH ./build_zeta_trig_tables > build_zeta_trig_tables.txt
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

  // zeta_even
  std::cout << "\nzeta_even\n";
  for (int i = 0; i <= 20; i += 2)
    {
      auto x = mpfr::mpreal(i, prec);
      std::cout << "  " << std::setw(w) << mpfr::zeta(x) << '\n';
    }

  // eta_even
  std::cout << "\neta_even\n";
  for (int i = 0; i <= 20; i += 2)
    {
      auto x = mpfr::mpreal(i, prec);
      auto pow2 = mpfr::pow(mpfr::mpreal(2, prec), 1 - x);
      std::cout << "  " << std::setw(w) << (1 - pow2) * mpfr::zeta(x) << '\n';
    }

  // lambda_even
  std::cout << "\nlambda_even\n";
  for (int i = 0; i <= 20; i += 2)
    {
      auto x = mpfr::mpreal(i, prec);
      auto pow2 = mpfr::pow(mpfr::mpreal(2, prec), -x);
      std::cout << "  " << std::setw(w) << (1 - pow2) * mpfr::zeta(x) << '\n';
    }

  // beta_odd
  std::cout << "\nbeta_odd\n";
  for (int i = 1; i <= 21; i += 2)
    {
      auto x = mpfr::mpreal(i, prec);
    }
}
