/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_bernoulli_2n_table build_bernoulli_2n_table.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./build_bernoulli_2n_table > build_bernoulli_2n_table.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o build_bernoulli_2n_table build_bernoulli_2n_table.cpp -lquadmath -lmpfr
./build_bernoulli_2n_table > build_bernoulli_2n_table.txt
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
  auto _S_2pi = 2 * mpfr::const_pi(prec);
  std::cout << "  2pi = " << std::setw(w) << _S_2pi << '\n';

  auto fact = mpfr::mpreal(2, prec);
  for (int n = 2; n <= 200; n += 2)
    {
      if ((n / 2) % 2 == 0)
	fact *= -1;
      fact *= mpfr::mpreal(n - 1, prec) / _S_2pi;
      fact *= mpfr::mpreal(n, prec) / _S_2pi;

      std::cout << "  " << std::setw(w) << fact * mpfr::zeta(mpfr::mpreal(n, prec)) << '\n';
    }
}

