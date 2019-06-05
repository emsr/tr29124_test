/*
$HOME/bin_specfun/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -Izeta -o test_dirichlet_eta test_dirichlet_eta.cpp -lquadmath -Lwrappers/debug -lwrap_cephes
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_dirichlet_eta > test_dirichlet_eta.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -Izeta -o test_dirichlet_eta test_dirichlet_eta.cpp -lquadmath -Lwrappers/debug -lwrap_cephes
PATH=wrappers/debug:$PATH ./test_dirichlet_eta > test_dirichlet_eta.txt
*/

#include <ext/cmath>
#include <iostream>
#include <iomanip>

int
main()
{
  std::cout.precision(__gnu_cxx::__digits10(1.0));
  std::cout << std::scientific;
  const auto w = 8 + std::cout.precision();

  for (int i = -500; i <= 500; ++i)
    {
      auto s = i * (0.01);
      std::cout << ' ' << s
		<< ' ' << std::setw(w) << std::riemann_zeta(s)
		<< ' ' << std::setw(w) << __gnu_cxx::dirichlet_eta(s)
		<< ' ' << std::setw(w) << __gnu_cxx::dirichlet_beta(s)
		<< ' ' << std::setw(w) << __gnu_cxx::dirichlet_lambda(s)
		<< '\n';
    }
}
