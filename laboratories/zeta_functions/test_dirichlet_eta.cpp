/**
 *
 */

#include <cmath>
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
