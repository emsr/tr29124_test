/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

int
main()
{
  std::cout.precision(emsr::digits10(1.0));
  std::cout << std::scientific;
  const auto w = 8 + std::cout.precision();

  for (int i = -500; i <= 500; ++i)
    {
      auto s = i * (0.01);
      std::cout << ' ' << s
		<< ' ' << std::setw(w) << emsr::riemann_zeta(s)
		<< ' ' << std::setw(w) << emsr::dirichlet_eta(s)
		<< ' ' << std::setw(w) << emsr::dirichlet_beta(s)
		<< ' ' << std::setw(w) << emsr::dirichlet_lambda(s)
		<< '\n';
    }
}
