/**
 *
 */

#include <cmath>
#include <iostream>

int
main()
{
  std::cout.precision(__gnu_cxx::__digits10(1.0));
  std::cout << std::scientific;

  double s = -1.0;
  std::cout << __gnu_cxx::dirichlet_eta(s) << '\n';
}
