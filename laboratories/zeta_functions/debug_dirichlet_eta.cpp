/**
 *
 */

#include <cmath>
#include <iostream>

#include <emsr/numeric_limits.h>
#include <emsr/sf_polylog.h>

int
main()
{
  std::cout.precision(emsr::digits10(1.0));
  std::cout << std::scientific;

  double s = -1.0;
  std::cout << emsr::dirichlet_eta(s) << '\n';
}
