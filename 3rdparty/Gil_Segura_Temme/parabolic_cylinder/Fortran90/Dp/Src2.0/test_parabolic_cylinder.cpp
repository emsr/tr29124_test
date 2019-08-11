
#include "parabolic_cylinder.h"
#include <iostream>
#include <iomanip>

void
test_parabolic_cylinder()
{
  std::cout.precision(16);
  auto w = 6 + std::cout.precision();

  int mode = 0, ierr = 0;
  for (int ia = 0; ia <= 5; ++ia)
    {
      double a = ia * 0.5;
      std::cout << " a = " << a << '\n';
      for (int ix = 0; ix <= 10; ++ix)
	{
	  double x = 0.1 * ix;
	  double uaxx[2], vaxx[2];
	  parab_cyl(a, x, mode, uaxx, vaxx, &ierr);
	  std::cout << ' ' << std::setw(w) << x
		    << ' ' << std::setw(w) << uaxx[0]
		    << ' ' << std::setw(w) << uaxx[1]
		    << ' ' << std::setw(w) << vaxx[0]
		    << ' ' << std::setw(w) << vaxx[1]
		    << '\n';
	}
    }
}

int
main()
{
  test_parabolic_cylinder();
}
