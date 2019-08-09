
#include "incomplete_gamma.h"
#include <iostream>
#include <iomanip>

void
test_incomplete_gamma()
{
  std::cout.precision(16);
  auto w = 6 + std::cout.precision();

  int ierr = 0;
  for (int ia = 1; ia <= 10; ++ia)
    {
      double a = ia * 0.5;
      std::cout << " a = " << a << '\n';
      for (int ix = 0; ix <= 10; ++ix)
	{
	  double x = 0.1 * ix;
	  double p, q;
	  inc_gamma(a, x, &p, &q, &ierr);
	  double xx = 0.0;
	  inv_inc_gamma (a, p, q, &xx, &ierr);
	  std::cout << ' ' << std::setw(w) << x
		    << ' ' << std::setw(w) << p
		    << ' ' << std::setw(w) << q
		    << ' ' << std::setw(w) << xx
		    << '\n';
	}
    }
}

int
main()
{
  test_incomplete_gamma();
}
