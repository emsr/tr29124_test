#include "toroidal_harmonic.h"
#include <iostream>
#include <iomanip>

void
test_tor_harmonic()
{
  std::cout.precision(16);
  auto w = 6 + std::cout.precision();
  int ifac = 1;
  double P[6], Q[6];
  int newn = 5;

  for (int l = 0; l <= 5; ++l)
    {
      std::cout << " l = " << l << '\n';
      for (int m = 0; m <= 5; ++m)
	{
	  std::cout << " m = " << m << '\n';
	  for (int i = 1; i <= 20; ++i)
	    {
	      double x = 1.0 + i * 0.1;
	      tor_harmonic(x, m, l, P, Q, &newn);
	      std::cout << ' ' << x
			<< ' ' << std::setw(w) << P[l]
			<< ' ' << std::setw(w) << Q[l] << '\n';
	    }
	}
    }
}

int
main()
{
  test_tor_harmonic();
}
