/**
 *
 */

//  This is to track the bomb in the bessel and neumann at x=8, nu=0
//  When this was 10*min J_0 and Y_0 would tank at x = 8.
//const _Tp fp_min = _Tp(10) * std::numeric_limits<_Tp>::min();
//  I bumped to 20*min but I should watch this!

#include <iostream>
#include <tr1/cmath>

int
main()
{
  double nu = 0.0;
  double xpre = 79.9;
  double xpost = 80.0;
  double jpre [[maybe_unused]] = std::tr1::cyl_bessel_j(nu, xpre);
  double jpost [[maybe_unused]] = std::tr1::cyl_bessel_j(nu, xpost);

  for (int i = 1; i <= 10000; ++i)
    {
      double x = xpost + i * 0.1;
      double jcurr = std::tr1::cyl_bessel_j(nu, x);
      std::cout << "x = " << x << "  j = " << jcurr << '\n';
    }

  return 0;
}
