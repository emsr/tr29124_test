/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <emsr/sf_legendre.h>

template<typename _Tp>
  void
  test_legendre_q()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    for (unsigned int l : {0, 1, 2, 3, 4, 5})
      {
	std::cout << "\n\n l = " << std::setw(2) << l << '\n';
	for (int i = -120; i <= +120; ++i)
	  {
	    auto x = _Tp(0.01L * i);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << emsr::legendre_q(l, x)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_legendre_q<float>();

  test_legendre_q<double>();

  test_legendre_q<long double>();

  //test_legendre_q<__float128>(1.0Q);

  return 0;
}
