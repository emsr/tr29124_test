/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -I../../quadrature/include -I../../cxx_summation/include -o test_legendre_q test_legendre_q.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:LD_LIBRARY_PATH$ ./test_legendre_q > test_legendre_q.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

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
		      << ' ' << std::setw(w) << __gnu_cxx::legendre_q(l, x)
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
