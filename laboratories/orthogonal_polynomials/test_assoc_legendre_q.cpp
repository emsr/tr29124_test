/*
$HOME/bin/bin/g++ -std=c++2a -g -Wall -Wextra -Wno-psabi -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -I../../quadrature/include -I../../cxx_summation/include -o test_assoc_legendre_q test_assoc_legendre_q.cpp -lquadmath
./test_assoc_legendre_q > test_assoc_legendre_q.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

/**
 * 
 */
template<typename _Tp>
  void
  test_assoc_legendre_q()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    for (auto l : {0, 1, 2, 3, 4 ,5})
      {
	for (auto m = 0; m <= l; ++m)
	  {
	    std::cout << "\n\n l = " << std::setw(2) << l
		      << "  m = " << std::setw(2) << m << '\n';
	    for (int i = -100; i <= 100; ++i)
	      {
		const auto x = _Tp(0.01 * i);
		//const auto Q = std::assoc_legendre_q(l, m, x);
		const auto Q = std::__detail::__assoc_legendre_q(l, m, x);
		std::cout << ' ' << std::setw(w) << x
			  << ' ' << std::setw(w) << Q.__Q_lm
			  << ' ' << std::setw(w) << Q.deriv()
			  << '\n';
	      }
	  }
      }
  }


int
main()
{
  test_assoc_legendre_q<float>();

  test_assoc_legendre_q<double>();

  test_assoc_legendre_q<long double>();

  //test_assoc_legendre_q<__float128>();

  return 0;
}