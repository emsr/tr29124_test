/**
 *
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
  test_assoc_legendre()
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
	    for (int i = -120; i <= 120; ++i)
	      {
		const auto x = _Tp(0.01 * i);
		//const auto P = std::assoc_legendre(l, m, x);
		const auto P = std::__detail::__assoc_legendre_p(l, m, x);
		std::cout << ' ' << std::setw(w) << x
			  << ' ' << std::setw(w) << P.__P_lm
			  << ' ' << std::setw(w) << P.deriv()
			  << '\n';
	      }
	  }
      }
  }


int
main()
{
  test_assoc_legendre<float>();

  test_assoc_legendre<double>();

  test_assoc_legendre<long double>();

  //test_assoc_legendre<__float128>();

  return 0;
}
