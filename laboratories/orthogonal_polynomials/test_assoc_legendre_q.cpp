/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <emsr/special_functions.h>

/**
 * 
 */
template<typename Tp>
  void
  test_assoc_legendre_q()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
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
		const auto x = Tp(0.01 * i);
		//const auto Q = emsr::assoc_legendre_q(l, m, x);
		const auto Q = emsr::detail::assoc_legendre_q(l, m, x);
		std::cout << ' ' << std::setw(w) << x
			  << ' ' << std::setw(w) << Q.Q_lm
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
