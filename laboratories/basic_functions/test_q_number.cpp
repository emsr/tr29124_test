/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

template<typename Real>
  void
  test_q_number()
  {
    std::cout.precision(std::numeric_limits<Real>::digits10);
    const auto w = 8 + std::cout.precision();

    for (auto i = -10; i <= 10; ++i)
      {
	auto a = i / Real{2};
	auto p = Real{1} / Real{2};
	std::cout << "\n  a = " << a << '\n';
	while (p > std::numeric_limits<Real>::epsilon())
	  {
	    auto q = Real{1} - p;
	    auto qk = std::pow(q, a);
	    std::cout << ' ' << std::setw(w) << q
		      << ' ' << std::setw(w) << (Real{1} - qk) / p
		      << '\n';
	    p /= 2;
	  }
      }
  }

int
main()
{
  test_q_number<double>();
}
