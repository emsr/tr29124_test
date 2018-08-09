/*
$HOME/bin/bin/g++ -std=c++2a -g -Wall -Wextra -o test_q_pochhammer test_q_pochhammer.cpp
./test_q_pochhammer > test_q_pochhammer.txt
*/

#include <cmath>
#include <iostream>
#include <iomanip>

template<typename Real>
  void
  test_q_pochhammer()
  {
    std::cout.precision(std::numeric_limits<Real>::digits10);
    const auto w = 8 + std::cout.precision();

    for (auto i = -10; i <= 10; ++i)
      {
	auto a = i / Real{2};
	std::cout << "\n  a = " << a << '\n';
	for (auto n = 1; n <= 10; ++n)
	  {
	    auto p = Real{1} / Real{2};
	    std::cout << "\n  n = " << n << '\n';
	    while (p > std::numeric_limits<Real>::epsilon())
	      {
		auto q = Real{1} - p;
		auto poch = Real{1};
		auto aqk = a;
		for (int k = 0; k < n; ++k)
		  {
		    poch *= (Real{1} - aqk);
		    aqk *= q;
		  }
		std::cout << ' ' << std::setw(w) << q
			  << ' ' << std::setw(w) << poch
			  << '\n';
		p /= 2;
	      }
	  }
      }
  }

int
main()
{
  test_q_pochhammer<double>();
}
