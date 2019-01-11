/*
$HOME/bin_specfun/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_chebyshev test_chebyshev.cpp -lquadmath
./test_chebyshev > test_chebyshev.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_chebyshev test_chebyshev.cpp -lquadmath
./test_chebyshev > test_chebyshev.txt
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

template<typename Tp>
  void
  test_chebyshev(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::vector<unsigned int> index{0, 1, 2, 3, 4, 5, 10, 20, 23, 50, 100};

    for (auto n : index)
      {
	std::cout << "\n\n n = " << std::setw(width) << n << '\n';
	std::cout << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "T"
		  << ' ' << std::setw(width) << "T'"
		  << ' ' << std::setw(width) << "U"
		  << ' ' << std::setw(width) << "U'"
		  << ' ' << std::setw(width) << "V"
		  << ' ' << std::setw(width) << "V'"
		  << ' ' << std::setw(width) << "W"
		  << ' ' << std::setw(width) << "W'"
		  << '\n';
	const auto del = Tp{1} / Tp{100};
	for (int i = -100; i <= 100; ++i)
	  {
	    auto x = del * i;
	    auto Ts = std::__detail::__chebyshev_t(n, x);
	    auto Us = std::__detail::__chebyshev_u(n, x);
	    auto Vs = std::__detail::__chebyshev_v(n, x);
	    auto Ws = std::__detail::__chebyshev_w(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << Ts.__T_n
		      << ' ' << std::setw(width) << Ts.deriv()
		      << ' ' << std::setw(width) << Us.__U_n
		      << ' ' << std::setw(width) << Us.deriv()
		      << ' ' << std::setw(width) << Vs.__V_n
		      << ' ' << std::setw(width) << Vs.deriv()
		      << ' ' << std::setw(width) << Ws.__W_n
		      << ' ' << std::setw(width) << Ws.deriv()
		      << '\n';
	  }
      }
  }
int
main()
{
  test_chebyshev(1.0);
}
