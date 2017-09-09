/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_charlier test_charlier.cpp
./test_charlier > test_charlier.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>

template<typename _Tp>
  _Tp
  __charlier_recur(int n, _Tp a, _Tp x)
  {
    auto Cnm2 = _Tp{1};
    if (n == 0)
      return Cnm2;

    auto Cnm1 = _Tp{1} - x / a;
    if (n == 1)
      return Cnm1;

    auto Cn = ((_Tp{1} + a - x) * Cnm1 - Cnm2) / a;
    for (int k = 3; k <= n; ++k)
      {
	Cnm2 = Cnm1;
	Cnm1 = Cn;
	Cn = ((_Tp(k - 1) + a - x) * Cnm1 - _Tp(k - 1) * Cnm2) / a;
      }

    return Cn;
  }

template<typename _Tp>
  void
  test_charlier()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    const auto N = 10;
    const auto a = _Tp{2};
    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto C = __charlier_recur(n, a, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << C
		      << '\n';
	  }
      }
  }

int
main()
{
  test_charlier<float>();
}
