/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_meixner test_meixner.cpp
./test_meixner > test_meixner.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>

template<typename _Tp>
  _Tp
  __meixner_recur(int n, _Tp beta, _Tp c, _Tp x)
  {
    auto Mnm2 = _Tp{1};
    if (n == 0)
      return Mnm2;

    auto Mnm1 = _Tp{1} + (_Tp{1} - _Tp{1} / c) * x / beta;
    if (n == 1)
      return Mnm1;

    auto nm1 = 1;
    auto cbnm1 = c * (beta + nm1);
    auto Mn = (((c - _Tp{1}) * x + nm1 + cbnm1) * Mnm1 - nm1 * Mnm2) / cbnm1;
    for (int k = 3; k <= n; ++k)
      {
	nm1 = _Tp(k - 1);
	cbnm1 = c * (beta + nm1);
	Mnm2 = Mnm1;
	Mnm1 = Mn;
	Mn = (((c - _Tp{1}) * x + nm1 + cbnm1) * Mnm1 - nm1 * Mnm2) / cbnm1;
      }

    return Mn;
  }

template<typename _Tp>
  void
  test_meixner()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    const auto N = 10;
    const auto beta = _Tp{0.5};
    const auto c = _Tp{2.0};
    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 100; ++i)
	  {
	    auto x = i * _Tp{0.1L};
	    auto M = __meixner_recur(n, beta, c, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << M
		      << '\n';
	  }
      }
  }

int
main()
{
  test_meixner<float>();
}
