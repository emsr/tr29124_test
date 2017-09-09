/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_krawtchouk test_krawtchouk.cpp
./test_krawtchouk > test_krawtchouk.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>

template<typename _Tp>
  _Tp
  __krawtchouk_recur(int n, _Tp p, int N, _Tp x)
  {
    auto Knm2 = _Tp{1};
    if (n == 0)
      return Knm2;

    auto Knm1 = _Tp{1} - x / p / _Tp(N);
    if (n == 1)
      return Knm1;

    auto pnn = p * _Tp(N - 1);
    auto nm1q = _Tp{1} * (_Tp{1} - p);
    auto Kn = ((pnn + nm1q - x) * Knm1 - nm1q * Knm2) / pnn;
    for (int k = 3; k <= n; ++k)
      {
	pnn = p * _Tp(N - k + 1);
	nm1q = _Tp(k - 1) * (_Tp{1} - p);
	Kn = ((pnn + nm1q - x) * Knm1 - nm1q * Knm2) / pnn;
      }

    return Kn;
  }

template<typename _Tp>
  void
  test_krawtchouk()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    const auto N = 10;
    const auto p = _Tp{0.5};
    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto K = __krawtchouk_recur(n, p, N, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << K
		      << '\n';
	  }
      }
  }

int
main()
{
  test_krawtchouk<float>();
}
