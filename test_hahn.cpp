/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_hahn test_hahn.cpp
./test_hahn > test_hahn.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>

template<typename _Tp>
  _Tp
  __hahn_recur(int n, _Tp alpha, _Tp beta, int N, _Tp x)
  {
    auto Qnm2 = _Tp{1};
    if (n == 0)
      return Qnm2;

    const auto ab = alpha + beta;

    auto Qnm1 = _Tp{1} - (ab + _Tp{2}) * x / (alpha = _Tp{1}) / _Tp(N);
    if (n == 1)
      return Qnm1;

    auto An = (ab + _Tp{2}) * (alpha + _Tp{2}) * _Tp(N - 1)
	    / (ab + _Tp{4}) / (ab + _Tp{5});
    auto Cn = (ab + _Tp(N + 2)) * (beta + _Tp{1})
	    / (ab + _Tp{4}) / (ab + _Tp{5});
    auto Qn = ((An + Cn - x) * Qnm1 - Cn * Qnm2) / An;

    for (int k = 3; k <= n; ++k)
      {
	auto An = (ab + _Tp(k)) * (alpha + _Tp(k)) * _Tp(N - k + 1)
		/ (ab + _Tp(2 * k)) / (ab + _Tp(2 * k + 1));
	auto Cn = _Tp(k - 1) * (ab + _Tp(N + k)) * (beta + _Tp(k - 1))
		/ (ab + _Tp(2 * k)) / (ab + _Tp(2 * k + 1));
	Qnm2 = Qnm1;
	Qnm1 = Qn;
        Qn = ((An + Cn - x) * Qnm1 - Cn * Qnm2) / An;
      }

    return Qn;
  }

template<typename _Tp>
  void
  test_hahn()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    const auto N = 10;
    const auto alpha = _Tp{0.5};
    const auto beta = _Tp{0.5};
    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto Q = __hahn_recur(n, alpha, beta, N, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << Q
		      << '\n';
	  }
      }
  }

int
main()
{
  test_hahn<float>();
}
