/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_dual_hahn test_dual_hahn.cpp
./test_dual_hahn > test_dual_hahn.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

/**
 * 
 */
template<typename _Tp>
  _Tp
  __dual_hahn_recur(int n, _Tp gamma, _Tp delta, int N, _Tp x)
  {
    auto Rnm2 = _Tp{1};
    if (n == 0)
      return Rnm2;

    const auto cd = gamma + delta;

    auto lambda = x * (x + cd + _Tp{1});

    auto Rnm1 = _Tp{1} - lambda / (gamma + _Tp{1}) / _Tp(N);
    if (n == 1)
      return Rnm1;

    auto An = (gamma + _Tp{2}) * _Tp(1 - N);
    auto Cn = -(delta + _Tp(N));
    auto Rn = ((An + Cn + lambda) * Rnm1 - Cn * Rnm2) / An;

    for (int k = 3; k <= n; ++k)
      {
	auto An = (gamma + _Tp(k)) * _Tp(k - N - 1);
	auto Cn = -_Tp(k - 1) * (delta + _Tp(N - k + 2));
	Rnm2 = Rnm1;
	Rnm1 = Rn;
        Rn = ((An + Cn + lambda) * Rnm1 - Cn * Rnm2) / An;
      }

    return Rn;
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  __dual_hahn(int n, _Tp gamma, _Tp delta, int N, _Tp x)
  {
    if (std::isnan(gamma))
      return gamma;
    else if (std::isnan(delta))
      return delta;
    else if (n > N)
      std::__throw_domain_error(__N("__dual_hahn: "
				"Degree must be less than or equal to big N."));
    else
      return __dual_hahn_recur(n, gamma, delta, N, x);
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_dual_hahn(_Tp gamma, _Tp delta, int N)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    auto lambda = [gamma, delta, N](_Tp x)
		  { return x * (x + gamma + delta + _Tp{1}); };

    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 200; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto R = __dual_hahn(n, gamma, delta, N, x);
	    std::cout << ' ' << std::setw(w) << lambda(x)
		      << ' ' << std::setw(w) << R
		      << '\n';
	  }
      }
  }

int
main()
{
  test_dual_hahn(1.0f/3.0f, 1.0f/2.0f, 10);
}
