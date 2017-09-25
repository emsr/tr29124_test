/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_continuous_dual_hahn test_continuous_dual_hahn.cpp
./test_continuous_dual_hahn > test_continuous_dual_hahn.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

template<typename _Tp>
  struct
  __continuous_dual_hahn_t
  {
    _Tp __value;
    _Tp __factor;
  };

/**
 * 
 */
template<typename _Tp>
  _Tp
  __continuous_dual_hahn_recur(int n, _Tp a, _Tp b, _Tp c, _Tp x)
  {
    auto Snm1 = _Tp{1};
    if (n == 0)
      return Snm1;

    const auto aa = a * a;
    const auto ab = a + b;
    const auto ac = a + c;
    const auto bc = b + c;
    const auto xx = x * x;

    auto fact = _Tp{1};

    auto fn = ab * ac;
    fact *= fn;
    auto Sn = _Tp{1} - n * (aa + xx) / fn;
    if (n == 1)
      return Sn;

    fn = (ab + _Tp{1}) * (ac + _Tp{1});
    fact *= fn;
    auto An = fn;
    auto Cn = bc;

    auto Snp1 = ((An + Cn - (aa + xx)) * Sn - Cn * Snm1) / An;

    for (int k = 2; k < n; ++k)
      {
	fn = (ab + _Tp(k)) * (ac + _Tp(k));
	fact *= fn;
	auto An = fn;
	auto Cn = _Tp(k) * (bc + _Tp(k - 1));

	Snm1 = Sn;
	Sn = Snp1;
        Snp1 = ((An + Cn - (aa + xx)) * Sn - Cn * Snm1) / An;
      }

    return fact * Snp1;
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  __continuous_dual_hahn(int n, _Tp a, _Tp b, _Tp c, _Tp x)
  {
    if (std::isnan(a))
      return a;
    else if (std::isnan(b))
      return b;
    else if (std::isnan(c))
      return c;
    else
      return __continuous_dual_hahn_recur(n, a, b, c, x);
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_continuous_dual_hahn(int n_max, _Tp a, _Tp b, _Tp c)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 220; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto S = __continuous_dual_hahn(n, a, b, c, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << S
		      << '\n';
	  }
      }
  }

int
main()
{
  test_continuous_dual_hahn(10, 2.0f, 2.0f, 2.0f);
}
