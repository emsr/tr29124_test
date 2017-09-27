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
 * Compute the dual Hahn polynomial by recursion:
 * @f[
 *    \lambda(x)R_n(x) = A_n R_{n+1} - (A_n + C_n) R_n + C_n R_{n-1}
 * @f]
 * where @f$ R_n(x) = R_n(\lambda(x); \gamma, \delta, N) @f$ and
 * @f[
 *    A_n = (n + \gamma + 1)(n - N)
 * @f]
 * and
 * @f[
 *    C_n = n(n - \delta - N - 1)
 * @f]
 */
template<typename _Tp, typename _TpX>
  _Tp
  __dual_hahn_recur(int n, _Tp gamma, _Tp delta, int N, _TpX x)
  {
    auto Rnm1 = _Tp{1};
    if (n == 0)
      return Rnm1;

    const auto cd = gamma + delta;

    auto lambda = _Tp(x) * (_Tp(x) + cd + _Tp{1});

    auto Rn = _Tp{1} - lambda / (gamma + _Tp{1}) / _Tp(N);
    if (n == 1)
      return Rn;

    auto An = (gamma + _Tp{2}) * _Tp(1 - N);
    auto Cn = -(delta + _Tp(N));
    auto Rnp1 = ((An + Cn + lambda) * Rn - Cn * Rnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	auto An = (gamma + _Tp(k + 1)) * _Tp(k - N);
	auto Cn = -_Tp(k) * (delta + _Tp(N - k + 1));
	Rnm1 = Rn;
	Rn = Rnp1;
        Rnp1 = ((An + Cn + lambda) * Rn - Cn * Rnm1) / An;
      }

    return Rnp1;
  }

/**
 * Compute the dual Hahn polynomial defined by
 * @f[
 *    R_n(\lambda(x); \gamma, \delta, N)
 *       = {}_3F_2(-n, -x, x + \gamma + \delta + 1; \gamma + 1, -N; 1)
 * @f]
 * where @f$ \lambda(x) = x(x + \gamma + \delta + 1) @f$.
 */
template<typename _Tp, typename _TpX>
  _Tp
  __dual_hahn(int n, _Tp gamma, _Tp delta, int N, _TpX x)
  {
    if (std::isnan(gamma))
      return gamma;
    else if (std::isnan(delta))
      return delta;
    else if (std::isnan(x))
      return x;
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
