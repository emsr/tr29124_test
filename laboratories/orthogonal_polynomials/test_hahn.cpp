/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_hahn test_hahn.cpp
./test_hahn > test_hahn.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

/**
 * Compute the Hahn polynomial by recursion:
 * @f[
 *    -xQ_n(x) = A_n Q_{n+1} - (A_n + C_n) Q_n + C_n Q_{n-1}
 * @f]
 * where @f$ Q_n(x) = Q_n(x;\alpha,\beta,N) @f$ and
 * @f[
 *    A_n = \frac{(n + \alpha + \beta + 1)(n + \alpha + 1)(N - n)}
 *               {(2n + \alpha + \beta + 1)(2n + \alpha + \beta + 2)}
 * @f]
 * and
 * @f[
 *    C_n = \frac{n(n + \alpha + \beta)(n + \beta)}
 *               {(2n + \alpha + \beta)(2n + \alpha + \beta + 1)}
 * @f]
 */
template<typename _Tp, typename _TpX>
  _Tp
  __hahn_recur(int n, _Tp alpha, _Tp beta, int N, _TpX x)
  {
    auto Qnm1 = _Tp{1};
    if (n == 0)
      return Qnm1;

    const auto ab = alpha + beta;

    auto Qn = _Tp{1} - _Tp(x) * (ab + _Tp{2}) / (alpha + _Tp{1}) / _Tp(N);
    if (n == 1)
      return Qn;

    auto An = (ab + _Tp{2}) * (alpha + _Tp{2}) * _Tp(N - 1)
	    / (ab + _Tp{4}) / (ab + _Tp{3});
    auto Cn = (ab + _Tp(N + 2)) * (beta + _Tp{1})
	    / (ab + _Tp{3}) / (ab + _Tp{2});
    auto Qnp1 = ((An + Cn - _Tp(x)) * Qn - Cn * Qnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	auto An = (ab + _Tp(k + 1)) * (alpha + _Tp(k + 1)) * _Tp(N - k)
		/ (ab + _Tp(2 * k + 2)) / (ab + _Tp(2 * k + 1));
	auto Cn = _Tp(k) * (ab + _Tp(N + k + 1)) * (beta + _Tp(k))
		/ (ab + _Tp(2 * k + 1)) / (ab + _Tp(2 * k));
	Qnm1 = Qn;
	Qn = Qnp1;
        Qnp1 = ((An + Cn - _Tp(x)) * Qn - Cn * Qnm1) / An;
      }

    return Qnp1;
  }

/**
 * Compute the Hahn polynomial defined by
 * @f[
 *    Q_n(x;\alpha,\beta,N)
 *     = {}_3F_2(-n, n + \alpha + \beta + 1, -x; \alpha + 1, -N; 1)
 * @f]
 */
template<typename _Tp, typename _TpX>
  _Tp
  __hahn(int n, _Tp alpha, _Tp beta, int N, _TpX x)
  {
    if (std::isnan(alpha))
      return alpha;
    else if (std::isnan(beta))
      return beta;
    else if (std::isnan(x))
      return x;
    else if (n > N)
      std::__throw_domain_error(__N("__hahn: "
				"Degree must be less than or equal to big N."));
    else
      return __hahn_recur(n, alpha, beta, N, x);
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_hahn(_Tp alpha, _Tp beta, int N)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto Q = __hahn(n, alpha, beta, N, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << Q
		      << '\n';
	  }
      }
  }

int
main()
{
  test_hahn(2.0f, 2.0f, 10);
}
