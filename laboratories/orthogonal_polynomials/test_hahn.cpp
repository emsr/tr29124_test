/**
 *
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
template<typename Tp, typename _TpX>
  Tp
  hahn_recur(int n, Tp alpha, Tp beta, int N, _TpX x)
  {
    auto Qnm1 = Tp{1};
    if (n == 0)
      return Qnm1;

    const auto ab = alpha + beta;

    auto Qn = Tp{1} - Tp(x) * (ab + Tp{2}) / (alpha + Tp{1}) / Tp(N);
    if (n == 1)
      return Qn;

    auto An = (ab + Tp{2}) * (alpha + Tp{2}) * Tp(N - 1)
	    / (ab + Tp{4}) / (ab + Tp{3});
    auto Cn = (ab + Tp(N + 2)) * (beta + Tp{1})
	    / (ab + Tp{3}) / (ab + Tp{2});
    auto Qnp1 = ((An + Cn - Tp(x)) * Qn - Cn * Qnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	auto An = (ab + Tp(k + 1)) * (alpha + Tp(k + 1)) * Tp(N - k)
		/ (ab + Tp(2 * k + 2)) / (ab + Tp(2 * k + 1));
	auto Cn = Tp(k) * (ab + Tp(N + k + 1)) * (beta + Tp(k))
		/ (ab + Tp(2 * k + 1)) / (ab + Tp(2 * k));
	Qnm1 = Qn;
	Qn = Qnp1;
        Qnp1 = ((An + Cn - Tp(x)) * Qn - Cn * Qnm1) / An;
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
template<typename Tp, typename _TpX>
  Tp
  hahn(int n, Tp alpha, Tp beta, int N, _TpX x)
  {
    if (std::isnan(alpha))
      return alpha;
    else if (std::isnan(beta))
      return beta;
    else if (std::isnan(x))
      return x;
    else if (n > N)
      throw std::domain_error("hahn: Degree must be less than or equal to big N.");
    else
      return hahn_recur(n, alpha, beta, N, x);
  }

/**
 * 
 */
template<typename Tp>
  void
  test_hahn(Tp alpha, Tp beta, int N)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * Tp{0.05L};
	    auto Q = hahn(n, alpha, beta, N, x);
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
