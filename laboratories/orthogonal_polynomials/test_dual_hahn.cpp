/**
 *
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
template<typename Tp, typename _TpX>
  Tp
  dual_hahn_recur(int n, Tp gamma, Tp delta, int N, _TpX x)
  {
    auto Rnm1 = Tp{1};
    if (n == 0)
      return Rnm1;

    const auto cd = gamma + delta;

    auto lambda = Tp(x) * (Tp(x) + cd + Tp{1});

    auto Rn = Tp{1} - lambda / (gamma + Tp{1}) / Tp(N);
    if (n == 1)
      return Rn;

    auto An = (gamma + Tp{2}) * Tp(1 - N);
    auto Cn = -(delta + Tp(N));
    auto Rnp1 = ((An + Cn + lambda) * Rn - Cn * Rnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	auto An = (gamma + Tp(k + 1)) * Tp(k - N);
	auto Cn = -Tp(k) * (delta + Tp(N - k + 1));
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
template<typename Tp, typename _TpX>
  Tp
  dual_hahn(int n, Tp gamma, Tp delta, int N, _TpX x)
  {
    if (std::isnan(gamma))
      return gamma;
    else if (std::isnan(delta))
      return delta;
    else if (std::isnan(x))
      return x;
    else if (n > N)
      throw std::domain_error("dual_hahn: Degree must be less than or equal to big N.");
    else
      return dual_hahn_recur(n, gamma, delta, N, x);
  }

/**
 * 
 */
template<typename Tp>
  void
  test_dual_hahn(Tp gamma, Tp delta, int N)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    auto lambda = [gamma, delta, N](Tp x)
		  { return x * (x + gamma + delta + Tp{1}); };

    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 200; ++i)
	  {
	    auto x = i * Tp{0.05L};
	    auto R = dual_hahn(n, gamma, delta, N, x);
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
