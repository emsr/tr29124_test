/**
 *
 */

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <wrap_burkhardt.h>

/**
 * Compute the Krawtchouk polynomial by recursion:
 * @f[
 *    -xK_n(x) = p(N - n) K_{n+1}(x)
 *             - [p(N - n) + n(1 - p)] K_n(x)
 *             + n(1 - p) K_{n-1}(x)
 * @f]
 * where @f$ K_n(x) = K_n(x; p, N) @f$.
 */
template<typename Tp, typename _TpX>
  Tp
  krawtchouk_recur(int n, Tp p, int N, _TpX x)
  {
    auto Knm1 = Tp{1};
    if (n == 0)
      return Knm1;

    auto Kn = Tp{1} - Tp(x) / p / Tp(N);
    if (n == 1)
      return Kn;

    const auto q = Tp{1} - p;
    auto pnn = p * Tp(N - 1);
    auto nq = q;
    auto Knp1 = ((pnn + nq - Tp(x)) * Kn - nq * Knm1) / pnn;
    for (int k = 2; k < n; ++k)
      {
	pnn = p * Tp(N - k);
	nq = Tp(k) * q;
	Knm1 = Kn;
	Kn = Knp1;
	Knp1 = ((pnn + nq - Tp(x)) * Kn - nq * Knm1) / pnn;
      }

    return Knp1;
  }

/**
 * Return the Krawtchouk polynomial defined by
 * @f[
 *    K_n(x; p, N) = {}_2F_1(-n, -x; -N; \frac{1}{p})
 * @f]
 */
template<typename Tp, typename _TpX>
  Tp
  krawtchouk(int n, Tp p, int N, _TpX x)
  {
    if (std::isnan(p))
      return p;
    else if (std::isnan(x))
      return x;
    else if (n > N)
      throw std::domain_error("krawtchouk: Degree must be less than or equal to big N.");
    else
      return krawtchouk_recur(n, p, N, x);
  }

/**
 * 
 */
template<typename Tp>
  void
  test_krawtchouk(Tp p, int N)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * Tp{0.05L};
	    auto K = krawtchouk(n, p, N, x);
	    auto K_test = burkhardt::krawtchouk(n, p, N, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << K
		      << ' ' << std::setw(w) << K_test
		      << '\n';
	  }
      }
  }

int
main()
{
  test_krawtchouk(0.5f, 10);
}
