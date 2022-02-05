/**
 *
 */

#include <iostream>
#include <iomanip>
#include <tr1/cmath>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

#include <emsr/quadrature_point.h>
#include <emsr/numeric_limits.h>

  /**
   * Compute the Laguerre polynomial ratio by continued fraction:
   * @f[
   *    \frac{L_n^{(\alpha)}(x)}{L'_n^{(\alpha)(x)} = \frac{x}{n-}
   *       \frac{(n+\alpha)n}{2n+\alpha-1-x-}
   *       \frac{(n-1+\alpha)(n-1)}{2n+\alpha-3-x-}
   *       \frac{(n-2+\alpha)(n-2)}{2n+\alpha-5-x-} ...
   *       \frac{1+\alpha}{1+\alpha-x}
   * @f]
   *
   * @see RICHARD J. MATHAR, GAUSS-LAGUERRE AND GAUSS-HERMITE QUADRATURE
   * ON 64, 96 AND 128 NODES
   *
   * @see T. S. Shao, T. C. Chen, and R. M. Frank,
   * Tables of zeros and Gaussian weights of certain associated Laguerre
   * polynomials and the related generalized Hermite polynomials,
   * Math. Comp. 18 (1964), no. 88, 598{616. MR 0166397 (29 #3674)
   */
  template<typename Tp>
    Tp
    laguerre_ratio(unsigned int n, Tp x)
    {
      auto frac = Tp{0};
      for (auto i = n - 1; i >= 0; --i)
	frac = Tp(n - i) * Tp(n - i)
		/ (Tp(2 * (n - i) - 1) - x - frac);
      frac = x / (n - frac);
      return frac;
    }

  /**
   * Compute the Laguerre polynomial ratio by continued fraction:
   * @f[
   *    \frac{L_n^{(\alpha)}(x)}{L'_n^{(\alpha)(x)} = \frac{x}{n-}
   *       \frac{(n+\alpha)n}{2n+\alpha-1-x-}
   *       \frac{(n-1+\alpha)(n-1)}{2n+\alpha-3-x-}
   *       \frac{(n-2+\alpha)(n-2)}{2n+\alpha-5-x-} \cdots
   *       \frac{1+\alpha}{1+\alpha-x}
   * @f]
   *
   * @see RICHARD J. MATHAR, GAUSS-LAGUERRE AND GAUSS-HERMITE QUADRATURE
   * ON 64, 96 AND 128 NODES
   *
   * @see T. S. Shao, T. C. Chen, and R. M. Frank,
   * Tables of zeros and Gaussian weights of certain associated Laguerre
   * polynomials and the related generalized Hermite polynomials,
   * Math. Comp. 18 (1964), no. 88, 598{616. MR 0166397 (29 #3674)
   */
  template<typename Tp, typename Ta>
    Tp
    laguerre_ratio(unsigned int n, Ta alpha, Tp x)
    {
      auto frac = Tp{0};
      for (auto i = n - 1; i >= 0; --i)
	frac = Tp(n - i + alpha) * Tp(n - i)
		/ (Tp(2 * (n - i) - 1 + alpha) - x - frac);
      frac = x / (n - frac);
      return frac;
    }

  /**
   * Return a vector of zeros of the Laguerre polynomial of degree n.
   */
  template<typename Tp>
    std::vector<emsr::QuadraturePoint<Tp>>
    laguerre_zeros(unsigned int n, Tp proto)
    {
      const auto s_eps = emsr::epsilon(proto);
      const unsigned int s_maxit = 1000u;

      std::vector<emsr::QuadraturePoint<Tp>> pt(n);

      auto z = Tp{0};
      auto w = Tp{0};
      for (auto i = 1u; i <= n; ++i)
	{
	  // Clever approximations for roots.
	  if (i == 1)
	    z += (3.0) / (1.0 + 2.4 * n);
	  else if (i == 2)
	    z += (15.0) / (1.0 + 2.5 * n);
	  else
	    {
	      auto ai = i - 2;
	      z += ((1.0 + 2.55 * ai) / (1.9 * ai))
		   * (z - pt[i - 3].point);
	    }
	  // Iterate TTRR for polynomial values
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto L2 = Tp{0};
	      auto L1 = Tp{1};
	      for (auto j = 1u; j <= n; ++j)
		{
		  auto L3 = L2;
		  L2 = L1;
		  L1 = ((Tp(2 * j - 1) - z) * L2
		       - (Tp(j - 1)) * L3) / Tp(j);
		}
	      // Derivative.
	      auto Lp = (Tp(n) * L1 - Tp(n) * L2) / z;
	      // Newton's rule for root.
	      auto z1 = z;
	      z = z1 - L1 / Lp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  w = Tp{-1} / (Lp * n * L2);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("laguerre_zeros: Too many iterations");
	   }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}
      return pt;
    }

  /**
   * Return an array of abscissae and weights for the Gauss-Laguerre rule.
   */
  template<typename Tp, typename Tn>
    std::vector<emsr::QuadraturePoint<Tp>>
    laguerre_zeros(unsigned int n, Tn alpha, Tp proto)
    {
      const auto s_eps = emsr::epsilon(proto);
      const unsigned int s_maxit = 1000;

      std::vector<emsr::QuadraturePoint<Tp>> pt(n);

      for (auto i = 1u; i <= n; ++i)
	{
	  auto z = Tp{0};
	  auto w = Tp{0};
	  // Clever approximations for roots.
	  if (i == 1)
	    z += (1.0 + alpha)
		 * (3.0 + 0.92 * alpha) / (1.0 + 2.4 * n + 1.8 * alpha);
	  else if (i == 2)
	    z += (15.0 + 6.25 * alpha) / (1.0 + 2.5 * n + 0.9 * alpha);
	  else
	    {
	      auto ai = i - 2;
	      z += ((1.0 + 2.55 * ai) / (1.9 * ai)
		     + 1.26 * ai * alpha / (1.0 + 3.5 * ai))
		   * (z - pt[i - 3].point) / (1.0 + 0.3 * alpha);
	    }
	  // Iterate TTRR for polynomial values
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto L2 = Tp{0};
	      auto L1 = Tp{1};
	      for (auto j = 1u; j <= n; ++j)
		{
		  auto L3 = L2;
		  L2 = L1;
		  L1 = ((Tp(2 * j - 1 + alpha) - z) * L2
			- (Tp(j - 1 + alpha)) * L3) / Tp(j);
		}
	      // Derivative.
	      auto Lp = (Tp(n) * L1 - Tp(n + alpha) * L2) / z;
	      // Newton's rule for root.
	      auto z1 = z;
	      z = z1 - L1 / Lp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  auto exparg = std::lgamma(Tp(alpha + n))
				- std::lgamma(Tp(n));
		  w = -std::exp(exparg) / (Lp * n * L2);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("laguerre_zeros: Too many iterations");
	   }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}
    return pt;
  }

template<typename Tp>
  void
  test_laguerre(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    for (int n = 0; n <= 50; ++n)
      {
	auto pt = laguerre_zeros(n, proto);
	std::cout << "\nl = " << std::setw(4) << n << ":\n";
	for (auto [z, w] : pt)
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << w
		    << '\n';
      }
  }

int
main()
{
  test_laguerre(1.0F);

  test_laguerre(1.0);

  test_laguerre(1.0L);

  //test_laguerre(1.0Q);

  return 0;
}

