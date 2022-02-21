
#include <limits>
#include <vector>
#include <iostream>
#include <iomanip>

#include <emsr/matrix.h>
#include <emsr/quadrature_point.h>

  /**
   * @brief  Return the Legendre polynomial by upward recursion
   * 	     on degree @f$ l @f$.
   *
   * The Legendre function of degree @f$ l @f$ and argument @f$ x @f$,
   * @f$ P_l(x) @f$, is defined by:
   * @f[
   *   P_l(x) = \frac{1}{2^l l!}\frac{d^l}{dx^l}(x^2 - 1)^{l}
   * @f]
   * This can be expressed as a series:
   * @f[
   *   P_l(x) = \frac{1}{2^l l!}\sum_{k=0}^{\lfloor l/2 \rfloor}
   *            \frac{(-1)^k(2l-2k)!}{k!(l-k)!(l-2k)!}x^{l-2k}
   * @f]
   *
   * @param  l  The degree of the Legendre polynomial.  @f$ l >= 0 @f$.
   * @param  x  The argument of the Legendre polynomial.
   */
  template<typename Tp>
    Tp
    legendre_p_norm(unsigned int l, Tp x)
    {
      const auto norm = std::sqrt(Tp(2 * l + 1) / Tp{2});

      if (x == Tp{+1})
	return norm * Tp{+1};
      else if (x == Tp{-1})
	return norm * (l % 2 == 1 ? Tp{-1} : Tp{+1});
      else
	{
	  auto P_lm2 = Tp{1};
	  if (l == 0)
	    return norm * P_lm2;

	  auto P_lm1 = x;
	  if (l == 1)
	    return norm * P_lm1;

	  auto P_l = Tp{2} * x * P_lm1 - P_lm2
		    - (x * P_lm1 - P_lm2) / Tp{2};
	  for (unsigned int ll = 3; ll <= l; ++ll)
	    {
	      P_lm2 = P_lm1;
	      P_lm1 = P_l;
	      // This arrangement is supposed to be better for roundoff
	      // protection, Arfken, 2nd Ed, Eq 12.17a.
	      P_l = Tp{2} * x * P_lm1 - P_lm2
		    - (x * P_lm1 - P_lm2) / Tp(ll);
	    }

	  return norm * P_l;
	}
    }

  /**
   * Return a Chebyshev polynomial of non-negative order @f$ n @f$
   * and real argument @f$ x @f$ by the recursion
   * @f[
   *    C_n(x) = 2xC_{n-1} - C_{n-2}
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   * @param C0 The value of the zeroth-order Chebyshev polynomial at @f$ x @f$
   * @param C1 The value of the first-order Chebyshev polynomial at @f$ x @f$
   */
  template<typename Tp>
    std::tuple<Tp, Tp, Tp>
    chebyshev_recur(unsigned int n, Tp x, Tp C0, Tp C1)
    {
      auto Ck = Tp{2} * x * C1 - C0;
      for (unsigned int j = 2; j < n; ++j)
      {
	C0 = C1;
	C1 = Ck;
	Ck = Tp{2} * x * C1 - C0;
      }
      return std::make_tuple(Ck, C1, C0);
    }

  /**
   * Return the Chebyshev polynomial of the first kind @f$ T_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the first kind is defined by:
   * @f[
   *    T_n(x) = \cos(n \theta)
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    Tp
    chebyshev_t(unsigned int n, Tp x)
    {
      auto T0 = Tp{1};
      if (n == 0)
	return T0;

      auto T1 = x;
      if (n == 1)
	return T1;

      auto Ts = chebyshev_recur(n, x, T0, T1);
      return std::get<0>(Ts);
    }

  /**
   * Return the Chebyshev polynomial of the second kind @f$ U_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the second kind is defined by:
   * @f[
   *    U_n(x) = \frac{\sin \left[(n + 1)\theta \right]}{\sin(\theta)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    Tp
    chebyshev_u(unsigned int n, Tp x)
    {
      auto U0 = Tp{1};
      if (n == 0)
	return U0;

      auto U1 = Tp{2} * x;
      if (n == 1)
	return U1;

      auto Us = chebyshev_recur(n, x, U0, U1);
      return std::get<0>(Us);
    }

  /**
   * *** stole this from tr29124/quadrature/build_*
   *
   * Return a vector of zeros of the Chebyshev function of the second kind
   * of order @f$ n @f$, @f$ U_n(x) @f$.
   * The zeros are given by:
   * @f[
   *   x_k = \cos\left(\frac{k\pi}{n + 1}\right), k \elem {1, ..., n}
   * @f]
   */
  template<typename Tp>
    std::vector<emsr::QuadraturePoint<Tp>>
    chebyshev_u_zeros(unsigned int n)
    {
      const auto S_pi = Tp{3.141592653589793238462643383279502884195Q};
      std::vector<emsr::QuadraturePoint<Tp>> pt(n);
      for (unsigned int k = 1; k <= n; ++k)
	{
	  auto arg = Tp(k) / Tp(n + 1);
	  auto half = (arg == Tp{0.5Q});
	  auto z = (half ? Tp{0} : std::cos(S_pi * arg));
	  auto w = S_pi * (Tp{1} - z) * (Tp{1} + z) / Tp(n + 1);
	  pt[k - 1].point = z;
	  pt[k - 1].weight = w;
	}
      return pt;
    }

/**
 * *** stole this from tr29124/quadrature/build_*
 *
 * @f[
 *    w_k = frac{c_k}{n}\left[1-\sum_{j=1}^{[n/2]}\frac{b_k}{4j^2-1}
 *            \cos\left(\frac{2jk\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 * 
 * @see Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules
 */
template<typename Tp>
  std::vector<emsr::QuadraturePoint<Tp>>
  build_clenshaw_curtis_sum(std::size_t n)
  {
    std::vector<emsr::QuadraturePoint<Tp>> out(n + 1);
    if (n == 1)
      {
	out[0].point = Tp{0};
	out[0].weight = Tp{2};
	return out;
      }
    else
      {
	const auto S_pi = Tp{3.141592653589793238462643383279502884195Q};
	auto uz = chebyshev_u_zeros<Tp>(n - 1);
	out[0].point = Tp{+1};
	out[0].weight = Tp{1} / (n * n - 1 + n % 2);
	for (auto k = 1u; k <= uz.size(); ++k)
	  {
	    out[k].point = uz[k - 1].point;

	    auto sum = Tp{0};
	    for (auto j = 1u; j <= n / 2; ++j)
	      {
		auto b = Tp(j == n / 2 ? 1 : 2);
		sum += b * std::cos(2 * j * k * S_pi / n)
		       / Tp(4 * j * j - 1);
	      }
	    auto w = Tp{2} * (Tp{1} - sum) / Tp(n);
	    out[k].weight = w;
	  }
	out[n].point = Tp{-1};
	out[n].weight = out[0].weight;
	return out;
      }
  }

template<typename Tp>
  void
  test_vandermonde()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nCQUAD Rules\n";
    for (const auto& n : {4, 8, 16, 32, 64})
      {
	std::cout << "\nClenshaw-Curtis " << n << "\n";
	auto ccvec = build_clenshaw_curtis_sum<Tp>(n);
	std::reverse(ccvec.begin(), ccvec.end());
	for (const auto& cc : ccvec)
	  {
	    std::cout << std::setw(w) << cc.point << ' '
		      << std::setw(w) << cc.weight << ' '
		      << '\n';
	  }

	// The Vandermonde-like matrices use polynomial values not monomial powers.
        std::vector<std::vector<Tp>> vdm(n + 1, std::vector<Tp>(n + 1));
	for (int i = 0; i <= n; ++i)
	  for (int j = 0; j <= n; ++j)
	    vdm[i][j] = legendre_p_norm(j, ccvec[i].point);
	using matrix_t = std::vector<std::vector<Tp>>;
	matrix_t inv(n + 1, std::vector<Tp>(n + 1));
	std::vector<int> index(n + 1);
	Tp parity;
	emsr::lu_decomp(n + 1, vdm, index, parity);
	emsr::lu_invert(n + 1, vdm, index, inv);
	std::cout << '\n';
	for (int i = 0; i <= n; ++i)
	  {
	    for (int j = 0; j <= n; ++j)
	      {
		if (std::abs(inv[i][j]) < 50 * std::numeric_limits<Tp>::epsilon())
		  inv[i][j] = Tp{};
		std::cout << ' ' << std::setw(w) << inv[i][j];
	      }
	    std::cout << '\n';
	  }
	std::cout << '\n';
      }
  }

int
main()
{
  test_vandermonde<long double>();
}

