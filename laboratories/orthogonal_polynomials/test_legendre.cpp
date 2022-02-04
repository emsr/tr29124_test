/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>
#include <emsr/quadrature_point.h>
#include <emsr/sf_legendre.h>
#include <emsr/sf_gamma.h> / factorial

  /**
   * Build a list of zeros and weights for the Gauss-Legendre integration rule
   * for the Legendre polynomial of degree @c l.
   */
  template<typename Tp>
    std::vector<emsr::QuadraturePoint<Tp>>
    legendre_zeros(unsigned int l, Tp proto = Tp{})
    {
      const auto _S_eps = emsr::epsilon(proto);
      const auto _S_pi = emsr::pi_v<Tp>;
      const unsigned int _S_maxit = 1000u;

      std::vector<emsr::QuadraturePoint<Tp>> pt(l);

      auto m = l / 2;

      // Treat the central zero for odd order specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large order.
      if (l & 1)
	{
	  if (l < emsr::detail::s_num_factorials<Tp>)
	    {
	      const auto lm = l - 1;
	      const auto lmfact = emsr::detail::factorial<Tp>(lm);
	      const auto mm = lm / 2;
	      const auto mmfact = emsr::detail::factorial<Tp>(mm);
	      auto Plm1 = ((lm & 1) ? -1 : 1)
			    * lmfact / mmfact / mmfact
			    / std::exp2(Tp(lm));
	      auto Ppl = l * Plm1;
	      pt[m].point = Tp{0};
	      pt[m].weight = Tp{2} / Ppl / Ppl;
	    }
	  else
	    {
	      const auto _S_ln2 = emsr::ln2_v<Tp>;
	      const auto lm = l - 1;
	      const auto lmfact = emsr::detail::log_factorial<Tp>(lm);
	      const auto mm = lm / 2;
	      const auto mmfact = emsr::detail::log_factorial<Tp>(mm);
	      auto Plm1 = (lm & 1 ? -1 : 1)
			  * std::exp(lmfact - 2 * mmfact - lm * _S_ln2);
	      auto Ppl = l * Plm1;
	      pt[m].point = Tp{0};
	      pt[m].weight = Tp{2} / Ppl / Ppl;
	    }
	}

      for (auto i = 1u; i <= m; ++i)
	{
	  // Clever approximation of root.
	  auto z = std::cos(_S_pi * (i - Tp{1} / Tp{4})
				    / (l + Tp{1} / Tp{2}));
	  auto z1 = z;
	  auto w = Tp{0};
	  for (auto its = 0u; its < _S_maxit; ++its)
	    {
	      // Compute P, P1, and P2 the Legendre polynomials of order
	      // l, l-1, l-2 respectively by iterating through the recursion
	      // relation for the Legendre polynomials.
	      // Compute Pp the derivative of the Legendre polynomial of order l.
	      auto P1 = Tp{0};
	      auto P = Tp{1};
	      for  (auto k = 1u; k <= l; ++k)
		{
		  auto P2 = P1;
		  P1 = P;
		  // Recursion for Legendre polynomials.
		  P = ((Tp{2} * k - Tp{1}) * z * P1
		      - (k - Tp{1}) * P2) / k;
		}
	      // Recursion for the derivative of The Legendre polynomial.
	      auto Pp = l * (z * P - P1) / (z * z - Tp{1});
	      z1 = z;
	      // Converge on root by Newton's method.
	      z = z1 - P / Pp;
	      if (std::abs(z - z1) < _S_eps)
		{
		  w = Tp{2} / ((Tp{1} - z * z) * Pp * Pp);
		  break;
		}
	      if (its > _S_maxit)
		throw std::logic_error("legendre_zeros: Too many iterations");
	    }

	  pt[i - 1].point = -z;
	  pt[l - i].point = z;
	  pt[i - 1].weight = w;
	  pt[l - i].weight = w;
	}

      return pt;
    }

/**
 * 
 */
template<typename Tp>
  void
  test_legendre(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    for (int l = 0; l <= 50; ++l)
      {
	std::cout << "\n\n";
	std::cout << ' ' << std::setw(w) << "x";
	std::cout << ' ' << std::setw(w) << "P_" << l << "(x)";
	std::cout << '\n';
	const auto del = Tp{1} / Tp{100};
	for (int i = -120; i <= 120; ++i)
	  {
	    auto x = i * del;
	    const auto P_l = emsr::detail::legendre_p(l, x);
	    const auto P_l0 = emsr::detail::assoc_legendre_p(l, 0, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << P_l.P_l
		      << ' ' << std::setw(w) << P_l0.P_lm
		      << ' ' << std::setw(w) << P_l.P_l - P_l0.P_lm
		      << ' ' << std::setw(w) << P_l.deriv()
		      << ' ' << std::setw(w) << P_l0.deriv()
		      << '\n';
	  }
      }

    for (unsigned int l = 0; l <= 3; ++l)
      {
	for (unsigned int m = l; m <= l; --m)
	  {
	    std::cout << "\n\n";
	    std::cout << ' ' << std::setw(w) << "x";
	    std::cout << ' ' << std::setw(w) << "P_" << l << "_" << m << "(x)";
	    std::cout << '\n';
	    const auto del = Tp{1} / Tp{100};
	    for (int i = -120; i <= 120; ++i)
	      {
		const auto x = i * del;
		const auto P = emsr::detail::assoc_legendre_p(l, m, x);
		std::cout << ' ' << std::setw(w) << x
			  << ' ' << std::setw(w) << P.P_lm
			  << ' ' << std::setw(w) << P.deriv()
			  << '\n';
	      }
	  }
      }

    for (int l = 0; l <= 1024; ++l)
      {
	auto pt = legendre_zeros(l, proto);
	std::cout << "\nl = " << std::setw(4) << l << ":\n";
	for (auto [z, w] : pt)
	  std::cout << ' ' << std::setw(w) << z
		    << ' ' << std::setw(w) << w
		    << '\n';
      }
  }

int
main()
{
  test_legendre(1.0F);

  test_legendre(1.0);

  test_legendre(1.0L);

  //test_legendre(1.0Q);

  return 0;
}

