/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

#include <legendre.tcc>
#include <emsr/sf_legendre.h>

  template<typename Tp>
    std::vector<emsr::QuadraturePoint<Tp>>
    lobatto_zeros(unsigned int l, Tp proto = Tp{})
    {
      const auto s_eps = emsr::epsilon(proto);
      const auto s_pi = emsr::pi_v<Tp>;
      const unsigned int s_maxit = 1000u;

      std::vector<emsr::QuadraturePoint<Tp>> pt(l);

      auto m = l / 2;

      // Treat the central zero for odd order specially.
      if (l & 1)
	{
	  auto lm = l - 1;
	  auto lmfact = emsr::detail::factorial<Tp>(lm);
	  auto mm = lm / 2;
	  auto mmfact = emsr::detail::factorial<Tp>(mm);
	  auto Plm1 = ((lm & 1) ? -1 : 1) * lmfact / mmfact / mmfact
			/ std::pow(Tp{2}, lm);
	  auto Ppl = l * Plm1;
	  pt[m].point = Tp{0};
	  pt[m].weight = Tp{2} / Ppl / Ppl;
	}

      for (auto i = 1u; i <= m; ++i)
	{
	  // Clever approximation of root.
	  auto z = std::cos(s_pi * (i - Tp{1} / Tp{4})
				    / (l + Tp{1} / Tp{2}));
	  auto z1 = z;
	  auto w = Tp{0};
	  for (auto its = 0u; its < s_maxit; ++its)
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
	      if (std::abs(z - z1) < s_eps)
		{
		  w = Tp{2} / ((Tp{1} - z * z) * Pp * Pp);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("lobatto_zeros: Too many iterations");
	    }

	  pt[i - 1].point = -z;
	  pt[l - i].point = z;
	  pt[i - 1].weight = w;
	  pt[l - i].weight = w;
	}

      return pt;
    }

template<typename Tp>
  void
  test_lobatto(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    for (int l = 0; l <= 50; ++l)
      {
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << "P_" << l << "(x)";
	std::cout << '\n';
	for (int i = -100; i <= 100; ++i)
	  {
	    auto x = i * Tp{0.01Q};
	    std::pair<Tp, Tp> Lo;
	    try
	      {
		Lo = lobatto(l, x);
	      }
	    catch (std::exception & err)
	      {
		std::cerr << err.what() << '\n';
	      }
	    std::cout << "  " << std::setw(width) << x;
	    std::cout << "  " << std::setw(width) << Lo.first;
	    std::cout << "  " << std::setw(width) << Lo.second;
	    std::cout << '\n';
	  }
      }
/*
    for (int l = 0; l <= 50; ++l)
      {
	auto pt = lobatto_zeros(l, proto);
	std::cout << "\nl = " << std::setw(4) << l << ":\n";
	for (auto [z, w] : pt)
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << w
		    << '\n';
      }
*/
  }

int
main()
{
  test_lobatto(1.0F);

  test_lobatto(1.0);

  test_lobatto(1.0L);

  //test_lobatto(1.0Q);

  return 0;
}

