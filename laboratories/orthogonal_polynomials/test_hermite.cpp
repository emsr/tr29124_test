/**
 *
 */

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <emsr/float128_io.h>
#include <new_hermite.tcc>
#include <emsr/continued_fractions.h>
#include <emsr/integration.h>

  /**
   * Compute the Hermite polynomial ratio by continued fraction:
   * @f[
   *    \frac{H_n(x)}{H_{n-1}(x)} = 2n\frac{H_n(x)}{H'_n(x)}
   *       = b_0 + K_{k=1}^{n-1}\left(\frac{-2(n-k)}{2x}\right)
   *      b_k = 2x, \mbox{ } k >= 0 \mbox{ } a_k = -2(n-k) \mbox{ } k >= 1
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
    hermite_ratio(unsigned int n, Tp x)
    {
      auto a = [n](std::size_t k, Tp) { return Tp(-2 * (n - k)); };
      using AFun = decltype(a);

      auto b = [x](std::size_t, Tp) { return Tp{2} * x; };
      using BFun = decltype(b);

      auto w = [](std::size_t, Tp) { return Tp{0}; };
      using TailFun = decltype(w);

      using CFrac = emsr::LentzContinuedFraction<Tp, AFun, BFun, TailFun>;
      CFrac Hrat(a, b, w);

      return Hrat();
    }

  /**
   * Build a vector of the Gauss-Hermite integration rule abscissae and weights.
   */
  template<typename Tp>
    std::vector<emsr::QuadraturePoint<Tp>>
    hermite_zeros(unsigned int n, Tp proto = Tp{})
    {
      const auto s_eps = emsr::epsilon(proto);
      const unsigned int s_maxit = 1000u;
      const auto s_pim4 = Tp{0.7511255444649424828587030047762276930510L};
      const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;

      std::vector<emsr::QuadraturePoint<Tp>> pt(n);

      const auto m = n / 2;

      // Treat the central zero for odd order specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large order.
      if (n & 1)
	{
	  if (n < emsr::detail::s_num_factorials<Tp>)
	    {
	      auto nm = n - 1;
	      auto nmfact = emsr::detail::factorial<Tp>(nm);
	      auto mm = nm / 2;
	      auto mmfact = emsr::detail::factorial<Tp>(mm);
	      auto Hnm1 = (mm & 1 ? Tp{-1} : Tp{1}) / mmfact;
	      pt[m].point = Tp{0};
	      pt[m].weight = s_sqrt_pi * std::pow(Tp{2}, Tp(n - 1))
				 / nmfact / Hnm1 / Hnm1 / n;
	    }
	  else
	    {
	      auto nm = n - 1;
	      auto nmfact = emsr::detail::log_factorial<Tp>(nm);
	      auto mm = nm / 2;
	      auto mmfact = emsr::detail::log_factorial<Tp>(mm);
	      pt[m].point = Tp{0};
	      pt[m].weight = s_sqrt_pi * std::pow(Tp{2}, Tp(n - 1))
				 *std::exp(-(nmfact - 2 * mmfact)) / n;
	    }
	}

      Tp z;
      Tp w = Tp{0};
      for (auto i = 1u; i <= m; ++i)
	{
	  if (i == 1)
	    z = std::sqrt(Tp(2 * n + 1))
		- 1.85575 * std::pow(Tp(2 * n + 1), -0.166667);
	  else if (i == 2)
	    z -= 1.14 * std::pow(Tp(n), 0.426) / z;
	  else if (i == 3)
	    z = 1.86 * z - 0.86 * pt[0].point;
	  else if (i == 4)
	    z = 1.91 * z - 0.91 * pt[1].point;
	  else
	    z = 2.0 * z - pt[i - 3].point;
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto H = s_pim4;
	      auto H1 = Tp{0};
	      for (auto k = 1u; k <= n; ++k)
		{
		  auto H2 = H1;
		  H1 = H;
		  H = z * std::sqrt(Tp{2} / k) * H1
		       - std::sqrt(Tp(k - 1) / Tp(k)) * H2;
		}
	      auto Hp = std::sqrt(Tp(2 * n)) * H1;
	      auto z1 = z;
	      z = z1 - H / Hp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  w = 2.0 / (Hp * Hp);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("hermite_zeros: Too many iterations");
	    }
	  pt[n - i].point = -z;
	  pt[n - i].weight = w;
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}

      return pt;
    }

template<typename Tp>
  void
  test_hermite(Tp proto = Tp{})
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
    const auto s_Ai0 = Tp{-2.3381074104597670384891972524467L};

    auto fname = [](std::string_view front, int n, std::string_view back)
		 {
		   std::ostringstream out;
		   out << front << n << back;
		   return out.str();
		 };

    std::cout.precision(emsr::digits10(proto));
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    std::cout << "\f\n\n  Table of integer sqrt\n";
    std::cout << "  =====================\n";
    for (int n = 0; n <= 100; ++n)
      std::cout << "  " << std::setw(4) << n
		<< "  " << std::setw(width) << std::sqrt(Tp(n)) << '\n';

    std::cout << "\f\n\n  Table of factorial sqrt\n";
    std::cout << "  =======================\n";
    auto fact = Tp{1};
    for (int n = 0; n <= 100; ++n)
      {
	std::cout << "  " << std::setw(4) << n
		  << "  " << std::setw(width) << fact << '\n';
	fact *= std::sqrt(Tp(n + 1));
      }

    std::cout << "\f\n\n  Examine asymptotic transition region\n";
    std::cout << "  ====================================\n";
    for (int n = 4; n <= 200; ++n)
      {
	auto xt = std::sqrt(Tp(2 * n + 1));
        auto delta = s_Ai0 / s_sqrt_2 / std::pow(n, Tp{1} / Tp{6});
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << "  " << std::setw(width) << fname("Ha2_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta2";
	std::cout << "  " << std::setw(width) << "delta2-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	const auto del = xt / Tp{2} / Tp{200};
	for (int i = 0; i <= 201; ++i)
          {
            auto x = Tp{3} * xt / Tp{4} + i * del;
            auto h = poly_hermite_recursion(n, x);
            // This sorta works... I think the old asymp is for He_n(x).
            //auto ht = std::exp2(n - 1) * poly_hermite_asymp_old(n, x / s_sqrt_2);
            auto ht = poly_hermite_asymp_old(n, x);
            auto h2 = poly_hermite_asymp(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - xt - delta) < del)
	      std::cout << "> ";
	    else if (std::abs(x - xt + delta) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
	    std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << "  " << std::setw(width) << h2
		      << "  " << std::setw(width) << h2 - h
		      << "  " << std::setw(width) << (h2 - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Examine asymptotic transition region\n";
    std::cout << "  ====================================\n";
    for (int n : {1000, 2000, 5000, 10000})
      {
	auto xt = std::sqrt(Tp(2 * n + 1));
        auto delta = s_Ai0 / s_sqrt_2 / std::pow(n, Tp{1} / Tp{6});
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << "  " << std::setw(width) << fname("Ha2_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta2";
	std::cout << "  " << std::setw(width) << "delta2-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	const auto del = xt / Tp{2} / Tp{200};
	for (int i = 1; i <= 201; ++i)
          {
            auto x = Tp{3} * xt / Tp{4} + i * del;
            auto h = poly_hermite_recursion(n, x);
            auto ht = poly_hermite_asymp_old(n, x);
            auto h2 = poly_hermite_asymp(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - xt - delta) < del)
	      std::cout << "> ";
	    else if (std::abs(x - xt + delta) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
	    std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << "  " << std::setw(width) << h2
		      << "  " << std::setw(width) << h2 - h
		      << "  " << std::setw(width) << (h2 - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Compare recursion and asymptotic\n";
    std::cout << "  ================================\n";
    for (int n = 0; n <= 50; ++n)
      {
	auto xt = std::sqrt(Tp(2 * n + 1));
        auto delta = s_Ai0 / s_sqrt_2 / std::pow(n, Tp{1} / Tp{6});
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';;
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';;
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	const auto del = Tp{1} / Tp{10};
	for (int i = 0; i <= 100; ++i)
          {
            auto x = i * del;
            auto h = poly_hermite_recursion(n, x);
            auto ht = poly_hermite_asymp(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - xt - delta) < del)
	      std::cout << "> ";
	    else if (std::abs(x - xt + delta) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
            std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Compare normalizations\n";
    std::cout << "  ======================\n";
    auto power = Tp{1};
    for (int n = 0; n <= 50; ++n)
      {
	// sqrt(factorial(n) * 2**n sqrt(pi))
	auto norm = std::exp(std::lgamma(n + 1) / Tp{2}) * std::sqrt(power * std::sqrt(s_pi));
	power *= Tp{2};
	std::cout << "  " << std::setw(width) << "n = " << n << '\n';
	std::cout << "  " << std::setw(width) << "norm = " << norm << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("H_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Hn_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("He_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Hen_", n, "(x)");
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	const auto del = Tp{1} / Tp{10};
	for (int i = 0; i <= 100; ++i)
          {
            auto x = i * del;
            auto h = poly_hermite_recursion(n, x);
            auto hn = poly_hermite_norm_recursion(n, x);
            auto he = poly_prob_hermite_recursion(n, x);
            auto hen = poly_prob_hermite_norm_recursion(n, x);
            std::cout << "  " << std::setw(width) << x;
            std::cout << "  " << std::setw(width) << h;
            std::cout << "  " << std::setw(width) << hn;
            std::cout << "  " << std::setw(width) << he;
            std::cout << "  " << std::setw(width) << hen;
            std::cout << '\n';
          }
      }

    for (int n = 0; n <= 50; ++n)
      {
	auto pt = hermite_zeros(n, proto);
	std::cout << "\nn = " << std::setw(4) << n << ":\n";
	for (auto [z, w] : pt)
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << w
		    << '\n';
      }

    return;
  }

template<typename Tp>
  void
  test_hermite_norm(Tp proto)
  {
    return;
  }


int
main()
{
  test_hermite(1.0F);

  test_hermite(1.0);

  test_hermite(1.0L);

  //test_hermite(1.0Q);

  return 0;
}
