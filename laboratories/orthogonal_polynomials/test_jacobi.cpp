/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <jacobi_small.hpp>

#include <emsr/float128_io.h>
#include <emsr/sf_jacobi.h>
#include <emsr/quadrature_point.h>

  template<typename Tp>
    std::vector<emsr::QuadraturePoint<Tp>>
    jacobi_zeros(unsigned int n, Tp alpha, Tp beta)
    {
      const auto s_eps = emsr::epsilon(alpha);
      const unsigned int s_maxit = 1000u;

      std::vector<emsr::QuadraturePoint<Tp>> pt(n);

      Tp z;
      Tp w = Tp{0};
      for (auto i = 1u; i <= n; ++i)
	{
	  if (i == 1)
	    {
	      auto an = alpha / n;
	      auto bn = beta / n;
	      auto r1 = (1.0 + alpha) * (2.78 / (4.0 + n * n)
			+ 0.768 * an / n);
	      auto r2 = 1.0 + 1.48 * an + 0.96 * bn
			+ 0.452 * an * an + 0.83 * an * bn;
	      z = 1.0 - r1 / r2;
	    }
	  else if (i == 2)
	    {
	      auto r1 = (4.1 + alpha)
			/ ((1.0 + alpha) * (1.0 + 0.156 * alpha));
	      auto r2 = 1.0
			+ 0.06 * (n - 8.0) * (1.0 + 0.12 * alpha) / n;
	      auto r3 = 1.0
		    + 0.012 * beta * (1.0 + 0.25 * std::abs(alpha)) / n;
	      z -= (1.0 - z) * r1 * r2 * r3;
	    }
	  else if (i == 3)
	    {
	      auto r1 = (1.67 + 0.28 * alpha) / (1.0 + 0.37 * alpha);
	      auto r2 = 1.0 + 0.22 * (n - 8.0) / n;
	      auto r3 = 1.0 + 8.0 * beta / ((6.28 + beta) * n * n);
	      z -= (pt[0].point - z) * r1 * r2 * r3;
	    }
	  else if (i == n - 1)
	    {
	      auto r1 = (1.0 + 0.235 * beta) / (0.766 + 0.119 * beta);
	      auto r2 = 1.0 / (1.0 + 0.639 * (n - 4.0)
						/ (1.0 + 0.71 * (n - 4.0)));
	      auto r3 = 1.0 / (1.0 + 20.0 * alpha
				/ ((7.5 + alpha) * n * n));
	      z += (z - pt[n - 4].point) * r1 * r2 * r3;
	    }
	  else if (i == n)
	    {
	      auto r1 = (1.0 + 0.37 * beta) / (1.67 + 0.28 * beta);
	      auto r2 = 1.0 / (1.0 + 0.22 * (n - 8.0) / n);
	      auto r3 = 1.0 / (1.0 + 8.0 * alpha
				 / ((6.28 + alpha) * n * n));
	      z += (z - pt[n - 3].point) * r1 * r2 * r3;
	    }
	  else
	    {
	      z = 3.0 * pt[i - 2].point
		  - 3.0 * pt[i - 3].point + pt[i - 4].point;
	    }

	  auto alphabeta = alpha + beta;
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto temp = Tp{2} + alphabeta;
	      auto P1 = (alpha - beta + temp * z) / Tp{2};
	      auto P2 = Tp{1};
	      for (auto j = 2u; j <= n; ++j)
		{
		  auto P3 = P2;
		  P2 = P1;
		  temp = Tp{2} * j + alphabeta;
		  auto a = Tp{2} * j * (j + alphabeta)
			   * (temp - Tp{2});
		  auto b = (temp - Tp{1})
			   * (alpha * alpha - beta * beta
				+ temp * (temp - Tp{2}) * z);
		  auto c = Tp{2} * (j - 1 + alpha)
			   * (j - 1 + beta) * temp;
		  P1 = (b * P2 - c * P3) / a;
		}
	      auto Pp = (n * (alpha - beta - temp * z) * P1
			   + Tp{2} * (n + alpha) * (n + beta) * P2)
			/ (temp * (Tp{1} - z * z));
	      auto z1 = z;
	      z = z1 - P1 / Pp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  w = std::exp(std::lgamma(alpha + Tp(n))
			       + std::lgamma(beta + Tp(n))
			       - std::lgamma(Tp(n + 1))
			       - std::lgamma(Tp(n + 1) + alphabeta))
		      * temp * std::pow(Tp{2}, alphabeta) / (Pp * P2);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("jacobi_zeros: Too many iterations");
	    }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}

      return pt;
    }

template<typename Tp>
  void
  test_jacobi(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 6;
    std::cout << "\njacobi\n";
    for (int n = 0; n <= 5; ++n)
      {
	for (int i = 0; i <= 3; ++i)
	  {
            auto alpha = i * Tp{1};
            for (int j = 0; j <= 3; ++j)
              {
        	auto beta = j * Tp{1};
        	std::cout << "n     = " << n << '\n';
        	std::cout << "alpha = " << alpha << '\n';
        	std::cout << "beta  = " << beta << '\n';
                Life::Jacobi<Tp> jac(n, alpha, beta);
		const auto del01 = Tp{1} / Tp{100};
		for (int k = 0; k <= 200; ++k)
        	  {
        	    auto x = (k - 100) * del01;
        	    std::cout << std::setw(width) << x
        	              << std::setw(width) << emsr::jacobi(n, alpha, beta, x)
        	              << std::setw(width) << jac(x)
        	              << '\n';
        	  }
        	std::cout << '\n';

		for (int n = 0; n <= 50; ++n)
		  {
		    auto pt = jacobi_zeros(n, alpha, beta);
		    std::cout << "\nn = " << std::setw(4) << n << ":\n";
		    for (auto [z, w] : pt)
		      std::cout << ' ' << std::setw(width) << z
				<< ' ' << std::setw(width) << w
				<< '\n';
		  }
	      }
            std::cout << '\n';
          }
      }
  }

int
main()
{
  test_jacobi(1.0F);

  test_jacobi(1.0);

  test_jacobi(1.0L);

  //test_jacobi(1.0Q);

  return 0;
}
