/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <jacobi_small.hpp>
#include <emsr/float128_io.h>

  template<typename _Tp>
    std::vector<emsr::QuadraturePoint<_Tp>>
    gegenbauer_zeros(unsigned int n, _Tp alpha)
    {
      const auto _S_eps = emsr::epsilon(alpha);
      const unsigned int _S_maxit = 1000u;
      std::vector<emsr::QuadraturePoint<_Tp>> pt(n);

      _Tp z;
      _Tp w = _Tp{0};
      for (auto i = 1u; i <= n; ++i)
	{
	  if (i == 1)
	    {
	      auto an = alpha / n;
	      auto an2 = an * an;
	      auto r1 = (1.0 + alpha) * (2.78 / (4.0 + n * n)
			+ 0.768 * an / n);
	      auto r2 = 1.0 + 1.48 * an + 0.96 * an + 1.282 * an2;
	      z = 1.0 - r1 / r2;
	    }
	  else if (i == 2)
	    {
	      auto r1 = (4.1 + alpha)
			/ ((1.0 + alpha) * (1.0 + 0.156 * alpha));
	      auto r2 = 1.0
			+ 0.06 * (n - 8.0) * (1.0 + 0.12 * alpha) / n;
	      auto r3 = 1.0
		    + 0.012 * alpha * (1.0 + 0.25 * std::abs(alpha)) / n;
	      z -= (1.0 - z) * r1 * r2 * r3;
	    }
	  else if (i == 3)
	    {
	      auto r1 = (1.67 + 0.28 * alpha) / (1.0 + 0.37 * alpha);
	      auto r2 = 1.0 + 0.22 * (n - 8.0) / n;
	      auto r3 = 1.0 + 8.0 * alpha / ((6.28 + alpha) * n * n);
	      z -= (pt[0].point - z) * r1 * r2 * r3;
	    }
	  else if (i == n - 1)
	    {
	      auto r1 = (1.0 + 0.235 * alpha) / (0.766 + 0.119 * alpha);
	      auto r2 = 1.0 / (1.0 + 0.639 * (n - 4.0)
						/ (1.0 + 0.71 * (n - 4.0)));
	      auto r3 = 1.0 / (1.0 + 20.0 * alpha
				/ ((7.5 + alpha) * n * n));
	      z += (z - pt[n - 4].point) * r1 * r2 * r3;
	    }
	  else if (i == n)
	    {
	      auto r1 = (1.0 + 0.37 * alpha) / (1.67 + 0.28 * alpha);
	      auto r2 = 1.0 / (1.0 + 0.22 * (n - 8.0) / n);
	      auto r3 = 1.0 / (1.0 + 8.0 * alpha
				 / ((6.28 + alpha) * n * n));
	      z += (z - pt[n - 3].point) * r1 * r2 * r3;
	    }
	  else
	    z = 3.0 * pt[i - 2].point
		- 3.0 * pt[i - 3].point + pt[i - 4].point;

	  auto __2alpha = _Tp{2} * alpha;
	  for (auto its = 1u; its <= _S_maxit; ++its)
	    {
	      auto temp = _Tp{2} + __2alpha;
	      auto C1 = (temp * z) / _Tp{2};
	      auto C2 = _Tp{1};
	      for (auto j = 2u; j <= n; ++j)
		{
		  auto C3 = C2;
		  C2 = C1;
		  temp = _Tp{2} * j + __2alpha;
		  auto a = _Tp{2} * j * (j + __2alpha)
			   * (temp - _Tp{2});
		  auto b = (temp - _Tp{1})
			   * temp * (temp - _Tp{2}) * z;
		  auto c = _Tp{2} * (j - 1 + alpha)
			   * (j - 1 + alpha) * temp;
		  C1 = (b * C2 - c * C3) / a;
		}
	      auto Cp = (n * (-temp * z) * C1
			   + _Tp{2} * (n + alpha) * (n + alpha) * C2)
			/ (temp * (_Tp{1} - z * z));
	      auto z1 = z;
	      z = z1 - C1 / Cp;
	      if (std::abs(z - z1) <= _S_eps)
		{
		  w = std::exp(std::lgamma(alpha + _Tp(n))
			       + std::lgamma(alpha + _Tp(n))
			       - std::lgamma(_Tp(n + 1))
			       - std::lgamma(_Tp(n + 1) + __2alpha))
		      * temp * std::pow(_Tp{2}, __2alpha) / (Cp * C2);
		  break;
		}
	      if (its > _S_maxit)
		throw std::logic_error("jacobi_zeros: Too many iterations");
	    }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}

      return pt;
    }

template<typename _Tp>
  void
  test_gegenbauer(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 6;

    std::cout << "\njacobi\n";
    for (int n = 0; n <= 5; ++n)
      {
	for (int i = 0; i <= 3; ++i)
	  {
            auto alpha = i * _Tp{1};
            auto jalpha = alpha - _Tp{1}/_Tp{2};
            auto jnorm = emsr::rising_factorial(_Tp{2} * alpha, n)
			/ emsr::rising_factorial(alpha + _Tp{1}/_Tp{2}, n);
            std::cout << "n     = " << n << '\n';
            std::cout << "alpha = " << alpha << '\n';
            Life::Jacobi<_Tp> jac(n, jalpha, jalpha);
	    const auto del = _Tp{1} / _Tp{100};
	    for (int k = 0; k <= 200; ++k)
              {
        	auto x = (k - 100) * del;
        	auto ggb = emsr::gegenbauer(n, alpha, x);
        	auto gj = jnorm * emsr::jacobi(n, jalpha, jalpha, x);
        	auto lj = jnorm * jac(x);
        	std::cout << std::setw(width) << x
        	          << std::setw(width) << ggb
        	          << std::setw(width) << gj
        	          << std::setw(width) << lj
        	          << std::setw(width) << ggb - lj
        	          << '\n';
              }
            std::cout << '\n';

	    for (int n = 0; n <= 50; ++n)
	      {
		auto pt = gegenbauer_zeros(n, alpha);
		std::cout << "\nn = " << std::setw(4) << n
			  << ", alpha = " << std::setw(width) << alpha << ":\n";
		for (auto [z, w] : pt)
		  std::cout << ' ' << std::setw(width) << z
		  	    << ' ' << std::setw(width) << w
		  	    << '\n';
	      }
          }
      }
  }

int
main()
{
  test_gegenbauer(1.0F);

  test_gegenbauer(1.0);

  test_gegenbauer(1.0L);

  //test_gegenbauer(1.0Q);

  return 0;
}
