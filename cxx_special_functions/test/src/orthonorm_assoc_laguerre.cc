
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_laguerre.h>

/* Orthonormality only works for integer alpha. */

// Try to manage the gamma ratio.
template<typename Tp>
  Tp
  gamma_ratio(int n, Tp alpha)
  {
    auto gaman1 = std::tgamma(Tp(1) + alpha);
    auto fact = gaman1;
    for (int k = 1; k <= n; ++k)
      fact *= (Tp(k) + alpha) / Tp(k);
    return fact;
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  norm_assoc_laguerre(int n1, int n2, Tp alpha, Tp x)
  {
    auto norm = gamma_ratio(n1, alpha);
    return std::pow(x, alpha) * std::exp(-x)
	 * emsr::assoc_laguerre(n1, alpha, x)
	 * emsr::assoc_laguerre(n2, alpha, x) / norm;
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_assoc_laguerre(Tp alpha)
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
	    auto func = [n1, n2, alpha](Tp x)
			-> Tp
			{ return norm_assoc_laguerre<Tp>(n1, n2, alpha, x); };

	    auto [result, error]
		= emsr::integrate_exp_sinh(func, Tp{0}, abs_prec, rel_prec);

            auto del = delta<Tp>(n1, n2) - result;
	    if (std::abs(del) > cmp_prec)
              {
		++num_errors;
        	std::cout << "n1 = " << n1 << "; n2 = " << n2
                          << "; res = " << result << "; del = " << del << '\n';
              }
	  }
      }
    return num_errors;
  }

int
main()
{
  int num_errors = 0;
  num_errors += test_assoc_laguerre<float>(0.0F);
  num_errors += test_assoc_laguerre<double>(0.0);
  num_errors += test_assoc_laguerre<long double>(0.0L);
  return num_errors;
}
