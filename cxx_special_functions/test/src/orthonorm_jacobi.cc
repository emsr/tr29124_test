
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_jacobi.h>

// Try to manage the four-gamma ratio.
// alpha > -1, beta > -1.
template<typename Tp>
  Tp
  gamma_ratio(int n, Tp alpha, Tp beta)
  {
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    if (std::abs(Tp{1} + alpha + beta) < s_eps)
      return Tp{0};
    else
      {
	auto gaman1 = std::tgamma(Tp{1} + alpha);
	auto gambn1 = std::tgamma(Tp{1} + beta);
	auto gamabn1 = std::tgamma(Tp{1} + alpha + beta);
	auto fact = gaman1 * gambn1 / gamabn1;
	for (int k = 1; k <= n; ++k)
	  fact *= (Tp(k) + alpha) * (Tp(k) + beta)
		/ (Tp(k) + alpha + beta) / Tp(k);
	return fact;
      }
  }

// Normalized Jacobi polynomial.
template<typename Tp>
  Tp
  normalized_jacobi(int n, Tp alpha, Tp beta, Tp x)
  {
    auto gam = gamma_ratio(n, alpha, beta);
    auto norm = std::sqrt(std::pow(Tp{2}, Tp{1} + alpha + beta)
	                * gam / (Tp(2 * n + 1) + alpha + beta));
    return emsr::jacobi(n, alpha, beta, x) / norm;
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int n2, Tp alpha, Tp beta, Tp x)
  {
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    if (std::abs(x - Tp{1}) < s_eps)
      return Tp{0};
    else if (std::abs(x + Tp{1}) < s_eps)
      return Tp{0};
    else
      {
	return std::pow(Tp{1} - x, alpha) * std::pow(Tp{1} + x, beta)
	     * normalized_jacobi(n1, alpha, beta, x)
	     * normalized_jacobi(n2, alpha, beta, x);
      }
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_jacobi(Tp alpha, Tp beta)
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    const bool singular = (alpha < Tp{0} || beta < Tp{0});

    const std::array<int, 10> degree{{0, 1, 2, 3, 5, 10, 15}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
            if (n2 > n1)
              continue; // No need to duplicate.

	    auto func = [n1, n2, alpha, beta](Tp x)
			-> Tp
			{ return integrand<Tp>(n1, n2, alpha, beta, x); };

	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func, Tp{-1}, Tp{1},
					             alpha, beta, 0, 0, abs_precision, rel_precision)
		: emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1}, abs_precision, rel_precision);

            auto del = delta<Tp>(n1, n2) - result;
	    if (std::abs(del) > cmp_precision)
              {
		++num_errors;
        	std::cout << "n1 = " << n1 << "; n2 = " << n2
                          << "; alpha = " << alpha << "; beta = " << beta
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
  num_errors += test_jacobi<float>(0.5F, 1.5F);
  num_errors += test_jacobi<double>(0.5, 1.5);
  num_errors += test_jacobi<long double>(0.5L, 1.5L);
  return num_errors;
}
