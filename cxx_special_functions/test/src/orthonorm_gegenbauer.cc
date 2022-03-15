
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_gegenbauer.h>
#include <emsr/math_constants.h>

// Normalized Gegenbauer polynomial.
template<typename Tp>
  Tp
  normalized_gegenbauer(int n, Tp lambda, Tp x)
  {
    constexpr auto s_pi = emsr::pi_v<Tp>;
    auto gama = std::tgamma(lambda);
    auto gamn2a = std::tgamma(n + Tp{2} * lambda);
    auto norm = std::sqrt(s_pi * std::pow(Tp{2}, Tp{1} - Tp{2} * lambda)
	                * gamn2a / emsr::factorial<Tp>(n) / (Tp(n) + lambda))
              / gama;
    return emsr::gegenbauer(n, lambda, x) / norm;
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int n2, Tp lambda, Tp x)
  {
    return std::pow(Tp{1} - x * x, lambda - Tp{0.5})
	 * normalized_gegenbauer(n1, lambda, x)
	 * normalized_gegenbauer(n2, lambda, x);
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_gegenbauer(Tp lambda)
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    const bool singular = (lambda < Tp{0.5});

    const std::array<int, 10> degree{{0, 1, 2, 3, 5, 10, 15}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
            if (n2 > n1)
              continue; // No need to duplicate.

	    auto func = [n1, n2, lambda](Tp x)
			-> Tp
			{ return integrand<Tp>(n1, n2, lambda, x); };

	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func,
					Tp{-1}, Tp{1},
					lambda - Tp{0.5}, lambda - Tp{0.5},
					0, 0,
					abs_precision, rel_precision)
		: emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1}, abs_precision, rel_precision);

            auto del = delta<Tp>(n1, n2) - result;
	    if (std::abs(del) > cmp_precision)
              {
		++num_errors;
        	std::cout << "n1 = " << n1 << "; n2 = " << n2 << "; lambda = " << lambda
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
  num_errors += test_gegenbauer<float>(0.5F);
  num_errors += test_gegenbauer<double>(0.5);
  num_errors += test_gegenbauer<long double>(0.5L);
  return num_errors;
}
