
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_hermite.h>
#include <emsr/math_constants.h>

// Normalized Hermite polynomial.
template<typename Tp>
  Tp
  normalized_hermite(int n, Tp x)
  {
    constexpr auto s_lnpi = emsr::lnpi_v<Tp>;
    constexpr auto s_ln2 = emsr::ln2_v<Tp>;
    const auto lnorm = Tp{0.5} * (Tp{0.5} * s_lnpi + Tp(n) * s_ln2 + emsr::lfactorial<Tp>(n));
    return emsr::hermite(n, x) * std::exp(-lnorm);
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int n2, Tp x)
  {
    return std::exp(-x * x)
	 * normalized_hermite(n1, x)
	 * normalized_hermite(n2, x);
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_hermite()
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
            if (n2 > n1)
              continue; // No need to duplicate.

	    auto func = [n1, n2](Tp x)
			-> Tp
			{ return integrand(n1, n2, x); };

	    auto [result, error]
		  = emsr::integrate_sinh_sinh(func, abs_precision, rel_precision);

            auto del = delta<Tp>(n1, n2) - result;
	    if (std::abs(del) > cmp_precision)
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
  num_errors += test_hermite<float>();
  num_errors += test_hermite<double>();
  num_errors += test_hermite<long double>();
  return num_errors;
}
