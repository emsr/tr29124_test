
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_laguerre.h>

// Laguerre polynomials are already normalized.

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int n2, Tp x)
  {
    return std::exp(-x)
	 * emsr::laguerre(n1, x)
	 * emsr::laguerre(n2, x);
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_laguerre()
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    const std::array<int, 10> degree{{0, 1, 2, 3, 5, 10, 15, 20}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
            if (n2 > n1)
              continue; // No need to duplicate.

	    auto func = [n1, n2](Tp x)
			-> Tp
			{ return integrand<Tp>(n1, n2, x); };

	    auto [result, error]
		= emsr::integrate_exp_sinh(func, Tp{0}, abs_precision, rel_precision);

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
  num_errors += test_laguerre<float>();
  num_errors += test_laguerre<double>();
  num_errors += test_laguerre<long double>();
  return num_errors;
}
