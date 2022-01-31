
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_chebyshev.h>
#include <emsr/math_constants.h>

// Normailized Chebyshev polynomial of the second kind.
template<typename Tp>
  Tp
  normalized_chebyshev_u(int n, Tp x)
  {
    const auto s_sqrtpid2 = emsr::sqrtpi_v<Tp> / emsr::sqrt2_v<Tp>;
    return emsr::chebyshev_u(n, x) / s_sqrtpid2;
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int n2, Tp x)
  {
    return normalized_chebyshev_u(n2, x)
	 * normalized_chebyshev_u(n1, x)
	 * std::sqrt(Tp{1} - x * x);
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_chebyshev_u()
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
			{return integrand(n1, n2, x);};

	    auto [result, error]
		= emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1}, abs_precision, rel_precision, 6);

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
  num_errors += test_chebyshev_u<float>();
  num_errors += test_chebyshev_u<double>();
  num_errors += test_chebyshev_u<long double>();
  return num_errors;
}
