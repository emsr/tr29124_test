
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_legendre.h>

// Normalized associated Legendre polynomial.
template<typename Tp>
  Tp
  normalized_assoc_legendre(int l, int m, Tp x)
  {
    return std::sqrt(Tp(2 * l + 1) / Tp{2})
	 * emsr::assoc_legendre(l, m, x);
  }

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int l1, int m1, int l2, int m2, Tp x)
  {
    return normalized_assoc_legendre(l1, m1, x)
         * normalized_assoc_legendre(l2, m2, x);
  }

template<typename Tp>
  Tp
  delta(int l1, int l2)
  { return l1 == l2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_assoc_legendre(int m1, int m2)
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    const std::array<int, 10> degree{{0, 1, 2, 3, 5, 10, 15, 20}};
    int num_errors = 0;

    for (const auto l1 : degree)
      {
	for (const auto l2 : degree)
	  {
            if (l2 > l1)
              continue; // No need to duplicate.

	    auto func = [l1, m1, l2, m2](Tp x)
			-> Tp
			{ return integrand(l1, m1, l2, m2, x); };

	    auto [result, error]
		= emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1}, abs_precision, rel_precision, 6);

            auto del = delta<Tp>(l1, l2) - result;
	    if (std::abs(del) > cmp_precision)
              {
		++num_errors;
        	std::cout << "l1 = " << l1 << "; m1 = " << m1 << "; l2 = " << l2 << "; m2 = " << m2
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
  num_errors += test_assoc_legendre<float>(0, 0);
  num_errors += test_assoc_legendre<double>(0, 0);
  num_errors += test_assoc_legendre<long double>(0, 0);
  return num_errors;
}
