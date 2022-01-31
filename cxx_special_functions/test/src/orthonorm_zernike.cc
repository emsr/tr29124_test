
#include <iostream>

#include <emsr/math_constants.h>
#include <emsr/integration.h>
#include <emsr/sf_jacobi.h>

#include "verify.h"

// Neumann's number
template<typename Tp>
  Tp
  epsilon(int m)
  { return m == 0 ? Tp{2} : Tp{1}; }

// Azimuthal integral of zernike product.

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  integrand(int n1, int m1, int n2, int m2, Tp rho)
  {
    const auto s_eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto s_eps = s_eps_factor * std::numeric_limits<Tp>::epsilon();
    const auto s_2pi = Tp{2} * emsr::pi_v<Tp>;

    // Normalized Zernike polynomials.
    auto nz1 = [n1, m1, rho](Tp phi)
	       -> Tp
	       { return std::sqrt(Tp(2 * n1 + 2)) * emsr::zernike(n1, m1, rho, phi); };
    auto nz2 = [n2, m2, rho](Tp phi)
	       -> Tp
	       { return std::sqrt(Tp(2 * n2 + 2)) * emsr::zernike(n2, m2, rho, phi); };

    auto fun = [rho, nz1, nz2](Tp phi)
		-> Tp
		{ return rho * nz1(phi) * nz2(phi); };

    auto val
	= emsr::integrate_tanh_sinh(fun, Tp{0}, Tp{s_2pi}, s_eps, s_eps, 8);

    return Tp{2} * val.result / s_2pi / epsilon<Tp>(m1);
  }

template<typename Tp>
  Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_zernike(int m1, int m2)
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = Tp{10} * rel_precision;

    const std::array<int, 10> degree{{0, 1, 2, 3, 4, 5, 8, 9, 16, 17}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	if ((n1 - m1) & 1)
	  continue;
	for (const auto n2 : degree)
	  {
            // No need to duplicate.
            if (n2 > n1)
              continue;
	    if (m2 > n2) // ?
	      continue;
	    if ((n2 - m2) & 1)
	      continue;

	    auto func = [n1, m1, n2, m2](Tp x)
			-> Tp
			{ return integrand(n1, m1, n2, m2, x); };

	    auto [result, error]
		= emsr::integrate_tanh_sinh(func, Tp{0}, Tp{1},
					abs_precision, rel_precision, 8);

            auto del = delta<Tp>(n1, m1, n2, m2) - result;
	    if (std::abs(del) > cmp_precision)
              {
		++num_errors;
        	std::cout << "n1 = " << n1 << "; m1 = " << m1 << "; n2 = " << n2 << "; m2 = " << m2
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
  num_errors += test_zernike<float>(1, 2);
  num_errors += test_zernike<double>(1, 2);
  num_errors += test_zernike<long double>(1, 2);
  return num_errors;
}
