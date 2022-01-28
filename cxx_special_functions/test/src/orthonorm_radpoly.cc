
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_jacobi.h>

#include "verify.h"

// Neumann's number
template<typename Tp>
  Tp
  epsilon(int m)
  { return m == 0 ? Tp{2} : Tp{1}; }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  norm_radpoly(int n1, int m1, int n2, int m2, Tp rho)
  {
    auto norm = Tp{1} / std::sqrt(Tp(2 * n1 + 2) * Tp(2 * n2 + 2));
    return rho
	 * emsr::radpoly(n1, m1, rho)
	 * emsr::radpoly(n2, m2, rho) / norm;
  }

template<typename Tp>
  Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_radpoly()
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 3, 4, 5, 8, 9, 16, 17}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
    for (int m1 = n1; m1 >= 0; m1 -= 2)
      {
	if ((n1 - m1) & 1)
	  continue;
	for (const auto n2 : degree)
	  {
	for (int m2 = n2; m2 >= 0; m2 -= 2)
	  {
	    //int m2 = m1;
	    if (m2 > n2)
	      continue;
	    if ((n2 - m2) & 1)
	      continue;
            // Orthonormality only works if the m's are the same.
            if (m2 != m1)
	      continue;

	    auto func = [n1, m1, n2, m2](Tp x)
			-> Tp
			{ return norm_radpoly(n1, m1, n2, m2, x); };

	    auto [result, error]
		= emsr::integrate_tanh_sinh(func, Tp{0}, Tp{1}, abs_prec, rel_prec, 8);

            auto del = delta<Tp>(n1, m1, n2, m2) - result;
	    if (std::abs(del) > cmp_prec)
              {
		++num_errors;
        	std::cout << "n1 = " << n1 << "; m1 = " << m1 << "; n2 = " << n2 << "; m2 = " << m2
                          << "; res = " << result << "; del = " << del << '\n';
              }
	  }
	  }
      }
      }
    return num_errors;
  }

int
main()
{
  int num_errors = 0;
  num_errors += test_radpoly<float>();
  num_errors += test_radpoly<double>();
  num_errors += test_radpoly<long double>();
  return num_errors;
}
