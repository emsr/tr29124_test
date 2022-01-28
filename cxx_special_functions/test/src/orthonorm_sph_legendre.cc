
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_legendre.h>

#include "verify.h"

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename Tp>
  Tp
  norm_sph_legendre(int l1, int m1, int l2, int m2, Tp theta)
  {
    const auto
    _S_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    return Tp{2} * _S_pi * std::sin(theta)
	 * emsr::sph_legendre(l1, m1, theta)
	 * emsr::sph_legendre(l2, m2, theta);
  }

template<typename Tp>
  Tp
  delta(int l1, int l2)
  { return l1 == l2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_sph_legendre(int m1, int m2)
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto
    _S_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int num_errors = 0;

    for (const auto l1 : degree)
      {
	if (m1 > l1)
	  continue;
	for (const auto l2 : degree)
	  {
	    if (m2 > l2)
	      continue;

	    auto func = [l1, m1, l2, m2](Tp theta)
			-> Tp
			{ return norm_sph_legendre(l1, m1, l2, m2, theta); };

	    auto [result, error]
		= emsr::integrate_tanh_sinh(func, Tp{0}, _S_pi, abs_prec, rel_prec, 6);

            auto del = delta<Tp>(l1, l2) - result;
	    if (std::abs(del) > cmp_prec)
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
  num_errors += test_sph_legendre<float>(0, 0);
  num_errors += test_sph_legendre<double>(0, 0);
  num_errors += test_sph_legendre<long double>(0, 0);
  return num_errors;
}
