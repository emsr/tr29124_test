
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_gegenbauer.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  norm_gegenbauer(int n1, int n2, Tp lambda, Tp x)
  {
    const auto
    _S_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    auto gama = std::tgamma(lambda);
    auto gamn2a = std::tgamma(n1 + Tp{2} * lambda);
    auto norm = _S_pi * std::pow(Tp{2}, Tp{1} - Tp{2} * lambda) * gamn2a
	      / emsr::factorial<Tp>(n1) / (Tp(n1) + lambda) / gama / gama;
    return std::pow(Tp{1} - x * x, lambda - Tp{0.5})
	 * emsr::gegenbauer(n1, lambda, x)
	 * emsr::gegenbauer(n2, lambda, x) / norm;
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
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = Tp{10} * rel_prec;

    const bool singular = (lambda < Tp{0.5});

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
	    auto func = [n1, n2, lambda](Tp x)
			-> Tp
			{ return norm_gegenbauer<Tp>(n1, n2, lambda, x); };

	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func,
					Tp{-1}, Tp{1},
					lambda - Tp{0.5}, lambda - Tp{0.5},
					0, 0,
					abs_prec, rel_prec)
		: emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1}, abs_prec, rel_prec);

            auto del = delta<Tp>(n1, n2) - result;
	    if (std::abs(del) > cmp_prec)
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
