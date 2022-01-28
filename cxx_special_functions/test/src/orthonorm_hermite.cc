
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_hermite.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  norm_hermite(int n1, int n2, Tp x)
  {
    const auto
    _S_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    auto lnorm = Tp{0.5} * (log(_S_pi) + Tp(n1 + n2) * log(Tp{2})
			  + emsr::lfactorial<Tp>(n1)
			  + emsr::lfactorial<Tp>(n2));
    return emsr::hermite(n2, x)
	 * std::exp(-x * x - lnorm)
	 * emsr::hermite(n1, x);
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
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
	    auto func = [n1, n2](Tp x)
			-> Tp
			{ return norm_hermite(n1, n2, x); };

	    auto [result, error]
		  = emsr::integrate_sinh_sinh(func, abs_prec, rel_prec);

            auto del = delta<Tp>(n1, n2) - result;
	    if (std::abs(del) > cmp_prec)
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
