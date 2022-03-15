
//  conf_hyperg
//  Test relations between special functions.

#include <limits>

#include <emsr/special_functions.h>
#include <specfun_stats.h>

#include <iostream>
#define VERIFY(A) \
  if (!(A)) \
    { \
      ++num_errors; \
      std::cout << "line " << __LINE__ \
	<< "  max_abs_frac = " << stats.max_abs_frac \
	<< '\n'; \
    } \
  else \
    { \
      std::cout << "Success: " \
	<< "  samples = " << stats.num \
	<< "  mean = " << stats.mean_diff() \
	<< "  variance = " << stats.variance_diff() \
	<< '\n'; \
    }

template<typename Tp>
  int
  test_cyl_bessel_j(Tp nu, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = std::pow(z / Tp{2}, nu)
		 * emsr::conf_hyperg_lim(nu + Tp{1}, -z * z / Tp{4})
		 / std::tgamma(nu + Tp{1});
	auto Jnu = emsr::cyl_bessel_j(nu, z);
	stats << std::make_pair(hyp, Jnu);
      }
    int num_errors = 0;
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    int num_errors = 0;

    num_errors += test_cyl_bessel_j(Tp{0.5L}, toler);
    num_errors += test_cyl_bessel_j(Tp{1.5L}, toler);
    num_errors += test_cyl_bessel_j(Tp{2.0L}, toler);
    num_errors += test_cyl_bessel_j(Tp{3.5L}, toler);
    num_errors += test_cyl_bessel_j(Tp{5.0L}, toler);

    return num_errors;
  }

int
main()
{
  return test<double>();
}
