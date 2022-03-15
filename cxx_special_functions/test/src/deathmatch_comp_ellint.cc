
//  comp_ellint_1/comp_ellint_2
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
  test_legendre(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    const auto _S_pi2 = emsr::pi_v<Tp> / Tp{2};
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto k = Tp{0.1L} * i;
	auto kp = std::sqrt(Tp{1} - k * k);
	auto K = emsr::comp_ellint_1(k);
	auto Kp = emsr::comp_ellint_1(kp);
	auto E = emsr::comp_ellint_2(k);
	auto Ep = emsr::comp_ellint_2(kp);
	stats << std::make_pair(K * Ep + E * Kp - K * Kp, _S_pi2);
      }
    int num_errors = 0;
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }


template<typename Tp>
  int
  test_theta_3(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
   const auto _S_pi2 = emsr::pi_v<Tp> / Tp{2};
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto k = Tp{0.1L} * i;
	auto q = emsr::ellnome(k);
	auto K = emsr::comp_ellint_1(k);
	auto tht3 = emsr::theta_3(Tp{0}, q);
	stats << std::make_pair(_S_pi2 * tht3 * tht3, K);
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

    num_errors += test_legendre(toler);
    num_errors += test_theta_3(toler);

    return num_errors;
  }

int
main()
{
  return test<double>();
}
