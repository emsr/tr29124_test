
//  hyperg
//  Test relations between special functions.

#include <limits>

#include <emsr/special_functions.h>

#include <specfun_testcase.h>
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
  test_log(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    int num_errors = 0;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = z * emsr::hyperg(Tp{1}, Tp{1}, Tp{2}, -z);
	auto ln1p = std::log1p(z);
	stats << std::make_pair(hyp, ln1p);
      }
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_pow(Tp a, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    int num_errors = 0;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = emsr::hyperg(a, Tp{1}, Tp{1}, z);
	auto binom = std::pow(Tp{1} - z, -a);
	stats << std::make_pair(hyp, binom);
      }
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_asin(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    int num_errors = 0;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = z * emsr::hyperg(Tp{0.5L}, Tp{0.5L}, Tp{1.5L}, z * z);
	auto arcsin = std::asin(z);
	stats << std::make_pair(hyp, arcsin);
      }
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_comp_ellint_1(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    const auto _S_pi2 = emsr::pi_v<Tp> / Tp{2};
    int num_errors = 0;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto k = Tp{0.1L} * i;
	auto hyp = _S_pi2 * emsr::hyperg(Tp{0.5L}, Tp{0.5L}, Tp{1}, k * k);
	auto K = std::comp_ellint_1(k);
	stats << std::make_pair(hyp, K);
      }
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_comp_ellint_2(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    const auto _S_pi2 = emsr::pi_v<Tp> / Tp{2};
    int num_errors = 0;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto k = Tp{0.1L} * i;
	auto hyp = _S_pi2 * emsr::hyperg(Tp{-0.5L}, Tp{0.5L}, Tp{1}, k * k);
	auto E = std::comp_ellint_2(k);
	stats << std::make_pair(hyp, E);
      }
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_jacobi(unsigned int n, Tp alpha, Tp beta,
	      Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    int num_errors = 0;
    Stats<Tp> stats(toler);
    for (auto i = -9; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = emsr::hyperg(-Tp(n), Tp(1 + n) + alpha + beta, Tp{1} + alpha, z);
	auto P = emsr::factorial<Tp>(n)
	       * emsr::jacobi(n, alpha, beta, Tp{1} - Tp{2} * z)
	       / emsr::rising_factorial(alpha + Tp{1}, int(n));
	stats << std::make_pair(hyp, P);
      }
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    int num_errors = 0;

    num_errors += test_log<Tp>(toler);

    num_errors += test_pow(Tp{0.5L}, toler);
    num_errors += test_pow(Tp{1.5L}, toler);
    num_errors += test_pow(Tp{5.0L}, toler);

    num_errors += test_asin<Tp>(toler);

    num_errors += test_comp_ellint_1<Tp>(toler);
    num_errors += test_comp_ellint_2<Tp>(toler);

    num_errors += test_jacobi(0, Tp{0.5L}, Tp{0.5L}, toler);
    num_errors += test_jacobi(3, Tp{0.5L}, Tp{0.5L}, toler);
    num_errors += test_jacobi(3, Tp{1.5L}, Tp{4.0L}, toler);

    return num_errors;
  }

int
main()
{
  return test<double>();
}
