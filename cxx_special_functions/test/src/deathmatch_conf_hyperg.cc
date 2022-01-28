
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
  test_erf(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    const auto _S_sqrt_pi = emsr::sqrtpi_v<Tp>;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = Tp{2} * z * emsr::conf_hyperg(Tp{0.5L}, Tp{1.5L}, -z * z) / _S_sqrt_pi;
	auto erf = std::erf(z);
	stats << std::make_pair(hyp, erf);
      }
    int num_errors = 0;
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_exp(Tp a, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = emsr::conf_hyperg(a, a, z);
	auto exp = std::exp(z);
	stats << std::make_pair(hyp, exp);
      }
    int num_errors = 0;
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_sinh(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = z * std::exp(-z) * emsr::conf_hyperg(Tp{1}, Tp{2}, Tp{2} * z);
	auto exp = std::sinh(z);
	stats << std::make_pair(hyp, exp);
      }
    int num_errors = 0;
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_kummer_xform_m(Tp a, Tp c, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto M = emsr::conf_hyperg(a, c, z);
	auto Mt = std::exp(z) * emsr::conf_hyperg(c - a, c, -z);
	stats << std::make_pair(M, Mt);
      }
    int num_errors = 0;
    VERIFY(stats.max_abs_frac < toler);
    return num_errors;
  }

template<typename Tp>
  int
  test_kummer_xform_u(Tp a, Tp c, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto U = emsr::tricomi_u(a, c, z);
	auto Ut = std::pow(z, Tp{1} - c)
		* emsr::tricomi_u(Tp{1} + a - c, Tp{1} - c, z);
	stats << std::make_pair(U, Ut);
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

    num_errors += test_erf<Tp>(toler);

    num_errors += test_sinh<Tp>(toler);

    num_errors += test_exp(Tp{1.5L}, toler);
    num_errors += test_exp(Tp{2.0L}, toler);
    num_errors += test_exp(Tp{3.5L}, toler);
    num_errors += test_exp(Tp{5.0L}, toler);

    num_errors += test_kummer_xform_m(Tp{1.5L}, Tp{0.5L}, toler);
    num_errors += test_kummer_xform_m(Tp{2.0L}, Tp{2.5L}, toler);
    num_errors += test_kummer_xform_m(Tp{3.5L}, Tp{0.5L}, toler);
    num_errors += test_kummer_xform_m(Tp{5.0L}, Tp{2.5L}, toler);

    num_errors += test_kummer_xform_u(Tp{1.5L}, Tp{0.5L}, toler);
    num_errors += test_kummer_xform_u(Tp{2.0L}, Tp{2.5L}, toler);
    num_errors += test_kummer_xform_u(Tp{3.5L}, Tp{0.5L}, toler);
    num_errors += test_kummer_xform_u(Tp{5.0L}, Tp{2.5L}, toler);

    return num_errors;
  }

int
main()
{
  return test<double>();
}
