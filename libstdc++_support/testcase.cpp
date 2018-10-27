#define STD !TR1

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#if STD
#  include <cmath>
#else
#  include <tr1/cmath>
#endif
#include <functional>
#include <utility>
#include <tuple>

#include "specfun_testcase.h"
#include "../wrappers/wrap_gsl.h"
#include "../wrappers/wrap_boost.h"
#include "../wrappers/wrap_burkhardt.h"

#include "testcase.tcc"

std::string
get_filename(const std::string & path,
	     const std::string & prefix,
	     const std::string & basename,
	     const std::string & extra,
	     const std::string & suffix)
{
  auto filename = path + "/" + prefix;
  filename += basename + extra + suffix;

  return filename;
}

/**
 * Return the Chebyshev polynomial of the first kind by trigonometric identity:
 * @f[
 *   T_n(x) = \cos(n\theta)
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename _Tp>
  _Tp
  __chebyshev_t_trig(unsigned int __n, _Tp __x)
  {
    auto __theta = std::acos(__x);
    return std::cos(__n * __theta);
  }

/**
 * Return the Chebyshev polynomial of the second kind by trigonometric identity:
 * @f[
 *   U_n(x) = \frac{\sin((n + 1)\theta)}{\sin(\theta)}
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename _Tp>
  _Tp
  __chebyshev_u_trig(unsigned int __n, _Tp __x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    if (std::abs(__x + _Tp{1}) < _S_eps)
      return (__n % 2 == 0 ? +1 : -1) * _Tp(__n + 1);
    else if (std::abs(__x - _Tp{1}) < _S_eps)
      return _Tp(__n + 1);
    else
      {
	auto __theta = std::acos(__x);
	return std::sin(_Tp(__n + 1) * __theta)
	     / std::sin(__theta);
      }
  }

/**
 * Return the Chebyshev polynomial of the third kind by trigonometric identity:
 * @f[
 *   V_n(x) = \frac{\cos((n + \frac{1}{2})\theta)}
 *                 {cos(\frac{1}{2}\theta)}
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename _Tp>
  _Tp
  __chebyshev_v_trig(unsigned int __n, _Tp __x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    if (std::abs(__x + _Tp{1}) < _S_eps)
      return (__n % 2 == 0 ? +1 : -1) * _Tp(2 * __n + 1);
    else
      {
	auto __theta = std::acos(__x);
	return std::cos(_Tp(__n + _Tp{1} / _Tp{2}) * __theta)
	     / std::cos(__theta / _Tp{2});
      }
  }

/**
 * Return the Chebyshev polynomial of the fourth kind by trigonometric identity:
 * @f[
 *   W_n(x) = \frac{\sin((n + \frac{1}{2})\theta)}
 *                 {\sin(\frac{1}{2}\theta)}
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename _Tp>
  _Tp
  __chebyshev_w_trig(unsigned int __n, _Tp __x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    if (std::abs(__x - _Tp{1}) < _S_eps)
      return _Tp(2 * __n + 1);
    else
      {
	auto __theta = std::acos(__x);
	return std::sin(_Tp(__n + _Tp{1} / _Tp{2}) * __theta)
	     / std::sin(__theta / _Tp{2});
      }
  }

double
__chebyshev_t(unsigned int __n, double __x)
{ return __chebyshev_t_trig<double>(__n, __x); }

double
__chebyshev_u(unsigned int __n, double __x)
{ return __chebyshev_u_trig<double>(__n, __x); }

double
__chebyshev_v(unsigned int __n, double __x)
{ return __chebyshev_v_trig<double>(__n, __x); }

double
__chebyshev_w(unsigned int __n, double __x)
{ return __chebyshev_w_trig<double>(__n, __x); }

template<typename Real>
  void
  harness()
  {
#if STD
    using __gnu_cxx::airy_ai;
    using __gnu_cxx::airy_bi;
    using       std::assoc_laguerre;
    using       std::assoc_legendre;
    using __gnu_cxx::bernoulli;
    using       std::beta;
    using __gnu_cxx::binomial;
    using __gnu_cxx::chebyshev_t;
    using __gnu_cxx::chebyshev_u;
    using __gnu_cxx::chebyshev_v;
    using __gnu_cxx::chebyshev_w;
    using __gnu_cxx::clausen;
    using __gnu_cxx::clausen_cl;
    using __gnu_cxx::clausen_sl;
    using __gnu_cxx::comp_ellint_d;
    using       std::comp_ellint_1;
    using       std::comp_ellint_2;
    using       std::comp_ellint_3;
    using __gnu_cxx::conf_hyperg;
    using __gnu_cxx::conf_hyperg_lim;
    using __gnu_cxx::coshint;
    using __gnu_cxx::cosint;
    using __gnu_cxx::cos_pi;
    using       std::cyl_bessel_i;
    using       std::cyl_bessel_j;
    using       std::cyl_bessel_k;
    using __gnu_cxx::cyl_hankel_1;
    using __gnu_cxx::cyl_hankel_2;
    using       std::cyl_neumann;
    using __gnu_cxx::dawson;
    using __gnu_cxx::debye;
    using __gnu_cxx::dilog;
    using __gnu_cxx::dirichlet_beta;
    using __gnu_cxx::dirichlet_eta;
    using __gnu_cxx::double_factorial;
    using       std::ellint_1;
    using       std::ellint_2;
    using       std::ellint_3;
    using __gnu_cxx::ellint_d;
    using __gnu_cxx::ellint_rc;
    using __gnu_cxx::ellint_rd;
    using __gnu_cxx::ellint_rf;
    using __gnu_cxx::ellint_rg;
    using __gnu_cxx::ellint_rj;
    using __gnu_cxx::euler;
    using __gnu_cxx::eulerian_1;
    using       std::expint;
    using __gnu_cxx::expint;
    using __gnu_cxx::factorial;
    using __gnu_cxx::falling_factorial;
    using __gnu_cxx::fresnel_c;
    using __gnu_cxx::fresnel_s;
    using __gnu_cxx::gegenbauer;
    using       std::hermite;
    using __gnu_cxx::heuman_lambda;
    using __gnu_cxx::hurwitz_zeta;
    using __gnu_cxx::hyperg;
    using __gnu_cxx::ibeta;
    using __gnu_cxx::ibetac;
    using __gnu_cxx::jacobi;
    using __gnu_cxx::jacobi_sn;
    using __gnu_cxx::jacobi_cn;
    using __gnu_cxx::jacobi_dn;
    using __gnu_cxx::jacobi_zeta;
    using       std::laguerre;
    using __gnu_cxx::lbinomial;
    using __gnu_cxx::ldouble_factorial;
    using       std::legendre;
    using __gnu_cxx::legendre_q;
    using __gnu_cxx::lfactorial;
    using __gnu_cxx::lfalling_factorial;
    using __gnu_cxx::lrising_factorial;
    using __gnu_cxx::owens_t;
    using __gnu_cxx::gamma_p;
    using __gnu_cxx::polygamma;
    using __gnu_cxx::digamma;
    using __gnu_cxx::gamma_q;
    using __gnu_cxx::radpoly;
    using       std::riemann_zeta;
    using __gnu_cxx::rising_factorial;
    using __gnu_cxx::sinc;
    using __gnu_cxx::sinc_pi;
    using __gnu_cxx::sinhc;
    using __gnu_cxx::sinhc_pi;
    using __gnu_cxx::sinhint;
    using __gnu_cxx::sinint;
    using __gnu_cxx::sin_pi;
    using       std::sph_bessel;
    using __gnu_cxx::sph_bessel_i;
    using __gnu_cxx::sph_bessel_k;
    using __gnu_cxx::sph_hankel_1;
    using __gnu_cxx::sph_hankel_2;
    using __gnu_cxx::sph_harmonic;
    using       std::sph_legendre;
    using       std::sph_neumann;
    using __gnu_cxx::stirling_1;
    using __gnu_cxx::stirling_2;
    using __gnu_cxx::tgamma_lower;
    using __gnu_cxx::tgamma;
    using __gnu_cxx::theta_1;
    using __gnu_cxx::theta_2;
    using __gnu_cxx::theta_3;
    using __gnu_cxx::theta_4;
    using __gnu_cxx::theta_s;
    using __gnu_cxx::theta_c;
    using __gnu_cxx::theta_d;
    using __gnu_cxx::theta_n;
    using __gnu_cxx::zernike;
#else
    using  std::tr1::assoc_laguerre;
    using  std::tr1::assoc_legendre;
    using  std::tr1::beta;
    using  std::tr1::comp_ellint_1;
    using  std::tr1::comp_ellint_2;
    using  std::tr1::comp_ellint_3;
    using  std::tr1::conf_hyperg;
    using  std::tr1::cyl_bessel_i;
    using  std::tr1::cyl_bessel_j;
    using  std::tr1::cyl_bessel_k;
    using  std::tr1::cyl_neumann;
    using  std::tr1::ellint_1;
    using  std::tr1::ellint_2;
    using  std::tr1::ellint_3;
    using  std::tr1::expint;
    using  std::tr1::hermite;
    using  std::tr1::hyperg;
    using  std::tr1::laguerre;
    using  std::tr1::legendre;
    using  std::tr1::riemann_zeta;
    using  std::tr1::sph_bessel;
    using  std::tr1::sph_legendre;
    using  std::tr1::sph_neumann;
#endif

    const auto _S_pi = __gnu_cxx::__math_constants<Real>::__pi;

    // Unsigned integer orders for various polynomials, harmonics.
    std::vector<unsigned int> vorder{0, 1, 2, 5, 10, 20, 50, 100};

    // Integer orders for various polynomials, harmonics, and spherical bessels.
    std::vector<int> iorder{0, 1, 2, 5, 10, 20, 50, 100};

    // ... corresponding "Real" integer orders for GSL.
    std::vector<Real> dvorder{0, 1, 2, 5, 10, 20, 50, 100};

    // Orders for cylindrical Bessel functions.
    std::vector<Real> cyl_neg_order{-5, -2, -1, -Real{2}/Real{3},
				    -Real{1}/Real{2}, -Real{1}/Real{3}};

    std::vector<Real> cyl_order{0, Real{1}/Real{3},
				Real{1}/Real{2}, Real{2}/Real{3},
				1, 2, 5, 10, 20, 50, 100};

    // Orders for spherical bessel functions.
    std::vector<unsigned int> sph_order{0, 1, 2, 5, 10, 20, 50, 100};

    const unsigned int num_phi = 10;
    std::vector<Real> vphid;
    for (unsigned int i = 0; i < num_phi - 1; ++i)
      vphid.push_back(Real{10} * i * _S_pi / Real{180});
    vphid.push_back(__gnu_cxx::__math_constants<Real>::__pi_half);

    std::vector<Real> vnopolesd;
    for (unsigned int i = 1; i < num_phi - 1; ++i)
      vnopolesd.push_back(Real{10} * i * _S_pi / Real{180});

    std::vector<Real> vab{0, Real{1}/Real{2}, 1, 2, 5, 10, 20};

    unsigned int test = 1;

    const std::string path = "check";

#if STD
    std::string nsname = "std";
    const std::string prefix = "/check_";
#else
    std::string nsname = "std::tr1";
    const std::string prefix = "/check_tr1_";
#endif

    std::string basename;
    std::string filename;

#if STD
    // Airy functions of the first kind.
    std::cout << "airy_ai\n" << std::flush;
    basename = "airy_ai";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_airy_ai(filename);
    maketest(airy_ai, gsl::airy_ai,
	     "testcase_airy_ai", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_airy_ai);

    // Airy functions of the second kind.
    std::cout << "airy_bi\n" << std::flush;
    basename = "airy_bi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_airy_bi(filename);
    maketest(airy_bi, gsl::airy_bi,
	     "testcase_airy_bi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_airy_bi);
#endif // STD

    // Associated Laguerre polynomials.
    std::cout << "assoc_laguerre\n" << std::flush;
    basename = "assoc_laguerre";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_assoc_laguerre(filename);
    maketest(assoc_laguerre, gsl::assoc_laguerre,
	     "testcase_assoc_laguerre", nsname, basename,
	     "n", vorder, "m", vorder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_assoc_laguerre);

    // Associated Legendre functions.
    std::cout << "assoc_legendre\n" << std::flush;
    basename = "assoc_legendre";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_assoc_legendre(filename);
    maketest(assoc_legendre, gsl::assoc_legendre,
	     "testcase_assoc_legendre", nsname, basename,
	     "l", vorder, "m", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_assoc_legendre);

    // Beta functions.
    std::cout << "beta\n" << std::flush;
    basename = "beta";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_beta(filename);
    maketest(beta, gsl::beta,
	     "testcase_beta", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 11),
	     "GSL",
	     file_beta);

    // Complete elliptic integrals of the first kind.
    // Avoid poles at |x| = 1.
    std::cout << "comp_ellint_1\n" << std::flush;
    basename = "comp_ellint_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_comp_ellint_1(filename);
    maketest(comp_ellint_1, beast::comp_ellint_1,
	     "testcase_comp_ellint_1", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "Boost",
	     file_comp_ellint_1);

    // Complete elliptic integrals of the second kind.
    // Avoid poles at |x| = 1.
    std::cout << "comp_ellint_2\n" << std::flush;
    basename = "comp_ellint_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_comp_ellint_2(filename);
    maketest(comp_ellint_2, beast::comp_ellint_2,
	     "testcase_comp_ellint_2", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "Boost",
	     file_comp_ellint_2);

    // Complete elliptic integrals of the third kind.
    // Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3\n" << std::flush;
    basename = "comp_ellint_3";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_comp_ellint_3(filename);
    maketest(comp_ellint_3, beast::comp_ellint_3,
	     "testcase_comp_ellint_3", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, false), 11),
	     "Boost",
	     file_comp_ellint_3);

#if STD
    // Confluent hypergeometric functions.
    // Skip the singularity at c = 0.
    std::cout << "conf_hyperg\n" << std::flush;
    basename = "conf_hyperg";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_conf_hyperg(filename);
    maketest(conf_hyperg, gsl::conf_hyperg,
	     "testcase_conf_hyperg", "__gnu_cxx", basename,
	     "a", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_conf_hyperg);
#else
    // Confluent hypergeometric functions.
    // Skip the singularity at c = 0.
    std::cout << "conf_hyperg\n" << std::flush;
    basename = "conf_hyperg";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_conf_hyperg(filename);
    maketest(conf_hyperg, gsl::conf_hyperg,
	     "testcase_conf_hyperg", nsname, basename,
	     "a", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_conf_hyperg);
#endif // STD

    // Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i\n" << std::flush;
    basename = "cyl_bessel_i";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_i(filename);
    test =
    maketest(cyl_bessel_i, beast::cyl_bessel_i,
	     "testcase_cyl_bessel_i", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_bessel_i, true, false);
    test =
    maketest(cyl_bessel_i, gsl::cyl_bessel_i,
	     "testcase_cyl_bessel_i", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_cyl_bessel_i, false, false, test);
    maketest(cyl_bessel_i, gsl::cyl_bessel_i,
	     "testcase_cyl_bessel_i", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_cyl_bessel_i, false, true, test);

    // Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j\n" << std::flush;
    basename = "cyl_bessel_j";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_j(filename);
    test =
    maketest(cyl_bessel_j, beast::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_bessel_j, true, false);
    test =
    maketest(cyl_bessel_j, gsl::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_cyl_bessel_j, false, false, test);
    test =
    maketest(cyl_bessel_j, gsl::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_cyl_bessel_j, false, false, test);
    maketest(cyl_bessel_j, gsl::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", {100},
	     "x", fill_argument(std::make_pair(Real{1000}, Real{2000}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_cyl_bessel_j, false, true, test);

    // Irregular modified cylindrical Bessel functions.
    // Skip the pole at the origin.
    std::cout << "cyl_bessel_k\n" << std::flush;
    basename = "cyl_bessel_k";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_k(filename);
    test =
    maketest(cyl_bessel_k, beast::cyl_bessel_k,
	     "testcase_cyl_bessel_k", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_bessel_k, true, false);
    test =
    maketest(cyl_bessel_k, gsl::cyl_bessel_k,
	     "testcase_cyl_bessel_k", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_bessel_k, false, false, test);
    maketest(cyl_bessel_k, gsl::cyl_bessel_k,
	     "testcase_cyl_bessel_k", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_bessel_k, false, true, test);

    // Cylindrical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "cyl_neumann\n" << std::flush;
    basename = "cyl_neumann";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_neumann(filename);
    test =
    maketest(cyl_neumann, beast::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_neumann, true, false);
    test =
    maketest(cyl_neumann, gsl::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_neumann, false, false, test);
    test =
    maketest(cyl_neumann, gsl::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_neumann, false, false, test);
    maketest(cyl_neumann, gsl::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", {100},
	     "x", fill_argument(std::make_pair(Real{1000}, Real{2000}),
				std::make_pair(false, true), 11),
	     "GSL",
	     file_cyl_neumann, false, true, test);

    // Elliptic integrals of the first kind.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_1\n" << std::flush;
    basename = "ellint_1";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_1(filename);
    maketest(ellint_1, beast::ellint_1,
	     "testcase_ellint_1", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     "Boost",
	     file_ellint_1);

    // Elliptic integrals of the second kind.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_2\n" << std::flush;
    basename = "ellint_2";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_2(filename);
    maketest(ellint_2, beast::ellint_2,
	     "testcase_ellint_2", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     "Boost",
	     file_ellint_2);

    // Elliptic integrals of the third kind.
    // Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3\n" << std::flush;
    basename = "ellint_3";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_ellint_3(filename);
    maketest(ellint_3, beast::ellint_3,
	     "testcase_ellint_3", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, false), 11),
	     "phi", vphid,
	     "Boost",
	     file_ellint_3);

    // Exponential integral.
    // Skip the pole at 0.
    std::cout << "expint\n" << std::flush;
    basename = "expint";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_expint(filename);
    test =
    maketest(expint, gsl::expint,
	     "testcase_expint", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-50}, Real{0}),
				std::make_pair(true, false), 51),
	     "GSL",
	     file_expint, true, false);
    maketest(expint, gsl::expint,
	     "testcase_expint", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{50}),
				std::make_pair(false, true), 51),
	     "GSL",
	     file_expint, false, true, test);

    // Hermite polynomials
    std::cout << "hermite\n" << std::flush;
    basename = "hermite";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hermite(filename);
    test =
    maketest(hermite, gsl::hermite,
	     "testcase_hermite", nsname, basename,
	     "n", vorder,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 201),
	     "GSL",
	     file_hermite, true, false);
    maketest(hermite, gsl::hermite,
	     "testcase_hermite", nsname, basename,
	     "n", {8, 18, 32, 50, 72, 128, 200, 1250, 5000},
	     "x", {Real{4}, Real{6}, Real{8}, Real{10}, Real{12},
		   Real{16}, Real{20}, Real{50}, Real{100}},
	     "GSL",
	     file_hermite, false, true, test);

#if STD
    // Hypergeometric functions.
    // Skip the singularity at c = 0.
    // Skip the singularity at x = -1.
    std::cout << "hyperg\n" << std::flush;
    basename = "hyperg";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hyperg(filename);
    maketest(hyperg, gsl::hyperg,
	     "testcase_hyperg", "__gnu_cxx", basename,
	     "a", vab,
	     "b", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 6),
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, false), 21),
	     "GSL",
	     file_hyperg);
#else
    // Hypergeometric functions.
    // Skip the singularity at c = 0.
    // Skip the singularity at x = -1.
    std::cout << "hyperg\n" << std::flush;
    basename = "hyperg";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hyperg(filename);
    maketest(hyperg, gsl::hyperg,
	     "testcase_hyperg", nsname, basename,
	     "a", vab,
	     "b", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 6),
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, false), 21),
	     "GSL",
	     file_hyperg);
#endif // STD

    // Laguerre polynomials.
    std::cout << "laguerre\n" << std::flush;
    basename = "laguerre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_laguerre(filename);
    maketest(laguerre, gsl::laguerre,
	     "testcase_laguerre", nsname, basename,
	     "n", vorder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_laguerre);

    // Legendre polynomials.
    std::cout << "legendre\n" << std::flush;
    basename = "legendre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_legendre(filename);
    maketest(legendre, gsl::legendre_p,
	     "testcase_legendre", nsname, basename,
	     "l", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_legendre);

    // Riemann zeta function.
    std::cout << "riemann_zeta\n" << std::flush;
    // Skip the pole at 1.
    basename = "riemann_zeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_riemann_zeta(filename);
    test =
    maketest(riemann_zeta, beast::riemann_zeta,
	     "testcase_riemann_zeta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{-10}, Real{1}),
				std::make_pair(true, false), 56),
	     "Boost",
	     file_riemann_zeta, true, false);
    maketest(riemann_zeta, beast::riemann_zeta,
	     "testcase_riemann_zeta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
				std::make_pair(false, true), 146),
	     "Boost",
	     file_riemann_zeta, false, true, test);

#if STD
    // Hurwitz zeta functions.
    std::cout << "hurwitz_zeta\n" << std::flush;
    // Skip the pole at 1.
    basename = "hurwitz_zeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hurwitz_zeta(filename);
    maketest(hurwitz_zeta, gsl::hurwitz_zeta,
	     "testcase_hurwitz_zeta", "__gnu_cxx", basename,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
				std::make_pair(false, true), 146),
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 26),
	     "GSL",
	     file_hurwitz_zeta);
#endif // STD

    // Spherical Bessel functions.
    std::cout << "sph_bessel\n" << std::flush;
    basename = "sph_bessel";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel(filename);
    test =
    maketest(sph_bessel, gsl::sph_bessel,
	     "testcase_sph_bessel", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_bessel, true, false);
    maketest(sph_bessel, gsl::sph_bessel,
	     "testcase_sph_bessel", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_bessel, false, true, test);

    // Spherical Legendre functions.
    std::cout << "sph_legendre\n" << std::flush;
    basename = "sph_legendre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_legendre(filename);
    maketest(sph_legendre, gsl::sph_legendre,
	     "testcase_sph_legendre", nsname, basename,
	     "l", vorder, "m", vorder,
	     "theta", fill_argument(std::make_pair(Real{0}, _S_pi),
				    std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_legendre);

    // Spherical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "sph_neumann\n" << std::flush;
    basename = "sph_neumann";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_neumann(filename);
    test =
    maketest(sph_neumann, gsl::sph_neumann,
	     "testcase_sph_neumann", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_sph_neumann, true, false);
    maketest(sph_neumann, gsl::sph_neumann,
	     "testcase_sph_neumann", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_sph_neumann, false, true, test);

#if STD
    // Carlson elliptic functions R_C.
    std::cout << "ellint_rc\n" << std::flush;
    basename = "ellint_rc";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rc(filename);
    maketest(ellint_rc, beast::ellint_rc,
	     "testcase_ellint_rc", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "Boost",
	     file_ellint_rc);

    // Carlson elliptic functions R_D.
    std::cout << "ellint_rd\n" << std::flush;
    basename = "ellint_rd";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rd(filename);
    maketest(ellint_rd, beast::ellint_rd,
	     "testcase_ellint_rd", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "Boost",
	     file_ellint_rd);

    // Carlson elliptic functions R_F.
    std::cout << "ellint_rf\n" << std::flush;
    basename = "ellint_rf";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rf(filename);
    maketest(ellint_rf, beast::ellint_rf,
	     "testcase_ellint_rf", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "Boost",
	     file_ellint_rf);

    // Carlson elliptic functions R_G.
    std::cout << "ellint_rg\n" << std::flush;
    basename = "ellint_rg";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rg(filename);
    maketest(ellint_rg, beast::ellint_rg,
	     "testcase_ellint_rg", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "Boost",
	     file_ellint_rg);

    // Carlson elliptic functions R_J.
    std::cout << "ellint_rj\n" << std::flush;
    basename = "ellint_rj";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rj(filename);
    maketest(ellint_rj, beast::ellint_rj,
	     "testcase_ellint_rj", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "p", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "Boost",
	     file_ellint_rj);

    // Dilogarithm functions.
    std::cout << "dilog\n" << std::flush;
    basename = "dilog";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_dilog(filename);
    maketest(dilog, gsl::dilog,
	     "testcase_dilog", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{1}),
				std::make_pair(true, true), 45),
	     "GSL",
	     file_dilog);

    // Log Gamma functions.
    std::cout << "lgamma\n" << std::flush;
    basename = "lgamma";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_lgamma(filename);
    maketest(lgamma, beast::lgamma,
	     "testcase_lgamma", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{-10}, Real{20}),
				std::make_pair(false, true), 151),
	     "Boost",
	     file_lgamma);

    // Gamma functions.
/*
    std::cout << "tgamma\n" << std::flush;
    basename = "tgamma";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_tgamma(filename);
    maketest(tgamma, beast::tgamma,
	     "testcase_tgamma", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "Boost",
	     file_tgamma);
*/
    // Upper incomplete Gamma functions.
    std::cout << "tgamma\n" << std::flush;
    basename = "tgamma";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_tgamma(filename);
    maketest(tgamma, beast::tgamma,
	     "testcase_tgamma", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     "Boost",
	     file_tgamma);

    // Lower incomplete Gamma functions.
    std::cout << "tgamma_lower\n" << std::flush;
    basename = "tgamma_lower";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_tgamma_lower(filename);
    maketest(tgamma_lower, beast::tgamma_lower,
	     "testcase_tgamma_lower", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     "Boost",
	     file_tgamma_lower);

    // Lower regularized incomplete Gamma functions.
    std::cout << "gamma_p\n" << std::flush;
    basename = "gamma_p";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_gamma_p(filename);
    maketest(gamma_p, gsl::gamma_p,
	     "testcase_gamma_p", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_gamma_p);

    // Upper regularized incomplete Gamma functions.
    std::cout << "gamma_q\n" << std::flush;
    basename = "gamma_q";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_gamma_q(filename);
    maketest(gamma_q, gsl::gamma_q,
	     "testcase_gamma_q", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_gamma_q);

    // Incomplete Beta functions.
    std::cout << "ibeta\n" << std::flush;
    basename = "ibeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ibeta(filename);
    maketest<Real, Real, Real, Real>(ibeta, gsl::ibeta,
	     "testcase_ibeta", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "b", fill_argument(std::make_pair(Real{5}, Real{0}),
				std::make_pair(true, false), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{1}),
				std::make_pair(false, false), 21),
	     "GSL",
	     file_ibeta);

    // Complementary incomplete Beta functions.
    std::cout << "ibetac\n" << std::flush;
    basename = "ibetac";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ibetac(filename);
    maketest<Real, Real, Real, Real>(ibetac, beast::ibetac,
	     "testcase_ibetac", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "b", fill_argument(std::make_pair(Real{5}, Real{0}),
				std::make_pair(true, false), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{1}),
				std::make_pair(false, false), 21),
	     "Boost",
	     file_ibetac);

    // Digamma or digamma functions.
    std::cout << "digamma\n" << std::flush;
    basename = "digamma";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_digamma(filename);
    const auto skip = Real{1} / Real{16};
    test =
    maketest(digamma, beast::digamma,
	     "testcase_digamma", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-10} + skip, Real{+10} + skip),
				std::make_pair(true, true), 201),
	     "Boost",
	     file_digamma, true, false);
    maketest(digamma, beast::digamma,
	     "testcase_digamma", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{1}, Real{100}),
				std::make_pair(true, true), 199),
	     "Boost",
	     file_digamma, false, true, test);

    // Psi or digamma functions.
    std::cout << "polygamma\n" << std::flush;
    basename = "polygamma";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_polygamma(filename);
    test =
    maketest(polygamma, beast::polygamma,
	     "testcase_polygamma", "__gnu_cxx", basename,
	     "m", {0U, 1U, 2U, 3U, 5U, 10U},
	     "x", fill_argument(std::make_pair(Real{-10} + skip, Real{+10} + skip),
				std::make_pair(true, true), 201),
	     "Boost",
	     file_polygamma, true, false);
    maketest(polygamma, beast::polygamma,
	     "testcase_polygamma", "__gnu_cxx", basename,
	     "m", {0U, 1U, 2U, 3U, 5U, 10U},
	     "x", fill_argument(std::make_pair(Real{1}, Real{100}),
				std::make_pair(true, true), 199),
	     "Boost",
	     file_polygamma, false, true, test);

    // Sine integral or Si functions.
    std::cout << "sinint\n" << std::flush;
    basename = "sinint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinint(filename);
    maketest(sinint, gsl::sinint,
	     "testcase_sinint", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_sinint);

    // Cosine integral or Ci functions.
    std::cout << "cosint\n" << std::flush;
    basename = "cosint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cosint(filename);
    maketest(cosint, gsl::cosint,
	     "testcase_cosint", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_cosint);

    // Hyperbolic sine integral or Shi functions.
    std::cout << "sinhint\n" << std::flush;
    basename = "sinhint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinhint(filename);
    maketest(sinhint, gsl::sinhint,
	     "testcase_sinhint", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_sinhint);

    // Hyperbolic cosine integral or Chi functions.
    std::cout << "coshint\n" << std::flush;
    basename = "coshint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_coshint(filename);
    maketest(coshint, gsl::coshint,
	     "testcase_coshint", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_coshint);

    // Dawson integral.
    std::cout << "dawson\n" << std::flush;
    basename = "dawson";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_dawson(filename);
    maketest(dawson, gsl::dawson,
	     "testcase_dawson", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+20}),
				std::make_pair(false, true), 201),
	     "GSL",
	     file_dawson);

    // Jacobian elliptic sine amplitude integrals.
    std::cout << "jacobi_sn\n" << std::flush;
    basename = "jacobi_sn";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi_sn(filename);
    maketest(jacobi_sn, beast::jacobi_sn,
	     "testcase_jacobi_sn", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101),
	     "Boost",
	     file_jacobi_sn);

    // Jacobian elliptic cosine amplitude integrals.
    std::cout << "jacobi_cn\n" << std::flush;
    basename = "jacobi_cn";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi_cn(filename);
    maketest(jacobi_cn, beast::jacobi_cn,
	     "testcase_jacobi_cn", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101),
	     "Boost",
	     file_jacobi_cn);

    // Jacobian elliptic delta amplitude integrals.
    std::cout << "jacobi_dn\n" << std::flush;
    basename = "jacobi_dn";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi_dn(filename);
    maketest(jacobi_dn, beast::jacobi_dn,
	     "testcase_jacobi_dn", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101),
	     "Boost",
	     file_jacobi_dn);

    // Exponential integral En.
    // Skip the pole at 0.
    std::cout << "expint_en\n" << std::flush;
    basename = "expint";
    filename = get_filename(path, prefix, basename, "_en", ".cc");
    std::ofstream file_expint_en(filename);
    maketest(expint, gsl::expint,
	     "testcase_expint_en", "__gnu_cxx", basename,
	     "n", {0, 1, 2, 3, 5},
	     "x", fill_argument(std::make_pair(Real{0}, Real{50}),
				std::make_pair(false, true), 51),
	     "GSL",
	     file_expint_en);

    // Fresnel cosine integral.
    std::cout << "fresnel_c\n" << std::flush;
    basename = "fresnel_c";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_fresnel_c(filename);
    maketest(fresnel_c, gsl::fresnel_c,
	     "testcase_fresnel_c", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(false, true), 401),
	     "GSL",
	     file_fresnel_c);

    // Fresnel sine integral.
    std::cout << "fresnel_s\n" << std::flush;
    basename = "fresnel_s";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_fresnel_s(filename);
    maketest(fresnel_s, gsl::fresnel_s,
	     "testcase_fresnel_s", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(false, true), 401),
	     "GSL",
	     file_fresnel_s);

    // Sinus cardinal function.
    std::cout << "sinc\n" << std::flush;
    basename = "sinc";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinc(filename);
    maketest(sinc, gsl::sinc,
	     "testcase_sinc", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(true, true), 401),
	     "GSL",
	     file_sinc);

    // Reperiodized sinus cardinal function.
    std::cout << "sinc_pi\n" << std::flush;
    basename = "sinc_pi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinc_pi(filename);
    maketest(sinc_pi, gsl::sinc_pi,
	     "testcase_sinc_pi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(true, true), 401),
	     "GSL",
	     file_sinc_pi);

    // Log rising factorials.
    std::cout << "lrising_factorial\n" << std::flush;
    basename = "lrising_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lrising_factorial(filename);
    maketest(lrising_factorial, beast::lrising_factorial,
	     "testcase_lrising_factorial", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}},
	     "Boost",
	     file_lrising_factorial, true, true);

    // Log falling factorials.
    std::cout << "lfalling_factorial\n" << std::flush;
    basename = "lfalling_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lfalling_factorial(filename);
    maketest(lfalling_factorial, beast::lfalling_factorial,
	     "testcase_lfalling_factorial", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {Real{0}, Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}},
	     "Boost",
	     file_lfalling_factorial, true, true);

    // Rising factorials.
    std::cout << "rising_factorial\n" << std::flush;
    basename = "rising_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_rising_factorial(filename);
    maketest(rising_factorial, beast::rising_factorial,
	     "testcase_rising_factorial", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", dvorder,
	     "Boost",
	     file_rising_factorial, true, true);

    // Falling factorials.
    std::cout << "falling_factorial\n" << std::flush;
    basename = "falling_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_falling_factorial(filename);
    maketest(falling_factorial, beast::falling_factorial,
	     "testcase_falling_factorial", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {Real{0}, Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}},
	     "Boost",
	     file_falling_factorial, true, true);

    // Regular modified spherical bessel functions.
    std::cout << "sph_bessel_i\n" << std::flush;
    basename = "sph_bessel_i";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel_i(filename);
    test =
    maketest(sph_bessel_i, gsl::sph_bessel_i,
	     "testcase_sph_bessel_i", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_bessel_i, true, false);
    maketest(sph_bessel_i, gsl::sph_bessel_i,
	     "testcase_sph_bessel_i", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_bessel_i, false, true, test);

    // Irregular modified spherical bessel functions.
    std::cout << "sph_bessel_k\n" << std::flush;
    basename = "sph_bessel_k";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel_k(filename);
    test =
    maketest(sph_bessel_k, gsl::sph_bessel_k,
	     "testcase_sph_bessel_k", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_sph_bessel_k, true, false);
    maketest(sph_bessel_k, gsl::sph_bessel_k,
	     "testcase_sph_bessel_k", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_sph_bessel_k, false, true, test);

    // Legendre functions of the second kind.
    std::cout << "legendre_q\n" << std::flush;
    basename = "legendre_q";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_legendre_q(filename);
    maketest(legendre_q, gsl::legendre_q,
	     "testcase_legendre_q", "__gnu_cxx", basename,
	     "l", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "GSL",
	     file_legendre_q);

    // Factorial.
    std::cout << "factorial\n" << std::flush;
    basename = "factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_factorial(filename);
    maketest(factorial<Real>, gsl::factorial,
	     "testcase_factorial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "GSL",
	     file_factorial);

    // Log factorial.
    std::cout << "lfactorial\n" << std::flush;
    basename = "lfactorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lfactorial(filename);
    maketest(lfactorial<Real>, gsl::lfactorial,
	     "testcase_lfactorial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 500U),
				std::make_pair(true, true), 501),
	     "GSL",
	     file_lfactorial);

    // Double factorial.
    std::cout << "double_factorial\n" << std::flush;
    basename = "double_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_double_factorial(filename);
    maketest(double_factorial<Real>, gsl::double_factorial,
	     "testcase_double_factorial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0, 50),
				std::make_pair(true, true), 51),
	     "GSL",
	     file_double_factorial);

    // Log double factorial.
    std::cout << "ldouble_factorial\n" << std::flush;
    basename = "ldouble_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_ldouble_factorial(filename);
    maketest(ldouble_factorial<Real>, gsl::ldouble_factorial,
	     "testcase_ldouble_factorial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0, 500),
				std::make_pair(true, true), 501),
	     "GSL",
	     file_ldouble_factorial);

    // Binomial coefficient.
    std::cout << "binomial\n" << std::flush;
    basename = "binomial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_binomial(filename);
    maketest(binomial<Real>, gsl::binomial,
	     "testcase_binomial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "k", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "GSL",
	     file_binomial);

    // Log binomial coefficient.
    std::cout << "lbinomial\n" << std::flush;
    basename = "lbinomial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lbinomial(filename);
    maketest(lbinomial<Real>, gsl::lbinomial,
	     "testcase_lbinomial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 200U),
				std::make_pair(true, true), 201),
	     "k", fill_argument(std::make_pair(0U, 200U),
				std::make_pair(true, true), 201),
	     "GSL",
	     file_lbinomial);

    // Gegenbauer polynomials.
    std::cout << "gegenbauer\n" << std::flush;
    basename = "gegenbauer";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_gegenbauer(filename);
    maketest(gegenbauer, gsl::gegenbauer,
	     "testcase_gegenbauer", "__gnu_cxx", basename,
	     "n", vorder,
	     "alpha", fill_argument(std::make_pair(Real{0}, Real{5}),
				    std::make_pair(true, true), 11),
             "x", fill_argument(std::make_pair(Real{0}, Real{20}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_gegenbauer);

    // Chebyshev polynomials of the first kind.
    std::cout << "chebyshev_t\n" << std::flush;
    basename = "chebyshev_t";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_chebyshev_t(filename);
    maketest(chebyshev_t, __chebyshev_t,
	     "testcase_chebyshev_t", "__gnu_cxx", basename,
	     "n", {0U, 1U, 5U, 8U, 10U, 20U, 40U, 100U},
	     "x", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_chebyshev_t);

    // Chebyshev polynomials of the second kind.
    std::cout << "chebyshev_u\n" << std::flush;
    basename = "chebyshev_u";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_chebyshev_u(filename);
    maketest(chebyshev_u, __chebyshev_u,
	     "testcase_chebyshev_u", "__gnu_cxx", basename,
	     "n", {0U, 1U, 5U, 8U, 10U, 20U, 40U, 100U},
	     "x", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_chebyshev_u);

    // Chebyshev polynomials of the third kind.
    std::cout << "chebyshev_v\n" << std::flush;
    basename = "chebyshev_v";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_chebyshev_v(filename);
    maketest(chebyshev_v, __chebyshev_v,
	     "testcase_chebyshev_v", "__gnu_cxx", basename,
	     "n", {0U, 1U, 5U, 8U, 10U, 20U, 40U, 100U},
	     "x", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_chebyshev_v);

    // Chebyshev polynomials of the fourth kind.
    std::cout << "chebyshev_w\n" << std::flush;
    basename = "chebyshev_w";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_chebyshev_w(filename);
    maketest(chebyshev_w, __chebyshev_w,
	     "testcase_chebyshev_w", "__gnu_cxx", basename,
	     "n", {0U, 1U, 5U, 8U, 10U, 20U, 40U, 100U},
	     "x", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_chebyshev_w);

    // Jacobi polynomials.
    std::cout << "jacobi\n" << std::flush;
    basename = "jacobi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi(filename);
    maketest(jacobi, gsl::jacobi,
	     "testcase_jacobi", "__gnu_cxx", basename,
	     "n", vorder,
	     "alpha", fill_argument(std::make_pair(Real{0}, Real{5}),
				    std::make_pair(true, true), 11),
             "beta", fill_argument(std::make_pair(Real{0}, Real{5}),
				   std::make_pair(true, true), 21),
             "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				   std::make_pair(true, true), 41),
	     "GSL",
	     file_jacobi);

    // Radial polynomials.
    std::cout << "radpoly\n" << std::flush;
    basename = "radpoly";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_radpoly(filename);
    maketest(radpoly, gsl::radpoly,
	     "testcase_radpoly", "__gnu_cxx", basename,
	     "n", vorder, "m", vorder,
             "rho", fill_argument(std::make_pair(Real{0}, Real{1}),
				  std::make_pair(true, true), 21),
	     "GSL",
	     file_radpoly);

    // Zernike polynomials.
    std::cout << "zernike\n" << std::flush;
    basename = "zernike";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_zernike(filename);
    maketest(zernike, gsl::zernike,
	     "testcase_zernike", "__gnu_cxx", basename,
	     "n", vorder, "m", iorder,
             "rho", fill_argument(std::make_pair(Real{0}, Real{1}),
				  std::make_pair(true, true), 21),
             "phi", vphid,
	     "GSL",
	     file_zernike);

    // Confluent hypergeometric limit functions.
    // Skip the singularity at c = 0.
    std::cout << "conf_hyperg_lim\n" << std::flush;
    basename = "conf_hyperg_lim";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_conf_hyperg_lim(filename);
    maketest(conf_hyperg_lim, gsl::conf_hyperg_lim,
	     "testcase_conf_hyperg_lim", "__gnu_cxx", basename,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_conf_hyperg_lim);

    // Heuman lambda functions.
    // Avoid poles at |x| = 1.
    std::cout << "heuman_lambda\n" << std::flush;
    basename = "heuman_lambda";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_heuman_lambda(filename);
    maketest(heuman_lambda, beast::heuman_lambda,
	     "testcase_heuman_lambda", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vnopolesd,
	     "Boost",
	     file_heuman_lambda);

    // Elliptic D integrals.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_d\n" << std::flush;
    basename = "ellint_d";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_d(filename);
    maketest(ellint_d, beast::ellint_d,
	     "testcase_ellint_d", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     "Boost",
	     file_ellint_d);

    // Complementary elliptic D integrals.
    // Avoid poles at |x| = 1.
    std::cout << "comp_ellint_d\n" << std::flush;
    basename = "comp_ellint_d";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_comp_ellint_d(filename);
    maketest(comp_ellint_d, beast::comp_ellint_d,
	     "testcase_comp_ellint_d", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "Boost",
	     file_comp_ellint_d);

    // Jacobi zeta functions.
    // Avoid poles at |x| = 1.
    std::cout << "jacobi_zeta\n" << std::flush;
    basename = "jacobi_zeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_jacobi_zeta(filename);
    maketest(jacobi_zeta, beast::jacobi_zeta,
	     "testcase_jacobi_zeta", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     "Boost",
	     file_jacobi_zeta);

    // Cylindrical Hankel functions of the first kind.
    std::cout << "cyl_hankel_1\n" << std::flush;
    basename = "cyl_hankel_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_hankel_1(filename);
    test =
    maketest(cyl_hankel_1, beast::cyl_hankel_1,
	     "testcase_cyl_hankel_1", "__gnu_cxx", basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_1, true, false);
    test =
    maketest(cyl_hankel_1, beast::cyl_hankel_1,
	     "testcase_cyl_hankel_1", "__gnu_cxx", basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_1, false, false, test);
    test =
    maketest(cyl_hankel_1, beast::cyl_hankel_1,
	     "testcase_cyl_hankel_1", "__gnu_cxx", basename,
	     "nu", {Real{-5}, Real{-2}, Real{0}, Real{2}, Real{5}},
	     "x", {Real{-10}, Real{-5}, Real{-2}, Real{-0.25}},
	     "Boost",
	     file_cyl_hankel_1, false, false, test);
    maketest(cyl_hankel_1, beast::cyl_hankel_1,
	     "testcase_cyl_hankel_1", "__gnu_cxx", basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_1, false, true, test);

    // Cylindrical Hankel functions of the second kind.
    std::cout << "cyl_hankel_2\n" << std::flush;
    basename = "cyl_hankel_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_hankel_2(filename);
    test =
    maketest(cyl_hankel_2, beast::cyl_hankel_2,
	     "testcase_cyl_hankel_2", "__gnu_cxx", basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_2, true, false);
    test =
    maketest(cyl_hankel_2, beast::cyl_hankel_2,
	     "testcase_cyl_hankel_2", "__gnu_cxx", basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_2, false, false, test);
    test =
    maketest(cyl_hankel_2, beast::cyl_hankel_2,
	     "testcase_cyl_hankel_2", "__gnu_cxx", basename,
	     "nu", {Real{-5}, Real{-2}, Real{0}, Real{2}, Real{5}},
	     "x", {Real{-10}, Real{-5}, Real{-2}, Real{-0.25}},
	     "Boost",
	     file_cyl_hankel_2, false, false, test);
    maketest(cyl_hankel_2, beast::cyl_hankel_2,
	     "testcase_cyl_hankel_2", "__gnu_cxx", basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_2, false, true, test);

    // Spherical Hankel functions of the first kind.
    std::cout << "sph_hankel_1\n" << std::flush;
    basename = "sph_hankel_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sph_hankel_1(filename);
    test =
    maketest(sph_hankel_1, beast::sph_hankel_1,
	     "testcase_sph_hankel_1", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_1, true, false);
    test =
    maketest(sph_hankel_1, beast::sph_hankel_1,
	     "testcase_sph_hankel_1", "__gnu_cxx", basename,
	     "nu", {0, 2, 5},
	     "x", {Real{-10}, Real{-5}, Real{-2}, Real{-0.25}},
	     "Boost",
	     file_sph_hankel_1, false, false, test);
    maketest(sph_hankel_1, beast::sph_hankel_1,
	     "testcase_sph_hankel_1", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_1, false, true, test);

    // Spherical Hankel functions of the second kind.
    std::cout << "sph_hankel_2\n" << std::flush;
    basename = "sph_hankel_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sph_hankel_2(filename);
    test =
    maketest(sph_hankel_2, beast::sph_hankel_2,
	     "testcase_sph_hankel_2", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_2, true, false);
    test =
    maketest(sph_hankel_2, beast::sph_hankel_2,
	     "testcase_sph_hankel_2", "__gnu_cxx", basename,
	     "nu", {0, 2, 5},
	     "x", {Real{-10}, Real{-5}, Real{-2}, Real{-0.25}},
	     "Boost",
	     file_sph_hankel_2, false, false, test);
    maketest(sph_hankel_2, beast::sph_hankel_2,
	     "testcase_sph_hankel_2", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_2, false, true, test);

    // Spherical harmonic functions.
    std::cout << "sph_harmonic\n" << std::flush;
    basename = "sph_harmonic";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_harmonic(filename);
    maketest(sph_harmonic, beast::sph_harmonic,
	     "testcase_sph_harmonic", "__gnu_cxx", basename,
	     "l", vorder, "m", iorder,
	     "theta", fill_argument(std::make_pair(Real{0}, _S_pi),
				    std::make_pair(true, true), 21),
	     "phi", vphid,
	     "Boost",
	     file_sph_harmonic);

    // Dirichlet eta function.
    std::cout << "dirichlet_eta\n" << std::flush;
    // Skip the pole at 1.
    basename = "dirichlet_eta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_dirichlet_eta(filename);
    test =
    maketest(dirichlet_eta, gsl::dirichlet_eta,
	     "testcase_dirichlet_eta", "__gnu_cxx", basename,
	     "s", fill_argument(std::make_pair(Real{-10}, Real{1}),
				std::make_pair(true, false), 56),
	     "GSL",
	     file_dirichlet_eta, true, false);
    maketest(dirichlet_eta, gsl::dirichlet_eta,
	     "testcase_dirichlet_eta", "__gnu_cxx", basename,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
				std::make_pair(false, true), 146),
	     "GSL",
	     file_dirichlet_eta, false, true, test);

    // Owens T functions.
    std::cout << "owens_t\n" << std::flush;
    basename = "owens_t";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_owens_t(filename);
    maketest(owens_t, beast::owens_t,
	     "testcase_owens_t", "__gnu_cxx", basename,
	     "h", fill_argument(std::make_pair(Real{-5}, Real{5}),
				std::make_pair(true, true), 41),
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_owens_t, true, true);

    // Clausen Cl function.
    std::cout << "clausen_cl\n" << std::flush;
    basename = "clausen_cl";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_clausen_cl(filename);
    maketest(clausen_cl<Real>, gsl::clausen_cl,
	     "testcase_clausen_cl", "__gnu_cxx", basename,
	     "m", fill_argument(std::make_pair(2U, 2U),
				std::make_pair(true, true), 1),
	     "w", fill_argument(std::make_pair(Real{-10}, Real{+10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_clausen_cl, true, true);

    // Bernoulli numbers.
    std::cout << "bernoulli\n" << std::flush;
    basename = "bernoulli";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_bernoulli(filename);
    maketest(bernoulli<Real>, beast::bernoulli,
	     "testcase_bernoulli", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 100U),
				std::make_pair(true, true), 101),
	     "Boost",
	     file_bernoulli);

    // Euler numbers.
    std::cout << "euler\n" << std::flush;
    basename = "euler";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_euler(filename);
    maketest(euler<Real>, burkhardt::euler,
	     "testcase_euler", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "Burkhardt",
	     file_euler);

    // Eulerian numbers of the first kind.
    std::cout << "eulerian_1\n" << std::flush;
    basename = "eulerian_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_eulerian_1(filename);
    maketest(eulerian_1<Real>, burkhardt::eulerian_1,
	     "testcase_eulerian_1", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 10U),
				std::make_pair(true, true), 11),
	     "m", fill_argument(std::make_pair(0U, 10U),
				std::make_pair(true, true), 11),
	     "Burkhardt",
	     file_eulerian_1);

    // Stirling numbers of the first kind.
    std::cout << "stirling_1\n" << std::flush;
    basename = "stirling_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_stirling_1(filename);
    maketest(stirling_1<Real>, burkhardt::stirling_1,
	     "testcase_stirling_1", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 10U),
				std::make_pair(true, true), 11),
	     "m", fill_argument(std::make_pair(0U, 10U),
				std::make_pair(true, true), 11),
	     "Burkhardt",
	     file_stirling_1);

    // Stirling numbers of the first kind.
    std::cout << "stirling_2\n" << std::flush;
    basename = "stirling_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_stirling_2(filename);
    maketest(stirling_2<Real>, burkhardt::stirling_2,
	     "testcase_stirling_2", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 10U),
				std::make_pair(true, true), 11),
	     "m", fill_argument(std::make_pair(0U, 10U),
				std::make_pair(true, true), 11),
	     "Burkhardt",
	     file_stirling_2);

    // Reperiodized sine function.
    std::cout << "sin_pi\n" << std::flush;
    basename = "sin_pi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sin_pi(filename);
    maketest(sin_pi, beast::sin_pi,
	     "testcase_sin_pi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+50}),
				std::make_pair(false, true), 701),
	     "Boost",
	     file_sin_pi);

    // Reperiodized cosine function.
    std::cout << "cos_pi\n" << std::flush;
    basename = "cos_pi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cos_pi(filename);
    maketest(cos_pi, beast::cos_pi,
	     "testcase_cos_pi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+50}),
				std::make_pair(false, true), 701),
	     "Boost",
	     file_cos_pi);

    // Debye functions.
    std::cout << "debye\n" << std::flush;
    basename = "debye";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_debye(filename);
    maketest(debye, gsl::debye,
	     "testcase_debye", "__gnu_cxx", basename,
	     "n", {1U, 2U, 3U, 4U, 5U, 6U},
	     "x", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_debye);

#endif // STD

  }


int
main()
{
  harness<double>();

  return 0;
}


