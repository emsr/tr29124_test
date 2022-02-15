
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <functional>
#include <utility>
#include <tuple>
#include <filesystem>

#include <emsr/special_functions.h>

#include <specfun_testcase.h>
#include <wrap_gsl.h>
#include <wrap_boost.h>
#include <wrap_burkhardt.h>

#include <testcase.tcc>

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
template<typename Tp>
  Tp
  chebyshev_t_trig(unsigned int n, Tp x)
  {
    auto theta = std::acos(x);
    return std::cos(n * theta);
  }

/**
 * Return the Chebyshev polynomial of the second kind by trigonometric identity:
 * @f[
 *   U_n(x) = \frac{\sin((n + 1)\theta)}{\sin(\theta)}
 * @f]
 * where @f$ \theta = \acos(x) @f$.
 */
template<typename Tp>
  Tp
  chebyshev_u_trig(unsigned int n, Tp x)
  {
    const auto _S_eps = emsr::epsilon(x);
    if (std::abs(x + Tp{1}) < _S_eps)
      return (n % 2 == 0 ? +1 : -1) * Tp(n + 1);
    else if (std::abs(x - Tp{1}) < _S_eps)
      return Tp(n + 1);
    else
      {
	auto theta = std::acos(x);
	return std::sin(Tp(n + 1) * theta)
	     / std::sin(theta);
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
template<typename Tp>
  Tp
  chebyshev_v_trig(unsigned int n, Tp x)
  {
    const auto _S_eps = emsr::epsilon(x);
    if (std::abs(x + Tp{1}) < _S_eps)
      return (n % 2 == 0 ? +1 : -1) * Tp(2 * n + 1);
    else
      {
	auto theta = std::acos(x);
	return std::cos(Tp(n + Tp{1} / Tp{2}) * theta)
	     / std::cos(theta / Tp{2});
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
template<typename Tp>
  Tp
  chebyshev_w_trig(unsigned int n, Tp x)
  {
    const auto _S_eps = emsr::epsilon(x);
    if (std::abs(x - Tp{1}) < _S_eps)
      return Tp(2 * n + 1);
    else
      {
	auto theta = std::acos(x);
	return std::sin(Tp(n + Tp{1} / Tp{2}) * theta)
	     / std::sin(theta / Tp{2});
      }
  }
/*
double
chebyshev_t(unsigned int n, double x)
{ return chebyshev_t_trig<double>(n, x); }

double
chebyshev_u(unsigned int n, double x)
{ return chebyshev_u_trig<double>(n, x); }

double
chebyshev_v(unsigned int n, double x)
{ return chebyshev_v_trig<double>(n, x); }

double
chebyshev_w(unsigned int n, double x)
{ return chebyshev_w_trig<double>(n, x); }
*/
template<typename Real>
  void
  harness()
  {
    namespace fsys = std::filesystem;

    const auto _S_pi = emsr::pi_v<Real>;

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
      vphid.push_back(Real{10} * i * emsr::rad_v<Real>);
    vphid.push_back(emsr::pi_v<Real> / Real{2});

    std::vector<Real> vnopolesd;
    for (unsigned int i = 1; i < num_phi - 1; ++i)
      vnopolesd.push_back(Real{10} * i * _S_pi / Real{180});

    std::vector<Real> vab{0, Real{1}/Real{2}, 1, 2, 5, 10, 20};

    unsigned int test = 1;

    const fsys::path path = "test/src/";
    fsys::create_directories(path);

    std::string nsname = "emsr";
    const std::string prefix = "/check_";

    std::string basename;
    std::string filename;

    // Airy functions of the first kind.
    std::cout << "airy_ai\n" << std::flush;
    basename = "airy_ai";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_airy_ai(filename);
    maketest(emsr::airy_ai, gsl::airy_ai,
	     "testcase_airy_ai", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_airy_ai);

    // Scaled Airy functions of the first kind.
    std::cout << "airy_ai_scaled\n" << std::flush;
    basename = "airy_ai_scaled";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_airy_ai_scaled(filename);
    maketest(emsr::airy_ai_scaled, gsl::airy_ai_scaled,
	     "testcase_airy_ai_scaled", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_airy_ai_scaled);

    // Airy functions of the second kind.
    std::cout << "airy_bi\n" << std::flush;
    basename = "airy_bi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_airy_bi(filename);
    maketest(emsr::airy_bi, gsl::airy_bi,
	     "testcase_airy_bi", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_airy_bi);

    // Scaled Airy functions of the second kind.
    std::cout << "airy_bi_scaled\n" << std::flush;
    basename = "airy_bi_scaled";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_airy_bi_scaled(filename);
    maketest(emsr::airy_bi_scaled, gsl::airy_bi_scaled,
	     "testcase_airy_bi_scaled", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_airy_bi_scaled);

    // Associated Laguerre polynomials.
    std::cout << "assoc_laguerre\n" << std::flush;
    basename = "assoc_laguerre";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_assoc_laguerre(filename);
    maketest(emsr::assoc_laguerre, gsl::assoc_laguerre,
	     "testcase_assoc_laguerre", nsname, basename,
	     "n", vorder, "alpha", vorder,
	     "x", fill_argument(std::make_pair(Real{-50}, Real{100}),
				std::make_pair(true, true), 16),
	     "GSL",
	     file_assoc_laguerre);

    // Associated Legendre functions.
    std::cout << "assoc_legendre\n" << std::flush;
    basename = "assoc_legendre";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_assoc_legendre(filename);
    maketest(emsr::assoc_legendre, beast::assoc_legendre,
	     "testcase_assoc_legendre", nsname, basename,
	     "l", vorder, "m", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_assoc_legendre);

    // Beta functions.
    std::cout << "beta\n" << std::flush;
    basename = "beta";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_beta(filename);
    maketest(emsr::beta, gsl::beta,
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
    maketest(emsr::comp_ellint_1, beast::comp_ellint_1,
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
    maketest(emsr::comp_ellint_2, beast::comp_ellint_2,
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
    maketest(emsr::comp_ellint_3, beast::comp_ellint_3,
	     "testcase_comp_ellint_3", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, false), 11),
	     "Boost",
	     file_comp_ellint_3);

    // Confluent hypergeometric functions.
    // Skip the singularity at c = 0.
    std::cout << "conf_hyperg\n" << std::flush;
    basename = "conf_hyperg";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_conf_hyperg(filename);
    maketest(emsr::conf_hyperg, gsl::conf_hyperg,
	     "testcase_conf_hyperg",
	     "emsr",
	      basename, "a", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_conf_hyperg);

    // Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i\n" << std::flush;
    basename = "cyl_bessel_i";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_i(filename);
    test =
    maketest(emsr::cyl_bessel_i, beast::cyl_bessel_i,
	     "testcase_cyl_bessel_i", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_bessel_i, false);
    test =
    maketest(emsr::cyl_bessel_i, gsl::cyl_bessel_i,
	     "testcase_cyl_bessel_i", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_cyl_bessel_i, false, test);
    maketest(emsr::cyl_bessel_i, gsl::cyl_bessel_i,
	     "testcase_cyl_bessel_i", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_cyl_bessel_i, true, test);

    // Scaled regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i_scaled\n" << std::flush;
    basename = "cyl_bessel_i_scaled";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_i_scaled(filename);
    test =
    maketest(emsr::cyl_bessel_i_scaled, gsl::cyl_bessel_i_scaled,
	     "testcase_cyl_bessel_i_scaled", nsname, basename,
	     "nu", {100},
	     "x", fill_argument(std::make_pair(Real{1000}, Real{2000}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_cyl_bessel_i_scaled);

    // Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j\n" << std::flush;
    basename = "cyl_bessel_j";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_j(filename);
    test =
    maketest(emsr::cyl_bessel_j, beast::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_bessel_j, false);
    test =
    maketest(emsr::cyl_bessel_j, beast::cyl_bessel_j,//gsl::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",//"GSL",
	     file_cyl_bessel_j, false, test);
    test =
    maketest(emsr::cyl_bessel_j, beast::cyl_bessel_j,//gsl::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",//"GSL",
	     file_cyl_bessel_j, false, test);
    maketest(emsr::cyl_bessel_j, beast::cyl_bessel_j,//gsl::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", {100},
	     "x", fill_argument(std::make_pair(Real{1000}, Real{2000}),
				std::make_pair(true, true), 11),
	     "Boost",//"GSL",
	     file_cyl_bessel_j, true, test);

    // Irregular modified cylindrical Bessel functions.
    // Skip the pole at the origin.
    std::cout << "cyl_bessel_k\n" << std::flush;
    basename = "cyl_bessel_k";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_k(filename);
    test =
    maketest(emsr::cyl_bessel_k, beast::cyl_bessel_k,
	     "testcase_cyl_bessel_k", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_bessel_k, false);
    test =
    maketest(emsr::cyl_bessel_k, gsl::cyl_bessel_k,
	     "testcase_cyl_bessel_k", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_bessel_k, false, test);
    maketest(emsr::cyl_bessel_k, gsl::cyl_bessel_k,
	     "testcase_cyl_bessel_k", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_bessel_k, true, test);

    // Scaled irregular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_k_scaled\n" << std::flush;
    basename = "cyl_bessel_k_scaled";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_k_scaled(filename);
    test =
    maketest(emsr::cyl_bessel_k_scaled, gsl::cyl_bessel_k_scaled,
	     "testcase_cyl_bessel_k_scaled", nsname, basename,
	     "nu", {100},
	     "x", fill_argument(std::make_pair(Real{1000}, Real{2000}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_cyl_bessel_k_scaled);

    // Cylindrical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "cyl_neumann\n" << std::flush;
    basename = "cyl_neumann";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_neumann(filename);
    test =
    maketest(emsr::cyl_neumann, beast::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_neumann, false);
    test =
    maketest(emsr::cyl_neumann, gsl::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_neumann, false, test);
    test =
    maketest(emsr::cyl_neumann, gsl::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_neumann, false, test);
    maketest(emsr::cyl_neumann, gsl::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", {100},
	     "x", fill_argument(std::make_pair(Real{1000}, Real{2000}),
				std::make_pair(false, true), 11),
	     "GSL",
	     file_cyl_neumann, true, test);

    // Elliptic integrals of the first kind.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_1\n" << std::flush;
    basename = "ellint_1";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_1(filename);
    maketest(emsr::ellint_1, beast::ellint_1,
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
    maketest(emsr::ellint_2, beast::ellint_2,
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
    maketest(emsr::ellint_3, beast::ellint_3,
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
    maketest(emsr::expint, gsl::expint,
	     "testcase_expint", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-50}, Real{0}),
				std::make_pair(true, false), 51),
	     "GSL",
	     file_expint, false);
    maketest(emsr::expint, gsl::expint,
	     "testcase_expint", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{50}),
				std::make_pair(false, true), 51),
	     "GSL",
	     file_expint, true, test);

    // Hermite polynomials
    std::cout << "hermite\n" << std::flush;
    basename = "hermite";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hermite(filename);
    test =
    maketest(emsr::hermite, gsl::hermite,
	     "testcase_hermite", nsname, basename,
	     "n", vorder,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 201),
	     "GSL",
	     file_hermite, false);
    maketest(emsr::hermite, gsl::hermite,
	     "testcase_hermite", nsname, basename,
	     "n", {8, 18, 32, 50, 72, 128, 200, 1250, 5000},
	     "x", {Real{4}, Real{6}, Real{8}, Real{10}, Real{12},
		   Real{16}, Real{20}, Real{50}, Real{100}},
	     "GSL",
	     file_hermite, true, test);

    // Propabilist Hermite polynomials
    std::cout << "hermite_he\n" << std::flush;
    basename = "hermite_he";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hermite_he(filename);
    test =
    maketest(emsr::hermite_he, gsl::hermite_he,
	     "testcase_hermite_he", nsname, basename,
	     "n", vorder,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 201),
	     "GSL",
	     file_hermite_he, false);
    maketest(emsr::hermite_he, gsl::hermite_he,
	     "testcase_hermite_he", nsname, basename,
	     "n", {8, 18, 32, 50, 72, 128, 200, 1250, 5000},
	     "x", {Real{4}, Real{6}, Real{8}, Real{10}, Real{12},
		   Real{16}, Real{20}, Real{50}, Real{100}},
	     "GSL",
	     file_hermite_he, true, test);

    // Hypergeometric functions.
    // Skip the singularity at c = 0.
    // Skip the singularity at x = -1.
    std::cout << "hyperg\n" << std::flush;
    basename = "hyperg";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hyperg(filename);
    maketest(emsr::hyperg, gsl::hyperg, "testcase_hyperg",
	     "emsr",
	     basename, "a", vab,
	     "b", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 6),
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, false), 21),
	     "GSL",
	     file_hyperg);

    // Laguerre polynomials.
    std::cout << "laguerre\n" << std::flush;
    basename = "laguerre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_laguerre(filename);
    maketest(emsr::laguerre, gsl::laguerre,
	     "testcase_laguerre", nsname, basename,
	     "n", vorder,
	     "x", fill_argument(std::make_pair(Real{-50}, Real{100}),
				std::make_pair(true, true), 31),
	     "GSL",
	     file_laguerre);

    // Legendre polynomials.
    std::cout << "legendre\n" << std::flush;
    basename = "legendre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_legendre(filename);
    maketest(emsr::legendre, beast::legendre_p,
	     "testcase_legendre", nsname, basename,
	     "l", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_legendre);

    // Riemann zeta function.
    std::cout << "riemann_zeta\n" << std::flush;
    // Skip the pole at 1.
    basename = "riemann_zeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_riemann_zeta(filename);
    test =
    maketest(emsr::riemann_zeta, beast::riemann_zeta,
	     "testcase_riemann_zeta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{-10}, Real{1}),
				std::make_pair(true, false), 56),
	     "Boost",
	     file_riemann_zeta, false);
    maketest(emsr::riemann_zeta, beast::riemann_zeta,
	     "testcase_riemann_zeta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
				std::make_pair(false, true), 146),
	     "Boost",
	     file_riemann_zeta, true, test);

    // Hurwitz zeta functions.
    std::cout << "hurwitz_zeta\n" << std::flush;
    // Skip the pole at 1.
    basename = "hurwitz_zeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hurwitz_zeta(filename);
    maketest(emsr::hurwitz_zeta, gsl::hurwitz_zeta,
	     "testcase_hurwitz_zeta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
				std::make_pair(false, true), 146),
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 26),
	     "GSL",
	     file_hurwitz_zeta);

    // Spherical Bessel functions.
    std::cout << "sph_bessel\n" << std::flush;
    basename = "sph_bessel";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel(filename);
    test =
    maketest(emsr::sph_bessel, gsl::sph_bessel,
	     "testcase_sph_bessel", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_bessel, false);
    maketest(emsr::sph_bessel, gsl::sph_bessel,
	     "testcase_sph_bessel", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_bessel, true, test);

    // Spherical Legendre functions.
    std::cout << "sph_legendre\n" << std::flush;
    basename = "sph_legendre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_legendre(filename);
    maketest(emsr::sph_legendre, beast::sph_legendre,
	     "testcase_sph_legendre", nsname, basename,
	     "l", vorder, "m", vorder,
	     "theta", fill_argument(std::make_pair(Real{0}, _S_pi),
				    std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_legendre);

    // Spherical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "sph_neumann\n" << std::flush;
    basename = "sph_neumann";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_neumann(filename);
    test =
    maketest(emsr::sph_neumann, gsl::sph_neumann,
	     "testcase_sph_neumann", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_sph_neumann, false);
    maketest(emsr::sph_neumann, gsl::sph_neumann,
	     "testcase_sph_neumann", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_sph_neumann, true, test);

    // Carlson elliptic functions R_C.
    std::cout << "ellint_rc\n" << std::flush;
    basename = "ellint_rc";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rc(filename);
    maketest(emsr::ellint_rc, beast::ellint_rc,
	     "testcase_ellint_rc", nsname, basename,
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
    maketest(emsr::ellint_rd, beast::ellint_rd,
	     "testcase_ellint_rd", nsname, basename,
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
    maketest(emsr::ellint_rf, beast::ellint_rf,
	     "testcase_ellint_rf", nsname, basename,
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
    maketest(emsr::ellint_rg, beast::ellint_rg,
	     "testcase_ellint_rg", nsname, basename,
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
    maketest(emsr::ellint_rj, beast::ellint_rj,
	     "testcase_ellint_rj", nsname, basename,
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
    maketest(emsr::dilog, gsl::dilog,
	     "testcase_dilog", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{1}),
				std::make_pair(true, true), 45),
	     "GSL",
	     file_dilog);

    // Log Gamma functions.
    std::cout << "lgamma\n" << std::flush;
    basename = "lgamma";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_lgamma(filename);
    maketest(emsr::lgamma, beast::lgamma,
	     "testcase_lgamma", nsname, basename,
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
    maketest(emsr::tgamma, beast::tgamma,
	     "testcase_tgamma", nsname, basename,
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
    maketest(emsr::tgamma, beast::tgamma,
	     "testcase_tgamma", nsname, basename,
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
    maketest(emsr::tgamma_lower, beast::tgamma_lower,
	     "testcase_tgamma_lower", nsname, basename,
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
    maketest(emsr::gamma_p, gsl::gamma_p,
	     "testcase_gamma_p", nsname, basename,
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
    maketest(emsr::gamma_q, gsl::gamma_q,
	     "testcase_gamma_q", nsname, basename,
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
    maketest<Real, Real, Real, Real>(emsr::ibeta, gsl::ibeta,
	     "testcase_ibeta", nsname, basename,
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
    maketest<Real, Real, Real, Real>(emsr::ibetac, beast::ibetac,
	     "testcase_ibetac", nsname, basename,
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
    maketest(emsr::digamma, beast::digamma,
	     "testcase_digamma", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-10} + skip, Real{+10} + skip),
				std::make_pair(true, true), 201),
	     "Boost",
	     file_digamma, false);
    maketest(emsr::digamma, beast::digamma,
	     "testcase_digamma", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{1}, Real{100}),
				std::make_pair(true, true), 199),
	     "Boost",
	     file_digamma, true, test);

    // Psi or digamma functions.
    std::cout << "polygamma\n" << std::flush;
    basename = "polygamma";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_polygamma(filename);
    test =
    maketest(emsr::polygamma, beast::polygamma,
	     "testcase_polygamma", nsname, basename,
	     "m", {0U, 1U, 2U, 3U, 5U, 10U},
	     "x", fill_argument(std::make_pair(Real{-10} + skip, Real{+10} + skip),
				std::make_pair(true, true), 201),
	     "Boost",
	     file_polygamma, false);
    maketest(emsr::polygamma, beast::polygamma,
	     "testcase_polygamma", nsname, basename,
	     "m", {0U, 1U, 2U, 3U, 5U, 10U},
	     "x", fill_argument(std::make_pair(Real{1}, Real{100}),
				std::make_pair(true, true), 199),
	     "Boost",
	     file_polygamma, true, test);

    // Sine integral or Si functions.
    std::cout << "sinint\n" << std::flush;
    basename = "sinint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinint(filename);
    maketest(emsr::sinint, gsl::sinint,
	     "testcase_sinint", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_sinint);

    // Cosine integral or Ci functions.
    std::cout << "cosint\n" << std::flush;
    basename = "cosint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cosint(filename);
    maketest(emsr::cosint, gsl::cosint,
	     "testcase_cosint", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_cosint);

    // Hyperbolic sine integral or Shi functions.
    std::cout << "sinhint\n" << std::flush;
    basename = "sinhint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinhint(filename);
    maketest(emsr::sinhint, gsl::sinhint,
	     "testcase_sinhint", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_sinhint);

    // Hyperbolic cosine integral or Chi functions.
    std::cout << "coshint\n" << std::flush;
    basename = "coshint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_coshint(filename);
    maketest(emsr::coshint, gsl::coshint,
	     "testcase_coshint", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_coshint);

    // Dawson integral.
    std::cout << "dawson\n" << std::flush;
    basename = "dawson";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_dawson(filename);
    maketest(emsr::dawson, gsl::dawson,
	     "testcase_dawson", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+20}),
				std::make_pair(false, true), 201),
	     "GSL",
	     file_dawson);

    // Jacobian elliptic sine amplitude integrals.
    std::cout << "jacobi_sn\n" << std::flush;
    basename = "jacobi_sn";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi_sn(filename);
    maketest(emsr::jacobi_sn, beast::jacobi_sn,
	     "testcase_jacobi_sn", nsname, basename,
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
    maketest(emsr::jacobi_cn, beast::jacobi_cn,
	     "testcase_jacobi_cn", nsname, basename,
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
    maketest(emsr::jacobi_dn, beast::jacobi_dn,
	     "testcase_jacobi_dn", nsname, basename,
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
    maketest(emsr::expint, gsl::expint,
	     "testcase_expint_en", nsname, basename,
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
    maketest(emsr::fresnel_c, gsl::fresnel_c,
	     "testcase_fresnel_c", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(false, true), 401),
	     "GSL",
	     file_fresnel_c);

    // Fresnel sine integral.
    std::cout << "fresnel_s\n" << std::flush;
    basename = "fresnel_s";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_fresnel_s(filename);
    maketest(emsr::fresnel_s, gsl::fresnel_s,
	     "testcase_fresnel_s", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(false, true), 401),
	     "GSL",
	     file_fresnel_s);

    // Sinus cardinal function.
    std::cout << "sinc\n" << std::flush;
    basename = "sinc";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinc(filename);
    maketest(emsr::sinc, gsl::sinc,
	     "testcase_sinc", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(true, true), 401),
	     "GSL",
	     file_sinc);

    // Reperiodized sinus cardinal function.
    std::cout << "sinc_pi\n" << std::flush;
    basename = "sinc_pi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinc_pi(filename);
    maketest(emsr::sinc_pi, gsl::sinc_pi,
	     "testcase_sinc_pi", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(true, true), 401),
	     "GSL",
	     file_sinc_pi);

    // Log rising factorials.
    std::cout << "lrising_factorial\n" << std::flush;
    basename = "lrising_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lrising_factorial(filename);
    maketest(emsr::lrising_factorial, beast::lrising_factorial,
	     "testcase_lrising_factorial", nsname, basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}},
	     "Boost",
	     file_lrising_factorial, true);

    // Log falling factorials.
    std::cout << "lfalling_factorial\n" << std::flush;
    basename = "lfalling_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lfalling_factorial(filename);
    maketest(emsr::lfalling_factorial, beast::lfalling_factorial,
	     "testcase_lfalling_factorial", nsname, basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {Real{0}, Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}},
	     "Boost",
	     file_lfalling_factorial, true);

    // Rising factorials.
    std::cout << "rising_factorial\n" << std::flush;
    basename = "rising_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_rising_factorial(filename);
    maketest(emsr::rising_factorial, beast::rising_factorial,
	     "testcase_rising_factorial", nsname, basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", dvorder,
	     "Boost",
	     file_rising_factorial, true);

    // Falling factorials.
    std::cout << "falling_factorial\n" << std::flush;
    basename = "falling_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_falling_factorial(filename);
    maketest(emsr::falling_factorial, beast::falling_factorial,
	     "testcase_falling_factorial", nsname, basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {Real{0}, Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}},
	     "Boost",
	     file_falling_factorial, true);

    // Regular modified spherical bessel functions.
    std::cout << "sph_bessel_i\n" << std::flush;
    basename = "sph_bessel_i";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel_i(filename);
    test =
    maketest(emsr::sph_bessel_i, gsl::sph_bessel_i,
	     "testcase_sph_bessel_i", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_bessel_i, false);
    maketest(emsr::sph_bessel_i, gsl::sph_bessel_i,
	     "testcase_sph_bessel_i", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_bessel_i, true, test);

    // Irregular modified spherical bessel functions.
    std::cout << "sph_bessel_k\n" << std::flush;
    basename = "sph_bessel_k";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel_k(filename);
    test =
    maketest(emsr::sph_bessel_k, gsl::sph_bessel_k,
	     "testcase_sph_bessel_k", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_sph_bessel_k, false);
    maketest(emsr::sph_bessel_k, gsl::sph_bessel_k,
	     "testcase_sph_bessel_k", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_sph_bessel_k, true, test);

    // Legendre functions of the second kind.
    std::cout << "legendre_q\n" << std::flush;
    basename = "legendre_q";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_legendre_q(filename);
    maketest(emsr::legendre_q, gsl::legendre_q,
	     "testcase_legendre_q", nsname, basename,
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
    maketest(emsr::factorial<Real>, gsl::factorial,
	     "testcase_factorial", nsname, basename,
	     "n", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "GSL",
	     file_factorial);

    // Log factorial.
    std::cout << "lfactorial\n" << std::flush;
    basename = "lfactorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lfactorial(filename);
    maketest(emsr::lfactorial<Real>, gsl::lfactorial,
	     "testcase_lfactorial", nsname, basename,
	     "n", fill_argument(std::make_pair(0U, 500U),
				std::make_pair(true, true), 501),
	     "GSL",
	     file_lfactorial);

    // Double factorial.
    std::cout << "double_factorial\n" << std::flush;
    basename = "double_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_double_factorial(filename);
    maketest(emsr::double_factorial<Real>, gsl::double_factorial,
	     "testcase_double_factorial", nsname, basename,
	     "n", fill_argument(std::make_pair(0, 50),
				std::make_pair(true, true), 51),
	     "GSL",
	     file_double_factorial);

    // Log double factorial.
    std::cout << "ldouble_factorial\n" << std::flush;
    basename = "ldouble_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_ldouble_factorial(filename);
    maketest(emsr::ldouble_factorial<Real>, gsl::ldouble_factorial,
	     "testcase_ldouble_factorial", nsname, basename,
	     "n", fill_argument(std::make_pair(0, 500),
				std::make_pair(true, true), 501),
	     "GSL",
	     file_ldouble_factorial);

    // Binomial coefficient.
    std::cout << "binomial\n" << std::flush;
    basename = "binomial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_binomial(filename);
    maketest(emsr::binomial<Real>, gsl::binomial,
	     "testcase_binomial", nsname, basename,
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
    maketest(emsr::lbinomial<Real>, gsl::lbinomial,
	     "testcase_lbinomial", nsname, basename,
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
    maketest(emsr::gegenbauer, gsl::gegenbauer,
	     "testcase_gegenbauer", nsname, basename,
	     "n", vorder,
	     "lambda", fill_argument(std::make_pair(Real{0}, Real{5}),
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
    maketest(emsr::chebyshev_t, chebyshev_t_trig<double>,
	     "testcase_chebyshev_t", nsname, basename,
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
    maketest(emsr::chebyshev_u, chebyshev_u_trig<double>,
	     "testcase_chebyshev_u", nsname, basename,
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
    maketest(emsr::chebyshev_v, chebyshev_v_trig<double>,
	     "testcase_chebyshev_v", nsname, basename,
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
    maketest(emsr::chebyshev_w, chebyshev_w_trig<double>,
	     "testcase_chebyshev_w", nsname, basename,
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
    maketest(emsr::jacobi, gsl::jacobi,
	     "testcase_jacobi", nsname, basename,
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
    maketest(emsr::radpoly, gsl::radpoly,
	     "testcase_radpoly", nsname, basename,
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
    maketest(emsr::zernike, gsl::zernike,
	     "testcase_zernike", nsname, basename,
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
    maketest(emsr::conf_hyperg_lim, gsl::conf_hyperg_lim,
	     "testcase_conf_hyperg_lim", nsname, basename,
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
    maketest(emsr::heuman_lambda, beast::heuman_lambda,
	     "testcase_heuman_lambda", nsname, basename,
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
    maketest(emsr::ellint_d, beast::ellint_d,
	     "testcase_ellint_d", nsname, basename,
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
    maketest(emsr::comp_ellint_d, beast::comp_ellint_d,
	     "testcase_comp_ellint_d", nsname, basename,
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
    maketest(emsr::jacobi_zeta, beast::jacobi_zeta,
	     "testcase_jacobi_zeta", nsname, basename,
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
    maketest(emsr::cyl_hankel_1, beast::cyl_hankel_1,
	     "testcase_cyl_hankel_1", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_1, false);
    test =
    maketest(emsr::cyl_hankel_1, beast::cyl_hankel_1,
	     "testcase_cyl_hankel_1", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_1, false, test);
    test =
    maketest(emsr::cyl_hankel_1, beast::cyl_hankel_1,
	     "testcase_cyl_hankel_1", nsname, basename,
	     "nu", {Real{-5}, Real{-2}, Real{0}, Real{2}, Real{5}},
	     "x", {Real{-10}, Real{-5}, Real{-2}, Real{-0.25}},
	     "Boost",
	     file_cyl_hankel_1, false, test);
    maketest(emsr::cyl_hankel_1, beast::cyl_hankel_1,
	     "testcase_cyl_hankel_1", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_1, true, test);

    // Cylindrical Hankel functions of the second kind.
    std::cout << "cyl_hankel_2\n" << std::flush;
    basename = "cyl_hankel_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_hankel_2(filename);
    test =
    maketest(emsr::cyl_hankel_2, beast::cyl_hankel_2,
	     "testcase_cyl_hankel_2", nsname, basename,
	     "nu", cyl_neg_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_2, false);
    test =
    maketest(emsr::cyl_hankel_2, beast::cyl_hankel_2,
	     "testcase_cyl_hankel_2", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_2, false, test);
    test =
    maketest(emsr::cyl_hankel_2, beast::cyl_hankel_2,
	     "testcase_cyl_hankel_2", nsname, basename,
	     "nu", {Real{-5}, Real{-2}, Real{0}, Real{2}, Real{5}},
	     "x", {Real{-10}, Real{-5}, Real{-2}, Real{-0.25}},
	     "Boost",
	     file_cyl_hankel_2, false, test);
    maketest(emsr::cyl_hankel_2, beast::cyl_hankel_2,
	     "testcase_cyl_hankel_2", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_cyl_hankel_2, true, test);

    // Spherical Hankel functions of the first kind.
    std::cout << "sph_hankel_1\n" << std::flush;
    basename = "sph_hankel_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sph_hankel_1(filename);
    test =
    maketest(emsr::sph_hankel_1, beast::sph_hankel_1,
	     "testcase_sph_hankel_1", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_1, false);
    test =
    maketest(emsr::sph_hankel_1, beast::sph_hankel_1,
	     "testcase_sph_hankel_1", nsname, basename,
	     "nu", {0, 2, 5},
	     "x", {Real{-10}, Real{-5}, Real{-2}, Real{-0.25}},
	     "Boost",
	     file_sph_hankel_1, false, test);
    maketest(emsr::sph_hankel_1, beast::sph_hankel_1,
	     "testcase_sph_hankel_1", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_1, true, test);

    // Spherical Hankel functions of the second kind.
    std::cout << "sph_hankel_2\n" << std::flush;
    basename = "sph_hankel_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sph_hankel_2(filename);
    test =
    maketest(emsr::sph_hankel_2, beast::sph_hankel_2,
	     "testcase_sph_hankel_2", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_2, false);
    test =
    maketest(emsr::sph_hankel_2, beast::sph_hankel_2,
	     "testcase_sph_hankel_2", nsname, basename,
	     "nu", {0, 2, 5},
	     "x", {Real{-10}, Real{-5}, Real{-2}, Real{-0.25}},
	     "Boost",
	     file_sph_hankel_2, false, test);
    maketest(emsr::sph_hankel_2, beast::sph_hankel_2,
	     "testcase_sph_hankel_2", nsname, basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_2, true, test);

    // Spherical harmonic functions.
    std::cout << "sph_harmonic\n" << std::flush;
    basename = "sph_harmonic";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_harmonic(filename);
    maketest(emsr::sph_harmonic, beast::sph_harmonic,
	     "testcase_sph_harmonic", nsname, basename,
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
    maketest(emsr::dirichlet_eta, gsl::dirichlet_eta,
	     "testcase_dirichlet_eta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{-10}, Real{1}),
				std::make_pair(true, false), 56),
	     "GSL",
	     file_dirichlet_eta, false);
    maketest(emsr::dirichlet_eta, gsl::dirichlet_eta,
	     "testcase_dirichlet_eta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
				std::make_pair(false, true), 146),
	     "GSL",
	     file_dirichlet_eta, true, test);

    // Owens T functions.
    std::cout << "owens_t\n" << std::flush;
    basename = "owens_t";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_owens_t(filename);
    maketest(emsr::owens_t, beast::owens_t,
	     "testcase_owens_t", nsname, basename,
	     "h", fill_argument(std::make_pair(Real{-5}, Real{5}),
				std::make_pair(true, true), 41),
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_owens_t, true);

    // Coulomb F functions.
    basename = "coulomb_f";
    std::cout << basename << '\n' << std::flush;
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_coulomb_f(filename);
    maketest(emsr::coulomb_f, gsl::coulomb_f,
             "testcase_coulomb_f", nsname, basename,
             "lambda", std::vector<Real>({Real{0}, Real{0.5}, Real{1}}),
             "eta", std::vector<Real>({Real{-2}, Real{0}, Real{2}, Real{10}}),
             "rho", fill_argument(std::make_pair(Real{0}, Real{20}),
				  std::make_pair(false, true), 41),
	     "GSL",
             file_coulomb_f, true);

    // Coulomb G functions.
    basename = "coulomb_g";
    std::cout << basename << '\n' << std::flush;
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_coulomb_g(filename);
    maketest(emsr::coulomb_g, gsl::coulomb_g,
             "testcase_coulomb_g", nsname, basename,
             "lambda", std::vector<Real>({Real{0}, Real{0.5}, Real{1}}),
             "eta", std::vector<Real>({Real{-2}, Real{0}, Real{2}, Real{10}}),
             "rho", fill_argument(std::make_pair(Real{0}, Real{20}),
				  std::make_pair(false, true), 41),
	     "GSL",
             file_coulomb_g, true);

    // Clausen Sl function?
    // Seems GSL only has Cl
    /*
    std::cout << "clausen_sl\n" << std::flush;
    basename = "clausen_sl";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_clausen_sl(filename);
    maketest(emsr::clausen_sl<Real>, gsl::clausen_sl,
	     "testcase_clausen_sl", nsname, basename,
	     "m", fill_argument(std::make_pair(2U, 2U),
				std::make_pair(true, true), 1),
	     "w", fill_argument(std::make_pair(Real{-10}, Real{+10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_clausen_sl, true);
    */

    // Clausen Cl function.
    std::cout << "clausen_cl\n" << std::flush;
    basename = "clausen_cl";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_clausen_cl(filename);
    maketest(emsr::clausen_cl<Real>, gsl::clausen_cl,
	     "testcase_clausen_cl", nsname, basename,
	     "m", fill_argument(std::make_pair(2U, 2U),
				std::make_pair(true, true), 1),
	     "w", fill_argument(std::make_pair(Real{-10}, Real{+10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_clausen_cl, true);

    // Bernoulli numbers.
    std::cout << "bernoulli\n" << std::flush;
    basename = "bernoulli";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_bernoulli(filename);
    maketest(emsr::bernoulli<Real>, beast::bernoulli,
	     "testcase_bernoulli", nsname, basename,
	     "n", fill_argument(std::make_pair(0U, 100U),
				std::make_pair(true, true), 101),
	     "Boost",
	     file_bernoulli);

    // Euler numbers.
    std::cout << "euler\n" << std::flush;
    basename = "euler";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_euler(filename);
    maketest(emsr::euler<Real>, burkhardt::euler,
	     "testcase_euler", nsname, basename,
	     "n", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "Burkhardt",
	     file_euler);

    // Eulerian numbers of the first kind.
    std::cout << "eulerian_1\n" << std::flush;
    basename = "eulerian_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_eulerian_1(filename);
    maketest(emsr::eulerian_1<Real>, burkhardt::eulerian_1,
	     "testcase_eulerian_1", nsname, basename,
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
    maketest(emsr::stirling_1<Real>, burkhardt::stirling_1,
	     "testcase_stirling_1", nsname, basename,
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
    maketest(emsr::stirling_2<Real>, burkhardt::stirling_2,
	     "testcase_stirling_2", nsname, basename,
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
    maketest(emsr::sin_pi, beast::sin_pi,
	     "testcase_sin_pi", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+50}),
				std::make_pair(false, true), 701),
	     "Boost",
	     file_sin_pi);

    // Reperiodized cosine function.
    std::cout << "cos_pi\n" << std::flush;
    basename = "cos_pi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cos_pi(filename);
    maketest(emsr::cos_pi, beast::cos_pi,
	     "testcase_cos_pi", nsname, basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+50}),
				std::make_pair(false, true), 701),
	     "Boost",
	     file_cos_pi);

    // Debye functions.
    std::cout << "debye\n" << std::flush;
    basename = "debye";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_debye(filename);
    maketest(emsr::debye, gsl::debye,
	     "testcase_debye", nsname, basename,
	     "n", {1U, 2U, 3U, 4U, 5U, 6U},
	     "x", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_debye);
  }


int
main()
{
  harness<double>();

  return 0;
}


