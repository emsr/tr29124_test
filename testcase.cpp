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

#define STD 1

#include "specfun_testcase.h"
#include "wrap_gsl.h"
#include "wrap_boost.h"
#include "wrap_burkhardt.h"

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


template<typename Real>
  void
  harness()
  {
#if STD
    std::string ns("std");
    using __gnu_cxx::airy_ai;
    using __gnu_cxx::airy_bi;
    using       std::assoc_laguerre;
    using       std::assoc_legendre;
    using __gnu_cxx::bernoulli;
    using       std::beta;
    using __gnu_cxx::bincoef;
    using __gnu_cxx::chebyshev_t;
    using __gnu_cxx::chebyshev_u;
    using __gnu_cxx::chebyshev_v;
    using __gnu_cxx::chebyshev_w;
    using __gnu_cxx::clausen;
    using __gnu_cxx::clausen_c;
    using __gnu_cxx::clausen_s;
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
    using       std::expint;
    using __gnu_cxx::expint;
    using __gnu_cxx::factorial;
    using __gnu_cxx::fresnel_c;
    using __gnu_cxx::fresnel_s;
    using __gnu_cxx::gegenbauer;
    using       std::hermite;
    using __gnu_cxx::heuman_lambda;
    using __gnu_cxx::hurwitz_zeta;
    using __gnu_cxx::hyperg;
    using __gnu_cxx::ibeta;
    using __gnu_cxx::jacobi;
    using __gnu_cxx::jacobi_sn;
    using __gnu_cxx::jacobi_cn;
    using __gnu_cxx::jacobi_dn;
    using __gnu_cxx::jacobi_zeta;
    using       std::laguerre;
    using __gnu_cxx::lbincoef;
    using __gnu_cxx::ldouble_factorial;
    using       std::legendre;
    using __gnu_cxx::legendre_q;
    using __gnu_cxx::lfactorial;
    using __gnu_cxx::lpochhammer_lower;
    using __gnu_cxx::lpochhammer;
    using __gnu_cxx::owens_t;
    using __gnu_cxx::pgamma;
    using __gnu_cxx::pochhammer_lower;
    using __gnu_cxx::pochhammer;
    using __gnu_cxx::psi;
    using __gnu_cxx::qgamma;
    using __gnu_cxx::radpoly;
    using       std::riemann_zeta;
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
    std::string ns("tr1");
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
    std::vector<Real> cyl_neg_order{-5, -2, -1, -Real{2.0L/3.0L},
				    -Real{0.5L}, -Real{1.0L/3.0L}};

    std::vector<Real> cyl_order{0, Real{1.0L/3.0L},
				Real{0.5L}, Real{2.0L/3.0L},
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

    std::vector<Real> vab{0, Real{0.5L}, 1, 2, 5, 10, 20};

    unsigned int test = 1;

    const std::string path = "check";
    const std::string prefix = "/check_";

    std::string nsname = "std";

    std::string basename;
    std::string filename;

#if STD
    // Airy functions of the first kind.
    std::cout << "airy_ai" << std::endl;
    basename = "airy_ai";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_airy_ai(filename.c_str());
    maketest(airy_ai, gsl::airy_ai,
	     "testcase_airy_ai", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_airy_ai);

    // Airy functions of the second kind.
    std::cout << "airy_bi" << std::endl;
    basename = "airy_bi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_airy_bi(filename.c_str());
    maketest(airy_bi, gsl::airy_bi,
	     "testcase_airy_bi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_airy_bi);
#endif // STD

    // Associated Laguerre polynomials.
    std::cout << "assoc_laguerre" << std::endl;
    basename = "assoc_laguerre";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_assoc_laguerre(filename.c_str());
    maketest(assoc_laguerre, gsl::assoc_laguerre,
	     "testcase_assoc_laguerre", nsname, basename,
	     "n", vorder, "m", vorder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_assoc_laguerre);

    // Associated Legendre functions.
    std::cout << "assoc_legendre" << std::endl;
    basename = "assoc_legendre";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_assoc_legendre(filename.c_str());
    maketest(assoc_legendre, gsl::assoc_legendre,
	     "testcase_assoc_legendre", nsname, basename,
	     "l", vorder, "m", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_assoc_legendre);

    // Beta functions.
    std::cout << "beta" << std::endl;
    basename = "beta";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_beta(filename.c_str());
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
    std::cout << "comp_ellint_1" << std::endl;
    basename = "comp_ellint_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_comp_ellint_1(filename.c_str());
    maketest(comp_ellint_1, gsl::comp_ellint_1,
	     "testcase_comp_ellint_1", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "GSL",
	     file_comp_ellint_1);

    // Complete elliptic integrals of the second kind.
    // Avoid poles at |x| = 1.
    std::cout << "comp_ellint_2" << std::endl;
    basename = "comp_ellint_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_comp_ellint_2(filename.c_str());
    maketest(comp_ellint_2, gsl::comp_ellint_2,
	     "testcase_comp_ellint_2", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "GSL",
	     file_comp_ellint_2);

    // Complete elliptic integrals of the third kind.
    // Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3" << std::endl;
    basename = "comp_ellint_3";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_comp_ellint_3(filename.c_str());
    maketest(comp_ellint_3, gsl::comp_ellint_3,
	     "testcase_comp_ellint_3", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, false), 11),
	     "GSL",
	     file_comp_ellint_3);

#if STD
    // Confluent hypergeometric functions.
    // Skip the singularity at c = 0.
    std::cout << "conf_hyperg" << std::endl;
    basename = "conf_hyperg";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_conf_hyperg(filename.c_str());
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
    std::cout << "conf_hyperg" << std::endl;
    basename = "conf_hyperg";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_conf_hyperg(filename.c_str());
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
    std::cout << "cyl_bessel_i" << std::endl;
    basename = "cyl_bessel_i";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_i(filename.c_str());
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
    std::cout << "cyl_bessel_j" << std::endl;
    basename = "cyl_bessel_j";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_j(filename.c_str());
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
    maketest(cyl_bessel_j, gsl::cyl_bessel_j,
	     "testcase_cyl_bessel_j", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_cyl_bessel_j, false, true, test);

    // Irregular modified cylindrical Bessel functions.
    // Skip the pole at the origin.
    std::cout << "cyl_bessel_k" << std::endl;
    basename = "cyl_bessel_k";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_bessel_k(filename.c_str());
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
    std::cout << "cyl_neumann" << std::endl;
    basename = "cyl_neumann";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_neumann(filename.c_str());
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
    maketest(cyl_neumann, gsl::cyl_neumann,
	     "testcase_cyl_neumann", nsname, basename,
	     "nu", cyl_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     "GSL",
	     file_cyl_neumann, false, true, test);

    // Elliptic integrals of the first kind.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_1" << std::endl;
    basename = "ellint_1";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_1(filename.c_str());
    maketest(ellint_1, gsl::ellint_1,
	     "testcase_ellint_1", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     "GSL",
	     file_ellint_1);

    // Elliptic integrals of the second kind.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_2" << std::endl;
    basename = "ellint_2";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_2(filename.c_str());
    maketest(ellint_2, gsl::ellint_2,
	     "testcase_ellint_2", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     "GSL",
	     file_ellint_2);

    // Elliptic integrals of the third kind.
    // Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3" << std::endl;
    basename = "ellint_3";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_ellint_3(filename.c_str());
    maketest(ellint_3, gsl::ellint_3,
	     "testcase_ellint_3", nsname, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, false), 11),
	     "phi", vphid,
	     "GSL",
	     file_ellint_3);

    // Exponential integral.
    // Skip the pole at 0.
    std::cout << "expint" << std::endl;
    basename = "expint";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_expint(filename.c_str());
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
    std::cout << "hermite" << std::endl;
    basename = "hermite";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hermite(filename.c_str());
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
	     "x", {4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 50.0, 100.0},
	     "GSL",
	     file_hermite, false, true, test);

#if STD
    // Hypergeometric functions.
    // Skip the singularity at c = 0.
    // Skip the singularity at x = -1.
    std::cout << "hyperg" << std::endl;
    basename = "hyperg";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hyperg(filename.c_str());
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
    std::cout << "hyperg" << std::endl;
    basename = "hyperg";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hyperg(filename.c_str());
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
    std::cout << "laguerre" << std::endl;
    basename = "laguerre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_laguerre(filename.c_str());
    maketest(laguerre, gsl::laguerre,
	     "testcase_laguerre", nsname, basename,
	     "n", vorder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_laguerre);

    // Legendre polynomials.
    std::cout << "legendre" << std::endl;
    basename = "legendre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_legendre(filename.c_str());
    maketest(legendre, gsl::legendre_p,
	     "testcase_legendre", nsname, basename,
	     "l", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_legendre);

    // Riemann zeta function.
    std::cout << "riemann_zeta" << std::endl;
    // Skip the pole at 1.
    basename = "riemann_zeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_riemann_zeta(filename.c_str());
    test =
    maketest(riemann_zeta, gsl::riemann_zeta,
	     "testcase_riemann_zeta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{-10}, Real{1}),
				std::make_pair(true, false), 56),
	     "GSL",
	     file_riemann_zeta, true, false);
    maketest(riemann_zeta, gsl::riemann_zeta,
	     "testcase_riemann_zeta", nsname, basename,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
				std::make_pair(false, true), 146),
	     "GSL",
	     file_riemann_zeta, false, true, test);

#if STD
    // Hurwitz zeta functions.
    std::cout << "hurwitz_zeta" << std::endl;
    // Skip the pole at 1.
    basename = "hurwitz_zeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_hurwitz_zeta(filename.c_str());
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
    std::cout << "sph_bessel" << std::endl;
    basename = "sph_bessel";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel(filename.c_str());
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
    std::cout << "sph_legendre" << std::endl;
    basename = "sph_legendre";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_legendre(filename.c_str());
    maketest(sph_legendre, gsl::legendre_sphPlm,
	     "testcase_sph_legendre", nsname, basename,
	     "l", vorder, "m", vorder,
	     "theta", fill_argument(std::make_pair(Real{0}, static_cast<Real>(M_PI)),
				    std::make_pair(true, true), 21),
	     "GSL",
	     file_sph_legendre);

    // Spherical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "sph_neumann" << std::endl;
    basename = "sph_neumann";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_neumann(filename.c_str());
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
    std::cout << "ellint_rc" << std::endl;
    basename = "ellint_rc";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rc(filename.c_str());
    maketest(ellint_rc, gsl::ellint_rc,
	     "testcase_ellint_rc", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "GSL",
	     file_ellint_rc);

    // Carlson elliptic functions R_D.
    std::cout << "ellint_rd" << std::endl;
    basename = "ellint_rd";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rd(filename.c_str());
    maketest(ellint_rd, gsl::ellint_rd,
	     "testcase_ellint_rd", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_ellint_rd);

    // Carlson elliptic functions R_F.
    std::cout << "ellint_rf" << std::endl;
    basename = "ellint_rf";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rf(filename.c_str());
    maketest(ellint_rf, gsl::ellint_rf,
	     "testcase_ellint_rf", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_ellint_rf);

    // Carlson elliptic functions R_G.
    std::cout << "ellint_rg" << std::endl;
    basename = "ellint_rg";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rg(filename.c_str());
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
    std::cout << "ellint_rj" << std::endl;
    basename = "ellint_rj";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_rj(filename.c_str());
    maketest(ellint_rj, gsl::ellint_rj,
	     "testcase_ellint_rj", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "p", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "GSL",
	     file_ellint_rj);

    // Dilogarithm functions.
    std::cout << "dilog" << std::endl;
    basename = "dilog";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_dilog(filename.c_str());
    maketest(dilog, gsl::dilog,
	     "testcase_dilog", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{1}),
				std::make_pair(true, true), 23),
	     "GSL",
	     file_dilog);

    // Log Gamma functions.
    std::cout << "lgamma" << std::endl;
    basename = "lgamma";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_lgamma(filename.c_str());
    maketest(lgamma, beast::lgamma,
	     "testcase_lgamma", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{-10}, Real{20}),
				std::make_pair(false, true), 151),
	     "Boost",
	     file_lgamma);

    // Upper incomplete Gamma functions.
    std::cout << "tgamma" << std::endl;
    basename = "tgamma";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_tgamma(filename.c_str());
    maketest(tgamma, beast::tgamma,
	     "testcase_tgamma", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     "Boost",
	     file_tgamma);

    // Lower incomplete Gamma functions.
    std::cout << "tgamma_lower" << std::endl;
    basename = "tgamma_lower";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_tgamma_lower(filename.c_str());
    maketest(tgamma_lower, beast::tgamma_lower,
	     "testcase_tgamma_lower", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     "Boost",
	     file_tgamma_lower);

    // Lower regularized incomplete Gamma functions.
    std::cout << "pgamma" << std::endl;
    basename = "pgamma";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_pgamma(filename.c_str());
    maketest(pgamma, gsl::pgamma,
	     "testcase_pgamma", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_pgamma);

    // Upper regularized incomplete Gamma functions.
    std::cout << "qgamma" << std::endl;
    basename = "qgamma";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_qgamma(filename.c_str());
    maketest(qgamma, gsl::qgamma,
	     "testcase_qgamma", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     "GSL",
	     file_qgamma);

    // Incomplete Beta functions.
    std::cout << "ibeta" << std::endl;
    basename = "ibeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ibeta(filename.c_str());
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

    // Digamma or psi functions.
    std::cout << "psi" << std::endl;
    basename = "psi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_psi(filename.c_str());
    test =
    maketest(psi, gsl::psi,
	     "testcase_psi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-9.9375}, +Real{10.0625}),
				std::make_pair(true, true), 801),
	     "GSL",
	     file_psi, true, false);
    maketest(psi, gsl::psi,
	     "testcase_psi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{1}, Real{100}),
				std::make_pair(true, true), 199),
	     "GSL",
	     file_psi, false, true, test);

    // Sine integral or Si functions.
    std::cout << "sinint" << std::endl;
    basename = "sinint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinint(filename.c_str());
    maketest(sinint, gsl::Si,
	     "testcase_sinint", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_sinint);

    // Cosine integral or Ci functions.
    std::cout << "cosint" << std::endl;
    basename = "cosint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cosint(filename.c_str());
    maketest(cosint, gsl::Ci,
	     "testcase_cosint", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_cosint);

    // Hyperbolic sine integral or Shi functions.
    std::cout << "sinhint" << std::endl;
    basename = "sinhint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinhint(filename.c_str());
    maketest(sinhint, gsl::Shi,
	     "testcase_sinhint", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_sinhint);

    // Hyperbolic cosine integral or Chi functions.
    std::cout << "coshint" << std::endl;
    basename = "coshint";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_coshint(filename.c_str());
    maketest(coshint, gsl::Chi,
	     "testcase_coshint", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
				std::make_pair(false, true), 101),
	     "GSL",
	     file_coshint);

    // Dawson integral.
    std::cout << "dawson" << std::endl;
    basename = "dawson";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_dawson(filename.c_str());
    maketest(dawson, gsl::dawson,
	     "testcase_dawson", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+20}),
				std::make_pair(false, true), 201),
	     "GSL",
	     file_dawson);

    // Jacobian elliptic sine amplitude integrals.
    std::cout << "jacobi_sn" << std::endl;
    basename = "jacobi_sn";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi_sn(filename.c_str());
    maketest(jacobi_sn, gsl::jacobi_sn,
	     "testcase_jacobi_sn", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101),
	     "GSL",
	     file_jacobi_sn);

    // Jacobian elliptic cosine amplitude integrals.
    std::cout << "jacobi_cn" << std::endl;
    basename = "jacobi_cn";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi_cn(filename.c_str());
    maketest(jacobi_cn, gsl::jacobi_cn,
	     "testcase_jacobi_cn", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101),
	     "GSL",
	     file_jacobi_cn);

    // Jacobian elliptic delta amplitude integrals.
    std::cout << "jacobi_dn" << std::endl;
    basename = "jacobi_dn";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi_dn(filename.c_str());
    maketest(jacobi_dn, gsl::jacobi_dn,
	     "testcase_jacobi_dn", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101),
	     "GSL",
	     file_jacobi_dn);

    // Exponential integral En.
    // Skip the pole at 0.
    std::cout << "expint_en" << std::endl;
    basename = "expint";
    filename = get_filename(path, prefix, basename, "_en", ".cc");
    std::ofstream file_expint_en(filename.c_str());
    maketest(expint, gsl::expint,
	     "testcase_expint_en", "__gnu_cxx", basename,
	     "n", {0, 1, 2, 3, 5},
	     "x", fill_argument(std::make_pair(Real{0}, Real{50}),
				std::make_pair(false, true), 51),
	     "GSL",
	     file_expint_en);

    // Fresnel cosine integral.
    std::cout << "fresnel_c" << std::endl;
    basename = "fresnel_c";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_fresnel_c(filename.c_str());
    maketest(fresnel_c, gsl::fresnel_c,
	     "testcase_fresnel_c", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(false, true), 401),
	     "GSL",
	     file_fresnel_c);

    // Fresnel sine integral.
    std::cout << "fresnel_s" << std::endl;
    basename = "fresnel_s";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_fresnel_s(filename.c_str());
    maketest(fresnel_s, gsl::fresnel_s,
	     "testcase_fresnel_s", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(false, true), 401),
	     "GSL",
	     file_fresnel_s);

    // Sinus cardinal function.
    std::cout << "sinc" << std::endl;
    basename = "sinc";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinc(filename.c_str());
    maketest(sinc, gsl::sinc,
	     "testcase_sinc", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(true, true), 401),
	     "GSL",
	     file_sinc);

    // Reperiodized sinus cardinal function.
    std::cout << "sinc_pi" << std::endl;
    basename = "sinc_pi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sinc_pi(filename.c_str());
    maketest(sinc_pi, gsl::sinc_pi,
	     "testcase_sinc_pi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
				std::make_pair(true, true), 401),
	     "GSL",
	     file_sinc_pi);

    // Log upper Pochhammer symbol.
    std::cout << "lpochhammer" << std::endl;
    basename = "lpochhammer";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lpochhammer(filename.c_str());
    maketest(lpochhammer, beast::lpochhammer,
	     "testcase_lpochhammer", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0},
	     "Boost",
	     file_lpochhammer, true, true);

    // Log lower Pochhammer symbol.
    std::cout << "lpochhammer_lower" << std::endl;
    basename = "lpochhammer_lower";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lpochhammer_lower(filename.c_str());
    maketest(lpochhammer_lower, beast::lpochhammer_lower,
	     "testcase_lpochhammer_lower", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0},
	     "Boost",
	     file_lpochhammer_lower, true, true);

    // Upper Pochhammer symbols (see boost::rising_factorial).
    std::cout << "pochhammer" << std::endl;
    basename = "pochhammer";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_pochhammer(filename.c_str());
    maketest(pochhammer, beast::pochhammer,
	     "testcase_pochhammer", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", dvorder,
	     "Boost",
	     file_pochhammer, true, true);

    // Lower Pochhammer symbols (see boost::falling_factorial).
    std::cout << "pochhammer_lower" << std::endl;
    basename = "pochhammer_lower";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_pochhammer_lower(filename.c_str());
    maketest(pochhammer_lower, beast::pochhammer_lower,
	     "testcase_pochhammer_lower", "__gnu_cxx", basename,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     "x", {0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0},
	     "Boost",
	     file_pochhammer_lower, true, true);

    // Regular modified spherical bessel functions.
    std::cout << "sph_bessel_i" << std::endl;
    basename = "sph_bessel_i";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel_i(filename.c_str());
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
    std::cout << "sph_bessel_k" << std::endl;
    basename = "sph_bessel_k";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_bessel_k(filename.c_str());
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
    std::cout << "legendre_q" << std::endl;
    basename = "legendre_q";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_legendre_q(filename.c_str());
    maketest(legendre_q, gsl::legendre_q,
	     "testcase_legendre_q", "__gnu_cxx", basename,
	     "l", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "GSL",
	     file_legendre_q);

    // Factorial.
    std::cout << "factorial" << std::endl;
    basename = "factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_factorial(filename.c_str());
    maketest(factorial<Real>, gsl::factorial,
	     "testcase_factorial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "GSL",
	     file_factorial);

    // Log factorial.
    std::cout << "lfactorial" << std::endl;
    basename = "lfactorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lfactorial(filename.c_str());
    maketest(lfactorial<Real>, gsl::lfactorial,
	     "testcase_lfactorial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 500U),
				std::make_pair(true, true), 501),
	     "GSL",
	     file_lfactorial);

    // Double factorial.
    std::cout << "double_factorial" << std::endl;
    basename = "double_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_double_factorial(filename.c_str());
    maketest(double_factorial<Real>, gsl::double_factorial,
	     "testcase_double_factorial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0, 50),
				std::make_pair(true, true), 51),
	     "GSL",
	     file_double_factorial);

    // Log double factorial.
    std::cout << "ldouble_factorial" << std::endl;
    basename = "ldouble_factorial";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_ldouble_factorial(filename.c_str());
    maketest(ldouble_factorial<Real>, gsl::ldouble_factorial,
	     "testcase_ldouble_factorial", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0, 500),
				std::make_pair(true, true), 501),
	     "GSL",
	     file_ldouble_factorial);

    // Binomial coefficient.
    std::cout << "bincoef" << std::endl;
    basename = "bincoef";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_bincoef(filename.c_str());
    maketest(bincoef<Real>, gsl::choose,
	     "testcase_bincoef", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "k", fill_argument(std::make_pair(0U, 50U),
				std::make_pair(true, true), 51),
	     "GSL",
	     file_bincoef);

    // Log binomial coefficient.
    std::cout << "lbincoef" << std::endl;
    basename = "lbincoef";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_lbincoef(filename.c_str());
    maketest(lbincoef<Real>, gsl::lnchoose,
	     "testcase_lbincoef", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 200U),
				std::make_pair(true, true), 201),
	     "k", fill_argument(std::make_pair(0U, 200U),
				std::make_pair(true, true), 201),
	     "GSL",
	     file_lbincoef);

    // Gegenbauer polynomials.
    std::cout << "gegenbauer" << std::endl;
    basename = "gegenbauer";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_gegenbauer(filename.c_str());
    maketest(gegenbauer, gsl::gegenpoly_n,
	     "testcase_gegenbauer", "__gnu_cxx", basename,
	     "n", vorder,
	     "alpha", fill_argument(std::make_pair(Real{0}, Real{5}),
				    std::make_pair(true, true), 11),
             "x", fill_argument(std::make_pair(Real{0}, Real{20}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_gegenbauer);

    // Chebyshev polynomials of the first kind.
    std::cout << "chebyshev_t - UNTESTED" << std::endl;
/*
    basename = "chebyshev_t";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_chebyshev_t(filename.c_str());
    maketest(chebyshev_t, gsl::chebyshev_t,
	     "testcase_chebyshev_t", "__gnu_cxx", basename,
	     "n", {0U, 1U, 5U, 8U, 10U, 20U, 40U, 100U},
	     "x", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_chebyshev_t);
*/
    // Chebyshev polynomials of the second kind.
    std::cout << "chebyshev_u - UNTESTED" << std::endl;

    // Chebyshev polynomials of the third kind.
    std::cout << "chebyshev_v - UNTESTED" << std::endl;

    // Chebyshev polynomials of the fourth kind.
    std::cout << "chebyshev_w - UNTESTED" << std::endl;

    // Jacobi polynomials.
    std::cout << "jacobi" << std::endl;
    basename = "jacobi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_jacobi(filename.c_str());
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
    std::cout << "radpoly" << std::endl;
    basename = "radpoly";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_radpoly(filename.c_str());
    maketest(radpoly, gsl::radpoly,
	     "testcase_radpoly", "__gnu_cxx", basename,
	     "n", vorder, "m", vorder,
             "rho", fill_argument(std::make_pair(Real{0}, Real{1}),
				  std::make_pair(true, true), 21),
	     "GSL",
	     file_radpoly);

    // Zernike polynomials.
    std::cout << "zernike" << std::endl;
    basename = "zernike";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_zernike(filename.c_str());
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
    std::cout << "conf_hyperg_lim" << std::endl;
    basename = "conf_hyperg_lim";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_conf_hyperg_lim(filename.c_str());
    maketest(conf_hyperg_lim, gsl::hyperg_0F1,
	     "testcase_conf_hyperg_lim", "__gnu_cxx", basename,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 21),
	     "GSL",
	     file_conf_hyperg_lim);

    // Heuman lambda functions.
    // Avoid poles at |x| = 1.
    std::cout << "heuman_lambda" << std::endl;
    basename = "heuman_lambda";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_heuman_lambda(filename.c_str());
    maketest(heuman_lambda, beast::heuman_lambda,
	     "testcase_heuman_lambda", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vnopolesd,
	     "Boost",
	     file_heuman_lambda);

    // Elliptic D integrals.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_d" << std::endl;
    basename = "ellint_d";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_ellint_d(filename.c_str());
    maketest(ellint_d, beast::ellint_d,
	     "testcase_ellint_d", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     "Boost",
	     file_ellint_d);

    // Complementary elliptic D integrals.
    // Avoid poles at |x| = 1.
    std::cout << "comp_ellint_d" << std::endl;
    basename = "comp_ellint_d";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_comp_ellint_d(filename.c_str());
    maketest(comp_ellint_d, beast::comp_ellint_d,
	     "testcase_comp_ellint_d", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "Boost",
	     file_comp_ellint_d);

    // Jacobi zeta functions.
    // Avoid poles at |x| = 1.
    std::cout << "jacobi_zeta" << std::endl;
    basename = "jacobi_zeta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_jacobi_zeta(filename.c_str());
    maketest(jacobi_zeta, beast::jacobi_zeta,
	     "testcase_jacobi_zeta", "__gnu_cxx", basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     "Boost",
	     file_jacobi_zeta);

    // Cylindrical Hankel functions of the first kind.
    std::cout << "cyl_hankel_1" << std::endl;
    basename = "cyl_hankel_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_hankel_1(filename.c_str());
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
    std::cout << "cyl_hankel_2" << std::endl;
    basename = "cyl_hankel_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cyl_hankel_2(filename.c_str());
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
    std::cout << "sph_hankel_1" << std::endl;
    basename = "sph_hankel_1";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sph_hankel_1(filename.c_str());
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
    std::cout << "sph_hankel_2" << std::endl;
    basename = "sph_hankel_2";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sph_hankel_2(filename.c_str());
    test =
    maketest(sph_hankel_2, beast::sph_hankel_2,
	     "testcase_sph_hankel_2", "__gnu_cxx", basename,
	     "n", sph_order,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_sph_hankel_2, true, false);
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
    std::cout << "sph_harmonic" << std::endl;
    basename = "sph_harmonic";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_sph_harmonic(filename.c_str());
    maketest(sph_harmonic, beast::sph_harmonic,
	     "testcase_sph_harmonic", "__gnu_cxx", basename,
	     "l", vorder, "m", iorder,
	     "theta", fill_argument(std::make_pair(Real{0}, static_cast<Real>(M_PI)),
				    std::make_pair(true, true), 21),
	     "phi", vphid,
	     "Boost",
	     file_sph_harmonic);

    // Dirichlet eta function.
    std::cout << "dirichlet_eta" << std::endl;
    // Skip the pole at 1.
    basename = "dirichlet_eta";
    filename = get_filename(path, prefix, basename, "", ".cc");
    std::ofstream file_dirichlet_eta(filename.c_str());
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
    std::cout << "owens_t" << std::endl;
    basename = "owens_t";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_owens_t(filename.c_str());
    maketest(owens_t, beast::owens_t,
	     "testcase_owens_t", "__gnu_cxx", basename,
	     "h", fill_argument(std::make_pair(Real{-5}, Real{5}),
				std::make_pair(true, true), 41),
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     "Boost",
	     file_owens_t, true, true);

    // Clausen C function.
    std::cout << "clausen_c" << std::endl;
    basename = "clausen_c";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_clausen_c(filename.c_str());
    maketest(clausen_c<Real>, gsl::clausen_c,
	     "testcase_clausen_c", "__gnu_cxx", basename,
	     "m", fill_argument(std::make_pair(2U, 2U),
				std::make_pair(true, true), 1),
	     "w", fill_argument(std::make_pair(Real{-10}, Real{+10}),
				std::make_pair(true, true), 41),
	     "GSL",
	     file_clausen_c, true, true);

    // Bernoulli numbers.
    std::cout << "bernoulli" << std::endl;
    basename = "bernoulli";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_bernoulli(filename.c_str());
    maketest(bernoulli<Real>, beast::bernoulli,
	     "testcase_bernoulli", "__gnu_cxx", basename,
	     "n", fill_argument(std::make_pair(0U, 100U),
				std::make_pair(true, true), 101),
	     "Boost",
	     file_bernoulli);

    // Reperiodized sine function.
    std::cout << "sin_pi" << std::endl;
    basename = "sin_pi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_sin_pi(filename.c_str());
    maketest(sin_pi, beast::sin_pi,
	     "testcase_sin_pi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+50}),
				std::make_pair(false, true), 701),
	     "Boost",
	     file_sin_pi);

    // Reperiodized cosine function.
    std::cout << "cos_pi" << std::endl;
    basename = "cos_pi";
    filename = get_filename(path, prefix, basename, "",  ".cc");
    std::ofstream file_cos_pi(filename.c_str());
    maketest(cos_pi, beast::cos_pi,
	     "testcase_cos_pi", "__gnu_cxx", basename,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+50}),
				std::make_pair(false, true), 701),
	     "Boost",
	     file_cos_pi);

#endif // STD

  }


int
main()
{
  harness<double>();

  return 0;
}


