/*
/home/ESmith-rowland/bin/bin/g++ -std=gnu++17 -fconcepts -g -I. -o testcase2 -I/usr/local/include -I/usr/local/include testcase2.cpp wrap_gsl.cpp wrap_boost.cpp lerchphi/Source/lerchphi.cpp /home/ESmith-rowland/tr29124_test/gslextras/Fresnel/fresnel.c -L/usr/local/lib -lgsl -lgslcblas -ljacobi
*/
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <filesystem>
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

#include "testcase2.tcc"

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

    using fun_t = Real (*)(Real);
    using fun2_t = Real (*)(Real, Real);
    using fun3_t = Real (*)(Real, Real, Real);
    using fun4_t = Real (*)(Real, Real, Real, Real);
    using funir_t = Real (*)(unsigned int, Real);
    using fun2ir_t = Real (*)(unsigned int, unsigned int, Real);

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
    Real phi[num_phi];
    for (unsigned int i = 0; i < num_phi; ++i)
      phi[i] = Real{10} * i * __gnu_cxx::__math_constants<Real>::__pi / Real{180};
    std::vector<Real> vphid(phi, phi + num_phi);

    std::vector<Real> vab{0, Real{0.5L}, 1, 2, 5, 10, 20};

    const std::string path = "check";
    const std::string prefix = "/check_";
    auto outfile = [path, prefix](const char * fun)
		   -> std::string
		   { return std::string(path + prefix + fun + ".cc"); };

#if STD
    // Airy functions of the first kind.
    auto test_airy_ai = 
    make_testcase2(make_test_function("__gnu_cxx::airy_ai", (fun_t)&airy_ai),
		   make_baseline_function("GSL", "gsl::airy_ai", (fun_t)&gsl::airy_ai),
		   [](Real){ return true; },
		   "testcase_airy",
		   make_argument("x", fill2(Real{-10}, Real{10}, 41)));
    std::ofstream file_airy_ai(outfile("airy_ai"));
    test_airy_ai(file_airy_ai);

    // Airy functions of the second kind.
    auto test_airy_bi = 
    make_testcase2(make_test_function("__gnu_cxx::airy_bi", (fun_t)&airy_bi),
		   make_baseline_function("GSL", "gsl::airy_bi", (fun_t)&gsl::airy_bi),
		   [](Real){ return true; },
		   "testcase_airy",
		   make_argument("x", fill2(Real{-10}, Real{10}, 41)));
    std::ofstream file_airy_bi(outfile("airy_bi"));
    test_airy_bi(file_airy_bi);
#endif // STD

    // Associated Laguerre polynomials.
    auto test_assoc_laguerre = 
    make_testcase2(make_test_function("std::assoc_laguerre", (fun2ir_t)&assoc_laguerre),
		   make_baseline_function("GSL", "gsl::assoc_laguerre", (fun2ir_t)&gsl::assoc_laguerre),
		   [](unsigned, unsigned, Real){ return true; },
		   "testcase_assoc_laguerre",
		   make_argument("n", vorder), make_argument("m", vorder),
		   make_argument("x", fill2(Real{-10}, Real{10}, 41)));
    std::ofstream file_assoc_laguerre(outfile("assoc_laguerre"));
    test_assoc_laguerre(file_assoc_laguerre);

    // Beta functions.
    auto test_beta = 
    make_testcase2(make_test_function("std::beta", (fun2_t)&beta),
		   make_baseline_function("GSL", "gsl::beta", (fun2_t)&gsl::beta),
		   [](Real){ return true; },
		   "testcase_beta",
		   make_argument("x", fill2(Real{0}, Real{100}, 11)),
		   make_argument("y", fill2(Real{0}, Real{100}, 11)));
    std::ofstream file_beta(outfile("beta"));
    test_beta(file_beta);

    // Complete elliptic integrals of the first kind.
    // Avoid poles at |x| = 1.
    auto test_comp_ellint_1 = 
    make_testcase2(make_test_function("std::comp_ellint_1", (fun_t)&comp_ellint_1),
		   make_baseline_function("GSL", "gsl::comp_ellint_1", (fun_t)&gsl::comp_ellint_1),
		   [](Real x){ return (std::abs(x) == 1 ? false : true); },
		   make_argument("k", fill2(Real{-1}, Real{1}, 21)));
    std::ofstream file_comp_ellint_1(outfile("comp_ellint_1"));
    test_comp_ellint_1(file_comp_ellint_1);

    // Complete elliptic integrals of the second kind.
    // Avoid poles at |x| = 1.
    auto test_comp_ellint_2 = 
    make_testcase2(make_test_function("std::comp_ellint_2", (fun_t)&comp_ellint_2),
		   make_baseline_function("GSL", "gsl::comp_ellint_2", (fun_t)&gsl::comp_ellint_2),
		   [](Real x){ return (std::abs(x) == 1 ? false : true); },
		   make_argument("k", fill2(Real{-1}, Real{1}, 21)));
    std::ofstream file_comp_ellint_2(outfile("comp_ellint_2"));
    test_comp_ellint_2(file_comp_ellint_2);

    // Complete elliptic integrals of the third kind.
    // Avoid poles at |x| = 1 and at nu = 1.
    auto test_comp_ellint_3 = 
    make_testcase2(make_test_function("std::comp_ellint_3", (fun2_t)&comp_ellint_3),
		   make_baseline_function("GSL", "gsl::comp_ellint_3", (fun2_t)&gsl::comp_ellint_3),
		   [](Real x){ return (std::abs(x) == 1 ? false : true); },
		   make_argument("k", fill2(Real{-1}, Real{1}, 21)),
		   make_argument("nu", fill2(Real{0}, Real{1}, 11)));
    std::ofstream file_comp_ellint_3(outfile("comp_ellint_3"));
    test_comp_ellint_3(file_comp_ellint_3);

#if STD
    // Confluent hypergeometric functions.
    // Skip the singularity at c = 0.
    auto test_conf_hyperg = 
    make_testcase2(make_test_function("__gnu_cxx::conf_hyperg", (fun3_t)&conf_hyperg),
		   make_baseline_function("GSL", "gsl::conf_hyperg", (fun3_t)&gsl::conf_hyperg),
		   [](Real, Real, Real c, Real){ return (c == 0 ? false : true); },
		   make_argument("a", vab),
		   make_argument("c", fill2(Real{0}, Real{10}, 11)),
		   make_argument("x", fill2(Real{-10}, Real{10}, 21)));
    std::ofstream file_conf_hyperg(outfile("conf_hyperg"));
    test_conf_hyperg(file_conf_hyperg);
#else
    // Confluent hypergeometric functions.
    // Skip the singularity at c = 0.
    auto test_conf_hyperg = 
    make_testcase2(make_test_function("std::conf_hyperg", (fun3_t)&conf_hyperg),
		   make_baseline_function("GSL", "gsl::conf_hyperg", (fun3_t)&gsl::conf_hyperg),
		   [](Real, Real, Real c, Real){ return (c == 0 ? false : true); },
		   make_argument("a", vab),
		   make_argument("c", fill2(Real{0}, Real{10}, 11)),
		   make_argument("x", fill2(Real{-10}, Real{10}, 21)));
    std::ofstream file_conf_hyperg(outfile("conf_hyperg"));
    test_conf_hyperg(file_conf_hyperg);
#endif // STD

    // Regular modified cylindrical Bessel functions.
    auto test_cyl_bessel_i = 
    make_testcase2(make_test_function("std::cyl_bessel_i", (fun2_t)&cyl_bessel_i),
		   make_baseline_function("Boost", "beast::cyl_bessel_i", (fun2_t)&beast::cyl_bessel_i),
		   [](Real, Real){ return true; },
		   make_argument("nu", cyl_neg_order),
		   make_argument("x", fill2(Real{0}, Real{5}, 21)));
    std::ofstream file_cyl_bessel_i(outfile("cyl_bessel_i"));
    test_cyl_bessel_i(file_cyl_bessel_i);

    // Cylindrical Bessel functions (of the first kind).
    auto test_cyl_bessel_j = 
    make_testcase2(make_test_function("std::cyl_bessel_j", (fun2_t)&cyl_bessel_j),
		   make_baseline_function("Boost", "beast::cyl_bessel_j", (fun2_t)&beast::cyl_bessel_j),
		   [](Real, Real){ return true; },
		   make_argument("nu", cyl_neg_order),
		   make_argument("x", fill2(Real{0}, Real{5}, 21)));
    std::ofstream file_cyl_bessel_j(outfile("cyl_bessel_j"));
    test_cyl_bessel_j(file_cyl_bessel_j);

    // Irregular modified cylindrical Bessel functions.
    // Skip the pole at the origin.
    auto test_cyl_bessel_k = 
    make_testcase2(make_test_function("std::cyl_bessel_k", (fun2_t)&cyl_bessel_k),
		   make_baseline_function("Boost", "beast::cyl_bessel_k", (fun2_t)&beast::cyl_bessel_k),
		   [](Real, Real x){ return (x == 0 ? false : true); },
		   make_argument("nu", cyl_neg_order),
		   make_argument("x", fill2(Real{0}, Real{5}, 21)));
    std::ofstream file_cyl_bessel_k(outfile("cyl_bessel_k"));
    test_cyl_bessel_k(file_cyl_bessel_k);

    // Cylindrical Neumann functions.
    // Skip the pole at the origin.
    auto test_cyl_neumann = 
    make_testcase2(make_test_function("std::cyl_neumann", (fun2_t)&cyl_neumann),
		   make_baseline_function("Boost", "beast::cyl_neumann", (fun2_t)&beast::cyl_neumann),
		   [](Real, Real x){ return (x == 0 ? false : true); },
		   make_argument("nu", cyl_neg_order),
		   make_argument("x", fill2(Real{0}, Real{5}, 21)));
    std::ofstream file_cyl_neumann(outfile("cyl_neumann"));
    test_cyl_neumann(file_cyl_neumann);

    // Elliptic integrals of the first kind.
    // Avoid poles at |x| = 1.
    auto test_ellint_1 = 
    make_testcase2(make_test_function("std::ellint_1", (fun2_t)&ellint_1),
		   make_baseline_function("GSL", "gsl::ellint_1", (fun2_t)&gsl::ellint_1),
		   [](Real x){ return (std::abs(x) == 1 ? false : true); },
		   make_argument("k", fill2(Real{-1}, Real{1}, 21)),
		   make_argument("phi", vphid));
    std::ofstream file_ellint_1(outfile("ellint_1"));
    test_ellint_1(file_ellint_1);

    // Elliptic integrals of the second kind.
    // Avoid poles at |x| = 1.
    auto test_ellint_2 = 
    make_testcase2(make_test_function("std::ellint_2", (fun2_t)&ellint_2),
		   make_baseline_function("GSL", "gsl::ellint_2", (fun2_t)&gsl::ellint_2),
		   [](Real x){ return (std::abs(x) == 1 ? false : true); },
		   make_argument("k", fill2(Real{-1}, Real{1}, 21)),
		   make_argument("phi", vphid));
    std::ofstream file_ellint_2(outfile("ellint_2"));
    test_ellint_2(file_ellint_2);

    // Elliptic integrals of the third kind.
    // Avoid poles at |x| = 1 and at nu = 1.
    auto test_ellint_3 = 
    make_testcase2(make_test_function("std::ellint_3", (fun3_t)&ellint_3),
		   make_baseline_function("GSL", "gsl::ellint_3", (fun3_t)&gsl::ellint_3),
		   [](Real x){ return (std::abs(x) == 1 ? false : true); },
		   make_argument("k", fill2(Real{-1}, Real{1}, 21)),
		   make_argument("nu", fill2(Real{0}, Real{1}, 11)),
		   make_argument("phi", vphid));
    std::ofstream file_ellint_3(outfile("ellint_3"));
    test_ellint_3(file_ellint_3);

    // Exponential integral.
    // Skip the pole at 0.
    auto test_expint = 
    make_testcase2(make_test_function("std::expint", (fun_t)&expint),
		   make_baseline_function("GSL", "gsl::expint", (fun_t)&gsl::expint),
		   [](Real x){ return (std::abs(x) == 1 ? false : true); },
		   make_argument("x", fill2(Real{-1}, Real{1}, 21)));
    std::ofstream file_expint(outfile("expint"));
    test_expint(file_expint);

    // Hermite polynomials
    auto test_hermite = 
    make_testcase2(make_test_function("std::hermite", (funir_t)&hermite),
		   make_baseline_function("GSL", "gsl::hermite", (funir_t)&gsl::hermite),
		   [](unsigned, Real){ return true; },
		   "testcase_hermite",
		   make_argument("n", vorder),
		   make_argument("x", fill2(Real{-10}, Real{10}, 201)));
    std::ofstream file_hermite(outfile("hermite"));
    test_hermite(file_hermite);

#if STD
    // Hypergeometric functions.
    // Skip the singularity at c = 0 and at x = -1.
    auto test_hyperg = 
    make_testcase2(make_test_function("__gnu_cxx::hyperg", (fun4_t)&hyperg),
		   make_baseline_function("GSL", "gsl::hyperg", (fun4_t)&gsl::hyperg),
		   [](Real, Real, Real c, Real x){ return ((c == 0 || x == -1) ? false : true); },
		   make_argument("a", vab), make_argument("b", vab),
		   make_argument("c", fill2(Real{0}, Real{10}, 11)),
		   make_argument("x", fill2(Real{-10}, Real{10}, 21)));
    std::ofstream file_hyperg(outfile("hyperg"));
    test_hyperg(file_hyperg);
#else
    // Hypergeometric functions.
    // Skip the singularity at c = 0 and at x = -1.
    auto test_hyperg = 
    make_testcase2(make_test_function("std::hyperg", (fun4_t)&hyperg),
		   make_baseline_function("GSL", "gsl::hyperg", (fun4_t)&gsl::hyperg),
		   [](Real, Real, Real c, Real x){ return ((c == 0 || x == -1) ? false : true); },
		   make_argument("a", vab), make_argument("b", vab),
		   make_argument("c", fill2(Real{0}, Real{10}, 11)),
		   make_argument("x", fill2(Real{-10}, Real{10}, 21)));
    std::ofstream file_hyperg(outfile("hyperg"));
    test_hyperg(file_hyperg);
#endif // STD

    // Laguerre polynomials.
    auto test_laguerre = 
    make_testcase2(make_test_function("std::laguerre", (funir_t)&laguerre),
		   make_baseline_function("GSL", "gsl::laguerre", (funir_t)&gsl::laguerre),
		   [](unsigned, Real){ return true; },
		   "testcase_laguerre",
		   make_argument("n", vorder),
		   make_argument("x", fill2(Real{0}, Real{100}, 21)));
    std::ofstream file_laguerre(outfile("laguerre"));
    test_laguerre(file_laguerre);

    // Legendre polynomials.
    auto test_legendre = 
    make_testcase2(make_test_function("std::legendre", (funir_t)&legendre),
		   make_baseline_function("GSL", "gsl::legendre", (funir_t)&gsl::legendre_p),
		   [](unsigned, Real){ return true; },
		   "testcase_legendre",
		   make_argument("l", vorder),
		   make_argument("x", fill2(Real{0}, Real{100}, 21)));
    std::ofstream file_legendre(outfile("legendre"));
    test_legendre(file_legendre);

    // Riemann zeta function.
    // Skip the pole at 1.
    auto test_riemann_zeta = 
    make_testcase2(make_test_function("std::riemann_zeta", (fun_t)&riemann_zeta),
		   make_baseline_function("GSL", "gsl::riemann_zeta", (fun_t)&gsl::riemann_zeta),
		   [](Real x){ return (x == 1 ? false : true); },
		   make_argument("s", fill2(Real{-10}, Real{+30}, 205)));
    std::ofstream file_riemann_zeta(outfile("riemann_zeta"));
    test_riemann_zeta(file_riemann_zeta);

#if STD
    // Hurwitz zeta functions.
    // Skip the pole at 1.
    auto test_hurwitz_zeta = 
    make_testcase2(make_test_function("std::hurwitz_zeta", (fun2_t)&hurwitz_zeta),
		   make_baseline_function("GSL", "gsl::hurwitz_zeta", (fun2_t)&gsl::hurwitz_zeta),
		   [](Real s, Real){ return (s == 1 ? false : true); },
		   make_argument("s", fill2(Real{-10}, Real{+30}, 205)),
		   make_argument("a", fill2(Real{0}, Real{5}, 26)));
    std::ofstream file_hurwitz_zeta(outfile("hurwitz_zeta"));
    test_hurwitz_zeta(file_hurwitz_zeta);
#endif // STD

    // Spherical Bessel functions.
    auto test_sph_bessel = 
    make_testcase2(make_test_function("std::sph_bessel", (funir_t)&sph_bessel),
		   make_baseline_function("GSL", "gsl::sph_bessel", (funir_t)&gsl::sph_bessel),
		   [](unsigned, Real){ return true; },
		   make_argument("n", sph_order),
		   make_argument("x", fill2(Real{0}, Real{5}, 21)));
    std::ofstream file_sph_bessel(outfile("sph_bessel"));
    test_sph_bessel(file_sph_bessel);

    // Spherical Legendre functions.
    auto test_sph_legendre = 
    make_testcase2(make_test_function("std::sph_legendre", sph_legendre),
		   make_baseline_function("GSL", "gsl::sph_legendre", gsl::legendre_sphPlm),
		   [](unsigned, unsigned, Real){ return true; },
		   make_argument("l", vorder), make_argument("m", vorder),
		   make_argument("x", fill2(Real{0}, static_cast<Real>(M_PI), 21)));
    std::ofstream file_sph_legendre(outfile("sph_legendre"));
    test_sph_legendre(file_sph_legendre);

    // Spherical Neumann functions.
    auto test_sph_neumann = 
    make_testcase2(make_test_function("std::sph_neumann", (funir_t)&sph_neumann),
		   make_baseline_function("GSL", "gsl::sph_neumann", (funir_t)&gsl::sph_neumann),
		   [](unsigned, Real){ return true; },
		   make_argument("n", sph_order),
		   make_argument("x", fill2(Real{0}, Real{5}, 21)));
    std::ofstream file_sph_neumann(outfile("sph_neumann"));
    test_sph_neumann(file_sph_neumann);

    // Spherical harmonic functions.
    auto test_sph_harmonic = 
    make_testcase2(make_test_function("std::sph_harmonic", (funir_t)&sph_harmonic),
		   make_baseline_function("GSL", "gsl::sph_harmonic", (funir_t)&gsl::sph_harmonic),
		   [](unsigned, unsigned, Real, Real){ return true; },
		   make_argument("l", vorder), make_argument("m", vorder),
		   make_argument("theta", fill2(Real{0}, static_cast<Real>(M_PI), 21)),
		   make_argument("phi", vphid));
    std::ofstream file_sph_harmonic(outfile("sph_harmonic"));
    test_sph_harmonic(file_sph_harmonic);

    // Dirichlet eta function.
    // Skip the pole at 1.
    auto test_dirichlet_eta = 
    make_testcase2(make_test_function("std::dirichlet_eta", (fun_t)&dirichlet_eta),
		   make_baseline_function("GSL", "gsl::dirichlet_eta", (fun_t)&gsl::dirichlet_eta),
		   [](Real x){ return (x == 1 ? false : true); },
		   make_argument("s", fill2(Real{-10}, Real{+30}, 205)));
    std::ofstream file_dirichlet_eta(outfile("dirichlet_eta"));
    test_dirichlet_eta(file_dirichlet_eta);

    // Owens T functions.
    auto test_owens_t = 
    make_testcase2(make_test_function("std::owens_t", (fun2_t)&owens_t),
		   make_baseline_function("GSL", "gsl::owens_t", (fun2_t)&gsl::owens_t),
		   [](Real, Real){ return true; },
		   make_argument("h", fill2(Real{-5}, Real{+5}, 41)),
		   make_argument("a", fill2(Real{0}, Real{5}, 21)));
    std::ofstream file_owens_t(outfile("owens_t"));
    test_owens_t(file_owens_t);

    // Clausen C function.
    auto test_clausen_c = 
    make_testcase2(make_test_function("std::clausen_c", (funir_t)&clausen_c),
		   make_baseline_function("GSL", "gsl::clausen_c", (funir_t)&gsl::clausen_c),
		   [](unsigned int, Real){ return true; },
		   make_argument("m", fill2(2U, 2U, 1)),
		   make_argument("w", fill2(Real{-10}, Real{+10}, 41)));
    std::ofstream file_clausen_c(outfile("clausen_c"));
    test_clausen_c(file_clausen_c);
  }


int
main()
{
  harness<double>();

  return 0;
}


