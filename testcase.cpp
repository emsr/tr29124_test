#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <experimental/string_view>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <functional>
#include <utility>
#include <tuple>

#define STD 1

#include "specfun_testcase.h"
#include "gsl_wrap.h"

#include "testcase.tcc"

std::string
get_filename(const std::string & path,
	     const std::string & prefix,
	     const std::string & funcname,
	     const std::string & extra,
	     const std::string & suffix)
{
  auto filename = path + "/" + prefix;
  filename += funcname + extra + suffix;

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
    using       std::beta;
    using __gnu_cxx::bincoef;
    using       std::comp_ellint_1;
    using       std::comp_ellint_2;
    using       std::comp_ellint_3;
    using __gnu_cxx::conf_hyperg;
    using __gnu_cxx::coshint;
    using __gnu_cxx::cosint;
    using       std::cyl_bessel_i;
    using       std::cyl_bessel_j;
    using       std::cyl_bessel_k;
    using       std::cyl_neumann;
    using __gnu_cxx::dawson;
    using __gnu_cxx::dilog;
    using __gnu_cxx::double_factorial;
    using       std::ellint_1;
    using       std::ellint_2;
    using       std::ellint_3;
    using __gnu_cxx::ellint_rc;
    using __gnu_cxx::ellint_rd;
    using __gnu_cxx::ellint_rf;
    using __gnu_cxx::ellint_rg;
    using __gnu_cxx::ellint_rj;
    using       std::expint;
    using __gnu_cxx::expint_e1;
    using __gnu_cxx::factorial;
    using __gnu_cxx::fresnel_c;
    using __gnu_cxx::fresnel_s;
    using __gnu_cxx::gamma_l;
    using __gnu_cxx::gamma_u;
    using __gnu_cxx::gegenbauer;
    using       std::hermite;
    using __gnu_cxx::hurwitz_zeta;
    using __gnu_cxx::hyperg;
    using __gnu_cxx::ibeta;
    using __gnu_cxx::jacobi_sn;
    using __gnu_cxx::jacobi_cn;
    using __gnu_cxx::jacobi_dn;
    using       std::laguerre;
    using __gnu_cxx::lbincoef;
    using __gnu_cxx::ldouble_factorial;
    using       std::legendre;
    using __gnu_cxx::legendre_q;
    using __gnu_cxx::lfactorial;
    using __gnu_cxx::lpochhammer_l;
    using __gnu_cxx::lpochhammer_u;
    using __gnu_cxx::pochhammer_l;
    using __gnu_cxx::pochhammer_u;
    using __gnu_cxx::psi;
    using       std::riemann_zeta;
    using __gnu_cxx::sinc;
    using __gnu_cxx::sinhint;
    using __gnu_cxx::sinint;
    using       std::sph_bessel;
    using __gnu_cxx::sph_bessel_i;
    using __gnu_cxx::sph_bessel_k;
    using       std::sph_legendre;
    using       std::sph_neumann;
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

    // Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
    std::vector<unsigned int> vorder{0, 1, 2, 5, 10, 20, 50, 100};

    // ... corresponding "Real" integer orders for GSL.
    std::vector<Real> dvorder{0, 1, 2, 5, 10, 20, 50, 100};

    // Orders for cylindrical Bessel functions.
    std::vector<Real> vborderd{0, Real{1.0L/3.0L},
			       Real{0.5L}, Real{2.0L/3.0L},
			       1, 2, 5, 10, 20, 50, 100};

    // Orders for spherical bessel functions.
    std::vector<unsigned int> sborder{ 0, 1, 2, 5, 10, 20, 50, 100 };

    const unsigned int num_phi = 10;
    Real phi[num_phi];
    for (unsigned int i = 0; i < num_phi; ++i)
      phi[i] = Real{10} * i * __gnu_cxx::__math_constants<Real>::__pi / Real{180};
    std::vector<Real> vphid(phi, phi + num_phi);

    std::vector<Real> vab{0, Real{0.5L}, 1, 2, 5, 10, 20};

    unsigned int test = 1;

    const std::string path = ".";
    const std::string prefix = "/check_";

    std::string nsname = "std";

    std::string funcname;
    std::string filename;

#if STD
    // Airy functions.
    std::cout << "airy_ai" << std::endl;
    funcname = "airy_ai";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_airy_ai(filename.c_str());
    maketest(airy_ai,
	     gsl::airy_ai,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
		    		std::make_pair(true, true), 41),
	     file_airy_ai);

    std::cout << "airy_bi" << std::endl;
    funcname = "airy_bi";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_airy_bi(filename.c_str());
    maketest(airy_bi, gsl::airy_bi,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
	    			std::make_pair(true, true), 41),
	     file_airy_bi);
#endif // STD

    // Associated Laguerre polynomials.
    std::cout << "assoc_laguerre" << std::endl;
    funcname = "assoc_laguerre";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_assoc_laguerre(filename.c_str());
    maketest(assoc_laguerre, gsl::laguerre_nm,
	     nsname, funcname,
	     "n", vorder, "m", vorder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 11),
	     file_assoc_laguerre);

    // Associated Legendre functions.
    std::cout << "assoc_legendre" << std::endl;
    funcname = "assoc_legendre";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_assoc_legendre(filename.c_str());
    maketest(assoc_legendre, gsl::legendre_Plm,
	     nsname, funcname,
	     "l", vorder, "m", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, true), 21),
	     file_assoc_legendre);

    // Beta function.
    std::cout << "beta" << std::endl;
    funcname = "beta";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_beta(filename.c_str());
    maketest(beta, gsl::beta,
	     nsname, funcname,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 11),
	     file_beta);

    // Complete elliptic integrals of the first kind.
    // Avoid poles at |x| = 1.
    std::cout << "comp_ellint_1" << std::endl;
    funcname = "comp_ellint_1";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_comp_ellint_1(filename.c_str());
    maketest(comp_ellint_1, gsl::ellint_Kcomp,
	     nsname, funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
	    			std::make_pair(false, false), 21),
	     file_comp_ellint_1);

    // Complete elliptic integrals of the second kind.
    // Avoid poles at |x| = 1.
    std::cout << "comp_ellint_2" << std::endl;
    funcname = "comp_ellint_2";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_comp_ellint_2(filename.c_str());
    maketest(comp_ellint_2, gsl::ellint_Ecomp,
	     nsname, funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
	    			std::make_pair(false, false), 21),
	     file_comp_ellint_2);

    // Complete elliptic integrals of the third kind.
    // Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3" << std::endl;
    funcname = "comp_ellint_3";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_comp_ellint_3(filename.c_str());
    maketest(comp_ellint_3, gsl::ellint_Pcomp,
	     nsname, funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, false), 11),
	     file_comp_ellint_3);

#if STD
    // Confluent hypergeometric functions.
    // Skip the singularity at c = 0.
    std::cout << "conf_hyperg" << std::endl;
    funcname = "conf_hyperg";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_conf_hyperg(filename.c_str());
    maketest(conf_hyperg, gsl::hyperg_1F1,
	     "__gnu_cxx", funcname,
	     "a", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 21),
	     file_conf_hyperg);
#else
    // Confluent hypergeometric functions.
    // Skip the singularity at c = 0.
    std::cout << "conf_hyperg" << std::endl;
    funcname = "conf_hyperg";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_conf_hyperg(filename.c_str());
    maketest(conf_hyperg, gsl::hyperg_1F1,
	     nsname, funcname,
	     "a", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 21),
	     file_conf_hyperg);
#endif // STD

    // Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i" << std::endl;
    funcname = "cyl_bessel_i";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_bessel_i(filename.c_str());
    test =
    maketest(cyl_bessel_i, gsl::bessel_Inu,
	     nsname, funcname,
	     "nu", vborderd,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     file_cyl_bessel_i, true, false);
    maketest(cyl_bessel_i, gsl::bessel_Inu,
	     nsname, funcname,
	     "nu", vborderd,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     file_cyl_bessel_i, false, true, test);

    // Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j" << std::endl;
    funcname = "cyl_bessel_j";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_bessel_j(filename.c_str());
    test =
    maketest(cyl_bessel_j, gsl::bessel_Jnu,
	     nsname, funcname,
	     "nu", vborderd,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     file_cyl_bessel_j, true, false);
    maketest(cyl_bessel_j, gsl::bessel_Jnu,
	     nsname, funcname,
	     "nu", vborderd,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     file_cyl_bessel_j, false, true, test);

    // Irregular modified cylindrical Bessel functions.
    // Skip the pole at the origin.
    std::cout << "cyl_bessel_k" << std::endl;
    funcname = "cyl_bessel_k";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_bessel_k(filename.c_str());
    test =
    maketest(cyl_bessel_k, gsl::bessel_Knu,
	     nsname, funcname,
	     "nu", vborderd,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     file_cyl_bessel_k, true, false);
    maketest(cyl_bessel_k, gsl::bessel_Knu,
	     nsname, funcname,
	     "nu", vborderd,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     file_cyl_bessel_k, false, true, test);

    // Cylindrical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "cyl_neumann" << std::endl;
    funcname = "cyl_neumann";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_neumann(filename.c_str());
    test =
    maketest(cyl_neumann, gsl::bessel_Ynu,
	     nsname, funcname,
	     "nu", vborderd,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
	     			std::make_pair(false, true), 21),
	     file_cyl_neumann, true, false);
    maketest(cyl_neumann, gsl::bessel_Ynu,
	     nsname, funcname,
	     "nu", vborderd,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     file_cyl_neumann, false, true, test);

    // Elliptic integrals of the first kind.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_1" << std::endl;
    funcname = "ellint_1";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_1(filename.c_str());
    maketest(ellint_1, gsl::ellint_F,
	     nsname, funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     file_ellint_1);

    // Elliptic integrals of the second kind.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_2" << std::endl;
    funcname = "ellint_2";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_2(filename.c_str());
    maketest(ellint_2, gsl::ellint_E,
	     nsname, funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "phi", vphid,
	     file_ellint_2);

    // Elliptic integrals of the third kind.
    // Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3" << std::endl;
    funcname = "ellint_3";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_ellint_3(filename.c_str());
    maketest(ellint_3, gsl::ellint_P,
	     nsname, funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(false, false), 21),
	     "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, false), 11),
	     "phi", vphid,
	     file_ellint_3);

    // Exponential integral.
    // Skip the pole at 0.
    std::cout << "expint" << std::endl;
    funcname = "expint";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_expint(filename.c_str());
    test =
    maketest(expint, gsl::expint_Ei,
	     nsname, funcname,
	     "x", fill_argument(std::make_pair(Real{-50}, Real{0}),
	    			std::make_pair(true, false), 51),
	     file_expint, true, false);
    maketest(expint, gsl::expint_Ei,
	     nsname, funcname,
	     "x", fill_argument(std::make_pair(Real{0}, Real{50}),
	    			std::make_pair(false, true), 51),
	     file_expint, false, true, test);

    // Hermite polynomials
    std::cout << "hermite" << std::endl;
    funcname = "hermite";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_hermite(filename.c_str());
    maketest(hermite, gsl::hermite,
	     nsname, funcname,
	     "n", vorder,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				std::make_pair(true, true), 201),
	     file_hermite);

#if STD
    // Hypergeometric functions.
    // Skip the singularity at c = 0.
    // Skip the singularity at x = -1.
    std::cout << "hyperg" << std::endl;
    funcname = "hyperg";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_hyperg(filename.c_str());
    maketest(hyperg, gsl::hyperg_2F1,
	     "__gnu_cxx", funcname,
	     "a", vab,
	     "b", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
		  		std::make_pair(false, true), 6),
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
		  		std::make_pair(true, false), 21),
	     file_hyperg);
#else
    // Hypergeometric functions.
    // Skip the singularity at c = 0.
    // Skip the singularity at x = -1.
    std::cout << "hyperg" << std::endl;
    funcname = "hyperg";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_hyperg(filename.c_str());
    maketest(hyperg, gsl::hyperg_2F1,
	     nsname, funcname,
	     "a", vab,
	     "b", vab,
	     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
		  		std::make_pair(false, true), 6),
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
		  		std::make_pair(true, false), 21),
	     file_hyperg);
#endif // STD

    // Laguerre polynomials.
    std::cout << "laguerre" << std::endl;
    funcname = "laguerre";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_laguerre(filename.c_str());
    maketest(laguerre, gsl::laguerre_n,
	     nsname, funcname,
	     "n", vorder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     file_laguerre);

    // Legendre polynomials.
    std::cout << "legendre" << std::endl;
    funcname = "legendre";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_legendre(filename.c_str());
    maketest(legendre, gsl::legendre_Pl,
	     nsname, funcname,
	     "l", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
				std::make_pair(true, true), 21),
	     file_legendre);

    // Riemann zeta function.
    std::cout << "riemann_zeta" << std::endl;
    // Skip the pole at 1.
    funcname = "riemann_zeta";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_riemann_zeta(filename.c_str());
    test =
    maketest(riemann_zeta, gsl::zeta,
	     nsname, funcname,
	     "s", fill_argument(std::make_pair(Real{-10}, Real{1}),
	    			std::make_pair(true, false), 56),
	     file_riemann_zeta, true, false);
    maketest(riemann_zeta, gsl::zeta,
	     nsname, funcname,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
	    			std::make_pair(false, true), 146),
	     file_riemann_zeta, false, true, test);

#if STD
    // Hurwitz zeta function.
    std::cout << "hurwitz_zeta" << std::endl;
    // Skip the pole at 1.
    funcname = "hurwitz_zeta";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_hurwitz_zeta(filename.c_str());
    maketest(hurwitz_zeta, gsl::hzeta,
	     "__gnu_cxx", funcname,
	     "s", fill_argument(std::make_pair(Real{1}, Real{30}),
				std::make_pair(false, true), 146),
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 26),
	     file_hurwitz_zeta);
#endif // STD

    // Spherical Bessel functions.
    std::cout << "sph_bessel" << std::endl;
    funcname = "sph_bessel";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_sph_bessel(filename.c_str());
    test =
    maketest(sph_bessel, gsl::bessel_jl,
	     nsname, funcname,
	     "n", sborder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     file_sph_bessel, true, false);
    maketest(sph_bessel, gsl::bessel_jl,
	     nsname, funcname,
	     "n", sborder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     file_sph_bessel, false, true, test);

    // Spherical Legendre functions.
    std::cout << "sph_legendre" << std::endl;
    funcname = "sph_legendre";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_sph_legendre(filename.c_str());
    maketest(sph_legendre, gsl::legendre_sphPlm,
	     nsname, funcname,
	     "l", vorder, "m", vorder,
	     "theta", fill_argument(std::make_pair(Real{0}, static_cast<Real>(M_PI)),
				    std::make_pair(true, true), 21),
	     file_sph_legendre);

    // Spherical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "sph_neumann" << std::endl;
    funcname = "sph_neumann";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_sph_neumann(filename.c_str());
    test =
    maketest(sph_neumann, gsl::bessel_yl,
	     nsname, funcname,
	     "n", sborder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     file_sph_neumann, true, false);
    maketest(sph_neumann, gsl::bessel_yl,
	     nsname, funcname,
	     "n", sborder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     file_sph_neumann, false, true, test);

#if STD
    // Carlson elliptic functions R_C.
    std::cout << "ellint_rc" << std::endl;
    funcname = "ellint_rc";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_rc(filename.c_str());
    maketest(ellint_rc, gsl::ellint_RC,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     file_ellint_rc);

    // Carlson elliptic functions R_D.
    std::cout << "ellint_rd" << std::endl;
    funcname = "ellint_rd";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_rd(filename.c_str());
    maketest(ellint_rd, gsl::ellint_RD,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     file_ellint_rd);

    // Carlson elliptic functions R_F.
    std::cout << "ellint_rf" << std::endl;
    funcname = "ellint_rf";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_rf(filename.c_str());
    maketest(ellint_rf, gsl::ellint_RF,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     file_ellint_rf);

    // Carlson elliptic functions R_J.
    std::cout << "ellint_rj" << std::endl;
    funcname = "ellint_rj";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_rj(filename.c_str());
    maketest(ellint_rj, gsl::ellint_RJ,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(true, true), 11),
	     "p", fill_argument(std::make_pair(Real{0}, +Real{5}),
				std::make_pair(false, true), 11),
	     file_ellint_rj);

    // Dilogarithm functions.
    std::cout << "dilog" << std::endl;
    funcname = "dilog";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_dilog(filename.c_str());
    maketest(dilog, gsl::dilog,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
	    			std::make_pair(true, true), 41),
	     file_dilog);

    // Upper incomplete Gamma functions.
    std::cout << "gamma_u" << std::endl;
    funcname = "gamma_u";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_gamma_u(filename.c_str());
    maketest(gamma_u, gsl::gamma_inc,
	     "__gnu_cxx", funcname,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 11),
	     file_gamma_u);

    // Incomplete Beta functions.
    std::cout << "ibeta" << std::endl;
    funcname = "ibeta";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ibeta(filename.c_str());
    maketest<Real, Real, Real, Real>(ibeta, gsl::beta_inc,
	     "__gnu_cxx", funcname,
	     "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 11),
	     "b", fill_argument(std::make_pair(Real{5}, Real{0}),
				std::make_pair(true, false), 11),
	     "x", fill_argument(std::make_pair(Real{0}, Real{1}),
				std::make_pair(false, false), 21),
	     file_ibeta);

    // Digamma or psi functions.
    std::cout << "psi" << std::endl;
    funcname = "psi";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_psi(filename.c_str());
    test =
    maketest(psi, gsl::psi,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{-9.9375}, +Real{10.0625}),
	    			std::make_pair(true, true), 801),
	     file_psi, true, false);
    maketest(psi, gsl::psi,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{1}, Real{100}),
	    			std::make_pair(true, true), 199),
	     file_psi, false, true, test);

    // Sine integral or Si functions.
    std::cout << "sinint" << std::endl;
    funcname = "sinint";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_sinint(filename.c_str());
    maketest(sinint, gsl::Si,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
	    			std::make_pair(false, true), 101),
	     file_sinint);

    // Cosine integral or Ci functions.
    std::cout << "cosint" << std::endl;
    funcname = "cosint";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cosint(filename.c_str());
    maketest(cosint, gsl::Ci,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
	    			std::make_pair(false, true), 101),
	     file_cosint);

    // Hyperbolic sine integral or Shi functions.
    std::cout << "sinhint" << std::endl;
    funcname = "sinhint";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_sinhint(filename.c_str());
    maketest(sinhint, gsl::Shi,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
	    			std::make_pair(false, true), 101),
	     file_sinhint);

    // Hyperbolic cosine integral or Chi functions.
    std::cout << "coshint" << std::endl;
    funcname = "coshint";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_coshint(filename.c_str());
    maketest(coshint, gsl::Chi,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
	    			std::make_pair(false, true), 101),
	     file_coshint);

    // Dawson integral.
    std::cout << "dawson" << std::endl;
    funcname = "dawson";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_dawson(filename.c_str());
    maketest(dawson, gsl::dawson,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+20}),
	    			std::make_pair(false, true), 201),
	     file_dawson);

    // Jacobian elliptic integrals.
    std::cout << "jacobi_sn" << std::endl;
    funcname = "jacobi_sn";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_jacobi_sn(filename.c_str());
    maketest(jacobi_sn, gsl::elljac_sn,
	     "__gnu_cxx", funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
		    		std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
		    		std::make_pair(true, true), 101),
	     file_jacobi_sn);

    // Jacobian elliptic integrals.
    std::cout << "jacobi_cn" << std::endl;
    funcname = "jacobi_cn";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_jacobi_cn(filename.c_str());
    maketest(jacobi_cn, gsl::elljac_cn,
	     "__gnu_cxx", funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
		    		std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
		    		std::make_pair(true, true), 101),
	     file_jacobi_cn);

    // Jacobian elliptic integrals.
    std::cout << "jacobi_dn" << std::endl;
    funcname = "jacobi_dn";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_jacobi_dn(filename.c_str());
    maketest(jacobi_dn, gsl::elljac_dn,
	     "__gnu_cxx", funcname,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
		    		std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
		    		std::make_pair(true, true), 101),
	     file_jacobi_dn);

    // Exponential integral E1.
    // Skip the pole at 0.
    std::cout << "expint_e1" << std::endl;
    funcname = "expint_e1";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_expint_e1(filename.c_str());
    test =
    maketest(expint_e1,
	     gsl::expint_E1,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{-50}, Real{0}),
	    			std::make_pair(true, false), 51),
	     file_expint_e1, true, false);
    maketest(expint_e1,
	     gsl::expint_E1,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{0}, Real{50}),
	    			std::make_pair(false, true), 51),
	     file_expint_e1, false, true, test);

    // Fresnel cosine integral.
    std::cout << "fresnel_c" << std::endl;
    funcname = "fresnel_c";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_fresnel_c(filename.c_str());
    maketest(fresnel_c, gsl::fresnel_c,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
	    			std::make_pair(false, true), 401),
	     file_fresnel_c);

    // Fresnel sine integral.
    std::cout << "fresnel_s" << std::endl;
    funcname = "fresnel_s";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_fresnel_s(filename.c_str());
    maketest(fresnel_s, gsl::fresnel_s,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
	    			std::make_pair(false, true), 401),
	     file_fresnel_s);

    // Sine cardinal function.
    std::cout << "sinc" << std::endl;
    funcname = "sinc";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_sinc(filename.c_str());
    maketest(sinc, gsl::sinc,
	     "__gnu_cxx", funcname,
	     "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
	    			std::make_pair(true, true), 401),
	     file_sinc);

    // Log upper Pochhammer symbol
    std::cout << "lpochhammer_u" << std::endl;
    funcname = "lpochhammer_u";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_lpochhammer_u(filename.c_str());
    maketest(lpochhammer_u, gsl::lnpoch,
	     "__gnu_cxx", funcname,
	     "a", {1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0},
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     file_lpochhammer_u, true, true);

    // Log lower Pochhammer symbol
    //std::cout << "lpochhammer_l" << std::endl;
    //funcname = "lpochhammer_l";
    //filename = get_filename(path, prefix, funcname, "",  ".cc");
    //std::ofstream file_lpochhammer_l(filename.c_str());

    // Upper Pochhammer symbols (see boost::rising_factorial)
    std::cout << "pochhammer_u" << std::endl;
    funcname = "pochhammer_u";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_pochhammer_u(filename.c_str());
    maketest(pochhammer_u, gsl::poch,
	     "__gnu_cxx", funcname,
	     "a", dvorder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     file_pochhammer_u, true, true);

    // Lower Pochhammer symbols (see boost::falling_factorial)
    //std::cout << "pochhammer_l" << std::endl;
    //funcname = "pochhammer_l";
    //filename = get_filename(path, prefix, funcname, "",  ".cc");
    //std::ofstream file_pochhammer_l(filename.c_str());

    // Regular modified spherical bessel functions.
    std::cout << "sph_bessel_i" << std::endl;
    funcname = "sph_bessel_i";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_sph_bessel_i(filename.c_str());
    test =
    maketest(sph_bessel_i, gsl::bessel_il,
	     "__gnu_cxx", funcname,
	     "n", sborder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(true, true), 21),
	     file_sph_bessel_i, true, false);
    maketest(sph_bessel_i, gsl::bessel_il,
	     "__gnu_cxx", funcname,
	     "n", sborder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(true, true), 21),
	     file_sph_bessel_i, false, true, test);

    // Irregular modified spherical bessel functions.
    std::cout << "sph_bessel_k" << std::endl;
    funcname = "sph_bessel_k";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_sph_bessel_k(filename.c_str());
    test =
    maketest(sph_bessel_k, gsl::bessel_kl,
	     "__gnu_cxx", funcname,
	     "n", sborder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				std::make_pair(false, true), 21),
	     file_sph_bessel_k, true, false);
    maketest(sph_bessel_k, gsl::bessel_kl,
	     "__gnu_cxx", funcname,
	     "n", sborder,
	     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
				std::make_pair(false, true), 21),
	     file_sph_bessel_k, false, true, test);

    // Legendre functions of the second kind.
    std::cout << "legendre_q" << std::endl;
    funcname = "legendre_q";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_legendre_q(filename.c_str());
    maketest(legendre_q, gsl::legendre_Ql,
	     nsname, funcname,
	     "l", vorder,
	     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
	        		std::make_pair(false, false), 21),
	     file_legendre_q);

    // Factorial.
    std::cout << "factorial" << std::endl;
    funcname = "factorial";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_factorial(filename.c_str());
    maketest(factorial<Real>, gsl::fact,
	     "__gnu_cxx", funcname,
	     "n", fill_argument(std::make_pair(0U, 50U),
	    			std::make_pair(true, true), 51),
	     file_factorial);

    // Log factorial.
    std::cout << "lfactorial" << std::endl;
    funcname = "lfactorial";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_lfactorial(filename.c_str());
    maketest(lfactorial<Real>, gsl::lnfact,
	     "__gnu_cxx", funcname,
	     "n", fill_argument(std::make_pair(0U, 500U),
	    			std::make_pair(true, true), 501),
	     file_lfactorial);

    // Double factorial.
    std::cout << "double_factorial" << std::endl;
    funcname = "double_factorial";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_double_factorial(filename.c_str());
    maketest(double_factorial<Real>, gsl::doublefact,
	     "__gnu_cxx", funcname,
	     "n", fill_argument(std::make_pair(0, 50),
	    			std::make_pair(true, true), 51),
	     file_double_factorial);

    // Log double factorial.
    std::cout << "ldouble_factorial" << std::endl;
    funcname = "ldouble_factorial";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_ldouble_factorial(filename.c_str());
    maketest(ldouble_factorial<Real>, gsl::lndoublefact,
	     "__gnu_cxx", funcname,
	     "n", fill_argument(std::make_pair(0, 500),
	    			std::make_pair(true, true), 501),
	     file_ldouble_factorial);

    // Binomial coefficient.
    std::cout << "bincoef" << std::endl;
    funcname = "bincoef";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_bincoef(filename.c_str());
    maketest(bincoef<Real>, gsl::choose,
	     "__gnu_cxx", funcname,
	     "n", fill_argument(std::make_pair(0U, 50U),
	    			std::make_pair(true, true), 51),
	     "k", fill_argument(std::make_pair(0U, 50U),
	    			std::make_pair(true, true), 51),
	     file_bincoef);

    // Log binomial coefficient.
    std::cout << "lbincoef" << std::endl;
    funcname = "lbincoef";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_lbincoef(filename.c_str());
    maketest(lbincoef<Real>, gsl::lnchoose,
	     "__gnu_cxx", funcname,
	     "n", fill_argument(std::make_pair(0U, 500U),
	    			std::make_pair(true, true), 501),
	     "k", fill_argument(std::make_pair(0U, 500U),
	    			std::make_pair(true, true), 501),
	     file_lbincoef);

    // Gegenbauer polynomials
    std::cout << "gegenbauer" << std::endl;
    funcname = "gegenbauer";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_gegenbauer(filename.c_str());
    maketest(gegenbauer, gsl::gegenpoly_n,
	     "__gnu_cxx", funcname,
	     "n", vorder,
	     "alpha", fill_argument(std::make_pair(Real{0}, Real{5}),
	    			    std::make_pair(true, true), 11),
             "x", fill_argument(std::make_pair(Real{0}, Real{20}),
	    			std::make_pair(true, true), 41),
	     file_gegenbauer);

#endif // STD
  }


int
main()
{
  harness<double>();

  return 0;
}


