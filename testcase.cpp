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
    //  Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
    std::vector<unsigned int> vorder{0, 1, 2, 5, 10, 20, 50, 100};

    //  ... corresponding "Real" integer orders for GSL.
    std::vector<Real> dvorder{0, 1, 2, 5, 10, 20, 50, 100};

    //  Orders for cylindrical Bessel functions.
    std::vector<Real> vborderd{0, Real{1.0L/3.0L},
			       Real{0.5L}, Real{2.0L/3.0L},
			       1, 2, 5, 10, 20, 50, 100};

    //  Orders for spherical bessel functions.
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

    //  Airy functions.
    std::cout << "airy_ai" << std::endl;
    funcname = "airy_ai";
    filename = get_filename(path, prefix, "airy", "",  ".cc");
    std::ofstream file_airy(filename.c_str());
    typedef Real airy(Real);
    test =
    maketest<Real, Real>((airy *)__gnu_cxx::airy_ai,
			 wrap_gsl_sf_airy_ai,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(-Real{10}, Real{10}),
					    std::make_pair(true, true), 41),
			 file_airy, true, false);

    std::cout << "airy_bi" << std::endl;
    funcname = "airy_bi";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    maketest<Real, Real>((airy *)__gnu_cxx::airy_bi,
			 wrap_gsl_sf_airy_bi,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(-Real{10}, Real{10}),
					    std::make_pair(true, true), 41),
			 file_airy, false, true, test);

    //  Associated Laguerre polynomials.
    std::cout << "assoc_laguerre" << std::endl;
    funcname = "assoc_laguerre";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_assoc_laguerre(filename.c_str());
    typedef Real assoc_laguerre(unsigned int, unsigned int, Real);
    maketest<Real, unsigned int, unsigned int, Real>((assoc_laguerre *)std::assoc_laguerre,
						     wrap_gsl_sf_laguerre_nm,
						     nsname, funcname,
						     "n", vorder, "m", vorder,
						     "x", fill_argument(std::make_pair(Real{0}, Real{100}),
									std::make_pair(true, true), 11),
						     file_assoc_laguerre);

    //  Associated Legendre functions.
    std::cout << "assoc_legendre" << std::endl;
    funcname = "assoc_legendre";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_assoc_legendre(filename.c_str());
    typedef Real assoc_legendre(unsigned int, unsigned int, Real);
    maketest<Real, unsigned int, unsigned int, Real>((assoc_legendre *)std::assoc_legendre,
						     wrap_gsl_sf_legendre_Plm,
						     nsname, funcname,
						     "l", vorder, "m", vorder,
						     "x", fill_argument(std::make_pair(-Real{1}, Real{1}),
									std::make_pair(true, true), 21),
						     file_assoc_legendre);

    //  Beta function.
    std::cout << "beta" << std::endl;
    funcname = "beta";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_beta(filename.c_str());
    typedef Real beta(Real, Real);
    maketest<Real, Real, Real>((beta *)std::beta,
			       wrap_gsl_sf_beta,
			       nsname, funcname,
			       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
						  std::make_pair(false, true), 11),
			       "y", fill_argument(std::make_pair(Real{0}, Real{100}),
						  std::make_pair(false, true), 11),
			       file_beta);

    //  Complete elliptic integrals of the first kind.
    //  Avoid poles at |x| = 1.
    std::cout << "comp_ellint_1" << std::endl;
    funcname = "comp_ellint_1";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_comp_ellint_1(filename.c_str());
    typedef Real comp_ellint_1(Real);
    maketest<Real, Real>((comp_ellint_1 *)std::comp_ellint_1,
			 wrap_gsl_sf_ellint_Kcomp,
			 nsname, funcname,
			 "k", fill_argument(std::make_pair(-Real{1}, Real{1}),
					    std::make_pair(false, false), 21),
			 file_comp_ellint_1);

    //  Complete elliptic integrals of the second kind.
    //  Avoid poles at |x| = 1.
    std::cout << "comp_ellint_2" << std::endl;
    funcname = "comp_ellint_2";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_comp_ellint_2(filename.c_str());
    typedef Real comp_ellint_2(Real);
    maketest<Real, Real>((comp_ellint_2 *)std::comp_ellint_2,
			 wrap_gsl_sf_ellint_Ecomp,
			 nsname, funcname,
			 "k", fill_argument(std::make_pair(-Real{1}, Real{1}),
					    std::make_pair(false, false), 21),
			 file_comp_ellint_2);

    //  Complete elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3" << std::endl;
    funcname = "comp_ellint_3";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_comp_ellint_3(filename.c_str());
    typedef Real comp_ellint_3(Real, Real);
    maketest<Real, Real, Real>((comp_ellint_3 *)std::comp_ellint_3,
			       wrap_gsl_sf_ellint_Pcomp,
			       nsname, funcname,
			       "k", fill_argument(std::make_pair(-Real{1}, Real{1}),
						  std::make_pair(false, false), 21),
			       "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
						   std::make_pair(true, false), 11),
			       file_comp_ellint_3);

    //  Confluent hypergeometric functions.
    //  Skip the singularity at c = 0.
    std::cout << "conf_hyperg" << std::endl;
    funcname = "conf_hyperg";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_conf_hyperg(filename.c_str());
    typedef Real conf_hyperg(Real, Real, Real);
    maketest<Real, Real, Real, Real>((conf_hyperg *)__gnu_cxx::conf_hyperg,
				     wrap_gsl_sf_hyperg_1F1,
				     "__gnu_cxx", funcname,
				     "a", vab,
				     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
							std::make_pair(false, true), 11),
				     "x", fill_argument(std::make_pair(-Real{10}, Real{10}),
							std::make_pair(true, true), 21),
				     file_conf_hyperg);

    //  Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i" << std::endl;
    funcname = "cyl_bessel_i";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_bessel_i(filename.c_str());
    typedef Real cyl_bessel_i(Real, Real);
    test =
    maketest<Real, Real, Real>((cyl_bessel_i*)std::cyl_bessel_i,
			       wrap_gsl_sf_bessel_Inu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(true, true), 21),
			       file_cyl_bessel_i, true, false);
    maketest<Real, Real, Real>((cyl_bessel_i*)std::cyl_bessel_i,
			       wrap_gsl_sf_bessel_Inu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
						  std::make_pair(true, true), 21),
			       file_cyl_bessel_i, false, true, test);

    //  Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j" << std::endl;
    funcname = "cyl_bessel_j";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_bessel_j(filename.c_str());
    typedef Real cyl_bessel_j(Real, Real);
    test =
    maketest<Real, Real, Real>((cyl_bessel_j*)std::cyl_bessel_j,
			       wrap_gsl_sf_bessel_Jnu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(true, true), 21),
			       file_cyl_bessel_j, true, false);
    maketest<Real, Real, Real>((cyl_bessel_j*)std::cyl_bessel_j,
			       wrap_gsl_sf_bessel_Jnu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
						  std::make_pair(true, true), 21),
			       file_cyl_bessel_j, false, true, test);

    //  Cylindrical Bessel functions (of the first kind) asymptotics.
  /*
    std::cout << "cyl_bessel_j asymptotics" << std::endl;
    funcname = "cyl_bessel_j_asymp";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_bessel_j_asymp(filename.c_str());
    maketest<Real, Real, Real>((cyl_bessel_j*)std::cyl_bessel_j,
			       wrap_gsl_sf_bessel_Jnu_asymp,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{1000}, Real{10000}),
						  std::make_pair(true, true), 11),
			       file_cyl_bessel_j_asymp);
  */
    //  Irregular modified cylindrical Bessel functions.
    // Skip the pole at the origin.
    std::cout << "cyl_bessel_k" << std::endl;
    funcname = "cyl_bessel_k";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_bessel_k(filename.c_str());
    typedef Real cyl_bessel_k(Real, Real);
    test =
    maketest<Real, Real, Real>((cyl_bessel_k*)std::cyl_bessel_k,
			       wrap_gsl_sf_bessel_Knu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(false, true), 21),
			       file_cyl_bessel_k, true, false);
    maketest<Real, Real, Real>((cyl_bessel_k*)std::cyl_bessel_k,
			       wrap_gsl_sf_bessel_Knu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
						  std::make_pair(false, true), 21),
			       file_cyl_bessel_k, false, true, test);

    //  Cylindrical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "cyl_neumann" << std::endl;
    funcname = "cyl_neumann";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_neumann(filename.c_str());
    typedef Real cyl_neumann(Real, Real);
    test =
    maketest<Real, Real, Real>((cyl_neumann*)std::cyl_neumann,
				  wrap_gsl_sf_bessel_Ynu,
				  nsname, funcname,
				  "nu", vborderd,
				  "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						     std::make_pair(false, true), 21),
				  file_cyl_neumann, true, false);
    maketest<Real, Real, Real>((cyl_neumann*)std::cyl_neumann,
			       wrap_gsl_sf_bessel_Ynu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
						  std::make_pair(false, true), 21),
			       file_cyl_neumann, false, true, test);

    //  Cylindrical Neumann functions asymptotics.
    // Skip the pole at the origin.
  /*
    std::cout << "cyl_neumann asymptotics" << std::endl;
    funcname = "cyl_neumann_asymp";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_neumann_asymp(filename.c_str());
    maketest<Real, Real, Real>((cyl_neumann*)std::cyl_neumann,
			       wrap_gsl_sf_bessel_Ynu_asymp,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{1000}, Real{10000}),
						  std::make_pair(false, true), 11),
			       file_cyl_neumann_asymp);
  */
    //  Elliptic integrals of the first kind.
    //  Avoid poles at |x| = 1.
    std::cout << "ellint_1" << std::endl;
    funcname = "ellint_1";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_1(filename.c_str());
    typedef Real ellint_1(Real, Real);
    maketest<Real, Real, Real>((ellint_1*)std::ellint_1,
			       wrap_gsl_sf_ellint_F,
			       nsname, funcname,
			       "k", fill_argument(std::make_pair(-Real{1}, Real{1}),
						  std::make_pair(false, false), 21),
			       "phi", vphid,
			       file_ellint_1);

    //  Elliptic integrals of the second kind.
    //  Avoid poles at |x| = 1.
    std::cout << "ellint_2" << std::endl;
    funcname = "ellint_2";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_2(filename.c_str());
    typedef Real ellint_2(Real, Real);
    maketest<Real, Real, Real>((ellint_2*)std::ellint_2,
			       wrap_gsl_sf_ellint_E,
			       nsname, funcname,
			       "k", fill_argument(std::make_pair(-Real{1}, Real{1}),
						  std::make_pair(false, false), 21),
			       "phi", vphid,
			       file_ellint_2);

    //  Elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3" << std::endl;
    funcname = "ellint_3";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_ellint_3(filename.c_str());
    typedef Real ellint_3(Real, Real, Real);
    maketest<Real, Real, Real, Real>((ellint_3*)std::ellint_3,
				     wrap_gsl_sf_ellint_P,
				     nsname, funcname,
				     "k", fill_argument(std::make_pair(-Real{1}, Real{1}),
							std::make_pair(false, false), 21),
				     "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
							 std::make_pair(true, false), 11),
				     "phi", vphid,
				     file_ellint_3);

    //  Exponential integral.
    //  Skip the pole at 0.
    std::cout << "expint" << std::endl;
    funcname = "expint";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_expint(filename.c_str());
    typedef Real expint(Real);
    test =
    maketest<Real, Real>((expint*)std::expint,
			 wrap_gsl_sf_expint_Ei,
			 nsname, funcname,
			 "x", fill_argument(std::make_pair(-Real{50}, Real{0}),
					    std::make_pair(true, false), 51),
			 file_expint, true, false);
    maketest<Real, Real>((expint*)std::expint,
			 wrap_gsl_sf_expint_Ei,
			 nsname, funcname,
			 "x", fill_argument(std::make_pair(Real{0}, Real{50}),
					    std::make_pair(false, true), 51),
			 file_expint, false, true, test);

    //  Hermite polynomials
    std::cout << "hermite  UNTESTED" << std::endl;

    //  Hypergeometric functions.
    //  Skip the singularity at c = 0.
    //  Skip the singularity at x = -1.
    std::cout << "hyperg" << std::endl;
    funcname = "hyperg";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_hyperg(filename.c_str());
    typedef Real hyperg(Real, Real, Real, Real);
    maketest<Real, Real, Real, Real, Real>((hyperg*)__gnu_cxx::hyperg,
					   wrap_gsl_sf_hyperg_2F1,
					   "__gnu_cxx", funcname,
					   "a", vab,
					   "b", vab,
					   "c", fill_argument(std::make_pair(Real{0}, Real{10}),
							      std::make_pair(false, true), 6),
					   "x", fill_argument(std::make_pair(-Real{1}, Real{1}),
							      std::make_pair(true, false), 21),
					   file_hyperg);

    //  Laguerre polynomials.
    std::cout << "laguerre" << std::endl;
    funcname = "laguerre";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_laguerre(filename.c_str());
    typedef Real laguerre(unsigned int, Real);
    maketest<Real, unsigned int, Real>((laguerre*)std::laguerre,
				       wrap_gsl_sf_laguerre_n,
				       nsname, funcname,
				       "n", vorder,
				       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
							  std::make_pair(true, true), 21),
				       file_laguerre);

    //  Legendre polynomials.
    std::cout << "legendre" << std::endl;
    funcname = "legendre";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_legendre(filename.c_str());
    typedef Real legendre(unsigned int, Real);
    maketest<Real, unsigned int, Real>((legendre*)std::legendre,
				       wrap_gsl_sf_legendre_Pl,
				       nsname, funcname,
				       "l", vorder,
				       "x", fill_argument(std::make_pair(-Real{1}, Real{1}),
							  std::make_pair(true, true), 21),
				       file_legendre);

    //  Riemann zeta function.
    std::cout << "riemann_zeta" << std::endl;
    //  Skip the pole at 1.
    funcname = "riemann_zeta";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_riemann_zeta(filename.c_str());
    typedef Real riemann_zeta(Real);
    test =
    maketest<Real, Real>((riemann_zeta*)std::riemann_zeta,
			 wrap_gsl_sf_zeta,
			 nsname, funcname,
			 "s", fill_argument(std::make_pair(-Real{10}, Real{1}),
					    std::make_pair(true, false), 56),
			 file_riemann_zeta, true, false);
    maketest<Real, Real>((riemann_zeta*)std::riemann_zeta,
			 wrap_gsl_sf_zeta,
			 nsname, funcname,
			 "s", fill_argument(std::make_pair(Real{1}, Real{30}),
					    std::make_pair(false, true), 146),
			 file_riemann_zeta, false, true, test);

    //  Hurwitz zeta function.
    std::cout << "hurwitz_zeta" << std::endl;
    //  Skip the pole at 1.
    funcname = "hurwitz_zeta";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_hurwitz_zeta(filename.c_str());
    typedef Real hurwitz_zeta(Real, Real);
    maketest<Real, Real, Real>((hurwitz_zeta*)__gnu_cxx::hurwitz_zeta,
			       wrap_gsl_sf_hzeta,
			       "__gnu_cxx", funcname,
			       "s", fill_argument(std::make_pair(Real{1}, Real{30}),
						  std::make_pair(false, true), 146),
			       "a", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(false, true), 26),
			       file_hurwitz_zeta);

    //  Spherical Bessel functions.
    std::cout << "sph_bessel" << std::endl;
    funcname = "sph_bessel";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_sph_bessel(filename.c_str());
    typedef Real sph_bessel(unsigned int, Real);
    test =
    maketest<Real, unsigned int, Real>((sph_bessel*)std::sph_bessel,
				       wrap_gsl_sf_bessel_jl,
				       nsname, funcname,
				       "n", sborder,
				       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
							  std::make_pair(true, true), 21),
				       file_sph_bessel, true, false);
    maketest<Real, unsigned int, Real>((sph_bessel*)std::sph_bessel,
				       wrap_gsl_sf_bessel_jl,
				       nsname, funcname,
				       "n", sborder,
				       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
							  std::make_pair(true, true), 21),
				       file_sph_bessel, false, true, test);

    //  Spherical Legendre functions.
    std::cout << "sph_legendre" << std::endl;
    funcname = "sph_legendre";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_sph_legendre(filename.c_str());
    typedef Real sph_legendre(unsigned int, unsigned int, Real);
    maketest<Real, unsigned int, unsigned int, Real>((sph_legendre*)std::sph_legendre, wrap_gsl_sf_legendre_sphPlm,
						     nsname, funcname,
						     "l", vorder, "m", vorder,
						     "theta", fill_argument(std::make_pair(Real{0}, static_cast<Real>(M_PI)),
									    std::make_pair(true, true), 21),
						     file_sph_legendre);

    //  Spherical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "sph_neumann" << std::endl;
    funcname = "sph_neumann";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_sph_neumann(filename.c_str());
    typedef Real sph_neumann(unsigned int, Real);
    test =
    maketest<Real, unsigned int, Real>((sph_neumann*)std::sph_neumann, wrap_gsl_sf_bessel_yl,
				       nsname, funcname,
				       "n", sborder,
				       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
							  std::make_pair(false, true), 21),
				       file_sph_neumann, true, false);
    maketest<Real, unsigned int, Real>((sph_neumann*)std::sph_neumann, wrap_gsl_sf_bessel_yl,
				       nsname, funcname,
				       "n", sborder,
				       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
							  std::make_pair(false, true), 21),
				       file_sph_neumann, false, true, test);

    //  Carlson elliptic functions R_C.
    std::cout << "ellint_rc" << std::endl;
    funcname = "ellint_rc";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_rc(filename.c_str());
    typedef Real ellint_rc(Real, Real);
    maketest<Real, Real, Real>((ellint_rc*)__gnu_cxx::ellint_rc, wrap_gsl_sf_ellint_RC,
			       "__gnu_cxx", funcname,
			       "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
						  std::make_pair(false, true), 11),
			       "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
						  std::make_pair(false, true), 11),
			       file_ellint_rc);

    //  Carlson elliptic functions R_D.
    std::cout << "ellint_rd" << std::endl;
    funcname = "ellint_rd";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_rd(filename.c_str());
    typedef Real ellint_rd(Real, Real, Real);
    maketest<Real, Real, Real, Real>((ellint_rd*)__gnu_cxx::ellint_rd, wrap_gsl_sf_ellint_RD,
			       "__gnu_cxx", funcname,
			       "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
						  std::make_pair(false, true), 11),
			       "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
						  std::make_pair(true, true), 11),
			       "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
						  std::make_pair(true, true), 11),
			       file_ellint_rd);

    //  Carlson elliptic functions R_F.
    std::cout << "ellint_rf" << std::endl;
    funcname = "ellint_rf";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_rf(filename.c_str());
    typedef Real ellint_rf(Real, Real, Real);
    maketest<Real, Real, Real, Real>((ellint_rf*)__gnu_cxx::ellint_rf, wrap_gsl_sf_ellint_RF,
			       "__gnu_cxx", funcname,
			       "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
						  std::make_pair(false, true), 11),
			       "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
						  std::make_pair(true, true), 11),
			       "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
						  std::make_pair(true, true), 11),
			       file_ellint_rf);

    //  Carlson elliptic functions R_J.
    std::cout << "ellint_rj" << std::endl;
    funcname = "ellint_rj";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_rj(filename.c_str());
    typedef Real ellint_rj(Real, Real, Real, Real);
    maketest<Real, Real, Real, Real, Real>((ellint_rj*)__gnu_cxx::ellint_rj, wrap_gsl_sf_ellint_RJ,
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

    //  Dilogarithm functions.
    std::cout << "dilog" << std::endl;
    funcname = "dilog";
    filename = get_filename(path, prefix, "dilog", "",  ".cc");
    std::ofstream file_dilog(filename.c_str());
    typedef Real dilog(Real);
    maketest<Real, Real>((dilog *)__gnu_cxx::dilog,
			 wrap_gsl_sf_dilog,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(-Real{10}, Real{10}),
					    std::make_pair(true, true), 41),
			 file_dilog);

    //  Upper incomplete Gamma functions.
    std::cout << "gamma_u" << std::endl;
    funcname = "gamma_u";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_gamma_u(filename.c_str());
    typedef Real gamma_u(Real, Real);
    maketest<Real, Real, Real>((gamma_u*)__gnu_cxx::gamma_u, wrap_gsl_sf_gamma_inc,
			       "__gnu_cxx", funcname,
			       "a", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(false, true), 11),
			       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(true, true), 11),
			       file_gamma_u);

    //  Incomplete Beta functions.
    std::cout << "ibeta" << std::endl;
    funcname = "ibeta";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ibeta(filename.c_str());
    typedef Real ibeta(Real, Real, Real);
    maketest<Real, Real, Real, Real>((ibeta*)__gnu_cxx::ibeta, wrap_gsl_sf_beta_inc,
			       "__gnu_cxx", funcname,
			       "a", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(false, true), 11),
			       "b", fill_argument(std::make_pair(Real{5}, Real{0}),
						  std::make_pair(true, false), 11),
			       "x", fill_argument(std::make_pair(Real{0}, Real{1}),
						  std::make_pair(false, false), 21),
			       file_ibeta);

    //  Digamma or psi functions.
    std::cout << "psi" << std::endl;
    funcname = "psi";
    filename = get_filename(path, prefix, "psi", "",  ".cc");
    std::ofstream file_psi(filename.c_str());
    typedef Real psi(Real);
    maketest<Real, Real>((psi *)__gnu_cxx::psi, wrap_gsl_sf_psi,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(-Real{9.875}, +Real{10.125}),
					    std::make_pair(true, true), 401),
			 file_psi);

  }


int
main()
{
  harness<double>();

  return 0;
}


