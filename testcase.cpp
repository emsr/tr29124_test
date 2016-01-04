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
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_airy_ai(filename.c_str());
    typedef Real airy(Real);
    maketest<Real, Real>((airy *)__gnu_cxx::airy_ai,
			 gsl::airy_ai,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
					    std::make_pair(true, true), 41),
			 file_airy_ai);

    std::cout << "airy_bi" << std::endl;
    funcname = "airy_bi";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_airy_bi(filename.c_str());
    maketest<Real, Real>((airy *)__gnu_cxx::airy_bi,
			 gsl::airy_bi,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
					    std::make_pair(true, true), 41),
			 file_airy_bi);

    //  Associated Laguerre polynomials.
    std::cout << "assoc_laguerre" << std::endl;
    funcname = "assoc_laguerre";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_assoc_laguerre(filename.c_str());
    typedef Real assoc_laguerre(unsigned int, unsigned int, Real);
    maketest<Real, unsigned int, unsigned int, Real>((assoc_laguerre *)std::assoc_laguerre,
						     gsl::laguerre_nm,
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
						     gsl::legendre_Plm,
						     nsname, funcname,
						     "l", vorder, "m", vorder,
						     "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
									std::make_pair(true, true), 21),
						     file_assoc_legendre);

    //  Beta function.
    std::cout << "beta" << std::endl;
    funcname = "beta";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_beta(filename.c_str());
    typedef Real beta(Real, Real);
    maketest<Real, Real, Real>((beta *)std::beta,
			       gsl::beta,
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
			 gsl::ellint_Kcomp,
			 nsname, funcname,
			 "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
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
			 gsl::ellint_Ecomp,
			 nsname, funcname,
			 "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
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
			       gsl::ellint_Pcomp,
			       nsname, funcname,
			       "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
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
				     gsl::hyperg_1F1,
				     "__gnu_cxx", funcname,
				     "a", vab,
				     "c", fill_argument(std::make_pair(Real{0}, Real{10}),
							std::make_pair(false, true), 11),
				     "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
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
			       gsl::bessel_Inu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(true, true), 21),
			       file_cyl_bessel_i, true, false);
    maketest<Real, Real, Real>((cyl_bessel_i*)std::cyl_bessel_i,
			       gsl::bessel_Inu,
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
			       gsl::bessel_Jnu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(true, true), 21),
			       file_cyl_bessel_j, true, false);
    maketest<Real, Real, Real>((cyl_bessel_j*)std::cyl_bessel_j,
			       gsl::bessel_Jnu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
						  std::make_pair(true, true), 21),
			       file_cyl_bessel_j, false, true, test);

    //  Irregular modified cylindrical Bessel functions.
    // Skip the pole at the origin.
    std::cout << "cyl_bessel_k" << std::endl;
    funcname = "cyl_bessel_k";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cyl_bessel_k(filename.c_str());
    typedef Real cyl_bessel_k(Real, Real);
    test =
    maketest<Real, Real, Real>((cyl_bessel_k*)std::cyl_bessel_k,
			       gsl::bessel_Knu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						  std::make_pair(false, true), 21),
			       file_cyl_bessel_k, true, false);
    maketest<Real, Real, Real>((cyl_bessel_k*)std::cyl_bessel_k,
			       gsl::bessel_Knu,
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
				  gsl::bessel_Ynu,
				  nsname, funcname,
				  "nu", vborderd,
				  "x", fill_argument(std::make_pair(Real{0}, Real{5}),
						     std::make_pair(false, true), 21),
				  file_cyl_neumann, true, false);
    maketest<Real, Real, Real>((cyl_neumann*)std::cyl_neumann,
			       gsl::bessel_Ynu,
			       nsname, funcname,
			       "nu", vborderd,
			       "x", fill_argument(std::make_pair(Real{0}, Real{100}),
						  std::make_pair(false, true), 21),
			       file_cyl_neumann, false, true, test);

    //  Elliptic integrals of the first kind.
    //  Avoid poles at |x| = 1.
    std::cout << "ellint_1" << std::endl;
    funcname = "ellint_1";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_ellint_1(filename.c_str());
    typedef Real ellint_1(Real, Real);
    maketest<Real, Real, Real>((ellint_1*)std::ellint_1,
			       gsl::ellint_F,
			       nsname, funcname,
			       "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
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
			       gsl::ellint_E,
			       nsname, funcname,
			       "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
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
				     gsl::ellint_P,
				     nsname, funcname,
				     "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
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
			 gsl::expint_Ei,
			 nsname, funcname,
			 "x", fill_argument(std::make_pair(Real{-50}, Real{0}),
					    std::make_pair(true, false), 51),
			 file_expint, true, false);
    maketest<Real, Real>((expint*)std::expint,
			 gsl::expint_Ei,
			 nsname, funcname,
			 "x", fill_argument(std::make_pair(Real{0}, Real{50}),
					    std::make_pair(false, true), 51),
			 file_expint, false, true, test);

    //  Hermite polynomials
    std::cout << "hermite" << std::endl;
    funcname = "hermite";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_hermite(filename.c_str());
    typedef Real hermite(unsigned int, Real);
    maketest<Real, unsigned int, Real>((hermite*)std::hermite,
				       gsl::hermite,
				       nsname, funcname,
				       "n", vorder,
				       "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
							  std::make_pair(true, true), 201),
				       file_hermite);

    //  Hypergeometric functions.
    //  Skip the singularity at c = 0.
    //  Skip the singularity at x = -1.
    std::cout << "hyperg" << std::endl;
    funcname = "hyperg";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_hyperg(filename.c_str());
    typedef Real hyperg(Real, Real, Real, Real);
    maketest<Real, Real, Real, Real, Real>((hyperg*)__gnu_cxx::hyperg,
					   gsl::hyperg_2F1,
					   "__gnu_cxx", funcname,
					   "a", vab,
					   "b", vab,
					   "c", fill_argument(std::make_pair(Real{0}, Real{10}),
							      std::make_pair(false, true), 6),
					   "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
							      std::make_pair(true, false), 21),
					   file_hyperg);

    //  Laguerre polynomials.
    std::cout << "laguerre" << std::endl;
    funcname = "laguerre";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_laguerre(filename.c_str());
    typedef Real laguerre(unsigned int, Real);
    maketest<Real, unsigned int, Real>((laguerre*)std::laguerre,
				       gsl::laguerre_n,
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
				       gsl::legendre_Pl,
				       nsname, funcname,
				       "l", vorder,
				       "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
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
			 gsl::zeta,
			 nsname, funcname,
			 "s", fill_argument(std::make_pair(Real{-10}, Real{1}),
					    std::make_pair(true, false), 56),
			 file_riemann_zeta, true, false);
    maketest<Real, Real>((riemann_zeta*)std::riemann_zeta,
			 gsl::zeta,
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
			       gsl::hzeta,
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
				       gsl::bessel_jl,
				       nsname, funcname,
				       "n", sborder,
				       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
							  std::make_pair(true, true), 21),
				       file_sph_bessel, true, false);
    maketest<Real, unsigned int, Real>((sph_bessel*)std::sph_bessel,
				       gsl::bessel_jl,
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
    maketest<Real, unsigned int, unsigned int, Real>((sph_legendre*)std::sph_legendre, gsl::legendre_sphPlm,
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
    maketest<Real, unsigned int, Real>((sph_neumann*)std::sph_neumann, gsl::bessel_yl,
				       nsname, funcname,
				       "n", sborder,
				       "x", fill_argument(std::make_pair(Real{0}, Real{5}),
							  std::make_pair(false, true), 21),
				       file_sph_neumann, true, false);
    maketest<Real, unsigned int, Real>((sph_neumann*)std::sph_neumann, gsl::bessel_yl,
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
    maketest<Real, Real, Real>((ellint_rc*)__gnu_cxx::ellint_rc, gsl::ellint_RC,
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
    maketest<Real, Real, Real, Real>((ellint_rd*)__gnu_cxx::ellint_rd, gsl::ellint_RD,
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
    maketest<Real, Real, Real, Real>((ellint_rf*)__gnu_cxx::ellint_rf, gsl::ellint_RF,
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
    maketest<Real, Real, Real, Real, Real>((ellint_rj*)__gnu_cxx::ellint_rj, gsl::ellint_RJ,
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
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_dilog(filename.c_str());
    typedef Real dilog(Real);
    maketest<Real, Real>((dilog *)__gnu_cxx::dilog,
			 gsl::dilog,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
					    std::make_pair(true, true), 41),
			 file_dilog);

    //  Upper incomplete Gamma functions.
    std::cout << "gamma_u" << std::endl;
    funcname = "gamma_u";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_gamma_u(filename.c_str());
    typedef Real gamma_u(Real, Real);
    maketest<Real, Real, Real>((gamma_u*)__gnu_cxx::gamma_u, gsl::gamma_inc,
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
    maketest<Real, Real, Real, Real>((ibeta*)__gnu_cxx::ibeta, gsl::beta_inc,
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
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_psi(filename.c_str());
    typedef Real psi(Real);
    test =
    maketest<Real, Real>((psi *)__gnu_cxx::psi, gsl::psi,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{-9.9375}, +Real{10.0625}),
					    std::make_pair(true, true), 801),
			 file_psi, true, false);
    maketest<Real, Real>((psi *)__gnu_cxx::psi, gsl::psi,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{1}, Real{100}),
					    std::make_pair(true, true), 199),
			 file_psi, false, true, test);

    //  Sine integral or Si functions.
    std::cout << "sinint" << std::endl;
    funcname = "sinint";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_sinint(filename.c_str());
    typedef Real sinint(Real);
    maketest<Real, Real>((sinint *)__gnu_cxx::sinint, gsl::Si,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
					    std::make_pair(false, true), 101),
			 file_sinint);

    //  Cosine integral or Ci functions.
    std::cout << "cosint" << std::endl;
    funcname = "cosint";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_cosint(filename.c_str());
    typedef Real cosint(Real);
    maketest<Real, Real>((cosint *)__gnu_cxx::cosint, gsl::Ci,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{0}, Real{+10}),
					    std::make_pair(false, true), 101),
			 file_cosint);

    //  Hyperbolic sine integral or Shi functions.
    std::cout << "sinhint" << std::endl;
    funcname = "sinhint";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_sinhint(filename.c_str());
    typedef Real sinhint(Real);
    maketest<Real, Real>((sinhint *)__gnu_cxx::sinhint, gsl::Shi,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
					    std::make_pair(false, true), 101),
			 file_sinhint);

    //  Hyperbolic cosine integral or Chi functions.
    std::cout << "coshint" << std::endl;
    funcname = "coshint";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_coshint(filename.c_str());
    typedef Real coshint(Real);
    maketest<Real, Real>((coshint *)__gnu_cxx::coshint, gsl::Chi,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{0}, Real{+5}),
					    std::make_pair(false, true), 101),
			 file_coshint);

    //  Dawson integral.
    std::cout << "dawson" << std::endl;
    funcname = "dawson";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_dawson(filename.c_str());
    typedef Real dawson(Real);
    maketest<Real, Real>((dawson *)__gnu_cxx::dawson, gsl::dawson,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{0}, Real{+20}),
					    std::make_pair(false, true), 201),
			 file_dawson);

    //  Jacobian elliptic integrals.
    std::cout << "jacobi_sn" << std::endl;
    funcname = "jacobi_sn";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_jacobi_sn(filename.c_str());
    typedef Real jacobi_sn(Real, Real);
    maketest<Real, Real, Real>((jacobi_sn *)__gnu_cxx::jacobi_sn, gsl::elljac_sn,
			 "__gnu_cxx", funcname,
			 "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
					    std::make_pair(false, true), 101),
			 "m", fill_argument(std::make_pair(Real{-1}, Real{+1}),
					    std::make_pair(false, true), 21),
			 file_jacobi_sn);

    //  Jacobian elliptic integrals.
    std::cout << "jacobi_cn" << std::endl;
    funcname = "jacobi_cn";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_jacobi_cn(filename.c_str());
    typedef Real jacobi_cn(Real, Real);
    maketest<Real, Real, Real>((jacobi_cn *)__gnu_cxx::jacobi_cn, gsl::elljac_cn,
			 "__gnu_cxx", funcname,
			 "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
					    std::make_pair(false, true), 101),
			 "m", fill_argument(std::make_pair(Real{-1}, Real{+1}),
					    std::make_pair(false, true), 21),
			 file_jacobi_cn);

    //  Jacobian elliptic integrals.
    std::cout << "jacobi_dn" << std::endl;
    funcname = "jacobi_dn";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_jacobi_dn(filename.c_str());
    typedef Real jacobi_dn(Real, Real);
    maketest<Real, Real, Real>((jacobi_dn *)__gnu_cxx::jacobi_dn, gsl::elljac_dn,
			 "__gnu_cxx", funcname,
			 "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
					    std::make_pair(false, true), 101),
			 "m", fill_argument(std::make_pair(Real{-1}, Real{+1}),
					    std::make_pair(false, true), 21),
			 file_jacobi_dn);

    //  Exponential integral E1.
    //  Skip the pole at 0.
    std::cout << "expint_e1" << std::endl;
    funcname = "expint_e1";
    filename = get_filename(path, prefix, funcname, "", ".cc");
    std::ofstream file_expint_e1(filename.c_str());
    typedef Real expint_e1(Real);
    test =
    maketest<Real, Real>((expint_e1*)__gnu_cxx::expint_e1,
			 gsl::expint_E1,
			 nsname, funcname,
			 "x", fill_argument(std::make_pair(Real{-50}, Real{0}),
					    std::make_pair(true, false), 51),
			 file_expint_e1, true, false);
    maketest<Real, Real>((expint_e1*)__gnu_cxx::expint_e1,
			 gsl::expint_E1,
			 nsname, funcname,
			 "x", fill_argument(std::make_pair(Real{0}, Real{50}),
					    std::make_pair(false, true), 51),
			 file_expint_e1, false, true, test);

    //  Fresnel cosine integral.
    std::cout << "fresnel_c" << std::endl;
    funcname = "fresnel_c";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_fresnel_c(filename.c_str());
    typedef Real fresnel_c(Real);
    maketest<Real, Real>((fresnel_c *)__gnu_cxx::fresnel_c, gsl::fresnel_c,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
					    std::make_pair(false, true), 401),
			 file_fresnel_c);

    //  Fresnel sine integral.
    std::cout << "fresnel_s" << std::endl;
    funcname = "fresnel_s";
    filename = get_filename(path, prefix, funcname, "",  ".cc");
    std::ofstream file_fresnel_s(filename.c_str());
    typedef Real fresnel_s(Real);
    maketest<Real, Real>((fresnel_s *)__gnu_cxx::fresnel_s, gsl::fresnel_s,
			 "__gnu_cxx", funcname,
			 "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
					    std::make_pair(false, true), 401),
			 file_fresnel_s);

  }


int
main()
{
  harness<double>();

  return 0;
}


