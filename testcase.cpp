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


std::string get_filename(const std::string & path,
                         const std::string & prefix,
                         const std::string & funcname,
                         const std::string & extra,
                         const std::string & suffix )
{
  std::string filename = path + "/" + prefix;
  filename += funcname + extra + suffix;

  return filename;
}

int
main()
{
  //  Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
  std::vector<unsigned int> vorder{ 0, 1, 2, 5, 10, 20, 50, 100 };

  //  ... corresponding "double" integer orders for GSL.
  std::vector<double> dvorder{ 0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0 };

  //  Orders for cylindrical Bessel functions.
  std::vector<double> vborderd{ 0.0, 1.0/3.0, 0.5, 2.0/3.0, 1.0, 2.0,
				5.0, 10.0, 20.0, 50.0, 100.0 };

  //  Orders for spherical bessel functions.
  std::vector<unsigned int> sborder{ 0, 1, 2, 5, 10, 20, 50, 100 };

  const unsigned long num_phi = 10;
  double phi[num_phi];
  for (unsigned int i = 0; i < num_phi; ++i)
    phi[i] = 10.0 * i * M_PI / 180.0;
  std::vector<double> vphid(phi, phi + num_phi);

  std::vector<double> vab{ 0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0 };

  unsigned int test = 1;

  const std::string path = ".";
  const std::string prefix = "/check_";

  std::string funcname;
  std::string filename;

  //  Airy functions.
  std::cout << "airy_ai" << std::endl;
  funcname = "airy_ai";
  filename = get_filename(path, prefix, "airy", "",  ".cc");
  std::ofstream file_airy(filename.c_str());
  typedef double airy(double);
  maketest<double, double>((airy *)__gnu_cxx::airy_ai,
                   	   wrap_gsl_sf_airy_ai,
                   	   "__gnu_cxx", funcname,
                   	   "x", fill_argument(std::make_pair(-10.0, 10.0),
                   			      std::make_pair(true, true), 41),
                   	   file_airy);

  std::cout << "airy_bi" << std::endl;
  funcname = "airy_bi";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  maketest<double, double>((airy *)__gnu_cxx::airy_bi,
                   	   wrap_gsl_sf_airy_bi,
                   	   "__gnu_cxx", funcname,
                   	   "x", fill_argument(std::make_pair(-10.0, 10.0),
                   			      std::make_pair(true, true), 41),
                   	   file_airy);

  //  5.2.1.1  Associated Laguerre polynomials.
  std::cout << "assoc_laguerre" << std::endl;
  funcname = "assoc_laguerre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_assoc_laguerre(filename.c_str());
  typedef double assoc_laguerre(unsigned int, unsigned int, double);
  maketest<double, unsigned int, unsigned int, double>((assoc_laguerre *)std::assoc_laguerre,
                                                       wrap_gsl_sf_laguerre_nm,
                                                       "std", funcname,
                                                       "n", vorder, "m", vorder,
                                                       "x", fill_argument(std::make_pair(0.0, 100.0),
                                                                          std::make_pair(true, true), 11),
                                                       file_assoc_laguerre);

  //  5.2.1.2  Associated Legendre functions.
  std::cout << "assoc_legendre" << std::endl;
  funcname = "assoc_legendre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_assoc_legendre(filename.c_str());
  typedef double assoc_legendre(unsigned int, unsigned int, double);
  maketest<double, unsigned int, unsigned int, double>((assoc_legendre *)std::assoc_legendre,
                                                       wrap_gsl_sf_legendre_Plm,
                                                       "std", funcname,
                                                       "l", vorder, "m", vorder,
                                                       "x", fill_argument(std::make_pair(-1.0, 1.0),
                                                                          std::make_pair(true, true), 21),
                                                       file_assoc_legendre);

  //  5.2.1.3  Beta function.
  std::cout << "beta" << std::endl;
  funcname = "beta";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_beta(filename.c_str());
  typedef double beta(double, double);
  maketest<double, double, double>((beta *)std::beta,
                                   wrap_gsl_sf_beta,
                                   "std", funcname,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(false, true), 11),
                                   "y", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(false, true), 11),
                                   file_beta);

  //  5.2.1.4  Complete elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "comp_ellint_1" << std::endl;
  funcname = "comp_ellint_1";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_comp_ellint_1(filename.c_str());
  typedef double comp_ellint_1(double);
  maketest<double, double>((comp_ellint_1 *)std::comp_ellint_1,
                           wrap_gsl_sf_ellint_Kcomp,
                           "std", funcname,
                           "k", fill_argument(std::make_pair(-1.0, 1.0),
                                              std::make_pair(false, false), 21),
                           file_comp_ellint_1);

  //  5.2.1.5  Complete elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "comp_ellint_2" << std::endl;
  funcname = "comp_ellint_2";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_comp_ellint_2(filename.c_str());
  typedef double comp_ellint_2(double);
  maketest<double, double>((comp_ellint_2 *)std::comp_ellint_2,
                           wrap_gsl_sf_ellint_Ecomp,
                           "std", funcname,
                           "k", fill_argument(std::make_pair(-1.0, 1.0),
                                              std::make_pair(false, false), 21),
                           file_comp_ellint_2);

  //  5.2.1.6  Complete elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "comp_ellint_3" << std::endl;
  funcname = "comp_ellint_3";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_comp_ellint_3(filename.c_str());
  typedef double comp_ellint_3(double, double);
  maketest<double, double, double>((comp_ellint_3 *)std::comp_ellint_3,
                                   wrap_gsl_sf_ellint_Pcomp,
                                   "std", funcname,
                                   "k", fill_argument(std::make_pair(-1.0, 1.0),
                                                      std::make_pair(false, false), 21),
                                   "nu", fill_argument(std::make_pair(0.0, 1.0),
                                                       std::make_pair(true, false), 11),
                                   file_comp_ellint_3);

  //  5.2.1.7  Confluent hypergeometric functions.
  //  Skip the singularity aat c = 0.
  std::cout << "conf_hyperg" << std::endl;
  funcname = "conf_hyperg";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_conf_hyperg(filename.c_str());
  typedef double conf_hyperg(double, double, double);
  maketest<double, double, double, double>((conf_hyperg *)__gnu_cxx::conf_hyperg,
                                           wrap_gsl_sf_hyperg_1F1,
                                           "__gnu_cxx", funcname,
                                           "a", vab,
                                           "c", fill_argument(std::make_pair(0.0, 10.0),
                                                              std::make_pair(false, true), 11),
                                           "x", fill_argument(std::make_pair(-10.0, 10.0),
                                                              std::make_pair(true, true), 21),
                                           file_conf_hyperg);

  //  5.2.1.8  Regular modified cylindrical Bessel functions.
  std::cout << "cyl_bessel_i" << std::endl;
  funcname = "cyl_bessel_i";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_bessel_i(filename.c_str());
  typedef double cyl_bessel_i(double, double);
  test = 
  maketest<double, double, double>((cyl_bessel_i*)std::cyl_bessel_i,
                                   wrap_gsl_sf_bessel_Inu,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 5.0),
                                                      std::make_pair(true, true), 21),
                                   file_cyl_bessel_i, true, false);
  maketest<double, double, double>((cyl_bessel_i*)std::cyl_bessel_i,
                                   wrap_gsl_sf_bessel_Inu,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(true, true), 21),
                                   file_cyl_bessel_i, false, true, test);

  //  5.2.1.9  Cylindrical Bessel functions (of the first kind).
  std::cout << "cyl_bessel_j" << std::endl;
  funcname = "cyl_bessel_j";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_bessel_j(filename.c_str());
  typedef double cyl_bessel_j(double, double);
  test = 
  maketest<double, double, double>((cyl_bessel_j*)std::cyl_bessel_j,
                                   wrap_gsl_sf_bessel_Jnu,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 5.0),
                                                      std::make_pair(true, true), 21),
                                   file_cyl_bessel_j, true, false);
  maketest<double, double, double>((cyl_bessel_j*)std::cyl_bessel_j,
                                   wrap_gsl_sf_bessel_Jnu,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(true, true), 21),
                                   file_cyl_bessel_j, false, true, test);

  //  5.2.1.9  Cylindrical Bessel functions (of the first kind) asymptotics.
/*
  std::cout << "cyl_bessel_j asymptotics" << std::endl;
  funcname = "cyl_bessel_j_asymp";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_bessel_j_asymp(filename.c_str());
  maketest<double, double, double>((cyl_bessel_j*)std::cyl_bessel_j,
                                   wrap_gsl_sf_bessel_Jnu_asymp,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(1000.0, 10000.0),
                                                      std::make_pair(true, true), 11),
                                   file_cyl_bessel_j_asymp);
*/
  //  5.2.1.10  Irregular modified cylindrical Bessel functions.
  // Skip the pole at the origin.
  std::cout << "cyl_bessel_k" << std::endl;
  funcname = "cyl_bessel_k";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_bessel_k(filename.c_str());
  typedef double cyl_bessel_k(double, double);
  test = 
  maketest<double, double, double>((cyl_bessel_k*)std::cyl_bessel_k,
                                   wrap_gsl_sf_bessel_Knu,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 5.0),
                                                      std::make_pair(false, true), 21),
                                   file_cyl_bessel_k, true, false);
  maketest<double, double, double>((cyl_bessel_k*)std::cyl_bessel_k,
                                   wrap_gsl_sf_bessel_Knu,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(false, true), 21),
                                   file_cyl_bessel_k, false, true, test);

  //  5.2.1.11  Cylindrical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "cyl_neumann" << std::endl;
  funcname = "cyl_neumann";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_neumann(filename.c_str());
  typedef double cyl_neumann(double, double);
  test = 
  maketest<double, double, double>((cyl_neumann*)std::cyl_neumann,
                                   wrap_gsl_sf_bessel_Ynu,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 5.0),
                                                      std::make_pair(false, true), 21),
                                   file_cyl_neumann, true, false);
  maketest<double, double, double>((cyl_neumann*)std::cyl_neumann,
                                   wrap_gsl_sf_bessel_Ynu,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(false, true), 21),
                                   file_cyl_neumann, false, true, test);

  //  5.2.1.11  Cylindrical Neumann functions asymptotics.
  // Skip the pole at the origin.
/*
  std::cout << "cyl_neumann asymptotics" << std::endl;
  funcname = "cyl_neumann_asymp";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_neumann_asymp(filename.c_str());
  maketest<double, double, double>((cyl_neumann*)std::cyl_neumann,
                                   wrap_gsl_sf_bessel_Ynu_asymp,
                                   "std", funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(1000.0, 10000.0),
                                                      std::make_pair(false, true), 11),
                                   file_cyl_neumann_asymp);
*/
  //  5.2.1.12  Elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "ellint_1" << std::endl;
  funcname = "ellint_1";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_ellint_1(filename.c_str());
  typedef double ellint_1(double, double);
  maketest<double, double, double>((ellint_1*)std::ellint_1,
                                   wrap_gsl_sf_ellint_F,
                                   "std", funcname,
                                   "k", fill_argument(std::make_pair(-1.0, 1.0),
                                                      std::make_pair(false, false), 21),
                                   "phi", vphid,
                                   file_ellint_1);

  //  5.2.1.13  Elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "ellint_2" << std::endl;
  funcname = "ellint_2";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_ellint_2(filename.c_str());
  typedef double ellint_2(double, double);
  maketest<double, double, double>((ellint_2*)std::ellint_2,
                                   wrap_gsl_sf_ellint_E,
                                   "std", funcname,
                                   "k", fill_argument(std::make_pair(-1.0, 1.0),
                                                      std::make_pair(false, false), 21),
                                   "phi", vphid,
                                   file_ellint_2);

  //  5.2.1.14  Elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "ellint_3" << std::endl;
  funcname = "ellint_3";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_ellint_3(filename.c_str());
  typedef double ellint_3(double, double, double);
  maketest<double, double, double, double>((ellint_3*)std::ellint_3,
                                           wrap_gsl_sf_ellint_P,
                                           "std", funcname,
                                           "k", fill_argument(std::make_pair(-1.0, 1.0),
                                                              std::make_pair(false, false), 21),
                                           "nu", fill_argument(std::make_pair(0.0, 1.0),
                                                               std::make_pair(true, false), 11),
                                           "phi", vphid,
                                           file_ellint_3);

  //  5.2.1.15  Exponential integral.
  //  Skip the pole at 0.
  std::cout << "expint" << std::endl;
  funcname = "expint";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_expint(filename.c_str());
  typedef double expint(double);
  test = 
  maketest<double, double>((expint*)std::expint,
                           wrap_gsl_sf_expint_Ei,
                           "std", funcname,
                           "x", fill_argument(std::make_pair(-50.0, 0.0),
                                              std::make_pair(true, false), 51),
                           file_expint, true, false);
  maketest<double, double>((expint*)std::expint,
                           wrap_gsl_sf_expint_Ei,
                           "std", funcname,
                           "x", fill_argument(std::make_pair(0.0, 50.0),
                                              std::make_pair(false, true), 51),
                           file_expint, false, true, test);

  //  5.2.1.16  Hermite polynomials
  std::cout << "hermite  UNTESTED" << std::endl;

  //  5.2.1.17  Hypergeometric functions.
  //  Skip the singularity at c = 0.
  //  Skip the singularity at x = -1.
  std::cout << "hyperg" << std::endl;
  funcname = "hyperg";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_hyperg(filename.c_str());
  typedef double hyperg(double, double, double, double);
  maketest<double, double, double, double, double>((hyperg*)__gnu_cxx::hyperg,
                                                   wrap_gsl_sf_hyperg_2F1,
                                                   "__gnu_cxx", funcname,
                                                   "a", vab,
                                                   "b", vab,
                                                   "c", fill_argument(std::make_pair(0.0, 10.0),
                                                                      std::make_pair(false, true), 6),
                                                   "x", fill_argument(std::make_pair(-1.0, 1.0),
                                                                      std::make_pair(false, true), 21),
                                                   file_hyperg);

  //  5.2.1.18  Laguerre polynomials.
  std::cout << "laguerre" << std::endl;
  funcname = "laguerre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_laguerre(filename.c_str());
  typedef double laguerre(unsigned int, double);
  maketest<double, unsigned int, double>((laguerre*)std::laguerre,
                                         wrap_gsl_sf_laguerre_n,
                                         "std", funcname,
                                         "n", vorder,
                                         "x", fill_argument(std::make_pair(0.0, 100.0),
                                                            std::make_pair(true, true), 21),
                                         file_laguerre);

  //  5.2.1.19  Legendre polynomials.
  std::cout << "legendre" << std::endl;
  funcname = "legendre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_legendre(filename.c_str());
  typedef double legendre(unsigned int, double);
  maketest<double, unsigned int, double>((legendre*)std::legendre,
                                         wrap_gsl_sf_legendre_Pl,
                                         "std", funcname,
                                         "l", vorder,
                                         "x", fill_argument(std::make_pair(-1.0, 1.0),
                                                            std::make_pair(true, true), 21),
                                         file_legendre);

  //  5.2.1.20  Riemann zeta function.
  std::cout << "riemann_zeta" << std::endl;
  //  Skip the pole at 1.
  funcname = "riemann_zeta";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_riemann_zeta(filename.c_str());
  typedef double riemann_zeta(double);
  test = 
  maketest<double, double>((riemann_zeta*)std::riemann_zeta,
                           wrap_gsl_sf_zeta,
                           "std", funcname,
                           "s", fill_argument(std::make_pair(-10.0, 1.0),
                                              std::make_pair(true, false), 56),
                           file_riemann_zeta, true, false);
  maketest<double, double>((riemann_zeta*)std::riemann_zeta,
                           wrap_gsl_sf_zeta,
                           "std", funcname,
                           "s", fill_argument(std::make_pair(1.0, 30.0),
                                              std::make_pair(false, true), 146),
                           file_riemann_zeta, false, true, test);

  //  Hurwitz zeta function.
  std::cout << "hurwitz_zeta" << std::endl;
  //  Skip the pole at 1.
  funcname = "hurwitz_zeta";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_hurwitz_zeta(filename.c_str());
  typedef double hurwitz_zeta(double, double);
  maketest<double, double, double>((hurwitz_zeta*)__gnu_cxx::hurwitz_zeta,
                        	   wrap_gsl_sf_hzeta,
                        	   "__gnu_cxx", funcname,
                        	   "s", fill_argument(std::make_pair(1.0, 30.0),
                                        	      std::make_pair(false, true), 146),
                        	   "a", fill_argument(std::make_pair(0.0, 5.0),
                                        	      std::make_pair(false, true), 26),
                        	   file_hurwitz_zeta);

  //  5.2.1.21  Spherical Bessel functions.
  std::cout << "sph_bessel" << std::endl;
  funcname = "sph_bessel";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_sph_bessel(filename.c_str());
  typedef double sph_bessel(unsigned int, double);
  test = 
  maketest<double, unsigned int, double>((sph_bessel*)std::sph_bessel,
                                         wrap_gsl_sf_bessel_jl,
                                         "std", funcname,
                                         "n", sborder,
                                         "x", fill_argument(std::make_pair(0.0, 5.0),
                                                            std::make_pair(true, true), 21),
                                         file_sph_bessel, true, false);
  maketest<double, unsigned int, double>((sph_bessel*)std::sph_bessel,
                                         wrap_gsl_sf_bessel_jl,
                                         "std", funcname,
                                         "n", sborder,
                                         "x", fill_argument(std::make_pair(0.0, 100.0),
                                                            std::make_pair(true, true), 21),
                                         file_sph_bessel, false, true, test);

  //  5.2.1.22  Spherical Legendre functions.
  std::cout << "sph_legendre" << std::endl;
  funcname = "sph_legendre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_sph_legendre(filename.c_str());
  typedef double sph_legendre(unsigned int, unsigned int, double);
  maketest<double, unsigned int, unsigned int, double>((sph_legendre*)std::sph_legendre, wrap_gsl_sf_legendre_sphPlm,
                                                       "std", funcname,
                                                       "l", vorder, "m", vorder,
                                                       "theta", fill_argument(std::make_pair(0.0, static_cast<double>(M_PI)),
                                                                              std::make_pair(true, true), 21),
                                                       file_sph_legendre);

  //  5.2.1.23  Spherical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "sph_neumann" << std::endl;
  funcname = "sph_neumann";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_sph_neumann(filename.c_str());
  typedef double sph_neumann(unsigned int, double);
  test = 
  maketest<double, unsigned int, double>((sph_neumann*)std::sph_neumann, wrap_gsl_sf_bessel_yl,
                                         "std", funcname,
                                         "n", sborder,
                                         "x", fill_argument(std::make_pair(0.0, 5.0),
                                                            std::make_pair(false, true), 21),
                                         file_sph_neumann, true, false);
  maketest<double, unsigned int, double>((sph_neumann*)std::sph_neumann, wrap_gsl_sf_bessel_yl,
                                         "std", funcname,
                                         "n", sborder,
                                         "x", fill_argument(std::make_pair(0.0, 100.0),
                                                            std::make_pair(false, true), 21),
                                         file_sph_neumann, false, true, test);

  return 0;
}


