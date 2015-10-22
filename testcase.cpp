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

#include "testcase.h"
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
  unsigned int order[] = { 0, 1, 2, 5, 10, 20, 50, 100 };
  const unsigned int num_order = sizeof(order) / sizeof(unsigned int);
  std::vector<unsigned int> vorder(order, order + num_order);

  //  ... corresponding "double" integer orders for GSL.
  double dorder[] = { 0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0 };
  const unsigned long num_dorder = sizeof(dorder) / sizeof(double);
  std::vector<double> dvorder(dorder, dorder + num_dorder);

  //  Orders for cylindrical Bessel functions.
  double border[] = { 0.0, 1.0/3.0, 0.5, 2.0/3.0, 1.0, 2.0,
                      5.0, 10.0, 20.0, 50.0, 100.0 };
  const unsigned long num_border = sizeof(border) / sizeof(double);
  std::vector<double> vborderd(border, border + num_border);

  //  Orders for spherical bessel functions.
  std::vector<unsigned int> sborder(order, order + num_order);

  const unsigned long num_phi = 10;
  double phi[num_phi];
  for (unsigned int i = 0; i < num_phi; ++i)
    phi[i] = 10.0 * i * M_PI / 180.0;
  std::vector<double> vphid(phi, phi + num_phi);

  double ab[] = { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0 };
  std::vector<double> vab( ab, ab + sizeof(ab) / sizeof(double) );

  const std::string path = ".";
  const std::string prefix = "/check_";

  std::string funcname;
  std::string filename;

  //  5.2.1.1  Associated Laguerre polynomials.
  std::cout << "assoc_laguerre\n";
  funcname = "assoc_laguerre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_assoc_laguerre(filename.c_str());
  typedef double assoc_laguerre(unsigned int, unsigned int, double);
  maketest<double, unsigned int, unsigned int, double>((assoc_laguerre *)std::assoc_laguerre,
                                                       wrap_gsl_sf_laguerre_nm,
                                                       funcname,
                                                       "n", vorder, "m", vorder,
                                                       "x", fill_argument(std::make_pair(0.0, 100.0),
                                                                          std::make_pair(true, true), 11),
                                                       file_assoc_laguerre);

  //  5.2.1.2  Associated Legendre functions.
  std::cout << "assoc_legendre\n";
  funcname = "assoc_legendre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_assoc_legendre(filename.c_str());
  typedef double assoc_legendre(unsigned int, unsigned int, double);
  maketest<double, unsigned int, unsigned int, double>((assoc_legendre *)std::assoc_legendre,
                                                       wrap_gsl_sf_legendre_Plm,
                                                       funcname,
                                                       "l", vorder, "m", vorder,
                                                       "x", fill_argument(std::make_pair(-1.0, 1.0),
                                                                          std::make_pair(true, true), 21),
                                                       file_assoc_legendre);

  //  5.2.1.3  Beta function.
  std::cout << "beta\n";
  funcname = "beta";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_beta(filename.c_str());
  typedef double beta(double, double);
  maketest<double, double, double>((beta *)std::beta,
                                   wrap_gsl_sf_beta,
                                   funcname,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(false, true), 11),
                                   "y", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(false, true), 11),
                                   file_beta);

  //  5.2.1.4  Complete elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "comp_ellint_1\n";
  funcname = "comp_ellint_1";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_comp_ellint_1(filename.c_str());
  typedef double comp_ellint_1(double);
  maketest<double, double>((comp_ellint_1 *)std::comp_ellint_1,
                           wrap_gsl_sf_ellint_Kcomp,
                           funcname,
                           "k", fill_argument(std::make_pair(-1.0, 1.0),
                                              std::make_pair(false, false), 21),
                           file_comp_ellint_1);

  //  5.2.1.5  Complete elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "comp_ellint_2\n";
  funcname = "comp_ellint_2";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_comp_ellint_2(filename.c_str());
  typedef double comp_ellint_2(double);
  maketest<double, double>((comp_ellint_2 *)std::comp_ellint_2,
                           wrap_gsl_sf_ellint_Ecomp,
                           funcname,
                           "k", fill_argument(std::make_pair(-1.0, 1.0),
                                              std::make_pair(false, false), 21),
                           file_comp_ellint_2);

  //  5.2.1.6  Complete elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "comp_ellint_3\n";
  funcname = "comp_ellint_3";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_comp_ellint_3(filename.c_str());
  typedef double comp_ellint_3(double, double);
  maketest<double, double, double>((comp_ellint_3 *)std::comp_ellint_3,
                                   wrap_gsl_sf_ellint_Pcomp,
                                   funcname,
                                   "k", fill_argument(std::make_pair(-1.0, 1.0),
                                                      std::make_pair(false, false), 21),
                                   "nu", fill_argument(std::make_pair(0.0, 1.0),
                                                       std::make_pair(true, false), 11),
                                   file_comp_ellint_3);

  //  5.2.1.7  Confluent hypergeometric functions.
  //  Skip the singularity aat c = 0.
/*
  std::cout << "conf_hyperg\n";
  funcname = "conf_hyperg";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_conf_hyperg(filename.c_str());
  typedef double conf_hyperg(double, double, double);
  maketest<double, double, double, double>((conf_hyperg *)std::conf_hyperg,
                                           wrap_gsl_sf_hyperg_1F1,
                                           funcname,
                                           "a", vab,
                                           "c", fill_argument(std::make_pair(0.0, 10.0),
                                                              std::make_pair(false, true), 11),
                                           "x", fill_argument(std::make_pair(-10.0, 10.0),
                                                              std::make_pair(true, true), 21),
                                           file_conf_hyperg);
*/

  //  5.2.1.8  Regular modified cylindrical Bessel functions.
  std::cout << "cyl_bessel_i\n";
  funcname = "cyl_bessel_i";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_bessel_i(filename.c_str());
  typedef double cyl_bessel_i(double, double);
  maketest<double, double, double>((cyl_bessel_i*)std::cyl_bessel_i,
                                   wrap_gsl_sf_bessel_Inu,
                                   funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(true, true), 21),
                                   file_cyl_bessel_i);

  //  5.2.1.9  Cylindrical Bessel functions (of the first kind).
  std::cout << "cyl_bessel_j\n";
  funcname = "cyl_bessel_j";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_bessel_j(filename.c_str());
  typedef double cyl_bessel_j(double, double);
  maketest<double, double, double>((cyl_bessel_j*)std::cyl_bessel_j,
                                   wrap_gsl_sf_bessel_Jnu,
                                   funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(true, true), 21),
                                   file_cyl_bessel_j);

  //  5.2.1.9  Cylindrical Bessel functions (of the first kind) asymptotics.
  std::cout << "cyl_bessel_j asymptotics\n";
  funcname = "cyl_bessel_j_asymp";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_bessel_j_asymp(filename.c_str());
  maketest<double, double, double>((cyl_bessel_j*)std::cyl_bessel_j,
                                   wrap_gsl_sf_bessel_Jnu_asymp,
                                   funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(1000.0, 10000.0),
                                                      std::make_pair(true, true), 11),
                                   file_cyl_bessel_j_asymp);

  //  5.2.1.10  Irregular modified cylindrical Bessel functions.
  // Skip the pole at the origin.
  std::cout << "cyl_bessel_k\n";
  funcname = "cyl_bessel_k";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_bessel_k(filename.c_str());
  typedef double cyl_bessel_k(double, double);
  maketest<double, double, double>((cyl_bessel_k*)std::cyl_bessel_k,
                                   wrap_gsl_sf_bessel_Knu,
                                   funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(false, true), 21),
                                   file_cyl_bessel_k);

  //  5.2.1.11  Cylindrical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "cyl_neumann\n";
  funcname = "cyl_neumann";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_neumann(filename.c_str());
  typedef double cyl_neumann(double, double);
  maketest<double, double, double>((cyl_neumann*)std::cyl_neumann,
                                   wrap_gsl_sf_bessel_Ynu,
                                   funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(false, true), 21),
                                   file_cyl_neumann);

  //  5.2.1.11  Cylindrical Neumann functions asymptotics.
  // Skip the pole at the origin.
  std::cout << "cyl_neumann asymptotics\n";
  funcname = "cyl_neumann_asymp";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_cyl_neumann_asymp(filename.c_str());
  maketest<double, double, double>((cyl_neumann*)std::cyl_neumann,
                                   wrap_gsl_sf_bessel_Ynu_asymp,
                                   funcname,
                                   "nu", vborderd,
                                   "x", fill_argument(std::make_pair(1000.0, 10000.0),
                                                      std::make_pair(false, true), 11),
                                   file_cyl_neumann_asymp);

  //  5.2.1.12  Elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "ellint_1\n";
  funcname = "ellint_1";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_ellint_1(filename.c_str());
  typedef double ellint_1(double, double);
  maketest<double, double, double>((ellint_1*)std::ellint_1,
                                   wrap_gsl_sf_ellint_F,
                                   funcname,
                                   "k", fill_argument(std::make_pair(-1.0, 1.0),
                                                      std::make_pair(false, false), 21),
                                   "phi", vphid,
                                   file_ellint_1);

  //  5.2.1.13  Elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "ellint_2\n";
  funcname = "ellint_2";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_ellint_2(filename.c_str());
  typedef double ellint_2(double, double);
  maketest<double, double, double>((ellint_2*)std::ellint_2,
                                   wrap_gsl_sf_ellint_E,
                                   funcname,
                                   "k", fill_argument(std::make_pair(-1.0, 1.0),
                                                      std::make_pair(false, false), 21),
                                   "phi", vphid,
                                   file_ellint_2);

  //  5.2.1.14  Elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "ellint_3\n";
  funcname = "ellint_3";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_ellint_3(filename.c_str());
  typedef double ellint_3(double, double, double);
  maketest<double, double, double, double>((ellint_3*)std::ellint_3,
                                           wrap_gsl_sf_ellint_P,
                                           funcname,
                                           "k", fill_argument(std::make_pair(-1.0, 1.0),
                                                              std::make_pair(false, false), 21),
                                           "nu", fill_argument(std::make_pair(0.0, 1.0),
                                                               std::make_pair(true, false), 11),
                                           "phi", vphid,
                                           file_ellint_3);

  //  5.2.1.15  Exponential integral.
  //  Skip the pole at 0.
  std::cout << "expint\n";
  funcname = "expint";
  filename = get_filename(path, prefix, funcname, "_neg",  ".cc");
  std::ofstream file_expint_neg(filename.c_str());
  typedef double expint(double);
  maketest<double, double>((expint*)std::expint,
                           wrap_gsl_sf_expint_Ei,
                           funcname,
                           "x", fill_argument(std::make_pair(-50.0, 0.0),
                                              std::make_pair(true, false), 51),
                           file_expint_neg);
  funcname = "expint";
  filename = get_filename(path, prefix, funcname, "_pos",  ".cc");
  std::ofstream file_expint_pos(filename.c_str());
  maketest<double, double>((expint*)std::expint,
                           wrap_gsl_sf_expint_Ei,
                           funcname,
                           "x", fill_argument(std::make_pair(0.0, 50.0),
                                              std::make_pair(false, true), 51),
                           file_expint_pos);

  //  5.2.1.16  Hermite polynomials
  std::cout << "hermite  UNTESTED\n";

/*
  //  5.2.1.17  Hypergeometric functions.
  //  Skip the singularity at c = 0.
  //  Skip the singularity at x = -1.
  std::cout << "5.2.1.17  hyperg\n";
  funcname = "hyperg";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_hyperg(filename.c_str());
  typedef double hyperg(double, double, double, double);
  maketest<double, double, double, double, double>((hyperg*)std::hyperg,
                                                   wrap_gsl_sf_hyperg_2F1,
                                                   funcname,
                                                   "a", vab,
                                                   "b", vab,
                                                   "c", fill_argument(std::make_pair(0.0, 10.0),
                                                                      std::make_pair(false, true), 6),
                                                   "x", fill_argument(std::make_pair(-1.0, 1.0),
                                                                      std::make_pair(false, true), 21),
                                                   file_hyperg);
*/

  //  5.2.1.18  Laguerre polynomials.
  std::cout << "laguerre\n";
  funcname = "laguerre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_laguerre(filename.c_str());
  typedef double laguerre(unsigned int, double);
  maketest<double, unsigned int, double>((laguerre*)std::laguerre,
                                         wrap_gsl_sf_laguerre_n,
                                         funcname,
                                         "n", vorder,
                                         "x", fill_argument(std::make_pair(0.0, 100.0),
                                                            std::make_pair(true, true), 21),
                                         file_laguerre);

  //  5.2.1.19  Legendre polynomials.
  std::cout << "legendre\n";
  funcname = "legendre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_legendre(filename.c_str());
  typedef double legendre(unsigned int, double);
  maketest<double, unsigned int, double>((legendre*)std::legendre,
                                         wrap_gsl_sf_legendre_Pl,
                                         funcname,
                                         "l", vorder,
                                         "x", fill_argument(std::make_pair(-1.0, 1.0),
                                                            std::make_pair(true, true), 21),
                                         file_legendre);

  //  5.2.1.20  Riemann zeta function.
  std::cout << "riemann_zeta\n";
  //  Skip the pole at 1.
  funcname = "riemann_zeta";
  filename = get_filename(path, prefix, funcname, "_neg",  ".cc");
  std::ofstream file_riemann_zeta_neg(filename.c_str());
  typedef double riemann_zeta(double);
  maketest<double, double>((riemann_zeta*)std::riemann_zeta,
                           wrap_gsl_sf_zeta,
                           funcname,
                           "x", fill_argument(std::make_pair(-10.0, 1.0),
                                              std::make_pair(true, false), 56),
                           file_riemann_zeta_neg);
  funcname = "riemann_zeta";
  filename = get_filename(path, prefix, funcname, "_pos",  ".cc");
  std::ofstream file_riemann_zeta_pos(filename.c_str());
  maketest<double, double>((riemann_zeta*)std::riemann_zeta,
                           wrap_gsl_sf_zeta,
                           funcname,
                           "x", fill_argument(std::make_pair(1.0, 30.0),
                                              std::make_pair(false, true), 146),
                           file_riemann_zeta_pos);

  //  5.2.1.21  Spherical Bessel functions.
  std::cout << "sph_bessel\n";
  funcname = "sph_bessel";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_sph_bessel(filename.c_str());
  typedef double sph_bessel(unsigned int, double);
  maketest<double, unsigned int, double>((sph_bessel*)std::sph_bessel,
                                         wrap_gsl_sf_bessel_jl,
                                         funcname,
                                         "n", sborder,
                                         "x", fill_argument(std::make_pair(0.0, 100.0),
                                                            std::make_pair(true, true), 21),
                                         file_sph_bessel);

  //  5.2.1.22  Spherical Legendre functions.
  std::cout << "sph_legendre\n";
  funcname = "sph_legendre";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_sph_legendre(filename.c_str());
  typedef double sph_legendre(unsigned int, unsigned int, double);
  maketest<double, unsigned int, unsigned int, double>((sph_legendre*)std::sph_legendre, wrap_gsl_sf_legendre_sphPlm,
                                                       funcname,
                                                       "l", vorder, "m", vorder,
                                                       "theta", fill_argument(std::make_pair(0.0, static_cast<double>(M_PI)),
                                                                              std::make_pair(true, true), 21),
                                                       file_sph_legendre);

  //  5.2.1.23  Spherical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "sph_neumann\n";
  funcname = "sph_neumann";
  filename = get_filename(path, prefix, funcname, "",  ".cc");
  std::ofstream file_sph_neumann(filename.c_str());
  typedef double sph_neumann(unsigned int, double);
  maketest<double, unsigned int, double>((sph_neumann*)std::sph_neumann, wrap_gsl_sf_bessel_yl,
                                         funcname,
                                         "n", sborder,
                                         "x", fill_argument(std::make_pair(0.0, 100.0),
                                                            std::make_pair(false, true), 21),
                                         file_sph_neumann);

  return 0;
}


