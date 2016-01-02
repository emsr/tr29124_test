
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <stdexcept>

#define STD 0

#if STD
#  include <cmath>
#  include <array>
#else
#  include <cmath>
#  include <tr1/cmath>
#  include <tr1/array>
#endif

#include "gsl_wrap.h"

#include "test_function.tcc"


///
///
///
int
main()
{

  //  Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
  const unsigned int num_order = 10;
  unsigned int order[num_order] = { 0, 1, 2, 3, 4, 5, 10, 20, 50, 100 };
  std::vector<unsigned int> vorder( order, order + num_order );

  //  ... corresponding "double" integer orders for GSL.
  double dorder[num_order] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 50.0, 100.0 };
  std::vector<double> dvorder( dorder, dorder + num_order );

  //  Orders for cylindrical Bessel functions.
  const unsigned long num_border = 12;
  double border[num_border] = { 0.0, 1.0/3.0, 0.5, 2.0/3.0, 1.0, 2.0, 3.0,
                                5.0, 10.0, 20.0, 50.0, 100.0 };
  std::vector<double> vborderd( border, border + num_border );

  //  Orders for spherical bessel functions.
  std::vector<unsigned int> sborder( order, order + num_order );

  const unsigned long num_nu = 10;
  double nu[num_nu];
  for ( unsigned int i = 0; i < num_nu; ++i )
    nu[i] = 10.0 * i * M_PI / 180.0;
  std::vector<double> vnud( nu, nu + num_nu );

  const unsigned long num_phi = 10;
  double phi[num_phi];
  for ( unsigned int i = 0; i < num_phi; ++i )
    phi[i] = 10.0 * i * M_PI / 180.0;
  std::vector<double> vphid( phi, phi + num_phi );


  std::string basename;

  //  Let's relearn C++ for the third time...
  std::tr1::function<float(unsigned int, unsigned int, float)> assoc_laguerref = std::tr1::assoc_laguerref;
  std::tr1::function<long double(unsigned int, unsigned int, long double)> assoc_laguerrel = std::tr1::assoc_laguerrel;
  //  You need a cast to pick an overload.
  std::tr1::function<double(unsigned int, unsigned int, double)> assoc_laguerred
    = (double(*)(unsigned int, unsigned int, double))std::tr1::assoc_laguerre;
  //  Typedefs work too.
  typedef double foo(unsigned int, unsigned int, double);
  std::tr1::function<foo> assoc_laguerre_foo = (foo *)std::tr1::assoc_laguerre;

  //  5.2.1.1  Associated Laguerre polynomials.
  std::cout << "5.2.1.1  assoc_laguerre" << std::endl;
  basename = "diff_assoc_laguerre";
  typedef double assoc_laguerre(unsigned int, unsigned int, double);
  rundiff<double, unsigned int, unsigned int, double>( (assoc_laguerre *)std::tr1::assoc_laguerre,
                                                       gsl::laguerre_nm, basename,
                                                       "n", vorder, "m", vorder,
                                                       "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                                           std::make_pair( true, true ) ) );


  //  5.2.1.2  Associated Legendre functions.
  std::cout << "5.2.1.2  assoc_legendre" << std::endl;
  basename = "diff_assoc_legendre";
  typedef double assoc_legendre(unsigned int, unsigned int, double);
  rundiff<double, unsigned int, unsigned int, double>( (assoc_legendre *)std::tr1::assoc_legendre,
                                                       gsl::legendre_Plm, basename,
                                                       "l", vorder, "m", vorder,
                                                       "x", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                                           std::make_pair( true, true ), 1001 ) );


  //  5.2.1.3  Beta function.
  std::cout << "5.2.1.3  beta" << std::endl;
  basename = "diff_beta";
  typedef double beta(double, double);
  rundiff<double, double, double>( (beta *)std::tr1::beta,
                                   gsl::beta, basename,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( false, true ) ),
                                   "y", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( false, true ) ) );


  //  5.2.1.4  Complete elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.4  comp_ellint_1" << std::endl;
  basename = "diff_comp_ellint_1";
  typedef double comp_ellint_1(double);
  rundiff<double, double>( (comp_ellint_1 *)std::tr1::comp_ellint_1,
                           gsl::ellint_Kcomp, basename,
                           "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                               std::make_pair( false, false ) ) );


  //  5.2.1.5  Complete elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.5  comp_ellint_2" << std::endl;
  basename = "diff_comp_ellint_2";
  typedef double comp_ellint_2(double);
  rundiff<double, double>( (comp_ellint_2 *)std::tr1::comp_ellint_2,
                           gsl::ellint_Ecomp, basename,
                           "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                               std::make_pair( false, false ) ) );


  //  5.2.1.6  Complete elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.6  comp_ellint_3" << std::endl;
  basename = "diff_comp_ellint_3";
  typedef double comp_ellint_3(double, double);
  rundiff<double, double, double>( (comp_ellint_3 *)std::tr1::comp_ellint_3,
                                   gsl::ellint_Pcomp, basename,
                                   "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                       std::make_pair( false, false ) ),
                                   "nu", vnud );


  //  5.2.1.7  Confluent hypergeometric functions.
  //  Skip the singularity at c = 0.
  std::cout << "5.2.1.7  conf_hyperg" << std::endl;
  basename = "diff_conf_hyperg";
  typedef double conf_hyperg(double, double, double);
  rundiff<double, double, double, double>( (conf_hyperg *)std::tr1::conf_hyperg,
                                           gsl::hyperg_1F1, basename,
                                           "a", fill_argument( std::make_pair( 0.0, 10.0 ),
                                                               std::make_pair( true, true ), 11 ),
                                           "c", fill_argument( std::make_pair( 0.0, 10.0 ),
                                                               std::make_pair( false, true ), 11 ),
                                           "x", fill_argument( std::make_pair( -10.0, 10.0 ),
                                                               std::make_pair( true, true ), 201 ) );


  //  5.2.1.8  Regular modified cylindrical Bessel functions.
  std::cout << "5.2.1.8  cyl_bessel_i" << std::endl;
  basename = "diff_cyl_bessel_i";
  typedef double cyl_bessel_i(double, double);
  rundiff<double, double, double>( (cyl_bessel_i*)std::tr1::cyl_bessel_i,
                                   gsl::bessel_Inu, basename,
                                   "nu", vborderd,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( true, true ), 1001 ) );


  //  5.2.1.9  Cylindrical Bessel functions (of the first kind).
  std::cout << "5.2.1.9  cyl_bessel_j" << std::endl;
  basename = "diff_cyl_bessel_j";
  typedef double cyl_bessel_j(double, double);
  rundiff<double, double, double>( (cyl_bessel_j*)std::tr1::cyl_bessel_j,
                                   gsl::bessel_Jnu, basename,
                                   "nu", vborderd,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( true, true ), 1001 ) );


  //  5.2.1.10  Irregular modified cylindrical Bessel functions.
  // Skip the pole at the origin.
  std::cout << "5.2.1.10  cyl_bessel_k" << std::endl;
  basename = "diff_cyl_bessel_k";
  typedef double cyl_bessel_k(double, double);
  rundiff<double, double, double>( (cyl_bessel_k*)std::tr1::cyl_bessel_k,
                                   gsl::bessel_Knu, basename,
                                   "nu", vborderd,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( false, true ), 1001 ) );


  //  5.2.1.11  Cylindrical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "5.2.1.11  cyl_neumann" << std::endl;
  basename = "diff_cyl_neumann";
  typedef double cyl_neumann(double, double);
  rundiff<double, double, double>( (cyl_neumann*)std::tr1::cyl_neumann,
                                   gsl::bessel_Ynu, basename,
                                   "nu", vborderd,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( false, true ), 1001 ) );


  //  5.2.1.12  Elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.12  ellint_1" << std::endl;
  basename = "diff_ellint_1";
  typedef double ellint_1(double, double);
  rundiff<double, double, double>( (ellint_1*)std::tr1::ellint_1,
                                   gsl::ellint_F, basename,
                                   "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                       std::make_pair( false, false ) ),
                                   "phi", vphid );


  //  5.2.1.13  Elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.13  ellint_2" << std::endl;
  basename = "diff_ellint_2";
  typedef double ellint_2(double, double);
  rundiff<double, double, double>( (ellint_2*)std::tr1::ellint_2,
                                   gsl::ellint_E, basename,
                                   "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                       std::make_pair( false, false ) ),
                                   "phi", vphid );


  //  5.2.1.14  Elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.14  ellint_3" << std::endl;
  basename = "diff_ellint_3";
  typedef double ellint_3(double, double, double);
  rundiff<double, double, double, double>( (ellint_3*)std::tr1::ellint_3,
                                           gsl::ellint_P, basename,
                                           "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                               std::make_pair( false, false ) ),
                                           "nu", vnud, "phi", vphid );


  //  5.2.1.15  Exponential integral.
  //  Skip the pole at 0.
  std::cout << "5.2.1.15  expint" << std::endl;
  basename = "diff_expint_neg";
  typedef double expint(double);
  rundiff<double, double>( (expint*)std::tr1::expint,
                           gsl::expint_Ei, basename,
                           "x", fill_argument( std::make_pair( -50.0, 0.0 ),
                                               std::make_pair( true, false ), 51 ) );
  basename = "diff_expint_pos";
  rundiff<double, double>( (expint*)std::tr1::expint,
                           gsl::expint_Ei, basename,
                           "x", fill_argument( std::make_pair( 0.0, 50.0 ),
                                               std::make_pair( false, true ), 51 ) );


  //  5.2.1.16  Hermite polynomials
  std::cout << "5.2.1.16  hermite  UNTESTED" << std::endl;
  //basename = "diff_hermite";
  //typedef double hermite(unsigned int, double);
  //rundiff<double, unsigned int, double>( (hermite*)std::tr1::hermite, gsl_xxx, basename, vorder,
  //                                       fill_argument( std::make_pair( -10.0, 10.0 ),
  //                                                      std::make_pair( true, true ) ) );


  //  5.2.1.17  Hypergeometric functions.
  //  Skip the singularity at c = 0.
  //  Skip the singularity at |x| = 1.
  std::cout << "5.2.1.17  hyperg" << std::endl;
  basename = "diff_hyperg";
  typedef double hyperg(double, double, double, double);
  rundiff<double, double, double, double, double>( (hyperg*)std::tr1::hyperg,
                                                   gsl::hyperg_2F1, basename,
                                                   "a", fill_argument( std::make_pair( 0.0, 10.0 ),
                                                                       std::make_pair( true, true ), 11 ),
                                                   "b", fill_argument( std::make_pair( 0.0, 10.0 ),
                                                                       std::make_pair( true, true ), 11 ),
                                                   "c", fill_argument( std::make_pair( 0.0, 10.0 ),
                                                                       std::make_pair( false, true ), 11 ),
                                                   "x", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                                       std::make_pair( false, false ), 21 ) );


  //  5.2.1.18  Laguerre polynomials.
  std::cout << "5.2.1.18  laguerre" << std::endl;
  basename = "diff_laguerre";
  typedef double laguerre(unsigned int, double);
  rundiff<double, unsigned int, double>( (laguerre*)std::tr1::laguerre,
                                         gsl::laguerre_n, basename,
                                         "n", vorder,
                                         "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                             std::make_pair( true, true ), 1001 ) );


  //  5.2.1.19  Legendre polynomials.
  std::cout << "5.2.1.19  legendre" << std::endl;
  basename = "diff_legendre";
  typedef double legendre(unsigned int, double);
  rundiff<double, unsigned int, double>( (legendre*)std::tr1::legendre,
                                         gsl::legendre_Pl, basename,
                                         "l", vorder,
                                         "x", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                             std::make_pair( true, true ), 1001 ) );


  //  5.2.1.20  Riemann zeta function.
  //  Skip the pole at 1.
  std::cout << "5.2.1.20  riemann_zeta" << std::endl;
  basename = "diff_riemann_zeta_neg";
  typedef double riemann_zeta(double);
  rundiff<double, double>( (riemann_zeta*)std::tr1::riemann_zeta,
                           gsl::zeta, basename,
                           "x", fill_argument( std::make_pair( -10.0, 1.0 ),
                                               std::make_pair( true, false ), 56 ) );
  basename = "diff_riemann_zeta_pos";
  rundiff<double, double>( (riemann_zeta*)std::tr1::riemann_zeta,
                           gsl::zeta, basename,
                           "x", fill_argument( std::make_pair( 1.0, 30.0 ),
                                               std::make_pair( false, true ), 146 ) );


  //  5.2.1.21  Spherical Bessel functions.
  std::cout << "5.2.1.21  sph_bessel" << std::endl;
  basename = "diff_sph_bessel";
  typedef double sph_bessel(unsigned int, double);
  rundiff<double, unsigned int, double>( (sph_bessel*)std::tr1::sph_bessel,
                                         gsl::bessel_jl, basename,
                                         "n", sborder,
                                         "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                             std::make_pair( true, true ), 1001 ) );


  //  5.2.1.22  Spherical Legendre functions.
  std::cout << "5.2.1.22  sph_legendre" << std::endl;
  basename = "diff_sph_legendre";
  typedef double sph_legendre(unsigned int, unsigned int, double);
  rundiff<double, unsigned int, unsigned int, double>( (sph_legendre*)std::tr1::sph_legendre,
                                                       gsl::legendre_sphPlm, basename,
                                                       "l", vorder, "m", vorder,
                                                       "theta", fill_argument( std::make_pair( 0.0, static_cast<double>(M_PI) ),
                                                                               std::make_pair( true, true ), 1001 ) );


  //  5.2.1.23  Spherical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "5.2.1.23  sph_neumann" << std::endl;
  basename = "diff_sph_neumann";
  typedef double sph_neumann(unsigned int, double);
  rundiff<double, unsigned int, double>( (sph_neumann*)std::tr1::sph_neumann,
                                         gsl::bessel_yl, basename,
                                         "n", sborder,
                                         "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                             std::make_pair( false, true ), 1001 ) );
                                                       


  return 0;
}

