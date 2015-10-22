
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <stdexcept>

#include <cmath>
#include <tr1/cmath>
#include <tr1/array>
#include <tr1/functional>
using namespace std::tr1::placeholders;

#include "gsl_wrap.h"

#include "test_func.tcc"


///
///
///
int main(int, char **)
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

  const unsigned long num_phi = 10;
  double phi[num_phi];
  for ( unsigned int i = 0; i < num_phi; ++i )
    phi[i] = 10.0 * i * M_PI / 180.0;
  std::vector<double> vphid( phi, phi + num_phi );

  double ab[] = { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0 };
  std::vector<double> vab( ab, ab + sizeof(ab) / sizeof(double) );

  std::string basename;


  //  5.2.1.1  Associated Laguerre polynomials.
  std::cout << "5.2.1.1  assoc_laguerre" << std::endl;
  basename = "diff_assoc_laguerre";
  rundiff<double, unsigned int, unsigned int, double>( std::tr1::assoc_laguerre,
                                                       wrap_gsl_sf_laguerre_nm, basename,
                                                       "n", vorder, "m", vorder,
                                                       "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                                           std::make_pair( true, true ) ) );


  //  5.2.1.2  Associated Legendre functions.
  std::cout << "5.2.1.2  assoc_legendre" << std::endl;
  basename = "diff_assoc_legendre";
  rundiff<double, unsigned int, unsigned int, double>( std::tr1::assoc_legendre,
                                                       wrap_gsl_sf_legendre_Plm, basename,
                                                       "l", vorder, "m", vorder,
                                                       "x", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                                           std::make_pair( true, true ), 1001 ) );


  //  5.2.1.3  Beta function.
  std::cout << "5.2.1.3  beta" << std::endl;
  basename = "diff_beta";
  rundiff<double, double, double>( std::tr1::beta,
                                   wrap_gsl_sf_beta, basename,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( false, true ) ),
                                   "y", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( false, true ) ) );


  //  5.2.1.4  Complete elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.4  comp_ellint_1" << std::endl;
  basename = "diff_comp_ellint_1";
  rundiff<double, double>( std::tr1::comp_ellint_1,
                           wrap_gsl_sf_ellint_Kcomp, basename,
                           "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                               std::make_pair( false, false ) ) );


  //  5.2.1.5  Complete elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.5  comp_ellint_2" << std::endl;
  basename = "diff_comp_ellint_2";
  rundiff<double, double>( std::tr1::comp_ellint_2,
                           wrap_gsl_sf_ellint_Ecomp, basename,
                           "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                               std::make_pair( false, false ) ) );


  //  5.2.1.6  Complete elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "5.2.1.6  comp_ellint_3" << std::endl;
  basename = "diff_comp_ellint_3";
  rundiff<double, double, double>( std::tr1::comp_ellint_3,
                                   wrap_gsl_sf_ellint_Pcomp, basename,
                                   "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                       std::make_pair( false, false ) ),
                                   "nu", fill_argument(std::make_pair(0.0, 1.0),
                                                       std::make_pair(true, false), 11) );


  //  5.2.1.7  Confluent hypergeometric functions.
  //  Skip the singularity at c = 0.
  std::cout << "5.2.1.7  conf_hyperg" << std::endl;
  basename = "diff_conf_hyperg";
  rundiff<double, double, double, double>( std::tr1::conf_hyperg,
                                           wrap_gsl_sf_hyperg_1F1, basename,
                                           "a", vab,//fill_argument( std::make_pair( 0.0, 10.0 ),
                                                //               std::make_pair( true, true ), 11 ),
                                           "c", fill_argument( std::make_pair( 0.0, 10.0 ),
                                                               std::make_pair( false, true ), 11 ),
                                           "x", fill_argument( std::make_pair( -10.0, 10.0 ),
                                                               std::make_pair( true, true ), 201 ) );


  //  5.2.1.8  Regular modified cylindrical Bessel functions.
  std::cout << "5.2.1.8  cyl_bessel_i" << std::endl;
  basename = "diff_cyl_bessel_i";
  rundiff<double, double, double>( std::tr1::cyl_bessel_i,
                                   wrap_gsl_sf_bessel_Inu, basename,
                                   "nu", vborderd,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( true, true ), 1001 ) );


  //  5.2.1.9  Cylindrical Bessel functions (of the first kind).
  std::cout << "5.2.1.9  cyl_bessel_j" << std::endl;
  basename = "diff_cyl_bessel_j";
  rundiff<double, double, double>( std::tr1::cyl_bessel_j,
                                   wrap_gsl_sf_bessel_Jnu, basename,
                                   "nu", vborderd,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( true, true ), 1001 ) );


  //  5.2.1.10  Irregular modified cylindrical Bessel functions.
  // Skip the pole at the origin.
  std::cout << "5.2.1.10  cyl_bessel_k" << std::endl;
  basename = "diff_cyl_bessel_k";
  rundiff<double, double, double>( std::tr1::cyl_bessel_k,
                                   wrap_gsl_sf_bessel_Knu, basename,
                                   "nu", vborderd,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( false, true ), 1001 ) );


  //  5.2.1.11  Cylindrical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "5.2.1.11  cyl_neumann" << std::endl;
  basename = "diff_cyl_neumann";
  rundiff<double, double, double>( std::tr1::cyl_neumann,
                                   wrap_gsl_sf_bessel_Ynu, basename,
                                   "nu", vborderd,
                                   "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                       std::make_pair( false, true ), 1001 ) );


  //  5.2.1.12  Elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.12  ellint_1" << std::endl;
  basename = "diff_ellint_1";
  rundiff<double, double, double>( std::tr1::ellint_1,
                                   wrap_gsl_sf_ellint_F, basename,
                                   "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                       std::make_pair( false, false ) ),
                                   "phi", vphid );


  //  5.2.1.13  Elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "5.2.1.13  ellint_2" << std::endl;
  basename = "diff_ellint_2";
  rundiff<double, double, double>( std::tr1::ellint_2,
                                   wrap_gsl_sf_ellint_E, basename,
                                   "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                       std::make_pair( false, false ) ),
                                   "phi", vphid );


  //  5.2.1.14  Elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "5.2.1.14  ellint_3" << std::endl;
  basename = "diff_ellint_3";
  rundiff<double, double, double, double>( std::tr1::ellint_3,
                                           wrap_gsl_sf_ellint_P, basename,
                                           "k", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                               std::make_pair( false, false ) ),
                                           "nu", fill_argument(std::make_pair(0.0, 1.0),
                                                               std::make_pair(true, false), 11),
                                           "phi", vphid );


  //  5.2.1.15  Exponential integral.
  //  Skip the pole at 0.
  std::cout << "5.2.1.15  expint" << std::endl;
  basename = "diff_expint_neg";
  rundiff<double, double>( std::tr1::expint,
                           wrap_gsl_sf_expint_Ei, basename,
                           "x", fill_argument( std::make_pair( -50.0, 0.0 ),
                                               std::make_pair( true, false ), 51 ) );
  basename = "diff_expint_pos";
  rundiff<double, double>( std::tr1::expint,
                           wrap_gsl_sf_expint_Ei, basename,
                           "x", fill_argument( std::make_pair( 0.0, 50.0 ),
                                               std::make_pair( false, true ), 51 ) );


  //  5.2.1.16  Hermite polynomials
  std::cout << "5.2.1.16  hermite  UNTESTED" << std::endl;
  //basename = "diff_hermite";
  //rundiff<double, unsigned int, double>( std::tr1::hermite, gsl_xxx, basename, vorder,
  //                                       fill_argument( std::make_pair( -10.0, 10.0 ),
  //                                                      std::make_pair( true, true ) ) );


  //  5.2.1.17  Hypergeometric functions.
  //  Skip the singularity at c = 0.
  //  Skip the singularities at x = -1.
  std::cout << "5.2.1.17  hyperg" << std::endl;
  basename = "diff_hyperg";
  rundiff<double, double, double, double, double>( std::tr1::hyperg,
                                                   wrap_gsl_sf_hyperg_2F1, basename,
                                                   "a", vab,//fill_argument( std::make_pair( 0.0, 10.0 ),
                                                        //               std::make_pair( true, true ), 11 ),
                                                   "b", vab,//fill_argument( std::make_pair( 0.0, 10.0 ),
                                                        //               std::make_pair( true, true ), 11 ),
                                                   "c", fill_argument( std::make_pair( 0.0, 10.0 ),
                                                                       std::make_pair( false, true ), 11 ),
                                                   "x", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                                       std::make_pair( false, true ), 21 ) );


  //  5.2.1.18  Laguerre polynomials.
  std::cout << "5.2.1.18  laguerre" << std::endl;
  basename = "diff_laguerre";
  rundiff<double, unsigned int, double>( std::tr1::laguerre,
                                         wrap_gsl_sf_laguerre_n, basename,
                                         "n", vorder,
                                         "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                             std::make_pair( true, true ), 1001 ) );


  //  5.2.1.19  Legendre polynomials.
  std::cout << "5.2.1.19  legendre" << std::endl;
  basename = "diff_legendre";
  rundiff<double, unsigned int, double>( std::tr1::legendre,
                                         wrap_gsl_sf_legendre_Pl, basename,
                                         "l", vorder,
                                         "x", fill_argument( std::make_pair( -1.0, 1.0 ),
                                                             std::make_pair( true, true ), 1001 ) );


  //  5.2.1.20  Riemann zeta function.
  //  Skip the pole at 1.
  std::cout << "5.2.1.20  riemann_zeta" << std::endl;
  basename = "diff_riemann_zeta_neg";
  rundiff<double, double>( std::tr1::riemann_zeta,
                           wrap_gsl_sf_zeta, basename,
                           "x", fill_argument( std::make_pair( -10.0, 1.0 ),
                                               std::make_pair( true, false ), 56 ) );
  basename = "diff_riemann_zeta_pos";
  rundiff<double, double>( std::tr1::riemann_zeta,
                           wrap_gsl_sf_zeta, basename,
                           "x", fill_argument( std::make_pair( 1.0, 30.0 ),
                                               std::make_pair( false, true ), 146 ) );


  //  5.2.1.21  Spherical Bessel functions.
  std::cout << "5.2.1.21  sph_bessel" << std::endl;
  basename = "diff_sph_bessel";
  rundiff<double, unsigned int, double>( std::tr1::sph_bessel,
                                         wrap_gsl_sf_bessel_jl, basename,
                                         "n", sborder,
                                         "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                             std::make_pair( true, true ), 1001 ) );


  //  5.2.1.22  Spherical Legendre functions.
  std::cout << "5.2.1.22  sph_legendre" << std::endl;
  basename = "diff_sph_legendre";
  rundiff<double, unsigned int, unsigned int, double>( std::tr1::sph_legendre,
                                                       wrap_gsl_sf_legendre_sphPlm, basename,
                                                       "l", vorder, "m", vorder,
                                                       "theta", fill_argument( std::make_pair( 0.0, static_cast<double>(M_PI) ),
                                                                               std::make_pair( true, true ), 1001 ) );


  //  5.2.1.23  Spherical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "5.2.1.23  sph_neumann" << std::endl;
  basename = "diff_sph_neumann";
  rundiff<double, unsigned int, double>( std::tr1::sph_neumann,
                                         wrap_gsl_sf_bessel_yl, basename,
                                         "n", sborder,
                                         "x", fill_argument( std::make_pair( 0.0, 100.0 ),
                                                             std::make_pair( false, true ), 1001 ) );
                                                       


  return 0;
}

