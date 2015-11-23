
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <stdexcept>

#define STD 1

#if STD
#  include <cmath>
#  include <array>
#else
#  include <cmath>
#  include <tr1/cmath>
#  include <tr1/array>
#endif

#include <gsl/gsl_sf.h>

#include "test_func.tcc"

#include "gsl_wrap.h"


///
///
///
int
main()
{
  //  Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
  std::vector<unsigned int> vorder{ 0, 1, 2, 3, 4, 5, 10, 20, 50, 100 };

  //  ... corresponding signed integer orders for GSL.
  std::vector<int> ivorder{ 0, 1, 2, 3, 4, 5, 10, 20, 50, 100 };

  //  ... corresponding "double" integer orders for GSL.
  std::vector<double> dvorder{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 50.0, 100.0 };

  //  Orders for cylindrical Bessel functions.
  //const unsigned long num_border = 12;
  //long double border[num_border] = { 0.0L, 1.0L/3.0L, 0.5L, 1.0L, 2.0L, 3.0L, 5.0L,
  //                          10.0L, 20.0L, 50.0L, 100.0L, 200.0L };
  const unsigned long num_border = 12;
  long double border[num_border] = { 0.0L, 1.0L/3.0L, 0.5L, 2.0L/3.0L, 1.0L, 2.0L, 3.0L,
                                     5.0L, 10.0L, 20.0L, 50.0L, 100.0L };
  std::vector<float> vborderf(border, border + num_border);
  std::vector<double> vborderd(border, border + num_border);
  std::vector<long double> vborderl(border, border + num_border);

  //  Orders for spherical bessel functions.
  std::vector<unsigned int> sborder{ 0, 1, 2, 3, 4, 5, 10, 20, 50, 100 };
  //  ... and for GSL.
  std::vector<int> isborder{ 0, 1, 2, 3, 4, 5, 10, 20, 50, 100 };

  const unsigned long num_phi = 10;
  long double phi[num_phi];
  for (unsigned int i = 0; i < num_phi; ++i)
    phi[i] = 10.0L * i * static_cast<long double>(M_PI) / 180.0L;
  std::vector<float> vphif(phi, phi + num_phi);
  std::vector<double> vphid(phi, phi + num_phi);
  std::vector<long double> vphil(phi, phi + num_phi);

#if 1
  std::string ns("std");
  using std::assoc_laguerre;
  using std::assoc_legendre;
  using std::beta;
  using std::comp_ellint_1;
  using std::comp_ellint_2;
  using std::comp_ellint_3;
  using __gnu_cxx::conf_hyperg;
  using std::cyl_bessel_i;
  using std::cyl_bessel_j;
  using std::cyl_bessel_k;
  using std::cyl_neumann;
  using std::ellint_1;
  using std::ellint_2;
  using std::ellint_3;
  using std::expint;
  using std::hermite;
  using __gnu_cxx::hyperg;
  using std::laguerre;
  using std::legendre;
  using std::riemann_zeta;
  using std::sph_bessel;
  using std::sph_legendre;
  using std::sph_neumann;
#else
  std::string ns("tr1");
  using std::tr1::assoc_laguerre;
  using std::tr1::assoc_legendre;
  using std::tr1::beta;
  using std::tr1::comp_ellint_1;
  using std::tr1::comp_ellint_2;
  using std::tr1::comp_ellint_3;
  using std::tr1::conf_hyperg;
  using std::tr1::cyl_bessel_i;
  using std::tr1::cyl_bessel_j;
  using std::tr1::cyl_bessel_k;
  using std::tr1::cyl_neumann;
  using std::tr1::ellint_1;
  using std::tr1::ellint_2;
  using std::tr1::ellint_3;
  using std::tr1::expint;
  using std::tr1::hermite;
  using std::tr1::hyperg;
  using std::tr1::laguerre;
  using std::tr1::legendre;
  using std::tr1::riemann_zeta;
  using std::tr1::sph_bessel;
  using std::tr1::sph_legendre;
  using std::tr1::sph_neumann;
#endif


  std::string basename;


  //  5.2.1.1  Associated Laguerre polynomials.
  //  double gsl_sf_laguerre_n(int n, double a, double x);
  std::cout << "5.2.1.1  assoc_laguerre" << std::endl;
  basename = "gsl_assoc_laguerre";
  runtest<double, int, double, double>(gsl_sf_laguerre_n, basename, ivorder, dvorder,
                                        fill_argument(std::make_pair(0.0, 100.0),
                                                       std::make_pair(true, true)));
  basename = ns + "_assoc_laguerre";
  runtest<float, unsigned int, unsigned int, float>(assoc_laguerre, basename, vorder, vorder,
                                                     fill_argument(std::make_pair(0.0F, 100.0F),
                                                                    std::make_pair(true, true)));
  runtest<double, unsigned int, unsigned int, double>(assoc_laguerre, basename, vorder, vorder,
                                                       fill_argument(std::make_pair(0.0, 100.0),
                                                                      std::make_pair(true, true)));
  runtest<long double, unsigned int, unsigned int, long double>(assoc_laguerre, basename, vorder, vorder,
                                                                 fill_argument(std::make_pair(0.0L, 100.0L),
                                                                                std::make_pair(true, true)));


  //  5.2.1.2  Associated Legendre functions.
  std::cout << "5.2.1.2  assoc_legendre" << std::endl;
  basename = "gsl_assoc_legendre";
  runtest<double, unsigned int, unsigned int, double>(wrap_gsl_sf_legendre_Plm,
                                                       basename, vorder, vorder,
                                                       fill_argument(std::make_pair(-1.0, 1.0),
                                                                      std::make_pair(true, true),
                                                                      1001));
  basename = ns + "_assoc_legendre";
  runtest<float, unsigned int, unsigned int, float>(assoc_legendre,
                                                     basename, vorder, vorder,
                                                     fill_argument(std::make_pair(-1.0F, 1.0F),
                                                                    std::make_pair(true, true),
                                                                    1001));
  runtest<double, unsigned int, unsigned int, double>(assoc_legendre,
                                                       basename, vorder, vorder,
                                                       fill_argument(std::make_pair(-1.0, 1.0),
                                                                      std::make_pair(true, true),
                                                                      1001));
  runtest<long double, unsigned int, unsigned int, long double>(assoc_legendre,
                                                                 basename, vorder, vorder,
                                                                 fill_argument(std::make_pair(-1.0L, 1.0L),
                                                                                std::make_pair(true, true),
                                                                                1001));


  //  5.2.1.3  Beta function.
  std::cout << "5.2.1.3  beta" << std::endl;
  basename = "gsl_beta";
  runtest<double, double, double>(wrap_gsl_sf_beta, basename,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(false, true)),
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(false, true)));
  basename = ns + "_beta";
  runtest<float, float, float>(beta, basename,
                                fill_argument(std::make_pair(0.0F, 100.0F),
                                               std::make_pair(false, true)),
                                fill_argument(std::make_pair(0.0F, 100.0F),
                                               std::make_pair(false, true)));
  runtest<double, double, double>(beta, basename,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(false, true)),
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(false, true)));
  runtest<long double, long double, long double>(beta, basename,
                                                  fill_argument(std::make_pair(0.0L, 100.0L),
                                                                 std::make_pair(false, true)),
                                                  fill_argument(std::make_pair(0.0L, 100.0L),
                                                                 std::make_pair(false, true)));


  //  5.2.1.4  Complete elliptic integrals of the first kind.
  std::cout << "5.2.1.4  comp_ellint_1" << std::endl;
  basename = ns + "_comp_ellint_1";
  basename = "gsl_comp_ellint_1";
  runtest<double, double>(wrap_gsl_sf_ellint_Kcomp, basename,
                           fill_argument(std::make_pair(-1.0, 1.0),
                                          std::make_pair(false, false)));  //  Avoid poles at |x| = 1.
  runtest<float, float>(comp_ellint_1, basename,
                         fill_argument(std::make_pair(-1.0F, 1.0F),
                                        std::make_pair(true, true)));
  runtest<double, double>(comp_ellint_1, basename,
                           fill_argument(std::make_pair(-1.0, 1.0),
                                          std::make_pair(true, true)));
  runtest<long double, long double>(comp_ellint_1, basename,
                                     fill_argument(std::make_pair(-1.0L, 1.0L),
                                                    std::make_pair(true, true)));


  //  5.2.1.5  Complete elliptic integrals of the second kind.
  std::cout << "5.2.1.5  comp_ellint_2" << std::endl;
  basename = "gsl_comp_ellint_2";
  runtest<double, double>(wrap_gsl_sf_ellint_Ecomp, basename,
                           fill_argument(std::make_pair(-1.0, 1.0),
                                          std::make_pair(false, false)));  //  Avoid poles at |x| = 1.
  basename = ns + "_comp_ellint_2";
  runtest<float, float>(comp_ellint_2, basename,
                         fill_argument(std::make_pair(-1.0F, 1.0F),
                                        std::make_pair(true, true)));
  runtest<double, double>(comp_ellint_2, basename,
                           fill_argument(std::make_pair(-1.0, 1.0),
                                          std::make_pair(true, true)));
  runtest<long double, long double>(comp_ellint_2, basename,
                                     fill_argument(std::make_pair(-1.0L, 1.0L),
                                                    std::make_pair(true, true)));


  //  5.2.1.6  Complete elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "5.2.1.6  comp_ellint_3" << std::endl;
  basename = "gsl_comp_ellint_3";
  runtest<double, double, double>(wrap_gsl_sf_ellint_Pcomp,
                                   basename,
                                   fill_argument(std::make_pair(-1.0, 1.0),
                                                  std::make_pair(false, false)),
                                   fill_argument(std::make_pair(0.0, 1.0),
                                                 std::make_pair(true, false), 11));
  basename = ns + "_comp_ellint_3";
  runtest<float, float, float>(comp_ellint_3,
                                basename,
                                fill_argument(std::make_pair(-1.0F, 1.0F),
                                               std::make_pair(true, true)),
                                fill_argument(std::make_pair(0.0F, 1.0F),
                                              std::make_pair(true, true), 11));
  runtest<double, double, double>(comp_ellint_3,
                                   basename,
                                   fill_argument(std::make_pair(-1.0, 1.0),
                                                  std::make_pair(true, true)),
                                   fill_argument(std::make_pair(0.0, 1.0),
                                                 std::make_pair(true, true), 11));
  runtest<long double, long double, long double>(comp_ellint_3,
                                     basename,
                                     fill_argument(std::make_pair(-1.0L, 1.0L),
                                                    std::make_pair(true, true)),
                                     fill_argument(std::make_pair(0.0L, 1.0L),
                                                   std::make_pair(true, true), 11));


  //  5.2.1.7  Confluent hypergeometric functions.
  std::cout << "5.2.1.7  conf_hyperg" << std::endl;
  basename = "gsl_conf_hyperg";
  runtest<double, double, double, double>(wrap_gsl_sf_hyperg_1F1,
                                           basename,
                                           fill_argument(std::make_pair(0.0, 10.0),
                                                          std::make_pair(true, true), 11),
                                           fill_argument(std::make_pair(0.0, 10.0),
                                                          std::make_pair(false, true), 11),  //  Skip the singularity
                                           fill_argument(std::make_pair(-10.0, 10.0),
                                                          std::make_pair(true, true), 201));
  basename = ns + "_conf_hyperg";
  runtest<float, float, float, float>(conf_hyperg,
                                       basename,
                                       fill_argument(std::make_pair(0.0F, 10.0F),
                                                      std::make_pair(true, true), 11),
                                       fill_argument(std::make_pair(0.0F, 10.0F),
                                                      std::make_pair(true, true), 11),
                                       fill_argument(std::make_pair(-10.0F, 10.0F),
                                                      std::make_pair(true, true), 201));
  runtest<double, double, double, double>(conf_hyperg,
                                           basename,
                                           fill_argument(std::make_pair(0.0, 10.0),
                                                          std::make_pair(true, true), 11),
                                           fill_argument(std::make_pair(0.0, 10.0),
                                                          std::make_pair(true, true), 11),
                                           fill_argument(std::make_pair(-10.0, 10.0),
                                                          std::make_pair(true, true), 201));
  runtest<long double, long double, long double, long double>(conf_hyperg,
                                                               basename,
                                                               fill_argument(std::make_pair(0.0L, 10.0L),
                                                                              std::make_pair(true, true), 11),
                                                               fill_argument(std::make_pair(0.0L, 10.0L),
                                                                              std::make_pair(true, true), 11),
                                                               fill_argument(std::make_pair(-10.0L, 10.0L),
                                                                              std::make_pair(true, true), 201));


  //  5.2.1.8  Regular modified cylindrical Bessel functions.
  std::cout << "5.2.1.8  cyl_bessel_i" << std::endl;
  basename = "gsl_cyl_bessel_i";
  runtest<double, double, double>(gsl_sf_bessel_Inu, basename, vborderd,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(true, true),
                                                  1001));
  basename = ns + "_cyl_bessel_i";
  runtest<float, float, float>(cyl_bessel_i, basename, vborderf,
                                fill_argument(std::make_pair(0.0F, 100.0F),
                                               std::make_pair(true, true),
                                               1001));
  runtest<double, double, double>(cyl_bessel_i, basename, vborderd,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(true, true),
                                                  1001));
  runtest<long double, long double, long double>(cyl_bessel_i, basename, vborderl,
                                                  fill_argument(std::make_pair(0.0L, 100.0L),
                                                                 std::make_pair(true, true),
                                                                 1001));


  //  5.2.1.9  Cylindrical Bessel functions (of the first kind).
  std::cout << "5.2.1.9  cyl_bessel_j" << std::endl;
  basename = "gsl_cyl_bessel_j";
  runtest<double, double, double>(gsl_sf_bessel_Jnu, basename, vborderd,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(true, true),
                                                  1001));
  basename = ns + "_cyl_bessel_j";
  runtest<float, float, float>(cyl_bessel_j, basename, vborderf,
                                fill_argument(std::make_pair(0.0F, 100.0F),
                                               std::make_pair(true, true),
                                               1001));
  runtest<double, double, double>(cyl_bessel_j, basename, vborderd,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(true, true),
                                                  1001));
  runtest<long double, long double, long double>(cyl_bessel_j, basename, vborderl,
                                                  fill_argument(std::make_pair(0.0L, 100.0L),
                                                                 std::make_pair(true, true),
                                                                 1001));


  //  5.2.1.10  Irregular modified cylindrical Bessel functions.
  std::cout << "5.2.1.10  cyl_bessel_k" << std::endl;
  basename = "gsl_cyl_bessel_k";
  runtest<double, double, double>(gsl_sf_bessel_Knu, basename, vborderd,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(false, true),  // Skip the pole at the origin.
                                                  1001));
  basename = ns + "_cyl_bessel_k";
  runtest<float, float, float>(cyl_bessel_k, basename, vborderf,
                                fill_argument(std::make_pair(0.0F, 100.0F),
                                               std::make_pair(true, true),
                                               1001));
  runtest<double, double, double>(cyl_bessel_k, basename, vborderd,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(true, true),
                                                  1001));
  runtest<long double, long double, long double>(cyl_bessel_k, basename, vborderl,
                                                  fill_argument(std::make_pair(0.0L, 100.0L),
                                                                 std::make_pair(true, true),
                                                                 1001));


  //  5.2.1.11  Cylindrical Neumann functions.
  std::cout << "5.2.1.11  cyl_neumann" << std::endl;
  basename = "gsl_cyl_neumann";
  runtest<double, double, double>(gsl_sf_bessel_Ynu, basename, vborderd,
                           fill_argument(std::make_pair(0.0, 100.0),
                                          std::make_pair(false, true),  // Skip the pole at the origin.
                                          1001));
  basename = ns + "_cyl_neumann";
  runtest<float, float, float>(cyl_neumann, basename, vborderf,
                                fill_argument(std::make_pair(0.0F, 100.0F),
                                               std::make_pair(true, true),
                                               1001));
  runtest<double, double, double>(cyl_neumann, basename, vborderd,
                                   fill_argument(std::make_pair(0.0, 100.0),
                                                  std::make_pair(true, true),
                                                  1001));
  runtest<long double, long double, long double>(cyl_neumann, basename, vborderl,
                                                  fill_argument(std::make_pair(0.0L, 100.0L),
                                                                 std::make_pair(true, true),
                                                                 1001));


  //  5.2.1.12  Elliptic integrals of the first kind.
  std::cout << "5.2.1.12  ellint_1" << std::endl;
  basename = "gsl_ellint_1";
  runtest<double, double, double>(wrap_gsl_sf_ellint_F,
                                   basename,
                                   fill_argument(std::make_pair(-1.0, 1.0),
                                                  std::make_pair(false, false)),  //  Avoid poles at |x| = 1.
                                   vphid);
  basename = ns + "_ellint_1";
  runtest<float, float, float>(ellint_1,
                                basename,
                                fill_argument(std::make_pair(-1.0F, 1.0F),
                                               std::make_pair(true, true)),
                                vphif);
  runtest<double, double, double>(ellint_1,
                                   basename,
                                   fill_argument(std::make_pair(-1.0, 1.0),
                                                  std::make_pair(true, true)),
                                   vphid);
  runtest<long double, long double, long double>(ellint_1,
                                     basename,
                                     fill_argument(std::make_pair(-1.0L, 1.0L),
                                                    std::make_pair(true, true)),
                                     vphil);


  //  5.2.1.13  Elliptic integrals of the second kind.
  std::cout << "5.2.1.13  ellint_2" << std::endl;
  basename = "gsl_ellint_2";
  runtest<double, double, double>(wrap_gsl_sf_ellint_E,
                                   basename,
                                   fill_argument(std::make_pair(-1.0, 1.0),
                                                  std::make_pair(false, false)),  //  Avoid poles at |x| = 1.
                                   vphid);
  basename = ns + "_ellint_2";
  runtest<float, float, float>(ellint_2,
                                basename,
                                fill_argument(std::make_pair(-1.0F, 1.0F),
                                               std::make_pair(true, true)),
                                vphif);
  runtest<double, double, double>(ellint_2,
                                   basename,
                                   fill_argument(std::make_pair(-1.0, 1.0),
                                                  std::make_pair(true, true)),
                                   vphid);
  runtest<long double, long double, long double>(ellint_2,
                                     basename,
                                     fill_argument(std::make_pair(-1.0L, 1.0L),
                                                    std::make_pair(true, true)),
                                     vphil);


  //  5.2.1.14  Elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "5.2.1.14  ellint_3" << std::endl;
  basename = "gsl_ellint_3";
  runtest<double, double, double, double>(wrap_gsl_sf_ellint_P,
                                           basename,
                                           fill_argument(std::make_pair(-1.0, 1.0),
                                                          std::make_pair(false, false)),
                                           fill_argument(std::make_pair(0.0, 1.0),
                                                         std::make_pair(true, false), 11),
                                           vphid);
  basename = ns + "_ellint_3";
  runtest<float, float, float, float>(ellint_3,
                                       basename,
                                       fill_argument(std::make_pair(-1.0F, 1.0F),
                                                      std::make_pair(true, true)),
                                       fill_argument(std::make_pair(0.0F, 1.0F),
                                                     std::make_pair(true, true), 11),
                                       vphif);
  runtest<double, double, double, double>(ellint_3,
                                           basename,
                                           fill_argument(std::make_pair(-1.0, 1.0),
                                                          std::make_pair(true, true)),
                                           fill_argument(std::make_pair(0.0, 1.0),
                                                         std::make_pair(true, true), 11),
                                           vphid);
  runtest<long double, long double, long double, long double>(ellint_3,
                                                               basename,
                                                               fill_argument(std::make_pair(-1.0L, 1.0L),
                                                               std::make_pair(true, true)),
                                                               fill_argument(std::make_pair(0.0L, 1.0L),
                                                                             std::make_pair(true, true), 11),
                                                               vphil);


  //  5.2.1.15  Exponential integrals.
  std::cout << "5.2.1.15  expint" << std::endl;
  //  Skip the pole at 0.
  basename = "gsl_expint_neg";
  runtest<double, double>(gsl_sf_expint_Ei, basename,
                           fill_argument(std::make_pair(-50.0, 0.0),
                                          std::make_pair(true, false), 51));
  basename = "gsl_expint_pos";
  runtest<double, double>(gsl_sf_expint_Ei, basename,
                           fill_argument(std::make_pair(0.0, 50.0),
                                          std::make_pair(false, true), 51));
  basename = ns + "_expint";
  runtest<float, float>(expint, basename,
                         fill_argument(std::make_pair(-50.0F, 50.0F),
                                        std::make_pair(true, true)));
  runtest<double, double>(expint, basename,
                           fill_argument(std::make_pair(-50.0, 50.0),
                                          std::make_pair(true, true)));
  runtest<long double, long double>(expint, basename,
                                     fill_argument(std::make_pair(-50.0L, 50.0L),
                                                    std::make_pair(true, true)));


  //  5.2.1.16  Hermite polynomials
  std::cout << "5.2.1.16  hermite" << std::endl;
  basename = "gsl_hermite";
  //runtest<double, unsigned int, double>(gsl_sf_hermite, basename, vorder,
  //                                       fill_argument(std::make_pair(-10.0, 10.0),
  //                                                      std::make_pair(true, true)));
  basename = ns + "_hermite";
  runtest<float, unsigned int, float>(hermite, basename, vorder,
                                       fill_argument(std::make_pair(-10.0F, 10.0F),
                                                      std::make_pair(true, true)));
  runtest<double, unsigned int, double>(hermite, basename, vorder,
                                         fill_argument(std::make_pair(-10.0, 10.0),
                                                        std::make_pair(true, true)));
  runtest<long double, unsigned int, long double>(hermite, basename, vorder,
                                                   fill_argument(std::make_pair(-10.0L, 10.0L),
                                                                  std::make_pair(true, true)));


  //  5.2.1.17  Hypergeometric functions.
  std::cout << "5.2.1.17  hyperg" << std::endl;
  basename = "gsl_hyperg";
  runtest<double, double, double, double, double>(wrap_gsl_sf_hyperg_2F1,
                                                   basename,
                                                   fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(true, true), 11),
                                                   fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(true, true), 11),
                                                   fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(false, true), 11),  //  Skip the singularity
                                                   fill_argument(std::make_pair(-1.0, 1.0),
                                                                  std::make_pair(false, true), 21));
  basename = ns + "_hyperg";
  runtest<float, float, float, float, float>(hyperg,
                                              basename,
                                              fill_argument(std::make_pair(0.0F, 10.0F),
                                                             std::make_pair(true, true), 11),
                                              fill_argument(std::make_pair(0.0F, 10.0F),
                                                             std::make_pair(true, true), 11),
                                              fill_argument(std::make_pair(0.0F, 10.0F),
                                                             std::make_pair(true, true), 11),
                                              fill_argument(std::make_pair(-1.0F, 1.0F),
                                                             std::make_pair(false, true), 21));
  runtest<double, double, double, double, double>(hyperg,
                                                   basename,
                                                   fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(true, true), 11),
                                                   fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(true, true), 11),
                                                   fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(true, true), 11),
                                                   fill_argument(std::make_pair(-1.0, 1.0),
                                                                  std::make_pair(false, true), 21));
  runtest<long double, long double, long double, long double, long double>(hyperg,
                                                                            basename,
                                                                            fill_argument(std::make_pair(0.0L, 10.0L),
                                                                                           std::make_pair(true, true), 11),
                                                                            fill_argument(std::make_pair(0.0L, 10.0L),
                                                                                           std::make_pair(true, true), 11),
                                                                            fill_argument(std::make_pair(0.0L, 10.0L),
                                                                                           std::make_pair(true, true), 11),
                                                                            fill_argument(std::make_pair(-1.0L, 1.0L),
                                                                                           std::make_pair(false, true), 21));


  //  5.2.1.18  Laguerre polynomials.
  std::cout << "5.2.1.18  laguerre" << std::endl;
  basename = "gsl_laguerre";
  runtest<double, unsigned int, double>(wrap_gsl_sf_laguerre_n, basename,
                                         vorder,
                                         fill_argument(std::make_pair(0.0, 100.0),
                                                        std::make_pair(true, true)));
  basename = ns + "_laguerre";
  runtest<float, unsigned int, float>(laguerre, basename, vorder,
                                       fill_argument(std::make_pair(0.0F, 100.0F),
                                                      std::make_pair(true, true)));
  runtest<double, unsigned int, double>(laguerre, basename, vorder,
                                         fill_argument(std::make_pair(0.0, 100.0),
                                                        std::make_pair(true, true)));
  runtest<long double, unsigned int, long double>(laguerre, basename, vorder,
                                                   fill_argument(std::make_pair(0.0L, 100.0L),
                                                                  std::make_pair(true, true)));


  ///  5.2.1.19  Legendre polynomials
  std::cout << "5.2.1.19  legendre" << std::endl;
  basename = "gsl_legendre";
  runtest<double, int, double>(gsl_sf_legendre_Pl, basename, ivorder,
                                fill_argument(std::make_pair(-1.0, 1.0),
                                               std::make_pair(true, true),
                                               1001));
  basename = ns + "_legendre";
  runtest<float, unsigned int, float>(legendre, basename, vorder,
                                       fill_argument(std::make_pair(-1.0F, 1.0F),
                                                      std::make_pair(true, true),
                                                      1001));
  runtest<double, unsigned int, double>(legendre, basename, vorder,
                                         fill_argument(std::make_pair(-1.0, 1.0),
                                                        std::make_pair(true, true),
                                                        1001));
  runtest<long double, unsigned int, long double>(legendre, basename, vorder,
                                                   fill_argument(std::make_pair(-1.0L, 1.0L),
                                                                  std::make_pair(true, true),
                                                                  1001));


  //  5.2.1.20  Riemann zeta function.
  std::cout << "5.2.1.20  riemann_zeta" << std::endl;
  //  Skip the pole at 1.
  basename = "gsl_riemann_zeta_neg";
  runtest<double, double>(wrap_gsl_sf_zeta, basename,
                           fill_argument(std::make_pair(-10.0, 1.0),
                                          std::make_pair(true, false), 56));
  basename = "gsl_riemann_zeta_pos";
  runtest<double, double>(wrap_gsl_sf_zeta, basename,
                           fill_argument(std::make_pair(1.0, 30.0),
                                          std::make_pair(false, true), 146));
  basename = ns + "_riemann_zeta";
  runtest<float, float>(riemann_zeta, basename,
                         fill_argument(std::make_pair(-10.0F, 30.0F),
                                        std::make_pair(true, true), 201));
  runtest<double, double>(riemann_zeta, basename,
                           fill_argument(std::make_pair(-10.0, 30.0),
                                          std::make_pair(true, true), 201));
  runtest<long double, long double>(riemann_zeta, basename,
                                     fill_argument(std::make_pair(-10.0L, 30.0L),
                                                    std::make_pair(true, true), 201));


  //  5.2.1.21  Spherical Bessel functions.
  std::cout << "5.2.1.21  sph_bessel" << std::endl;
  basename = "gsl_sph_bessel";
  runtest<double, int, double>(gsl_sf_bessel_jl, basename, isborder,
                                fill_argument(std::make_pair(0.0, 100.0),
                                               std::make_pair(true, true),
                                               1001));
  basename = ns + "_sph_bessel";
  runtest<float, unsigned int, float>(sph_bessel, basename, sborder,
                                       fill_argument(std::make_pair(0.0F, 100.0F),
                                                      std::make_pair(true, true),
                                                      1001));
  runtest<double, unsigned int, double>(sph_bessel, basename, sborder,
                                         fill_argument(std::make_pair(0.0, 100.0),
                                                        std::make_pair(true, true),
                                                        1001));
  runtest<long double, unsigned int, long double>(sph_bessel, basename, sborder,
                                                   fill_argument(std::make_pair(0.0L, 100.0L),
                                                                  std::make_pair(true, true),
                                                                  1001));

  //  5.2.1.21  Spherical Legendre functions.
  std::cout << "5.2.1.22  sph_legendre" << std::endl;
  basename = "gsl_sph_legendre";
  runtest<double, unsigned int, unsigned int, double>(wrap_gsl_sf_legendre_sphPlm, basename, vorder, vorder,
                                                       fill_argument(std::make_pair(0.0, static_cast<double>(M_PI)),
                                                                      std::make_pair(true, true),
                                                                      1001));
  basename = ns + "_sph_legendre";
  runtest<float, unsigned int, unsigned int, float>(sph_legendre, basename, vorder, vorder,
                                                     fill_argument(std::make_pair(0.0F, 100.0F),
                                                                    std::make_pair(true, true),
                                                                    1001));
  runtest<double, unsigned int, unsigned int, double>(sph_legendre, basename, vorder, vorder,
                                                       fill_argument(std::make_pair(0.0, 100.0),
                                                                      std::make_pair(true, true),
                                                                      1001));
  runtest<long double, unsigned int, unsigned int, long double>(sph_legendre, basename, vorder, vorder,
                                                                 fill_argument(std::make_pair(0.0L, 100.0L),
                                                                                std::make_pair(true, true),
                                                                                1001));



  //  5.2.1.23  Spherical Neumann functions.
  std::cout << "5.2.1.23  sph_neumann" << std::endl;
  basename = "gsl_sph_neumann";
  runtest<double, int, double>(gsl_sf_bessel_yl, basename, isborder,
                                fill_argument(std::make_pair(0.0, 100.0),
                                               std::make_pair(false, true),  // Skip the pole at the origin.
                                               1001));
  basename = ns + "_sph_neumann";
  runtest<float, unsigned int, float>(sph_neumann, basename, sborder,
                                       fill_argument(std::make_pair(0.0F, 100.0F),
                                                      std::make_pair(true, true),
                                                      1001));
  runtest<double, unsigned int, double>(sph_neumann, basename, sborder,
                                         fill_argument(std::make_pair(0.0, 100.0),
                                                        std::make_pair(true, true),
                                                        1001));
  runtest<long double, unsigned int, long double>(sph_neumann, basename, sborder,
                                                   fill_argument(std::make_pair(0.0L, 100.0L),
                                                                  std::make_pair(true, true),
                                                                  1001));


  return 0;
}

