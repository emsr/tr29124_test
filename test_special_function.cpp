
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
template<typename _Tp>
int
do_test<_Tp>()
{
  //  Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
  std::vector<unsigned int> uiorder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

  //  ... corresponding signed integer orders for GSL.
  std::vector<int> iorder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

  //  ... corresponding double "integer" orders for GSL.
  std::vector<double> duiorder{0.0, 1.0, 2.0, 3.0, 4.0,
			       5.0, 10.0, 20.0, 50.0, 100.0};

  //  Orders for cylindrical Bessel functions.
  std::vector<double> border{_Tp{0}, _Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1},
                              _Tp{2}, _Tp{3}, _Tp{5}, _Tp{10}, _Tp{20}, _Tp{50}, _Tp{100}};

  //  Orders for spherical Bessel functions.
  std::vector<unsigned int> uisborder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};
  //  ... and for GSL.
  std::vector<int> isborder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

  const unsigned long num_phi = 10;
  long double phi[num_phi];
  for (unsigned int i = 0; i < num_phi; ++i)
    phi[i] = _Tp{10} * i * static_cast<_Tp>(M_PI) / _Tp{180};
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

  //  Associated Laguerre polynomials.
  //  double gsl_sf_laguerre_n(int n, double a, double x);
  std::cout << "assoc_laguerre" << std::endl;
  basename = "gsl_assoc_laguerre";
  runtest<double, int, double, double>(gsl_sf_laguerre_n, basename, iorder, duiorder,
                                       fill_argument(std::make_pair(0.0, 100.0),
                                                     std::make_pair(true, true)));
  basename = ns + "_assoc_laguerre";
  runtest<_Tp, unsigned int, unsigned int, _Tp>(assoc_laguerre, basename, uiorder, uiorder,
                                                fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                                              std::make_pair(true, true)));


  //  Associated Legendre functions.
  std::cout << "assoc_legendre" << std::endl;
  basename = "gsl_assoc_legendre";
  runtest<double, unsigned int, unsigned int, double>(wrap_gsl_sf_legendre_Plm, basename, uiorder, uiorder,
                                                      fill_argument(std::make_pair(-1.0, 1.0),
                                                                    std::make_pair(true, true),
                                                      1001));
  basename = ns + "_assoc_legendre";
  runtest<_Tp, unsigned int, unsigned int, _Tp>(assoc_legendre, basename, uiorder, uiorder,
                                                fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                                              std::make_pair(true, true),
                                                1001));


  //  Beta function.
  std::cout << "beta" << std::endl;
  basename = "gsl_beta";
  runtest<double, double, double>(wrap_gsl_sf_beta, basename,
                                  fill_argument(std::make_pair(0.0, 100.0),
                                                std::make_pair(false, true)),
                                  fill_argument(std::make_pair(0.0, 100.0),
                                                std::make_pair(false, true)));
  basename = ns + "_beta";
  runtest<_Tp, _Tp, _Tp>(beta, basename,
                         fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                       std::make_pair(false, true)),
                         fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                       std::make_pair(false, true)));


  //  Complete elliptic integrals of the first kind.
  std::cout << "comp_ellint_1" << std::endl;
  basename = ns + "_comp_ellint_1";
  basename = "gsl_comp_ellint_1";
  runtest<double, double>(wrap_gsl_sf_ellint_Kcomp, basename,
                          fill_argument(std::make_pair(-1.0, 1.0),
                                        std::make_pair(false, false)));  //  Avoid poles at |x| = 1.
  runtest<_Tp, _Tp>(comp_ellint_1, basename,
                    fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                  std::make_pair(true, true)));


  //  Complete elliptic integrals of the second kind.
  std::cout << "comp_ellint_2" << std::endl;
  basename = "gsl_comp_ellint_2";
  runtest<double, double>(wrap_gsl_sf_ellint_Ecomp, basename,
                          fill_argument(std::make_pair(-1.0, 1.0),
                                        std::make_pair(false, false)));  //  Avoid poles at |x| = 1.
  basename = ns + "_comp_ellint_2";
  runtest<_Tp, _Tp>(comp_ellint_2, basename,
                    fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                  std::make_pair(true, true)));


  //  Complete elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "comp_ellint_3" << std::endl;
  basename = "gsl_comp_ellint_3";
  runtest<double, double, double>(wrap_gsl_sf_ellint_Pcomp, basename,
                                  fill_argument(std::make_pair(-1.0, 1.0),
                                                std::make_pair(false, false)),
                                  fill_argument(std::make_pair(0.0, 1.0),
                                                std::make_pair(true, false), 11));
  basename = ns + "_comp_ellint_3";
  runtest<_Tp, _Tp, _Tp>(comp_ellint_3, basename,
                         fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                       std::make_pair(true, true)),
                         fill_argument(std::make_pair(_Tp{0}, _Tp{1}),
                                       std::make_pair(true, true), 11));


  //  Confluent hypergeometric functions.
  std::cout << "conf_hyperg" << std::endl;
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
  runtest<_Tp, _Tp, _Tp, _Tp>(conf_hyperg, basename,
                              fill_argument(std::make_pair(_Tp{0}, _Tp{10}),
                                            std::make_pair(true, true), 11),
                              fill_argument(std::make_pair(_Tp{0}, _Tp{10}),
                                            std::make_pair(true, true), 11),
                              fill_argument(std::make_pair(-_Tp{10}, _Tp{10}),
                                            std::make_pair(true, true), 201));


  //  Regular modified cylindrical Bessel functions.
  std::cout << "cyl_bessel_i" << std::endl;
  basename = "gsl_cyl_bessel_i";
  runtest<double, double, double>(gsl_sf_bessel_Inu, basename, border,
                                  fill_argument(std::make_pair(0.0, 100.0),
                                                std::make_pair(true, true),
                                  1001));
  basename = ns + "_cyl_bessel_i";
  runtest<_Tp, _Tp, _Tp>(cyl_bessel_i, basename, border,
                         fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                       std::make_pair(true, true),
                         1001));


  //  Cylindrical Bessel functions (of the first kind).
  std::cout << "cyl_bessel_j" << std::endl;
  basename = "gsl_cyl_bessel_j";
  runtest<double, double, double>(gsl_sf_bessel_Jnu, basename, border,
                                  fill_argument(std::make_pair(0.0, 100.0),
                                                std::make_pair(true, true),
                                  1001));
  basename = ns + "_cyl_bessel_j";
  runtest<_Tp, _Tp, _Tp>(cyl_bessel_j, basename, border,
                         fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                       std::make_pair(true, true),
                         1001));


  //  Irregular modified cylindrical Bessel functions.
  std::cout << "cyl_bessel_k" << std::endl;
  basename = "gsl_cyl_bessel_k";
  runtest<double, double, double>(gsl_sf_bessel_Knu, basename, border,
                                  fill_argument(std::make_pair(0.0, 100.0),
                                                std::make_pair(false, true),  // Skip the pole at the origin.
                                  1001));
  basename = ns + "_cyl_bessel_k";
  runtest<_Tp, _Tp, _Tp>(cyl_bessel_k, basename, border,
                         fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                       std::make_pair(true, true),
                         1001));


  //  Cylindrical Neumann functions.
  std::cout << "cyl_neumann" << std::endl;
  basename = "gsl_cyl_neumann";
  runtest<double, double, double>(gsl_sf_bessel_Ynu, basename, border,
                                  fill_argument(std::make_pair(0.0, 100.0),
                                                std::make_pair(false, true),  // Skip the pole at the origin.
                                   1001));
  basename = ns + "_cyl_neumann";
  runtest<_Tp, _Tp, _Tp>(cyl_neumann, basename, border,
                         fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                       std::make_pair(true, true),
                         1001));


  //  Elliptic integrals of the first kind.
  std::cout << "ellint_1" << std::endl;
  basename = "gsl_ellint_1";
  runtest<double, double, double>(wrap_gsl_sf_ellint_F,
                                  basename,
                                  fill_argument(std::make_pair(-1.0, 1.0),
                                                std::make_pair(false, false)),  //  Avoid poles at |x| = 1.
                                  vphid);
  basename = ns + "_ellint_1";
  runtest<_Tp, _Tp, _Tp>(ellint_1,
                         basename,
                         fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                       std::make_pair(true, true)),
                         vphid);


  //  Elliptic integrals of the second kind.
  std::cout << "ellint_2" << std::endl;
  basename = "gsl_ellint_2";
  runtest<double, double, double>(wrap_gsl_sf_ellint_E, basename,
                                  fill_argument(std::make_pair(-1.0, 1.0),
                                                std::make_pair(false, false)),  //  Avoid poles at |x| = 1.
                                  vphid);
  basename = ns + "_ellint_2";
  runtest<_Tp, _Tp, _Tp>(ellint_2, basename,
                         fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                                std::make_pair(true, true)),
                         vphid);


  //  Elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "ellint_3" << std::endl;
  basename = "gsl_ellint_3";
  runtest<double, double, double, double>(wrap_gsl_sf_ellint_P,
                                          basename,
                                          fill_argument(std::make_pair(-1.0, 1.0),
                                                        std::make_pair(false, false)),
                                          fill_argument(std::make_pair(0.0, 1.0),
                                                        std::make_pair(true, false), 11),
                                          vphid);
  basename = ns + "_ellint_3";
  runtest<_Tp, _Tp, _Tp, _Tp>(ellint_3, basename,
                              fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                            std::make_pair(true, true)),
                              fill_argument(std::make_pair(_Tp{0}, _Tp{1}),
                                            std::make_pair(true, true), 11),
                              vphid);


  //  Exponential integrals.
  std::cout << "expint" << std::endl;
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
  runtest<_Tp, _Tp>(expint, basename,
                    fill_argument(std::make_pair(-_Tp{50}, _Tp{50}),
                                  std::make_pair(true, true)));


  //  Hermite polynomials
  std::cout << "hermite" << std::endl;
  basename = "gsl_hermite";
  //runtest<double, unsigned int, double>(gsl_sf_hermite, basename, uiorder,
  //                                      fill_argument(std::make_pair(-10.0, 10.0),
  //                                                    std::make_pair(true, true)));
  basename = ns + "_hermite";
  runtest<_Tp, unsigned int, _Tp>(hermite, basename, uiorder,
                                  fill_argument(std::make_pair(-_Tp{10}, _Tp{10}),
                                                std::make_pair(true, true)));


  //  Hypergeometric functions.
  std::cout << "hyperg" << std::endl;
  basename = "gsl_hyperg";
  runtest<double, double, double, double, double>(wrap_gsl_sf_hyperg_2F1, basename,
                                                  fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(true, true), 11),
                                                  fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(true, true), 11),
                                                  fill_argument(std::make_pair(0.0, 10.0),
                                                                  std::make_pair(false, true), 11),  //  Skip the singularity
                                                  fill_argument(std::make_pair(-1.0, 1.0),
                                                                  std::make_pair(false, true), 21));
  basename = ns + "_hyperg";
  runtest<_Tp, _Tp, _Tp, _Tp, _Tp>(hyperg, basename,
                                   fill_argument(std::make_pair(_Tp{0}, _Tp{10}),
                                                 std::make_pair(true, true), 11),
                                   fill_argument(std::make_pair(_Tp{0}, _Tp{10}),
                                                 std::make_pair(true, true), 11),
                                   fill_argument(std::make_pair(_Tp{0}, _Tp{10}),
                                                 std::make_pair(true, true), 11),
                                   fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                                 std::make_pair(false, true), 21));


  //  Laguerre polynomials.
  std::cout << "laguerre" << std::endl;
  basename = "gsl_laguerre";
  runtest<double, unsigned int, double>(wrap_gsl_sf_laguerre_n, basename,
                                        uiorder,
                                        fill_argument(std::make_pair(0.0, 100.0),
                                                      std::make_pair(true, true)));
  basename = ns + "_laguerre";
  runtest<_Tp, unsigned int, _Tp>(laguerre, basename, uiorder,
                                  fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                                std::make_pair(true, true)));


  ///  Legendre polynomials
  std::cout << "legendre" << std::endl;
  basename = "gsl_legendre";
  runtest<double, int, double>(gsl_sf_legendre_Pl, basename, iorder,
                               fill_argument(std::make_pair(-1.0, 1.0),
                                             std::make_pair(true, true),
                               1001));
  basename = ns + "_legendre";
  runtest<_Tp, unsigned int, _Tp>(legendre, basename, uiorder,
                                  fill_argument(std::make_pair(-_Tp{1}, _Tp{1}),
                                                std::make_pair(true, true),
                                  1001));


  //  Riemann zeta function.
  std::cout << "riemann_zeta" << std::endl;
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
  runtest<_Tp, _Tp>(riemann_zeta, basename,
                    fill_argument(std::make_pair(-_Tp{10}, _Tp{30}),
                                  std::make_pair(true, true), 201));


  //  Spherical Bessel functions.
  std::cout << "sph_bessel" << std::endl;
  basename = "gsl_sph_bessel";
  runtest<double, int, double>(gsl_sf_bessel_jl, basename, isborder,
                               fill_argument(std::make_pair(0.0, 100.0),
                                             std::make_pair(true, true),
                               1001));
  basename = ns + "_sph_bessel";
  runtest<_Tp, unsigned int, _Tp>(sph_bessel, basename, uisborder,
                                  fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                                std::make_pair(true, true),
                                  1001));

  //  Spherical Legendre functions.
  std::cout << "sph_legendre" << std::endl;
  basename = "gsl_sph_legendre";
  runtest<double, unsigned int, unsigned int, double>(wrap_gsl_sf_legendre_sphPlm, basename, uiorder, uiorder,
                                                      fill_argument(std::make_pair(0.0, static_cast<double>(M_PI)),
                                                                    std::make_pair(true, true),
                                                      1001));
  basename = ns + "_sph_legendre";
  runtest<_Tp, unsigned int, unsigned int, _Tp>(sph_legendre, basename, uiorder, uiorder,
                                                fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                                              std::make_pair(true, true),
                                                1001));



  //  Spherical Neumann functions.
  std::cout << "sph_neumann" << std::endl;
  basename = "gsl_sph_neumann";
  runtest<double, int, double>(gsl_sf_bessel_yl, basename, isborder,
                               fill_argument(std::make_pair(0.0, 100.0),
                                             std::make_pair(false, true),  // Skip the pole at the origin.
                                             1001));
  basename = ns + "_sph_neumann";
  runtest<_Tp, unsigned int, _Tp>(sph_neumann, basename, uisborder,
                                  fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
                                                std::make_pair(true, true),
                                  1001));


  return 0;
}

int
main()
{
  do_test<long double>();
  do_test<double>();
  do_test<float>();
}
