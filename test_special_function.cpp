
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

#include "float128.h"
#include "test_func.tcc"
#include "gsl_wrap.h"


///
///
///
template<typename _Tp>
  void
  do_test()
  {
    //  Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
    std::vector<unsigned int> uiorder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

    //  ... corresponding signed integer orders for GSL.
    std::vector<int> iorder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

    //  ... corresponding double "integer" orders for GSL.
    std::vector<double> dorder{0.0, 1.0, 2.0, 3.0, 4.0,
				 5.0, 10.0, 20.0, 50.0, 100.0};

    //  Orders for cylindrical Bessel functions.
    std::vector<_Tp> border{_Tp{0}, _Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}, _Tp{2}/_Tp{3}, _Tp{1},
			    _Tp{2}, _Tp{3}, _Tp{5}, _Tp{10}, _Tp{20}, _Tp{50}, _Tp{100}};
    //  ... and for GSL.
    std::vector<double> dborder{0.0, 1.0/3.0, 1.0/2.0, 2.0/3.0, 1.0,
				2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0};

    //  Orders for spherical Bessel functions.
    std::vector<unsigned int> uisborder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};
    //  ... and for GSL.
    std::vector<int> isborder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

    const unsigned long num_phi = 19; // 0 - 180 degrees.
    _Tp phi[num_phi];
    for (unsigned int i = 0; i < num_phi; ++i)
      phi[i] = _Tp{10} * i * static_cast<_Tp>(M_PI) / _Tp{180};
    std::vector<double> vphid(phi, phi + num_phi);
    std::vector<_Tp> vphi(phi, phi + num_phi);

#if STD
    std::string ns("std");
    using __gnu_cxx::airy_ai;
    using __gnu_cxx::airy_bi;
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
    using __gnu_cxx::hurwitz_zeta;
    using std::sph_bessel;
    using std::sph_legendre;
    using std::sph_neumann;
    using __gnu_cxx::ellint_rc;
    using __gnu_cxx::ellint_rd;
    using __gnu_cxx::ellint_rf;
    using __gnu_cxx::ellint_rj;
    using __gnu_cxx::dilog;
    using __gnu_cxx::gamma_u;
    using __gnu_cxx::ibeta;
    using __gnu_cxx::psi;
    using __gnu_cxx::sinint;
    using __gnu_cxx::cosint;
    using __gnu_cxx::sinhint;
    using __gnu_cxx::coshint;
    using __gnu_cxx::jacobi_sn;
    using __gnu_cxx::jacobi_cn;
    using __gnu_cxx::jacobi_dn;
    using __gnu_cxx::fresnel_s;
    using __gnu_cxx::fresnel_c;
    using __gnu_cxx::dawson;
    using __gnu_cxx::sinc;
    using __gnu_cxx::expint_e1;
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

#if STD
    //  Airy Ai functions.
    std::cout << "airy_ai" << std::endl;
    basename = "gsl_airy_ai";
    runtest<double, double>(gsl::airy_ai, basename,
			    fill_argument(std::make_pair(-10.0, +10.0),
					  std::make_pair(true, true), 41));
    basename = ns + "_airy_ai";
    runtest<_Tp, _Tp>(airy_ai, basename,
		      fill_argument(std::make_pair(_Tp{-10}, +_Tp{10}),
				    std::make_pair(true, true), 41));

    //  Airy Bi functions.
    std::cout << "airy_bi" << std::endl;
    basename = "gsl_airy_bi";
    runtest<double, double>(gsl::airy_bi, basename,
			    fill_argument(std::make_pair(-10.0, +10.0),
					  std::make_pair(true, true), 41));
    basename = ns + "_airy_bi";
    runtest<_Tp, _Tp>(airy_bi, basename,
		      fill_argument(std::make_pair(_Tp{-10}, +_Tp{10}),
				    std::make_pair(true, true), 41));
#endif // STD

    //  Associated Laguerre polynomials.
    //  double gsl_sf_laguerre_n(int n, double a, double x);
    std::cout << "assoc_laguerre" << std::endl;
    basename = "gsl_assoc_laguerre";
    runtest<double, int, double, double>(gsl_sf_laguerre_n, basename, iorder, dorder,
					 fill_argument(std::make_pair(0.0, 100.0),
						       std::make_pair(true, true)));
    basename = ns + "_assoc_laguerre";
    runtest<_Tp, unsigned int, unsigned int, _Tp>(assoc_laguerre, basename, uiorder, uiorder,
						  fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
								std::make_pair(true, true)));


    //  Associated Legendre functions.
    std::cout << "assoc_legendre" << std::endl;
    basename = "gsl_assoc_legendre";
    runtest<double, unsigned int, unsigned int, double>(gsl::legendre_Plm, basename, uiorder, uiorder,
							fill_argument(std::make_pair(-1.0, 1.0),
								      std::make_pair(true, true),
							1001));
    basename = ns + "_assoc_legendre";
    runtest<_Tp, unsigned int, unsigned int, _Tp>(assoc_legendre, basename, uiorder, uiorder,
						  fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
								std::make_pair(true, true),
						  1001));


    //  Beta function.
    std::cout << "beta" << std::endl;
    basename = "gsl_beta";
    runtest<double, double, double>(gsl::beta, basename,
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
    runtest<double, double>(gsl::ellint_Kcomp, basename,
			    fill_argument(std::make_pair(-1.0, 1.0),
					  std::make_pair(false, false)));  //  Avoid poles at |x| = 1.
    basename = ns + "_comp_ellint_1";
    runtest<_Tp, _Tp>(comp_ellint_1, basename,
		      fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
				    std::make_pair(true, true)));


    //  Complete elliptic integrals of the second kind.
    std::cout << "comp_ellint_2" << std::endl;
    basename = "gsl_comp_ellint_2";
    runtest<double, double>(gsl::ellint_Ecomp, basename,
			    fill_argument(std::make_pair(-1.0, 1.0),
					  std::make_pair(false, false)));  //  Avoid poles at |x| = 1.
    basename = ns + "_comp_ellint_2";
    runtest<_Tp, _Tp>(comp_ellint_2, basename,
		      fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
				    std::make_pair(true, true)));


    //  Complete elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3" << std::endl;
    basename = "gsl_comp_ellint_3";
    runtest<double, double, double>(gsl::ellint_Pcomp, basename,
				    fill_argument(std::make_pair(-1.0, 1.0),
						  std::make_pair(false, false)),
				    fill_argument(std::make_pair(0.0, 1.0),
						  std::make_pair(true, false), 11));
    basename = ns + "_comp_ellint_3";
    runtest<_Tp, _Tp, _Tp>(comp_ellint_3, basename,
			   fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
					 std::make_pair(true, true)),
			   fill_argument(std::make_pair(_Tp{0}, _Tp{1}),
					 std::make_pair(true, true), 11));


    //  Confluent hypergeometric functions.
    std::cout << "conf_hyperg" << std::endl;
    basename = "gsl_conf_hyperg";
    runtest<double, double, double, double>(gsl_sf_hyperg_1F1,
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
				fill_argument(std::make_pair(_Tp{-10}, _Tp{10}),
					      std::make_pair(true, true), 201));


    //  Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i" << std::endl;
    basename = "gsl_cyl_bessel_i";
    runtest<double, double, double>(gsl::bessel_Inu, basename, dborder,
				    fill_argument(std::make_pair(0.0, 100.0),
						  std::make_pair(true, true), 1001));
    basename = ns + "_cyl_bessel_i";
    runtest<_Tp, _Tp, _Tp>(cyl_bessel_i, basename, border,
			   fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
					 std::make_pair(true, true), 1001));


    //  Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j" << std::endl;
    basename = "gsl_cyl_bessel_j";
    runtest<double, double, double>(gsl::bessel_Jnu, basename, dborder,
				    fill_argument(std::make_pair(0.0, 100.0),
						  std::make_pair(true, true), 1001));
    basename = ns + "_cyl_bessel_j";
    runtest<_Tp, _Tp, _Tp>(cyl_bessel_j, basename, border,
			   fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
					 std::make_pair(true, true), 1001));


    //  Irregular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_k" << std::endl;
    basename = "gsl_cyl_bessel_k";
    runtest<double, double, double>(gsl::bessel_Knu, basename, dborder,
				    fill_argument(std::make_pair(0.0, 100.0),
						  std::make_pair(false, true),  // Skip the pole at the origin.
						  1001));
    basename = ns + "_cyl_bessel_k";
    runtest<_Tp, _Tp, _Tp>(cyl_bessel_k, basename, border,
			   fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
					 std::make_pair(true, true), 1001));


    //  Cylindrical Neumann functions.
    std::cout << "cyl_neumann" << std::endl;
    basename = "gsl_cyl_neumann";
    runtest<double, double, double>(gsl::bessel_Ynu, basename, dborder,
				    fill_argument(std::make_pair(0.0, 100.0),
						  std::make_pair(false, true),  // Skip the pole at the origin.
						  1001));
    basename = ns + "_cyl_neumann";
    runtest<_Tp, _Tp, _Tp>(cyl_neumann, basename, border,
			   fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
					 std::make_pair(true, true), 1001));


    //  Elliptic integrals of the first kind.
    std::cout << "ellint_1" << std::endl;
    basename = "gsl_ellint_1";
    runtest<double, double, double>(gsl::ellint_F,
				    basename,
				    fill_argument(std::make_pair(-1.0, 1.0),
						  std::make_pair(false, false)),  //  Avoid poles at |x| = 1.
						  vphid);
    basename = ns + "_ellint_1";
    runtest<_Tp, _Tp, _Tp>(ellint_1,
			   basename,
			   fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
					 std::make_pair(true, true)), vphi);


    //  Elliptic integrals of the second kind.
    std::cout << "ellint_2" << std::endl;
    basename = "gsl_ellint_2";
    runtest<double, double, double>(gsl::ellint_E, basename,
				    fill_argument(std::make_pair(-1.0, 1.0),
						  std::make_pair(false, false)),  //  Avoid poles at |x| = 1.
						  vphid);
    basename = ns + "_ellint_2";
    runtest<_Tp, _Tp, _Tp>(ellint_2, basename,
			   fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
						  std::make_pair(true, true)), vphi);


    //  Elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3" << std::endl;
    basename = "gsl_ellint_3";
    runtest<double, double, double, double>(gsl::ellint_P,
					    basename,
					    fill_argument(std::make_pair(-1.0, 1.0),
							  std::make_pair(false, false)),
					    fill_argument(std::make_pair(0.0, 1.0),
							  std::make_pair(true, false), 11), vphid);
    basename = ns + "_ellint_3";
    runtest<_Tp, _Tp, _Tp, _Tp>(ellint_3, basename,
				fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
					      std::make_pair(true, true)),
				fill_argument(std::make_pair(_Tp{0}, _Tp{1}),
					      std::make_pair(true, true), 11), vphi);


    //  Exponential integrals.
    std::cout << "expint" << std::endl;
    //  Skip the pole at 0.
    basename = "gsl_expint_neg";
    runtest<double, double>(gsl::expint_Ei, basename,
			    fill_argument(std::make_pair(-50.0, 0.0),
					  std::make_pair(true, false), 51));
    basename = "gsl_expint_pos";
    runtest<double, double>(gsl::expint_Ei, basename,
			    fill_argument(std::make_pair(0.0, 50.0),
					  std::make_pair(false, true), 51));
    basename = ns + "_expint";
    runtest<_Tp, _Tp>(expint, basename,
		      fill_argument(std::make_pair(_Tp{-50}, _Tp{50}),
				    std::make_pair(true, true)));


    //  Hermite polynomials
    std::cout << "hermite" << std::endl;
    basename = "gsl_hermite";
    runtest<double, unsigned int, double>(gsl::hermite, basename, uiorder,
    				    fill_argument(std::make_pair(-10.0, 10.0),
    						  std::make_pair(true, true)));
    basename = ns + "_hermite";
    runtest<_Tp, unsigned int, _Tp>(hermite, basename, uiorder,
				    fill_argument(std::make_pair(_Tp{-10}, _Tp{10}),
						  std::make_pair(true, true)));


    //  Hypergeometric functions.
    std::cout << "hyperg" << std::endl;
    basename = "gsl_hyperg";
    runtest<double, double, double, double, double>(gsl_sf_hyperg_2F1, basename,
						    fill_argument(std::make_pair(0.0, 10.0),
								  std::make_pair(true, true), 11),
						    fill_argument(std::make_pair(0.0, 10.0),
								  std::make_pair(true, true), 11),
						    fill_argument(std::make_pair(0.0, 10.0),
								  std::make_pair(false, true), 11),  //  Skip the singularity
						    fill_argument(std::make_pair(-1.0, 1.0),
								  std::make_pair(true, false), 21));
    basename = ns + "_hyperg";
    runtest<_Tp, _Tp, _Tp, _Tp, _Tp>(hyperg, basename,
				     fill_argument(std::make_pair(_Tp{0}, _Tp{10}),
						   std::make_pair(true, true), 11),
				     fill_argument(std::make_pair(_Tp{0}, _Tp{10}),
						   std::make_pair(true, true), 11),
				     fill_argument(std::make_pair(_Tp{0}, _Tp{10}),
						   std::make_pair(true, true), 11),
				     fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
						   std::make_pair(false, true), 21));


    //  Laguerre polynomials.
    std::cout << "laguerre" << std::endl;
    basename = "gsl_laguerre";
    runtest<double, unsigned int, double>(gsl::laguerre_n, basename,
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
					       std::make_pair(true, true), 1001));
    basename = ns + "_legendre";
    runtest<_Tp, unsigned int, _Tp>(legendre, basename, uiorder,
				    fill_argument(std::make_pair(_Tp{-1}, _Tp{1}),
						  std::make_pair(true, true), 1001));


    //  Riemann zeta function.
    std::cout << "riemann_zeta" << std::endl;
    //  Skip the pole at 1.
    basename = "gsl_riemann_zeta_neg";
    runtest<double, double>(gsl::zeta, basename,
			    fill_argument(std::make_pair(-10.0, 1.0),
					  std::make_pair(true, false), 56));
    basename = "gsl_riemann_zeta_pos";
    runtest<double, double>(gsl::zeta, basename,
			    fill_argument(std::make_pair(1.0, 30.0),
					  std::make_pair(false, true), 146));
    basename = ns + "_riemann_zeta";
    runtest<_Tp, _Tp>(riemann_zeta, basename,
		      fill_argument(std::make_pair(_Tp{-10}, _Tp{30}),
				    std::make_pair(true, true), 201));

#if STD
    //  Hurwitz zeta function.
    std::cout << "hurwitz_zeta" << std::endl;
    //  Skip the pole at 1.
    basename = "gsl_hurwitz_zeta";
    runtest<double, double, double>(gsl::hzeta, basename,
				    fill_argument(std::make_pair(1.0, 30.0),
						  std::make_pair(false, true), 146),
				    fill_argument(std::make_pair(0.0, 5.0),
						  std::make_pair(false, true), 26));
    basename = ns + "_hurwitz_zeta";
    runtest<_Tp, _Tp, _Tp>(hurwitz_zeta, basename,
			   fill_argument(std::make_pair(_Tp{1}, _Tp{30}),
					 std::make_pair(true, true), 146),
			   fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					 std::make_pair(true, true), 26));
#endif // STD


    //  Spherical Bessel functions.
    std::cout << "sph_bessel" << std::endl;
    basename = "gsl_sph_bessel";
    runtest<double, int, double>(gsl_sf_bessel_jl, basename, isborder,
				 fill_argument(std::make_pair(0.0, 100.0),
					       std::make_pair(true, true), 1001));
    basename = ns + "_sph_bessel";
    runtest<_Tp, unsigned int, _Tp>(sph_bessel, basename, uisborder,
				    fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
						  std::make_pair(true, true), 1001));

    //  Spherical Legendre functions.
    std::cout << "sph_legendre" << std::endl;
    basename = "gsl_sph_legendre";
    runtest<double, unsigned int, unsigned int, double>(gsl::legendre_sphPlm, basename, uiorder, uiorder,
							fill_argument(std::make_pair(0.0, static_cast<double>(M_PI)),
								      std::make_pair(true, true), 1001));
    basename = ns + "_sph_legendre";
    runtest<_Tp, unsigned int, unsigned int, _Tp>(sph_legendre, basename, uiorder, uiorder,
						  fill_argument(std::make_pair(_Tp{0}, _Tp{100}),
								std::make_pair(true, true), 1001));



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
						  std::make_pair(true, true), 1001));

#if STD
    //  Carlson elliptic functions R_C.
    std::cout << "ellint_rc" << std::endl;
    basename = "gsl_ellint_rc";
    runtest<double, double, double>(gsl::ellint_RC, basename,
				    fill_argument(std::make_pair(0.0, 5.0),
						  std::make_pair(false, true), 11),
				    fill_argument(std::make_pair(0.0, 5.0),
						  std::make_pair(false, true), 11));
    basename = ns + "_ellint_rc";
    runtest<_Tp, _Tp, _Tp>(ellint_rc, basename,
			   fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					 std::make_pair(false, true), 11),
			   fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					 std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_D.
    std::cout << "ellint_rd" << std::endl;
    basename = "gsl_ellint_rd";
    runtest<double, double, double, double>(gsl::ellint_RD, basename,
					    fill_argument(std::make_pair(0.0, 5.0),
							  std::make_pair(false, true), 11),
					    fill_argument(std::make_pair(0.0, 5.0),
							  std::make_pair(true, true), 11),
					    fill_argument(std::make_pair(0.0, 5.0),
							  std::make_pair(false, true), 11));
    basename = ns + "_ellint_rd";
    runtest<_Tp, _Tp, _Tp, _Tp>(ellint_rd, basename,
				fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					      std::make_pair(false, true), 11),
				fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					      std::make_pair(true, true), 11),
				fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					      std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_F.
    std::cout << "ellint_rf" << std::endl;
    basename = "gsl_ellint_rf";
    runtest<double, double, double, double>(gsl::ellint_RF, basename,
					    fill_argument(std::make_pair(0.0, 5.0),
							  std::make_pair(false, true), 11),
					    fill_argument(std::make_pair(0.0, 5.0),
							  std::make_pair(true, true), 11),
					    fill_argument(std::make_pair(0.0, 5.0),
							  std::make_pair(false, true), 11));
    basename = ns + "_ellint_rf";
    runtest<_Tp, _Tp, _Tp, _Tp>(ellint_rf, basename,
				fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					      std::make_pair(false, true), 11),
				fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					      std::make_pair(true, true), 11),
				fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					      std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_J.
    std::cout << "ellint_rj" << std::endl;
    basename = "gsl_ellint_rj";
    runtest<double, double, double, double, double>(gsl::ellint_RJ, basename,
						    fill_argument(std::make_pair(0.0, 5.0),
								  std::make_pair(false, true), 11),
						    fill_argument(std::make_pair(0.0, 5.0),
								  std::make_pair(true, true), 11),
						    fill_argument(std::make_pair(0.0, 5.0),
								  std::make_pair(false, true), 11),
						    fill_argument(std::make_pair(0.0, 5.0),
								  std::make_pair(false, true), 11));
    basename = ns + "_ellint_rj";
    runtest<_Tp, _Tp, _Tp, _Tp, _Tp>(ellint_rj, basename,
				     fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
						   std::make_pair(false, true), 11),
				     fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
						   std::make_pair(true, true), 11),
				     fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
						   std::make_pair(false, true), 11),
				     fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
						   std::make_pair(false, true), 11));

    //  Dilogarithm functions.
    std::cout << "dilog" << std::endl;
    basename = "gsl_dilog";
    runtest<double, double>(gsl::dilog, basename,
			    fill_argument(std::make_pair(-10.0, 1.0),
					  std::make_pair(true, true), 23));

    basename = ns + "_dilog";
    runtest<_Tp, _Tp>(dilog, basename,
		      fill_argument(std::make_pair(_Tp{-10}, _Tp{1}),
				    std::make_pair(true, true), 23));

    //  Upper incomplete Gamma functions.
    std::cout << "gamma_u" << std::endl;
    basename = "gsl_gamma_u";
    runtest<double, double, double>(gsl::gamma_inc, basename,
				    fill_argument(std::make_pair(0.0, 5.0),
						  std::make_pair(false, true), 11),
				    fill_argument(std::make_pair(0.0, 5.0),
						  std::make_pair(true, true), 11));

    basename = ns + "_gamma_u";
    runtest<_Tp, _Tp, _Tp>(gamma_u, basename,
			   fill_argument(std::make_pair(_Tp{0}, +_Tp{5}),
					 std::make_pair(false, true), 11),
			   fill_argument(std::make_pair(_Tp{0}, +_Tp{5}),
					 std::make_pair(true, true), 11));

    //  Incomplete Beta functions.
    std::cout << "ibeta" << std::endl;
    basename = "gsl_ibeta";
    runtest<double, double, double, double>(gsl::beta_inc, basename,
					    fill_argument(std::make_pair(0.0, 5.0),
							  std::make_pair(false, true), 11),
					    fill_argument(std::make_pair(5.0, 0.0),
							  std::make_pair(true, false), 11),
					    fill_argument(std::make_pair(0.0, 1.0),
							  std::make_pair(false, false), 21));
    basename = ns + "_ibeta";
    runtest<_Tp, _Tp, _Tp, _Tp>(ibeta, basename,
				fill_argument(std::make_pair(_Tp{0}, _Tp{5}),
					      std::make_pair(false, true), 11),
				fill_argument(std::make_pair(_Tp{5}, _Tp{0}),
					      std::make_pair(false, true), 11),
				fill_argument(std::make_pair(_Tp{0}, _Tp{1}),
					      std::make_pair(false, false), 21));

    //  Digamma or psi functions.
    std::cout << "psi" << std::endl;
    basename = "gsl_psi";
    runtest<double, double>(gsl::psi, basename,
			    fill_argument(std::make_pair(-9.875, 10.125),
					  std::make_pair(true, true), 41));
    basename = ns + "_psi";
    runtest<_Tp, _Tp>(psi, basename,
		      fill_argument(std::make_pair(_Tp{-9.9375}, _Tp{10.0625}),
				    std::make_pair(true, true), 801));

    //  Sine integral or Si functions.
    std::cout << "gsl_sinint" << std::endl;
    basename = "sinint";
    runtest<double, double>(gsl_sf_Si, basename,
			    fill_argument(std::make_pair(0.0, +10.0),
					  std::make_pair(false, true), 101));
    basename = ns + "_sinint";
    runtest<_Tp, _Tp>(sinint, basename,
		      fill_argument(std::make_pair(_Tp{0}, _Tp{+10}),
				    std::make_pair(false, true), 101));

    //  Cosine integral or Ci functions.
    std::cout << "gsl_cosint" << std::endl;
    basename = "cosint";
    runtest<double, double>(gsl_sf_Ci, basename,
			    fill_argument(std::make_pair(0.0, +10.0),
					  std::make_pair(false, true), 101));
    basename = ns + "_cosint";
    runtest<_Tp, _Tp>(cosint, basename,
		      fill_argument(std::make_pair(_Tp{0}, _Tp{+10}),
				    std::make_pair(false, true), 101));

    //  Hyperbolic sine integral or Shi functions.
    std::cout << "gsl_sinhint" << std::endl;
    basename = "sinhint";
    runtest<double, double>(gsl_sf_Shi, basename,
			    fill_argument(std::make_pair(0.0, +5.0),
					  std::make_pair(false, true), 101));
    basename = ns + "_sinhint";
    runtest<_Tp, _Tp>(sinhint, basename,
		      fill_argument(std::make_pair(_Tp{0}, _Tp{+5}),
				    std::make_pair(false, true), 101));

    //  Hyperbolic cosine integral or Chi functions.
    std::cout << "gsl_coshint" << std::endl;
    basename = "coshint";
    runtest<double, double>(gsl_sf_Chi, basename,
			    fill_argument(std::make_pair(0.0, +5.0),
					  std::make_pair(false, true), 101));
    basename = ns + "_coshint";
    runtest<_Tp, _Tp>(coshint, basename,
		      fill_argument(std::make_pair(_Tp{0}, _Tp{+5}),
				    std::make_pair(false, true), 101));

    //  Dawson integral.
    std::cout << "gsl_dawson" << std::endl;
    basename = "dawson";
    runtest<double, double>(gsl::dawson, basename,
			    fill_argument(std::make_pair(0.0, +5.0),
					  std::make_pair(false, true), 101));
    basename = ns + "_dawson";
    runtest<_Tp, _Tp>(dawson, basename,
		      fill_argument(std::make_pair(_Tp{0}, _Tp{+5}),
				    std::make_pair(false, true), 101));

    //  Exponential integral E1.
    std::cout << "gsl_expint_e1" << std::endl;
    basename = "expint_e1";
    runtest<double, double>(gsl::expint_e1, basename,
			    fill_argument(std::make_pair(0.0, +5.0),
					  std::make_pair(false, true), 101));
    basename = ns + "_expint_e1";
    runtest<_Tp, _Tp>(expint_e1, basename,
		      fill_argument(std::make_pair(_Tp{0}, _Tp{+5}),
				    std::make_pair(false, true), 101));

    //  Sine cardinal function.
    std::cout << "gsl_sinc" << std::endl;
    basename = "sinc";
    runtest<double, double>(gsl::sinc, basename,
			    fill_argument(std::make_pair(-20.0, +20.0),
					  std::make_pair(false, true), 401));
    basename = ns + "_sinc";
    runtest<_Tp, _Tp>(sinc, basename,
		      fill_argument(std::make_pair(_Tp{20}, _Tp{+20}),
				    std::make_pair(false, true), 401));

#endif // STD

    return;
  }

int
main()
{
#if LOCAL
  do_test<__float128>();
#endif
  do_test<long double>();
  do_test<double>();
  do_test<float>();
}
