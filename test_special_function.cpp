
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
template<typename Real>
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
    std::vector<Real> border{0, Real{1}/Real{3}, Real{1}/Real{2}, Real{2}/Real{3}, 1,
			    2, 3, 5, 10, 20, 50, 100};
    //  ... and for GSL.
    std::vector<double> dborder{0.0, 1.0/3.0, 1.0/2.0, 2.0/3.0, 1.0,
				2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0};

    //  Orders for spherical Bessel functions.
    std::vector<unsigned int> uisborder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};
    //  ... and for GSL.
    std::vector<int> isborder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

    const unsigned long num_phi = 19; // 0 - 180 degrees.
    Real phi[num_phi];
    for (unsigned int i = 0; i < num_phi; ++i)
      phi[i] = Real{10} * i * static_cast<Real>(M_PI) / Real{180};
    std::vector<double> vphid(phi, phi + num_phi);
    std::vector<Real> vphi(phi, phi + num_phi);

#if STD
    std::string ns("std");
    using __gnu_cxx::airy_ai;
    using __gnu_cxx::airy_bi;
    using       std::assoc_laguerre;
    using       std::assoc_legendre;
    using       std::beta;
    using __gnu_cxx::bincoef;
    using __gnu_cxx::chebyshev_t;
    using __gnu_cxx::chebyshev_u;
    using __gnu_cxx::chebyshev_v;
    using __gnu_cxx::chebyshev_w;
    using __gnu_cxx::clausen;
    using __gnu_cxx::clausen_c;
    using __gnu_cxx::clausen_s;
    using __gnu_cxx::comp_ellint_d;
    using       std::comp_ellint_1;
    using       std::comp_ellint_2;
    using       std::comp_ellint_3;
    using __gnu_cxx::conf_hyperg;
    using __gnu_cxx::conf_hyperg_lim;
    using __gnu_cxx::coshint;
    using __gnu_cxx::cosint;
    using       std::cyl_bessel_i;
    using       std::cyl_bessel_j;
    using       std::cyl_bessel_k;
    using __gnu_cxx::cyl_hankel_1;
    using __gnu_cxx::cyl_hankel_2;
    using       std::cyl_neumann;
    using __gnu_cxx::dawson;
    using __gnu_cxx::dilog;
    using __gnu_cxx::dirichlet_beta;
    using __gnu_cxx::dirichlet_eta;
    using __gnu_cxx::double_factorial;
    using       std::ellint_1;
    using       std::ellint_2;
    using       std::ellint_3;
    using __gnu_cxx::ellint_d;
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
    using __gnu_cxx::heuman_lambda;
    using __gnu_cxx::hurwitz_zeta;
    using __gnu_cxx::hyperg;
    using __gnu_cxx::ibeta;
    using __gnu_cxx::jacobi;
    using __gnu_cxx::jacobi_sn;
    using __gnu_cxx::jacobi_cn;
    using __gnu_cxx::jacobi_dn;
    using __gnu_cxx::jacobi_zeta;
    using       std::laguerre;
    using __gnu_cxx::lbincoef;
    using __gnu_cxx::ldouble_factorial;
    using       std::legendre;
    using __gnu_cxx::legendre_q;
    using __gnu_cxx::lfactorial;
    using __gnu_cxx::lpochhammer_l;
    using __gnu_cxx::lpochhammer_u;
    using __gnu_cxx::owens_t;
    using __gnu_cxx::pgamma;
    using __gnu_cxx::pochhammer_l;
    using __gnu_cxx::pochhammer_u;
    using __gnu_cxx::psi;
    using __gnu_cxx::qgamma;
    using __gnu_cxx::radpoly;
    using       std::riemann_zeta;
    using __gnu_cxx::sinhc;
    using __gnu_cxx::sinhc_pi;
    using __gnu_cxx::sinc;
    using __gnu_cxx::sinc_pi;
    using __gnu_cxx::sinhint;
    using __gnu_cxx::sinint;
    using       std::sph_bessel;
    using __gnu_cxx::sph_bessel_i;
    using __gnu_cxx::sph_bessel_k;
    using __gnu_cxx::sph_hankel_1;
    using __gnu_cxx::sph_hankel_2;
    using __gnu_cxx::sph_harmonic;
    using       std::sph_legendre;
    using       std::sph_neumann;
    using __gnu_cxx::zernike;
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

    std::string basename;

#if STD
    //  Airy Ai functions.
    std::cout << "airy_ai" << std::endl;
    basename = "gsl_airy_ai";
    runtest(gsl::airy_ai, basename,
	    fill_argument(std::make_pair(-10.0, +10.0),
			  std::make_pair(true, true), 41));
    basename = ns + "_airy_ai";
    runtest(airy_ai, basename,
	    fill_argument(std::make_pair(Real{-10}, Real{+10}),
	        	  std::make_pair(true, true), 41));

    //  Airy Bi functions.
    std::cout << "airy_bi" << std::endl;
    basename = "gsl_airy_bi";
    runtest(gsl::airy_bi, basename,
	    fill_argument(std::make_pair(-10.0, +10.0),
			  std::make_pair(true, true), 41));
    basename = ns + "_airy_bi";
    runtest(airy_bi, basename,
	    fill_argument(std::make_pair(Real{-10}, Real{+10}),
	        	  std::make_pair(true, true), 41));
#endif // STD

    //  Associated Laguerre polynomials.
    //  double gsl_sf_laguerre_n(int n, double a, double x);
    std::cout << "assoc_laguerre" << std::endl;
    basename = "gsl_assoc_laguerre";
    runtest(gsl_sf_laguerre_n, basename, iorder, dorder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(true, true)));
    basename = ns + "_assoc_laguerre";
    runtest(assoc_laguerre, basename, uiorder, uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		  	  std::make_pair(true, true)));


    //  Associated Legendre functions.
    std::cout << "assoc_legendre" << std::endl;
    basename = "gsl_assoc_legendre";
    runtest(gsl::legendre_Plm, basename, uiorder, uiorder,
	    fill_argument(std::make_pair(-1.0, 1.0),
	    		  std::make_pair(true, true), 1001));
    basename = ns + "_assoc_legendre";
    runtest(assoc_legendre, basename, uiorder, uiorder,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		  	  std::make_pair(true, true), 1001));


    //  Beta function.
    std::cout << "beta" << std::endl;
    basename = "gsl_beta";
    runtest(gsl::beta, basename,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true)),
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true)));
    basename = ns + "_beta";
    runtest(beta, basename,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(false, true), 101),
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(false, true), 101));


    //  Complete elliptic integrals of the first kind.
    std::cout << "comp_ellint_1" << std::endl;
    basename = ns + "_comp_ellint_1";
    basename = "gsl_comp_ellint_1";
    runtest(gsl::comp_ellint_1, basename,
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101));  //  Avoid poles at |x| = 1.
    basename = ns + "_comp_ellint_1";
    runtest(comp_ellint_1, basename,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	        	  std::make_pair(true, true), 101));


    //  Complete elliptic integrals of the second kind.
    std::cout << "comp_ellint_2" << std::endl;
    basename = "gsl_comp_ellint_2";
    runtest(gsl::comp_ellint_2, basename,
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101));  //  Avoid poles at |x| = 1.
    basename = ns + "_comp_ellint_2";
    runtest(comp_ellint_2, basename,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	        	  std::make_pair(true, true), 101));


    //  Complete elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3" << std::endl;
    basename = "gsl_comp_ellint_3";
    runtest(gsl::comp_ellint_3, basename,
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false)),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(true, false), 11));
    basename = ns + "_comp_ellint_3";
    runtest(comp_ellint_3, basename,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 	  std::make_pair(true, true)),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
		 	  std::make_pair(true, true), 11));


    //  Confluent hypergeometric functions.
    std::cout << "conf_hyperg" << std::endl;
    basename = "gsl_conf_hyperg";
    runtest(gsl_sf_hyperg_1F1,
	    basename,
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(false, true), 11),  //  Skip the singularity
	    fill_argument(std::make_pair(-10.0, 10.0),
			  std::make_pair(true, true), 201));
    basename = ns + "_conf_hyperg";
    runtest(conf_hyperg, basename,
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{-10}, Real{10}),
	    		  std::make_pair(true, true), 201));


    //  Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i" << std::endl;
    basename = "gsl_cyl_bessel_i";
    runtest(gsl::cyl_bessel_i, basename, dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(true, true), 1001));
    basename = ns + "_cyl_bessel_i";
    runtest(cyl_bessel_i, basename, border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j" << std::endl;
    basename = "gsl_cyl_bessel_j";
    runtest(gsl::cyl_bessel_j, basename, dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(true, true), 1001));
    basename = ns + "_cyl_bessel_j";
    runtest(cyl_bessel_j, basename, border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Irregular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_k" << std::endl;
    basename = "gsl_cyl_bessel_k";
    runtest(gsl::cyl_bessel_k, basename, dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true),  // Skip the pole at the origin.
			  1001));
    basename = ns + "_cyl_bessel_k";
    runtest(cyl_bessel_k, basename, border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Cylindrical Neumann functions.
    std::cout << "cyl_neumann" << std::endl;
    basename = "gsl_cyl_neumann";
    runtest(gsl::cyl_neumann, basename, dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true),  // Skip the pole at the origin.
			  1001));
    basename = ns + "_cyl_neumann";
    runtest(cyl_neumann, basename, border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Elliptic integrals of the first kind.
    std::cout << "ellint_1" << std::endl;
    basename = "gsl_ellint_1";
    runtest(gsl::ellint_1,
	    basename,
	    fill_argument(std::make_pair(-1.0, 1.0),
	        	  std::make_pair(false, false), 101),  //  Avoid poles at |x| = 1.
	        	  vphid);
    basename = ns + "_ellint_1";
    runtest(ellint_1,
	    basename,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 	  std::make_pair(true, true), 101), vphi);


    //  Elliptic integrals of the second kind.
    std::cout << "ellint_2" << std::endl;
    basename = "gsl_ellint_2";
    runtest(gsl::ellint_2, basename,
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101),  //  Avoid poles at |x| = 1.
			  vphid);
    basename = ns + "_ellint_2";
    runtest(ellint_2, basename,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 		   std::make_pair(true, true), 101), vphi);


    //  Elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3" << std::endl;
    basename = "gsl_ellint_3";
    runtest(gsl::ellint_3,
	    basename,
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false)),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(true, false), 11), vphid);
    basename = ns + "_ellint_3";
    runtest(ellint_3, basename,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	    		  std::make_pair(true, true)),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
	    		  std::make_pair(true, true), 11), vphi);


    //  Exponential integrals.
    std::cout << "expint" << std::endl;
    //  Skip the pole at 0.
    basename = "gsl_expint_neg";
    runtest(gsl::expint_Ei, basename,
	    fill_argument(std::make_pair(-50.0, 0.0),
			  std::make_pair(true, false), 51));
    basename = "gsl_expint_pos";
    runtest(gsl::expint_Ei, basename,
	    fill_argument(std::make_pair(0.0, 50.0),
			  std::make_pair(false, true), 51));
    basename = ns + "_expint";
    runtest(expint, basename,
	    fill_argument(std::make_pair(Real{-50}, Real{50}),
	        	  std::make_pair(true, true), 101));


    //  Hermite polynomials
    std::cout << "hermite" << std::endl;
    basename = "gsl_hermite";
    runtest(gsl::hermite, basename, uiorder,
    	    fill_argument(std::make_pair(-10.0, 10.0),
    			  std::make_pair(true, true), 101));
    basename = ns + "_hermite";
    runtest(hermite, basename, uiorder,
	    fill_argument(std::make_pair(Real{-10}, Real{10}),
			  std::make_pair(true, true), 101));


    //  Hypergeometric functions.
    std::cout << "hyperg" << std::endl;
    basename = "gsl_hyperg";
    runtest(gsl_sf_hyperg_2F1, basename,
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(false, true), 11),  //  Skip the singularity
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(true, false), 21));
    basename = ns + "_hyperg";
    runtest(hyperg, basename,
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	        	  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	        	  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	        	  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	        	  std::make_pair(false, true), 21));


    //  Laguerre polynomials.
    std::cout << "laguerre" << std::endl;
    basename = "gsl_laguerre";
    runtest(gsl::laguerre, basename,
	    uiorder,
	    fill_argument(std::make_pair(0.0, 100.0),
		  	  std::make_pair(true, true), 101));
    basename = ns + "_laguerre";
    runtest(laguerre, basename, uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 101));


    ///  Legendre polynomials
    std::cout << "legendre" << std::endl;
    basename = "gsl_legendre";
    runtest(gsl_sf_legendre_Pl, basename, iorder,
	    fill_argument(std::make_pair(-1.0, 1.0),
	   		  std::make_pair(true, true), 1001));
    basename = ns + "_legendre";
    runtest(legendre, basename, uiorder,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
			  std::make_pair(true, true), 1001));


    //  Riemann zeta function.
    std::cout << "riemann_zeta" << std::endl;
    //  Skip the pole at 1.
    basename = "gsl_riemann_zeta_neg";
    runtest(gsl::riemann_zeta, basename,
	    fill_argument(std::make_pair(-10.0, 1.0),
			  std::make_pair(true, false), 56));
    basename = "gsl_riemann_zeta_pos";
    runtest(gsl::riemann_zeta, basename,
	    fill_argument(std::make_pair(1.0, 30.0),
			  std::make_pair(false, true), 146));
    basename = ns + "_riemann_zeta";
    runtest(riemann_zeta, basename,
	    fill_argument(std::make_pair(Real{-10}, Real{30}),
	        	  std::make_pair(true, true), 201));

#if STD
    //  Hurwitz zeta function.
    std::cout << "hurwitz_zeta" << std::endl;
    //  Skip the pole at 1.
    basename = "gsl_hurwitz_zeta";
    runtest(gsl::hurwitz_zeta, basename,
	    fill_argument(std::make_pair(1.0, 30.0),
			  std::make_pair(false, true), 146),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 26));
    basename = ns + "_hurwitz_zeta";
    runtest(hurwitz_zeta, basename,
	    fill_argument(std::make_pair(Real{1}, Real{30}),
		 	  std::make_pair(true, true), 146),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(true, true), 26));
#endif // STD


    //  Spherical Bessel functions.
    std::cout << "sph_bessel" << std::endl;
    basename = "gsl_sph_bessel";
    runtest(gsl_sf_bessel_jl, basename, isborder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(true, true), 1001));
    basename = ns + "_sph_bessel";
    runtest(sph_bessel, basename, uisborder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 1001));

    //  Spherical Legendre functions.
    std::cout << "sph_legendre" << std::endl;
    basename = "gsl_sph_legendre";
    runtest(gsl::legendre_sphPlm, basename, uiorder, uiorder,
	    fill_argument(std::make_pair(0.0, static_cast<double>(M_PI)),
	    		  std::make_pair(true, true), 1001));
    basename = ns + "_sph_legendre";
    runtest(sph_legendre, basename, uiorder, uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		  	  std::make_pair(true, true), 1001));



    //  Spherical Neumann functions.
    std::cout << "sph_neumann" << std::endl;
    basename = "gsl_sph_neumann";
    runtest(gsl_sf_bessel_yl, basename, isborder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(false, true),  // Skip the pole at the origin.
	   		  1001));
    basename = ns + "_sph_neumann";
    runtest(sph_neumann, basename, uisborder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 1001));

#if STD
    //  Carlson elliptic functions R_C.
    std::cout << "ellint_rc" << std::endl;
    basename = "gsl_ellint_rc";
    runtest(gsl::ellint_rc, basename,
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    basename = ns + "_ellint_rc";
    runtest(ellint_rc, basename,
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_D.
    std::cout << "ellint_rd" << std::endl;
    basename = "gsl_ellint_rd";
    runtest(gsl::ellint_rd, basename,
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    basename = ns + "_ellint_rd";
    runtest(ellint_rd, basename,
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_F.
    std::cout << "ellint_rf" << std::endl;
    basename = "gsl_ellint_rf";
    runtest(gsl::ellint_rf, basename,
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    basename = ns + "_ellint_rf";
    runtest(ellint_rf, basename,
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_J.
    std::cout << "ellint_rj" << std::endl;
    basename = "gsl_ellint_rj";
    runtest(gsl::ellint_rj, basename,
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    basename = ns + "_ellint_rj";
    runtest(ellint_rj, basename,
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	        	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	        	  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	        	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	        	  std::make_pair(false, true), 11));

    //  Dilogarithm functions.
    std::cout << "dilog" << std::endl;
    basename = "gsl_dilog";
    runtest(gsl::dilog, basename,
	    fill_argument(std::make_pair(-10.0, 1.0),
			  std::make_pair(true, true), 23));

    basename = ns + "_dilog";
    runtest<Real, Real>(dilog, basename,
		      fill_argument(std::make_pair(Real{-10}, Real{1}),
				    std::make_pair(true, true), 23));

    //  Upper incomplete Gamma functions.
    std::cout << "gamma_u" << std::endl;
    basename = "gsl_gamma_u";
    runtest(gsl::gamma_u, basename,
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11));

    basename = ns + "_gamma_u";
    runtest(gamma_u, basename,
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
		 	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
		 	  std::make_pair(true, true), 11));

    //  Incomplete Beta functions.
    std::cout << "ibeta" << std::endl;
    basename = "gsl_ibeta";
    runtest(gsl::ibeta, basename,
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(5.0, 0.0),
			  std::make_pair(true, false), 11),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(false, false), 21));
    basename = ns + "_ibeta";
    runtest(ibeta, basename,
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{5}, Real{0}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
	    		  std::make_pair(false, false), 21));

    //  Digamma or psi functions.
    std::cout << "psi" << std::endl;
    basename = "gsl_psi";
    runtest(gsl::psi, basename,
	    fill_argument(std::make_pair(-9.875, 10.125),
			  std::make_pair(true, true), 41));
    basename = ns + "_psi";
    runtest(psi, basename,
	    fill_argument(std::make_pair(Real{-9.9375}, Real{10.0625}),
	        	  std::make_pair(true, true), 801));

    //  Sine integral or Si functions.
    std::cout << "gsl_sinint" << std::endl;
    basename = "sinint";
    runtest(gsl_sf_Si, basename,
	    fill_argument(std::make_pair(0.0, +10.0),
			  std::make_pair(false, true), 101));
    basename = ns + "_sinint";
    runtest(sinint, basename,
	    fill_argument(std::make_pair(Real{0}, Real{+10}),
	        	  std::make_pair(false, true), 101));

    //  Cosine integral or Ci functions.
    std::cout << "gsl_cosint" << std::endl;
    basename = "cosint";
    runtest(gsl_sf_Ci, basename,
	    fill_argument(std::make_pair(0.0, +10.0),
			  std::make_pair(false, true), 101));
    basename = ns + "_cosint";
    runtest(cosint, basename,
	    fill_argument(std::make_pair(Real{0}, Real{+10}),
	        	  std::make_pair(false, true), 101));

    //  Hyperbolic sine integral or Shi functions.
    std::cout << "gsl_sinhint" << std::endl;
    basename = "sinhint";
    runtest(gsl_sf_Shi, basename,
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    basename = ns + "_sinhint";
    runtest(sinhint, basename,
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Hyperbolic cosine integral or Chi functions.
    std::cout << "gsl_coshint" << std::endl;
    basename = "coshint";
    runtest(gsl_sf_Chi, basename,
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    basename = ns + "_coshint";
    runtest(coshint, basename,
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Dawson integral.
    std::cout << "gsl_dawson" << std::endl;
    basename = "dawson";
    runtest(gsl::dawson, basename,
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    basename = ns + "_dawson";
    runtest(dawson, basename,
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Exponential integral E1.
    std::cout << "gsl_expint_e1" << std::endl;
    basename = "expint_e1";
    runtest(gsl::expint_E1, basename,
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    basename = ns + "_expint_e1";
    runtest(expint_e1, basename,
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Normalized sinus cardinal function.
    std::cout << "gsl_sinc" << std::endl;
    basename = "sinc";
    runtest(gsl::sinc, basename,
	    fill_argument(std::make_pair(-20.0, +20.0),
			  std::make_pair(false, true), 401));
    basename = ns + "_sinc";
    runtest(sinc, basename,
	    fill_argument(std::make_pair(Real{20}, Real{+20}),
	        	  std::make_pair(false, true), 401));

    //  Sinus cardinal function.
    std::cout << "gsl_sinc_pi" << std::endl;
    basename = "sinc_pi";
    runtest(gsl::sinc_pi, basename,
	    fill_argument(std::make_pair(-20.0, +20.0),
			  std::make_pair(false, true), 401));
    basename = ns + "_sinc_pi";
    runtest(sinc_pi, basename,
	    fill_argument(std::make_pair(Real{20}, Real{+20}),
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
