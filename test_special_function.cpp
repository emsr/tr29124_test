
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

#include <bits/float128.h>
#include "test_func.tcc"
#include "wrap_gsl.h"
#include "wrap_burkhardt.h"


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
    using __gnu_cxx::bernoulli;
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
    using __gnu_cxx::expint;
    using __gnu_cxx::factorial;
    using __gnu_cxx::fresnel_c;
    using __gnu_cxx::fresnel_s;
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
    using __gnu_cxx::lpochhammer_lower;
    using __gnu_cxx::lpochhammer;
    using __gnu_cxx::owens_t;
    using __gnu_cxx::pgamma;
    using __gnu_cxx::pochhammer_lower;
    using __gnu_cxx::pochhammer;
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
    using __gnu_cxx::tgamma_lower;
    using __gnu_cxx::tgamma;
    using __gnu_cxx::theta_1;
    using __gnu_cxx::theta_2;
    using __gnu_cxx::theta_3;
    using __gnu_cxx::theta_4;
    using __gnu_cxx::theta_s;
    using __gnu_cxx::theta_c;
    using __gnu_cxx::theta_d;
    using __gnu_cxx::theta_n;
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
    runtest(gsl::airy_ai, "gsl_airy_ai",
	    fill_argument(std::make_pair(-10.0, +10.0),
			  std::make_pair(true, true), 41));
    runtest(airy_ai, ns + "_airy_ai",
	    fill_argument(std::make_pair(Real{-10}, Real{+10}),
	        	  std::make_pair(true, true), 41));

    //  Airy Bi functions.
    std::cout << "airy_bi" << std::endl;
    runtest(gsl::airy_bi, "gsl_airy_bi",
	    fill_argument(std::make_pair(-10.0, +10.0),
			  std::make_pair(true, true), 41));
    runtest(airy_bi, ns + "_airy_bi",
	    fill_argument(std::make_pair(Real{-10}, Real{+10}),
	        	  std::make_pair(true, true), 41));
#endif // STD

    //  Associated Laguerre polynomials.
    //  double gsl_sf_laguerre_n(int n, double a, double x);
    std::cout << "assoc_laguerre" << std::endl;
    runtest(gsl_sf_laguerre_n, "gsl_assoc_laguerre", iorder, dorder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(true, true)));
    runtest(assoc_laguerre, ns + "_assoc_laguerre", uiorder, uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		  	  std::make_pair(true, true)));


    //  Associated Legendre functions.
    std::cout << "assoc_legendre" << std::endl;
    runtest(gsl::assoc_legendre, "gsl_assoc_legendre", uiorder, uiorder,
	    fill_argument(std::make_pair(-1.0, 1.0),
	    		  std::make_pair(true, true), 1001));
    runtest(assoc_legendre, ns + "_assoc_legendre", uiorder, uiorder,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		  	  std::make_pair(true, true), 1001));


    //  Beta function.
    std::cout << "beta" << std::endl;
    runtest(gsl::beta, "gsl_beta",
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true)),
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true)));
    runtest(beta, ns + "_beta",
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(false, true), 101),
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(false, true), 101));


    //  Complete elliptic integrals of the first kind.
    std::cout << "comp_ellint_1" << std::endl;
    runtest(gsl::comp_ellint_1, "gsl_comp_ellint_1",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101));  //  Avoid poles at |x| = 1.
    runtest(comp_ellint_1, ns + "_comp_ellint_1",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	        	  std::make_pair(true, true), 101));


    //  Complete elliptic integrals of the second kind.
    std::cout << "comp_ellint_2" << std::endl;
    runtest(gsl::comp_ellint_2, "gsl_comp_ellint_2",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101));  //  Avoid poles at |x| = 1.
    runtest(comp_ellint_2, ns + "_comp_ellint_2",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	        	  std::make_pair(true, true), 101));


    //  Complete elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3" << std::endl;
    runtest(gsl::comp_ellint_3, "gsl_comp_ellint_3",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false)),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(true, false), 11));
    runtest(comp_ellint_3, ns + "_comp_ellint_3",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 	  std::make_pair(true, true)),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
		 	  std::make_pair(true, true), 11));


    //  Confluent hypergeometric functions.
    std::cout << "conf_hyperg" << std::endl;
    runtest(gsl_sf_hyperg_1F1,
	    "gsl_conf_hyperg",
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(false, true), 11),  //  Skip the singularity
	    fill_argument(std::make_pair(-10.0, 10.0),
			  std::make_pair(true, true), 201));
    runtest(conf_hyperg, ns + "_conf_hyperg",
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{-10}, Real{10}),
	    		  std::make_pair(true, true), 201));


    //  Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i" << std::endl;
    runtest(gsl::cyl_bessel_i, "gsl_cyl_bessel_i", dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(true, true), 1001));
    runtest(cyl_bessel_i, ns + "_cyl_bessel_i", border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j" << std::endl;
    runtest(gsl::cyl_bessel_j, "gsl_cyl_bessel_j", dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(true, true), 1001));
    runtest(cyl_bessel_j, ns + "_cyl_bessel_j", border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Irregular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_k" << std::endl;
    runtest(gsl::cyl_bessel_k, "gsl_cyl_bessel_k", dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true),  // Skip the pole at the origin.
			  1001));
    runtest(cyl_bessel_k, ns + "_cyl_bessel_k", border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Cylindrical Neumann functions.
    std::cout << "cyl_neumann" << std::endl;
    runtest(gsl::cyl_neumann, "gsl_cyl_neumann", dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true),  // Skip the pole at the origin.
			  1001));
    runtest(cyl_neumann, ns + "_cyl_neumann", border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Elliptic integrals of the first kind.
    std::cout << "ellint_1" << std::endl;
    runtest(gsl::ellint_1, "gsl_ellint_1",
	    fill_argument(std::make_pair(-1.0, 1.0),
	        	  std::make_pair(false, false), 101),  //  Avoid poles at |x| = 1.
	        	  vphid);
    runtest(ellint_1, ns + "_ellint_1",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 	  std::make_pair(true, true), 101), vphi);


    //  Elliptic integrals of the second kind.
    std::cout << "ellint_2" << std::endl;
    runtest(gsl::ellint_2, "gsl_ellint_2",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101),  //  Avoid poles at |x| = 1.
			  vphid);
    runtest(ellint_2, ns + "_ellint_2",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 		   std::make_pair(true, true), 101), vphi);


    //  Elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3" << std::endl;
    runtest(gsl::ellint_3, "gsl_ellint_3",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false)),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(true, false), 11), vphid);
    runtest(ellint_3, ns + "_ellint_3",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	    		  std::make_pair(true, true)),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
	    		  std::make_pair(true, true), 11), vphi);


    //  Exponential integrals.
    std::cout << "expint" << std::endl;
    //  Skip the pole at 0.
    runtest(gsl::expint, "gsl_expint_neg",
	    fill_argument(std::make_pair(-50.0, 0.0),
			  std::make_pair(true, false), 51));
    runtest(gsl::expint, "gsl_expint_pos",
	    fill_argument(std::make_pair(0.0, 50.0),
			  std::make_pair(false, true), 51));
    runtest(expint, ns + "_expint",
	    fill_argument(std::make_pair(Real{-50}, Real{50}),
	        	  std::make_pair(true, true), 101));


    //  Hermite polynomials
    std::cout << "hermite" << std::endl;
    runtest(gsl::hermite, "gsl_hermite", uiorder,
    	    fill_argument(std::make_pair(-10.0, 10.0),
    			  std::make_pair(true, true), 101));
    runtest(hermite, ns + "_hermite", uiorder,
	    fill_argument(std::make_pair(Real{-10}, Real{10}),
			  std::make_pair(true, true), 101));


    //  Hypergeometric functions.
    std::cout << "hyperg" << std::endl;
    runtest(gsl_sf_hyperg_2F1, "gsl_hyperg",
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(false, true), 11),  //  Skip the singularity
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(true, false), 21));
    runtest(hyperg, ns + "_hyperg",
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
    runtest(gsl::laguerre, "gsl_laguerre",
	    uiorder,
	    fill_argument(std::make_pair(0.0, 100.0),
		  	  std::make_pair(true, true), 101));
    runtest(laguerre, ns + "_laguerre", uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 101));


    ///  Legendre polynomials
    std::cout << "legendre" << std::endl;
    runtest(gsl_sf_legendre_Pl, "gsl_legendre", iorder,
	    fill_argument(std::make_pair(-1.0, 1.0),
	   		  std::make_pair(true, true), 1001));
    runtest(legendre, ns + "_legendre", uiorder,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
			  std::make_pair(true, true), 1001));


    //  Riemann zeta function.
    std::cout << "riemann_zeta" << std::endl;
    //  Skip the pole at 1.
    runtest(gsl::riemann_zeta, "gsl_riemann_zeta_neg",
	    fill_argument(std::make_pair(-10.0, 1.0),
			  std::make_pair(true, false), 56));
    runtest(gsl::riemann_zeta, "gsl_riemann_zeta_pos",
	    fill_argument(std::make_pair(1.0, 30.0),
			  std::make_pair(false, true), 146));
    runtest(riemann_zeta, ns + "_riemann_zeta",
	    fill_argument(std::make_pair(Real{-10}, Real{30}),
	        	  std::make_pair(true, true), 201));

#if STD
    //  Hurwitz zeta function.
    std::cout << "hurwitz_zeta" << std::endl;
    //  Skip the pole at 1.
    runtest(gsl::hurwitz_zeta, "gsl_hurwitz_zeta",
	    fill_argument(std::make_pair(1.0, 30.0),
			  std::make_pair(false, true), 146),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 26));
    runtest(hurwitz_zeta, ns + "_hurwitz_zeta",
	    fill_argument(std::make_pair(Real{1}, Real{30}),
		 	  std::make_pair(true, true), 146),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(true, true), 26));
#endif // STD


    //  Spherical Bessel functions.
    std::cout << "sph_bessel" << std::endl;
    runtest(gsl_sf_bessel_jl, "gsl_sph_bessel", isborder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(true, true), 1001));
    runtest(sph_bessel, ns + "_sph_bessel", uisborder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 1001));

    //  Spherical Legendre functions.
    std::cout << "sph_legendre" << std::endl;
    runtest(gsl::legendre_sphPlm, "gsl_sph_legendre", uiorder, uiorder,
	    fill_argument(std::make_pair(0.0, static_cast<double>(M_PI)),
	    		  std::make_pair(true, true), 1001));
    runtest(sph_legendre, ns + "_sph_legendre", uiorder, uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		  	  std::make_pair(true, true), 1001));



    //  Spherical Neumann functions.
    std::cout << "sph_neumann" << std::endl;
    runtest(gsl_sf_bessel_yl, "gsl_sph_neumann", isborder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(false, true),  // Skip the pole at the origin.
	   		  1001));
    runtest(sph_neumann, ns + "_sph_neumann", uisborder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 1001));

#if STD
    //  Carlson elliptic functions R_C.
    std::cout << "ellint_rc" << std::endl;
    runtest(gsl::ellint_rc, "gsl_ellint_rc",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    runtest(ellint_rc, ns + "_ellint_rc",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_D.
    std::cout << "ellint_rd" << std::endl;
    runtest(gsl::ellint_rd, "gsl_ellint_rd",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    runtest(ellint_rd, ns + "_ellint_rd",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_F.
    std::cout << "ellint_rf" << std::endl;
    runtest(gsl::ellint_rf, "gsl_ellint_rf",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    runtest(ellint_rf, ns + "_ellint_rf",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_J.
    std::cout << "ellint_rj" << std::endl;
    runtest(gsl::ellint_rj, "gsl_ellint_rj",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    runtest(ellint_rj, ns + "_ellint_rj",
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
    runtest(gsl::dilog, "gsl_dilog",
	    fill_argument(std::make_pair(-10.0, 1.0),
			  std::make_pair(true, true), 23));

    runtest<Real, Real>(dilog, ns + "_dilog",
		      fill_argument(std::make_pair(Real{-10}, Real{1}),
				    std::make_pair(true, true), 23));

    //  Upper incomplete Gamma functions.
    std::cout << "tgamma" << std::endl;
    runtest(gsl::tgamma, "gsl_tgamma",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11));

    runtest(tgamma, ns + "_tgamma",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
		 	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
		 	  std::make_pair(true, true), 11));

    //  Incomplete Beta functions.
    std::cout << "ibeta" << std::endl;
    runtest(gsl::ibeta, "gsl_ibeta",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(5.0, 0.0),
			  std::make_pair(true, false), 11),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(false, false), 21));
    runtest(ibeta, ns + "_ibeta",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{5}, Real{0}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
	    		  std::make_pair(false, false), 21));

    //  Digamma or psi functions.
    std::cout << "psi" << std::endl;
    runtest(gsl::psi, "gsl_psi",
	    fill_argument(std::make_pair(-9.875, 10.125),
			  std::make_pair(true, true), 41));
    runtest(psi, ns + "_psi",
	    fill_argument(std::make_pair(Real{-9.9375}, Real{10.0625}),
	        	  std::make_pair(true, true), 801));

    //  Sine integral or Si functions.
    std::cout << "gsl_sinint" << std::endl;
    runtest(gsl_sf_Si, "sinint",
	    fill_argument(std::make_pair(0.0, +10.0),
			  std::make_pair(false, true), 101));
    runtest(sinint, ns + "_sinint",
	    fill_argument(std::make_pair(Real{0}, Real{+10}),
	        	  std::make_pair(false, true), 101));

    //  Cosine integral or Ci functions.
    std::cout << "gsl_cosint" << std::endl;
    runtest(gsl_sf_Ci, "cosint",
	    fill_argument(std::make_pair(0.0, +10.0),
			  std::make_pair(false, true), 101));
    runtest(cosint, ns + "_cosint",
	    fill_argument(std::make_pair(Real{0}, Real{+10}),
	        	  std::make_pair(false, true), 101));

    //  Hyperbolic sine integral or Shi functions.
    std::cout << "gsl_sinhint" << std::endl;
    runtest(gsl_sf_Shi, "sinhint",
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    runtest(sinhint, ns + "_sinhint",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Hyperbolic cosine integral or Chi functions.
    std::cout << "gsl_coshint" << std::endl;
    runtest(gsl_sf_Chi, "coshint",
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    runtest(coshint, ns + "_coshint",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Dawson integral.
    std::cout << "gsl_dawson" << std::endl;
    runtest(gsl::dawson, "dawson",
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    runtest(dawson, ns + "_dawson",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Exponential integral E_n.
    std::cout << "gsl_expint_en" << std::endl;
    runtest(gsl::expint, "expint_en",
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    runtest(expint, ns + "_expint_en",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Normalized sinus cardinal function.
    std::cout << "gsl_sinc" << std::endl;
    runtest(gsl::sinc, "sinc",
	    fill_argument(std::make_pair(-20.0, +20.0),
			  std::make_pair(false, true), 401));
    runtest(sinc, ns + "_sinc",
	    fill_argument(std::make_pair(Real{20}, Real{+20}),
	        	  std::make_pair(false, true), 401));

    //  Sinus cardinal function.
    std::cout << "gsl_sinc_pi" << std::endl;
    runtest(gsl::sinc_pi, "sinc_pi",
	    fill_argument(std::make_pair(-20.0, +20.0),
			  std::make_pair(false, true), 401));
    runtest(sinc_pi, ns + "_sinc_pi",
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
  //do_test<long double>();
  do_test<double>();
  //do_test<float>();
}
