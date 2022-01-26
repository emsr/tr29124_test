/**
 *
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>

#include <cmath>
#include <array>

#include <gsl/gsl_sf.h>

#include <emsr/float128_io.h>
#include <test_func.tcc>
#include <wrap_gsl.h>
#include <wrap_burkhardt.h>


///
///
///
template<typename Real>
  void
  do_test(Real proto = Real{})
  {
    const auto s_pi = emsr::pi_v<Real>;

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
      phi[i] = Real{10} * i * s_pi / Real{180};
    std::vector<double> vphid(phi, phi + num_phi);
    std::vector<Real> vphi(phi, phi + num_phi);

    std::string ns("std");
    using emsr::airy_ai;
    using emsr::airy_bi;
    using emsr::assoc_laguerre;
    using emsr::assoc_legendre;
    using emsr::bernoulli;
    using emsr::beta;
    using emsr::binomial;
    using emsr::chebyshev_t;
    using emsr::chebyshev_u;
    using emsr::chebyshev_v;
    using emsr::chebyshev_w;
    using emsr::clausen;
    using emsr::clausen_cl;
    using emsr::clausen_sl;
    using emsr::comp_ellint_d;
    using emsr::comp_ellint_1;
    using emsr::comp_ellint_2;
    using emsr::comp_ellint_3;
    using emsr::conf_hyperg;
    using emsr::conf_hyperg_lim;
    using emsr::coshint;
    using emsr::cosint;
    using emsr::cyl_bessel_i;
    using emsr::cyl_bessel_j;
    using emsr::cyl_bessel_k;
    using emsr::cyl_hankel_1;
    using emsr::cyl_hankel_2;
    using emsr::cyl_neumann;
    using emsr::dawson;
    using emsr::debye;
    using emsr::digamma;
    using emsr::dilog;
    using emsr::dirichlet_beta;
    using emsr::dirichlet_eta;
    using emsr::double_factorial;
    using emsr::ellint_1;
    using emsr::ellint_2;
    using emsr::ellint_3;
    using emsr::ellint_d;
    using emsr::ellint_rc;
    using emsr::ellint_rd;
    using emsr::ellint_rf;
    using emsr::ellint_rg;
    using emsr::ellint_rj;
    using emsr::expint;
    using emsr::expint;
    using emsr::factorial;
    using emsr::fresnel_c;
    using emsr::fresnel_s;
    using emsr::gegenbauer;
    using emsr::hermite;
    using emsr::heuman_lambda;
    using emsr::hurwitz_zeta;
    using emsr::hyperg;
    using emsr::ibeta;
    using emsr::jacobi;
    using emsr::jacobi_sn;
    using emsr::jacobi_cn;
    using emsr::jacobi_dn;
    using emsr::jacobi_zeta;
    using emsr::laguerre;
    using emsr::lbinomial;
    using emsr::ldouble_factorial;
    using emsr::legendre;
    using emsr::legendre_q;
    using emsr::lfactorial;
    using emsr::lfalling_factorial;
    using emsr::lrising_factorial;
    using emsr::owens_t;
    using emsr::gamma_p;
    using emsr::falling_factorial;
    using emsr::rising_factorial;
    using emsr::gamma_q;
    using emsr::radpoly;
    using emsr::riemann_zeta;
    using emsr::sinhc;
    using emsr::sinhc_pi;
    using emsr::sinc;
    using emsr::sinc_pi;
    using emsr::sinhc;
    using emsr::sinhc_pi;
    using emsr::sinhint;
    using emsr::sinint;
    using emsr::sph_bessel;
    using emsr::sph_bessel_i;
    using emsr::sph_bessel_k;
    using emsr::sph_hankel_1;
    using emsr::sph_hankel_2;
    using emsr::sph_harmonic;
    using emsr::sph_legendre;
    using emsr::sph_neumann;
    using emsr::tgamma_lower;
    using emsr::tgamma;
    using emsr::theta_1;
    using emsr::theta_2;
    using emsr::theta_3;
    using emsr::theta_4;
    using emsr::theta_s;
    using emsr::theta_c;
    using emsr::theta_d;
    using emsr::theta_n;
    using emsr::zernike;

    std::string basename;

    //  Airy Ai functions.
    std::cout << "airy_ai\n" << std::flush;
    runtest(gsl::airy_ai, "gsl_airy_ai",
	    fill_argument(std::make_pair(-10.0, +10.0),
			  std::make_pair(true, true), 41));
    runtest<Real>(airy_ai, ns + "_airy_ai",
	    fill_argument(std::make_pair(Real{-10}, Real{+10}),
	        	  std::make_pair(true, true), 41));

    //  Airy Bi functions.
    std::cout << "airy_bi\n" << std::flush;
    runtest(gsl::airy_bi, "gsl_airy_bi",
	    fill_argument(std::make_pair(-10.0, +10.0),
			  std::make_pair(true, true), 41));
    runtest<Real>(airy_bi, ns + "_airy_bi",
	    fill_argument(std::make_pair(Real{-10}, Real{+10}),
	        	  std::make_pair(true, true), 41));

    //  Associated Laguerre polynomials.
    std::cout << "assoc_laguerre\n" << std::flush;
    runtest(gsl::assoc_laguerre, "gsl_assoc_laguerre", uiorder, uiorder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(true, true)));
    runtest<Real>(assoc_laguerre, ns + "_assoc_laguerre", uiorder, uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		  	  std::make_pair(true, true)));


    //  Associated Legendre functions.
    std::cout << "assoc_legendre\n" << std::flush;
    runtest(gsl::assoc_legendre, "gsl_assoc_legendre", uiorder, uiorder,
	    fill_argument(std::make_pair(-1.0, 1.0),
	    		  std::make_pair(true, true), 1001));
    runtest<Real>(assoc_legendre, ns + "_assoc_legendre", uiorder, uiorder,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		  	  std::make_pair(true, true), 1001));


    //  Beta function.
    std::cout << "beta\n" << std::flush;
    runtest(gsl::beta, "gsl_beta",
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true)),
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true)));
    runtest<Real>(beta, ns + "_beta",
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(false, true), 101),
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(false, true), 101));


    //  Complete elliptic integrals of the first kind.
    std::cout << "comp_ellint_1\n" << std::flush;
    runtest(gsl::comp_ellint_1, "gsl_comp_ellint_1",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101));  //  Avoid poles at |x| = 1.
    runtest<Real>(comp_ellint_1, ns + "_comp_ellint_1",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	        	  std::make_pair(true, true), 101));


    //  Complete elliptic integrals of the second kind.
    std::cout << "comp_ellint_2\n" << std::flush;
    runtest(gsl::comp_ellint_2, "gsl_comp_ellint_2",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101));  //  Avoid poles at |x| = 1.
    runtest<Real>(comp_ellint_2, ns + "_comp_ellint_2",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	        	  std::make_pair(true, true), 101));


    //  Complete elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3\n" << std::flush;
    runtest(gsl::comp_ellint_3, "gsl_comp_ellint_3",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false)),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(true, false), 11));
    runtest<Real>(comp_ellint_3, ns + "_comp_ellint_3",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 	  std::make_pair(true, true)),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
		 	  std::make_pair(true, true), 11));


    //  Confluent hypergeometric functions.
    std::cout << "conf_hyperg\n" << std::flush;
    runtest(gsl::conf_hyperg,
	    "gsl_conf_hyperg",
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(false, true), 11),  //  Skip the singularity
	    fill_argument(std::make_pair(-10.0, 10.0),
			  std::make_pair(true, true), 201));
    runtest<Real>(conf_hyperg, ns + "_conf_hyperg",
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{-10}, Real{10}),
	    		  std::make_pair(true, true), 201));


    //  Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i\n" << std::flush;
    runtest(gsl::cyl_bessel_i, "gsl_cyl_bessel_i", dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(true, true), 1001));
    runtest<Real>(cyl_bessel_i, ns + "_cyl_bessel_i", border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j\n" << std::flush;
    runtest(gsl::cyl_bessel_j, "gsl_cyl_bessel_j", dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(true, true), 1001));
    runtest<Real>(cyl_bessel_j, ns + "_cyl_bessel_j", border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Irregular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_k\n" << std::flush;
    runtest(gsl::cyl_bessel_k, "gsl_cyl_bessel_k", dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true),  // Skip the pole at the origin.
			  1001));
    runtest<Real>(cyl_bessel_k, ns + "_cyl_bessel_k", border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Cylindrical Neumann functions.
    std::cout << "cyl_neumann\n" << std::flush;
    runtest(gsl::cyl_neumann, "gsl_cyl_neumann", dborder,
	    fill_argument(std::make_pair(0.0, 100.0),
			  std::make_pair(false, true),  // Skip the pole at the origin.
			  1001));
    runtest<Real>(cyl_neumann, ns + "_cyl_neumann", border,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		 	  std::make_pair(true, true), 1001));


    //  Elliptic integrals of the first kind.
    std::cout << "ellint_1\n" << std::flush;
    runtest(gsl::ellint_1, "gsl_ellint_1",
	    fill_argument(std::make_pair(-1.0, 1.0),
	        	  std::make_pair(false, false), 101),  //  Avoid poles at |x| = 1.
	        	  vphid);
    runtest<Real>(ellint_1, ns + "_ellint_1",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 	  std::make_pair(true, true), 101), vphi);


    //  Elliptic integrals of the second kind.
    std::cout << "ellint_2\n" << std::flush;
    runtest(gsl::ellint_2, "gsl_ellint_2",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false), 101),  //  Avoid poles at |x| = 1.
			  vphid);
    runtest<Real>(ellint_2, ns + "_ellint_2",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
		 		   std::make_pair(true, true), 101), vphi);


    //  Elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3\n" << std::flush;
    runtest(gsl::ellint_3, "gsl_ellint_3",
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(false, false)),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(true, false), 11), vphid);
    runtest<Real>(ellint_3, ns + "_ellint_3",
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	    		  std::make_pair(true, true)),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
	    		  std::make_pair(true, true), 11), vphi);


    //  Exponential integrals.
    std::cout << "expint\n" << std::flush;
    //  Skip the pole at 0.
    runtest(gsl::expint, "gsl_expint_neg",
	    fill_argument(std::make_pair(-50.0, 0.0),
			  std::make_pair(true, false), 51));
    runtest(gsl::expint, "gsl_expint_pos",
	    fill_argument(std::make_pair(0.0, 50.0),
			  std::make_pair(false, true), 51));
    runtest<Real>(expint, ns + "_expint",
	    fill_argument(std::make_pair(Real{-50}, Real{50}),
	        	  std::make_pair(true, true), 101));


    //  Hermite polynomials
    std::cout << "hermite\n" << std::flush;
    runtest(gsl::hermite, "gsl_hermite", uiorder,
    	    fill_argument(std::make_pair(-10.0, 10.0),
    			  std::make_pair(true, true), 101));
    runtest<Real>(hermite, ns + "_hermite", uiorder,
	    fill_argument(std::make_pair(Real{-10}, Real{10}),
			  std::make_pair(true, true), 101));


    //  Hypergeometric functions.
    std::cout << "hyperg\n" << std::flush;
    runtest(gsl::hyperg, "gsl_hyperg",
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 10.0),
			  std::make_pair(false, true), 11),  //  Skip the singularity
	    fill_argument(std::make_pair(-1.0, 1.0),
			  std::make_pair(true, false), 21));
    runtest<Real>(hyperg, ns + "_hyperg",
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	        	  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	        	  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{10}),
	        	  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
	        	  std::make_pair(false, true), 21));


    //  Laguerre polynomials.
    std::cout << "laguerre\n" << std::flush;
    runtest(gsl::laguerre, "gsl_laguerre",
	    uiorder,
	    fill_argument(std::make_pair(0.0, 100.0),
		  	  std::make_pair(true, true), 101));
    runtest<Real>(laguerre, ns + "_laguerre", uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 101));


    ///  Legendre polynomials
    std::cout << "legendre\n" << std::flush;
    runtest(gsl::legendre_p, "gsl_legendre", uiorder,
	    fill_argument(std::make_pair(-1.0, 1.0),
	   		  std::make_pair(true, true), 1001));
    runtest<Real>(legendre, ns + "_legendre", uiorder,
	    fill_argument(std::make_pair(Real{-1}, Real{1}),
			  std::make_pair(true, true), 1001));


    //  Riemann zeta function.
    std::cout << "riemann_zeta\n" << std::flush;
    //  Skip the pole at 1.
    runtest(gsl::riemann_zeta, "gsl_riemann_zeta_neg",
	    fill_argument(std::make_pair(-10.0, 1.0),
			  std::make_pair(true, false), 56));
    runtest(gsl::riemann_zeta, "gsl_riemann_zeta_pos",
	    fill_argument(std::make_pair(1.0, 30.0),
			  std::make_pair(false, true), 146));
    runtest<Real>(riemann_zeta, ns + "_riemann_zeta",
	    fill_argument(std::make_pair(Real{-10}, Real{30}),
	        	  std::make_pair(true, true), 201));

    //  Hurwitz zeta function.
    std::cout << "hurwitz_zeta\n" << std::flush;
    //  Skip the pole at 1.
    runtest(gsl::hurwitz_zeta, "gsl_hurwitz_zeta",
	    fill_argument(std::make_pair(1.0, 30.0),
			  std::make_pair(false, true), 146),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 26));
    runtest<Real>(hurwitz_zeta, ns + "_hurwitz_zeta",
	    fill_argument(std::make_pair(Real{1}, Real{30}),
		 	  std::make_pair(true, true), 146),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(true, true), 26));


    //  Spherical Bessel functions.
    std::cout << "sph_bessel\n" << std::flush;
    runtest(gsl::sph_bessel, "gsl_sph_bessel", uisborder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(true, true), 1001));
    runtest<Real>(sph_bessel, ns + "_sph_bessel", uisborder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 1001));

    //  Spherical Legendre functions.
    std::cout << "sph_legendre\n" << std::flush;
    runtest(gsl::sph_legendre, "gsl_sph_legendre", uiorder, uiorder,
	    fill_argument(std::make_pair(0.0, 100.0),
	    		  std::make_pair(true, true), 1001));
    runtest<Real>(sph_legendre, ns + "_sph_legendre", uiorder, uiorder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
		  	  std::make_pair(true, true), 1001));



    //  Spherical Neumann functions.
    std::cout << "sph_neumann\n" << std::flush;
    runtest(gsl::sph_neumann, "gsl_sph_neumann", uisborder,
	    fill_argument(std::make_pair(0.0, 100.0),
	   		  std::make_pair(false, true),  // Skip the pole at the origin.
	   		  1001));
    runtest<Real>(sph_neumann, ns + "_sph_neumann", uisborder,
	    fill_argument(std::make_pair(Real{0}, Real{100}),
			  std::make_pair(true, true), 1001));

    //  Carlson elliptic functions R_C.
    std::cout << "ellint_rc\n" << std::flush;
    runtest(gsl::ellint_rc, "gsl_ellint_rc",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    runtest<Real>(ellint_rc, ns + "_ellint_rc",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
		 	  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_D.
    std::cout << "ellint_rd\n" << std::flush;
    runtest(gsl::ellint_rd, "gsl_ellint_rd",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    runtest<Real>(ellint_rd, ns + "_ellint_rd",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_F.
    std::cout << "ellint_rf\n" << std::flush;
    runtest(gsl::ellint_rf, "gsl_ellint_rf",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    runtest<Real>(ellint_rf, ns + "_ellint_rf",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_J.
    std::cout << "ellint_rj\n" << std::flush;
    runtest(gsl::ellint_rj, "gsl_ellint_rj",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11));
    runtest<Real>(ellint_rj, ns + "_ellint_rj",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	        	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	        	  std::make_pair(true, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	        	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	        	  std::make_pair(false, true), 11));

    //  Dilogarithm functions.
    std::cout << "dilog\n" << std::flush;
    runtest(gsl::dilog, "gsl_dilog",
	    fill_argument(std::make_pair(-10.0, 1.0),
			  std::make_pair(true, true), 23));

    runtest<Real, Real>(dilog, ns + "_dilog",
		      fill_argument(std::make_pair(Real{-10}, Real{1}),
				    std::make_pair(true, true), 23));

    //  Upper incomplete Gamma functions.
    std::cout << "tgamma\n" << std::flush;
    runtest(gsl::tgamma, "gsl_tgamma",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(true, true), 11));

    runtest<Real>(tgamma, ns + "_tgamma",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
		 	  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
		 	  std::make_pair(true, true), 11));

    //  Incomplete Beta functions.
    std::cout << "ibeta\n" << std::flush;
    runtest(gsl::ibeta, "gsl_ibeta",
	    fill_argument(std::make_pair(0.0, 5.0),
			  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(5.0, 0.0),
			  std::make_pair(true, false), 11),
	    fill_argument(std::make_pair(0.0, 1.0),
			  std::make_pair(false, false), 21));
    runtest<Real>(ibeta, ns + "_ibeta",
	    fill_argument(std::make_pair(Real{0}, Real{5}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{5}, Real{0}),
	    		  std::make_pair(false, true), 11),
	    fill_argument(std::make_pair(Real{0}, Real{1}),
	    		  std::make_pair(false, false), 21));

    //  Digamma or psi functions.
    std::cout << "digamma\n" << std::flush;
    runtest(gsl::digamma, "gsl_digamma",
	    fill_argument(std::make_pair(-9.875, 10.125),
			  std::make_pair(true, true), 41));
    runtest<Real>(digamma, ns + "_digamma",
	    fill_argument(std::make_pair(Real{-9.9375}, Real{10.0625}),
	        	  std::make_pair(true, true), 801));

    //  Sine integral or Si functions.
    std::cout << "gsl_sinint\n" << std::flush;
    runtest(gsl::sinint, "sinint",
	    fill_argument(std::make_pair(0.0, +10.0),
			  std::make_pair(false, true), 101));
    runtest<Real>(sinint, ns + "_sinint",
	    fill_argument(std::make_pair(Real{0}, Real{+10}),
	        	  std::make_pair(false, true), 101));

    //  Cosine integral or Ci functions.
    std::cout << "gsl_cosint\n" << std::flush;
    runtest(gsl::cosint, "cosint",
	    fill_argument(std::make_pair(0.0, +10.0),
			  std::make_pair(false, true), 101));
    runtest<Real>(cosint, ns + "_cosint",
	    fill_argument(std::make_pair(Real{0}, Real{+10}),
	        	  std::make_pair(false, true), 101));

    //  Hyperbolic sine integral or Shi functions.
    std::cout << "gsl_sinhint\n" << std::flush;
    runtest(gsl::sinhint, "sinhint",
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    runtest<Real>(sinhint, ns + "_sinhint",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Hyperbolic cosine integral or Chi functions.
    std::cout << "gsl_coshint\n" << std::flush;
    runtest(gsl::coshint, "coshint",
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    runtest<Real>(coshint, ns + "_coshint",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Dawson integral.
    std::cout << "gsl_dawson\n" << std::flush;
    runtest(gsl::dawson, "dawson",
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    runtest<Real>(dawson, ns + "_dawson",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Exponential integral E_n.
    std::cout << "gsl_expint_en\n" << std::flush;
    runtest(gsl::expint, "expint_en",
	    fill_argument(std::make_pair(0.0, +5.0),
			  std::make_pair(false, true), 101));
    runtest<Real>(expint, ns + "_expint_en",
	    fill_argument(std::make_pair(Real{0}, Real{+5}),
	        	  std::make_pair(false, true), 101));

    //  Normalized sinus cardinal function.
    std::cout << "gsl_sinc\n" << std::flush;
    runtest(gsl::sinc, "sinc",
	    fill_argument(std::make_pair(-20.0, +20.0),
			  std::make_pair(false, true), 401));
    runtest<Real>(sinc, ns + "_sinc",
	    fill_argument(std::make_pair(Real{20}, Real{+20}),
	        	  std::make_pair(false, true), 401));

    //  Sinus cardinal function.
    std::cout << "gsl_sinc_pi\n" << std::flush;
    runtest(gsl::sinc_pi, "sinc_pi",
	    fill_argument(std::make_pair(-20.0, +20.0),
			  std::make_pair(false, true), 401));
    runtest<Real>(sinc_pi, ns + "_sinc_pi",
	    fill_argument(std::make_pair(Real{20}, Real{+20}),
	        	  std::make_pair(false, true), 401));

    return;
  }

int
main()
{
  do_test<__float128>();
  do_test<long double>();
  do_test<double>();
  do_test<float>();
}
