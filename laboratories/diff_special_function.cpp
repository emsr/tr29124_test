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
#include <emsr/special_functions.h>
#include <wrap_gsl.h>
#include <wrap_boost.h>
#include <wrap_burkhardt.h>

#include <test_func.tcc>

///
///
///
int
main()
{
    using _TpGSL = double;
    using Real = double;

    const auto s_pi = emsr::pi_v<Real>;

    //  Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
    std::vector<unsigned int> vorder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

    // Integer orders for various polynomials, harmonics, and spherical bessels.
    std::vector<int> iorder{0, 1, 2, 5, 10, 20, 50, 100};

    //  ... corresponding "_TpGSL" integer orders for GSL.
    std::vector<_TpGSL> dvorder{std::begin(vorder), std::end(vorder)};

    //  Orders for spherical bessel functions.
    std::vector<unsigned int> sborder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

    // Orders for cylindrical Bessel functions.
    std::vector<Real> cyl_neg_order{-5, -2, -1, -Real{2.0Q/3.0Q},
				    -Real{0.5Q}, -Real{1.0Q/3.0Q}};

    //  Orders for cylindrical Bessel functions.
    std::vector<_TpGSL> cyl_order{0, _TpGSL{1}/_TpGSL{3},
				 _TpGSL{0.5Q}, _TpGSL{2}/_TpGSL{3},
				 1, 2, 3, 5, 10, 20, 50, 100};

    // Orders for spherical bessel functions.
    std::vector<unsigned int> sph_order{0, 1, 2, 5, 10, 20, 50, 100};

    const unsigned int num_phi = 19; // 0 - 180 degrees.
    _TpGSL phi[num_phi];
    for (unsigned int i = 0; i < num_phi; ++i)
      phi[i] = _TpGSL{10} * i * s_pi / _TpGSL{180};
    std::vector<_TpGSL> vphid(phi, phi + num_phi);

    std::vector<_TpGSL> vab{0, 0.5, 1, 2, 5, 10, 20};

    std::string basename;

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
    using emsr::cos_pi;
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
    using emsr::falling_factorial;
    using emsr::fresnel_c;
    using emsr::fresnel_s;
    using emsr::gegenbauer;
    using emsr::hermite;
    using emsr::heuman_lambda;
    using emsr::hurwitz_zeta;
    using emsr::hyperg;
    using emsr::ibeta;
    using emsr::ibetac;
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
    using emsr::gamma_q;
    using emsr::radpoly;
    using emsr::riemann_zeta;
    using emsr::rising_factorial;
    using emsr::sinhc;
    using emsr::sinhc_pi;
    using emsr::sinc;
    using emsr::sinc_pi;
    using emsr::sinhc;
    using emsr::sinhc_pi;
    using emsr::sinhint;
    using emsr::sinint;
    using emsr::sin_pi;
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

    //  Airy Ai functions.
    std::cout << "airy_ai" << '\n';
    basename = "diff_airy_ai";
    rundiff(airy_ai, gsl::airy_ai, basename,
	    "x", fill_argument(std::make_pair(Real{-10}, Real{+10}),
			       std::make_pair(true, true), 41));

    //  Airy Bi functions.
    std::cout << "airy_bi" << '\n';
    basename = "diff_airy_bi";
    rundiff(airy_bi, gsl::airy_bi, basename,
	    "x", fill_argument(std::make_pair(Real{-10}, Real{+10}),
			       std::make_pair(true, true), 41));

    //  Associated Laguerre polynomials.
    std::cout << "assoc_laguerre" << '\n';
    basename = "diff_assoc_laguerre";
    rundiff(assoc_laguerre, gsl::assoc_laguerre, basename,
	    "n", vorder, "m", vorder,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
	    		       std::make_pair(true, true), 101));


    //  Associated Legendre functions.
    std::cout << "assoc_legendre" << '\n';
    basename = "diff_assoc_legendre";
    rundiff(assoc_legendre, gsl::assoc_legendre, basename,
	    "l", vorder, "m", vorder,
	    "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
	    		       std::make_pair(true, true), 1001));


    //  Beta function.
    std::cout << "beta" << '\n';
    basename = "diff_beta";
    rundiff(beta, gsl::beta, basename,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
			       std::make_pair(false, true), 101),
	    "y", fill_argument(std::make_pair(Real{0}, Real{100}),
			       std::make_pair(false, true), 101));


    // Binomial coefficient.
    std::cout << "binomial" << '\n';
    basename = "diff_binomial";
    rundiff(binomial<Real>, gsl::binomial, basename,
	    "n", fill_argument(std::make_pair(0U, 50U),
	    		       std::make_pair(true, true), 51),
	    "k", fill_argument(std::make_pair(0U, 50U),
	    		       std::make_pair(true, true), 51));

    // Log binomial coefficient.
    std::cout << "lbinomial" << '\n';
    basename = "diff_lbinomial";
    rundiff(lbinomial<Real>, gsl::lbinomial, basename,
	    "n", fill_argument(std::make_pair(0U, 200U),
			       std::make_pair(true, true), 201),
	    "k", fill_argument(std::make_pair(0U, 200U),
			       std::make_pair(true, true), 201));


    //  Complete elliptic integrals of the first kind.
    //  Avoid poles at |x| = 1.
    std::cout << "comp_ellint_1" << '\n';
    basename = "diff_comp_ellint_1";
    rundiff(comp_ellint_1, gsl::comp_ellint_1, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 101));


    //  Complete elliptic integrals of the second kind.
    //  Avoid poles at |x| = 1.
    std::cout << "comp_ellint_2" << '\n';
    basename = "diff_comp_ellint_2";
    rundiff(comp_ellint_2, gsl::comp_ellint_2, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 101));


    //  Complete elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "comp_ellint_3" << '\n';
    basename = "diff_comp_ellint_3";
    rundiff(comp_ellint_3, gsl::comp_ellint_3, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 101),
	    "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				std::make_pair(true, false), 11));


    //  Confluent hypergeometric functions.
    //  Skip the singularity at c = 0.
    std::cout << "conf_hyperg" << '\n';
    basename = "diff_conf_hyperg";
    rundiff(conf_hyperg, gsl::conf_hyperg, basename,
	    "a", vab,
	    "c", fill_argument(std::make_pair(Real{0}, Real{10}),
			       std::make_pair(false, true), 11),
	    "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
			       std::make_pair(true, true), 201));


    //  Regular modified cylindrical Bessel functions.
    std::cout << "cyl_bessel_i" << '\n';
    basename = "diff_cyl_bessel_i";
    rundiff(cyl_bessel_i, gsl::cyl_bessel_i, basename,
	    "nu", cyl_order,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
			       std::make_pair(true, true), 1001));


    //  Cylindrical Bessel functions (of the first kind).
    std::cout << "cyl_bessel_j" << '\n';
    basename = "diff_cyl_bessel_j";
    rundiff(cyl_bessel_j, gsl::cyl_bessel_j, basename,
	    "nu", cyl_order,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
			       std::make_pair(true, true), 1001));


    //  Irregular modified cylindrical Bessel functions.
    // Skip the pole at the origin.
    std::cout << "cyl_bessel_k" << '\n';
    basename = "diff_cyl_bessel_k";
    rundiff(cyl_bessel_k, gsl::cyl_bessel_k, basename,
	    "nu", cyl_order,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
			       std::make_pair(false, true), 1001));


    // Cylindrical Hankel functions of the first kind.
    std::cout << "cyl_hankel_1" << '\n';
    basename = "cyl_hankel_1";
    rundiff(cyl_hankel_1, beast::cyl_hankel_1, basename,
	    "nu", cyl_order,
	    "x", fill_argument(fill_argument(std::make_pair(Real{0}, Real{5}),
					     std::make_pair(false, true), 21),
			       fill_argument(std::make_pair(Real{0}, Real{100}),
					     std::make_pair(false, true), 21)));

    // Cylindrical Hankel functions of the second kind.
    std::cout << "cyl_hankel_2" << '\n';
    basename = "cyl_hankel_2";
    rundiff(cyl_hankel_2, beast::cyl_hankel_2, basename,
	    "nu", cyl_order,
	    "x", fill_argument(fill_argument(std::make_pair(Real{0}, Real{5}),
					     std::make_pair(false, true), 21),
			       fill_argument(std::make_pair(Real{0}, Real{100}),
					     std::make_pair(false, true), 21)));


    //  Cylindrical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "cyl_neumann" << '\n';
    basename = "diff_cyl_neumann";
    rundiff(cyl_neumann, gsl::cyl_neumann, basename,
	    "nu", cyl_order,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
			       std::make_pair(false, true), 1001), 101);


    //  Elliptic integrals of the first kind.
    //  Avoid poles at |x| = 1.
    std::cout << "ellint_1" << '\n';
    basename = "diff_ellint_1";
    rundiff(ellint_1, gsl::ellint_1, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 101),
	    "phi", vphid);


    //  Elliptic integrals of the second kind.
    //  Avoid poles at |x| = 1.
    std::cout << "ellint_2" << '\n';
    basename = "diff_ellint_2";
    rundiff(ellint_2, gsl::ellint_2, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 101),
	    "phi", vphid);


    //  Elliptic integrals of the third kind.
    //  Avoid poles at |x| = 1 and at nu = 1.
    std::cout << "ellint_3" << '\n';
    basename = "diff_ellint_3";
    rundiff(ellint_3, gsl::ellint_3, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 101),
	    "nu", fill_argument(std::make_pair(Real{0}, Real{1}),
				std::make_pair(true, false), 11),
	    "phi", vphid);


    //  Exponential integral.
    //  Skip the pole at 0.
    std::cout << "expint" << '\n';
    basename = "diff_expint_neg";
    rundiff(expint, gsl::expint, basename,
	    "x", fill_argument(std::make_pair(Real{-50}, Real{0}),
			       std::make_pair(true, false), 51));
    basename = "diff_expint_pos";
    rundiff(expint, gsl::expint, basename,
	    "x", fill_argument(std::make_pair(Real{0}, Real{50}),
			       std::make_pair(false, true), 51));

    // Dawson integral.
    std::cout << "dawson" << '\n';
    basename = "diff_dawson";
    rundiff(dawson, gsl::dawson, basename,
	    "x", fill_argument(std::make_pair(Real{0}, Real{+20}),
			       std::make_pair(false, true), 201));


    //  Hermite polynomials
    std::cout << "hermite" << '\n';
    basename = "diff_hermite";
    rundiff(hermite, gsl::hermite, basename,
	    "n", vorder,
	    "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
	  		       std::make_pair(true, true), 101));


    //  Hypergeometric functions.
    //  Skip the singularity at c = 0.
    //  Skip the singularities at x = -1.
    std::cout << "hyperg" << '\n';
    basename = "diff_hyperg";
    rundiff(hyperg, gsl::hyperg, basename,
	    "a", vab, "b", vab,
	    "c", fill_argument(std::make_pair(Real{0}, Real{10}),
			       std::make_pair(false, true), 11),
	    "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(true, false), 21));


    //  Laguerre polynomials.
    std::cout << "laguerre" << '\n';
    basename = "diff_laguerre";
    rundiff(laguerre, gsl::laguerre, basename,
	    "n", vorder,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
	  		       std::make_pair(true, true), 1001));


    //  Legendre polynomials.
    std::cout << "legendre" << '\n';
    basename = "diff_legendre";
    rundiff(legendre, gsl::legendre_p, basename,
	    "l", vorder,
	    "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
	  		       std::make_pair(true, true), 1001));


    //  Riemann zeta function.
    //  Skip the pole at 1.
    std::cout << "riemann_zeta" << '\n';
    basename = "diff_riemann_zeta_neg";
    rundiff(riemann_zeta, gsl::riemann_zeta, basename,
	    "x", fill_argument(std::make_pair(Real{-10}, Real{1}),
			       std::make_pair(true, false), 56));
    basename = "diff_riemann_zeta_pos";
    rundiff(riemann_zeta, gsl::riemann_zeta, basename,
	    "x", fill_argument(std::make_pair(Real{1}, Real{30}),
			       std::make_pair(false, true), 146));


    //  Hurwitz zeta function.
    std::cout << "hurwitz_zeta" << '\n';
    //  Skip the pole at 1.
    basename = "diff_hurwitz_zeta";
    rundiff(hurwitz_zeta, gsl::hurwitz_zeta, basename,
	    "s", fill_argument(std::make_pair(Real{1}, Real{30}),
			       std::make_pair(false, true), 146),
	    "a", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 26));


    //  Spherical Bessel functions.
    std::cout << "sph_bessel" << '\n';
    basename = "diff_sph_bessel";
    rundiff(sph_bessel, gsl::sph_bessel, basename,
	    "n", sborder,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
	  		       std::make_pair(true, true), 1001));


    //  Spherical Legendre functions.
    std::cout << "sph_legendre" << '\n';
    basename = "diff_sph_legendre";
    rundiff(sph_legendre, gsl::sph_legendre, basename,
	    "l", vorder, "m", vorder,
	    "theta", fill_argument(std::make_pair(Real{0}, s_pi),
	    			   std::make_pair(true, true), 1001));


    //  Spherical Neumann functions.
    // Skip the pole at the origin.
    std::cout << "sph_neumann" << '\n';
    basename = "diff_sph_neumann";
    rundiff(sph_neumann, gsl::sph_neumann, basename,
	    "n", sborder,
	    "x", fill_argument(std::make_pair(Real{0}, Real{100}),
	  		       std::make_pair(false, true), 1001));


    //  Carlson elliptic functions R_C.
    std::cout << "ellint_rc" << '\n';
    basename = "diff_ellint_rc";
    rundiff(ellint_rc, gsl::ellint_rc, basename,
	    "x", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 11),
	    "y", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_D.
    std::cout << "ellint_rd" << '\n';
    basename = "diff_ellint_rd";
    rundiff(ellint_rd, gsl::ellint_rd, basename,
	    "x", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 11),
	    "y", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(true, true), 11),
	    "z", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 11));

    //  Carlson elliptic functions R_F.
    std::cout << "ellint_rf" << '\n';
    basename = "diff_ellint_rf";
    rundiff(ellint_rf, gsl::ellint_rf, basename,
	    "x", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 11),
	    "y", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(true, true), 11),
	    "z", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 11));

      // Carlson elliptic functions R_G.
      std::cout << "ellint_rg" << '\n';
      basename = "diff_ellint_rg";
      rundiff(ellint_rg, beast::ellint_rg, basename,
	      "x", fill_argument(std::make_pair(Real{0}, +Real{5}),
				 std::make_pair(false, true), 11),
	      "y", fill_argument(std::make_pair(Real{0}, +Real{5}),
				 std::make_pair(true, true), 11),
	      "z", fill_argument(std::make_pair(Real{0}, +Real{5}),
				 std::make_pair(true, true), 11));

    //  Carlson elliptic functions R_J.
    std::cout << "ellint_rj" << '\n';
    basename = "diff_ellint_rj";
    rundiff(ellint_rj, gsl::ellint_rj, basename,
	    "x", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 11),
	    "y", fill_argument(std::make_pair(Real{0}, Real{5}),
		 	       std::make_pair(true, true), 11),
	    "z", fill_argument(std::make_pair(Real{0}, Real{5}),
		 	       std::make_pair(false, true), 11),
	    "p", fill_argument(std::make_pair(Real{0}, Real{5}),
		 	       std::make_pair(false, true), 11));

    //  Dilogarithm functions.
    std::cout << "dilog" << '\n';
    basename = "diff_dilog";
    rundiff(dilog, gsl::dilog, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{-10}, _TpGSL{1}),
			       std::make_pair(true, true), 23));

    //  Upper incomplete Gamma functions.
    std::cout << "tgamma" << '\n';
    basename = "diff_tgamma";
    rundiff(tgamma, gsl::tgamma, basename,
	    "a", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{5}),
			       std::make_pair(false, true), 11),
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{5}),
			       std::make_pair(true, true), 11));

      // Lower incomplete Gamma functions.
      std::cout << "tgamma_lower" << '\n';
      basename = "diff_tgamma_lower";
      rundiff(tgamma_lower, gsl::tgamma_lower, basename,
	      "a", fill_argument(std::make_pair(Real{0}, Real{5}),
				 std::make_pair(false, true), 11),
	      "x", fill_argument(std::make_pair(Real{0}, Real{5}),
				 std::make_pair(true, true), 11));

    //  Incomplete Beta functions.
    std::cout << "ibeta" << '\n';
    basename = "diff_ibeta";
    rundiff(ibeta, gsl::ibeta, basename,
	    "a", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{5}),
			       std::make_pair(false, true), 11),
	    "b", fill_argument(std::make_pair(_TpGSL{5}, _TpGSL{0}),
			       std::make_pair(true, false), 11),
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{1}),
			       std::make_pair(false, false), 21));

    //  Complementary incomplete Beta functions.
    std::cout << "ibetac" << '\n';
    basename = "diff_ibetac";
    rundiff(ibetac, beast::ibetac, basename,
	    "a", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{5}),
			       std::make_pair(false, true), 11),
	    "b", fill_argument(std::make_pair(_TpGSL{5}, _TpGSL{0}),
			       std::make_pair(true, false), 11),
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{1}),
			       std::make_pair(false, false), 21));

    //  Digamma or psi functions.
    std::cout << "digamma" << '\n';
    basename = "diff_digamma";
    rundiff(digamma, gsl::digamma, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{-9.9375Q}, _TpGSL{10.0625Q}),
			       std::make_pair(true, true), 801));

    //  Sine integral or Si functions.
    std::cout << "sinint" << '\n';
    basename = "diff_sinint";
    rundiff(sinint, gsl::sinint, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+10}),
			       std::make_pair(false, true), 101));

    //  Cosine integral or Ci functions.
    std::cout << "cosint" << '\n';
    basename = "diff_cosint";
    rundiff(cosint, gsl::cosint, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+10}),
			       std::make_pair(false, true), 101));

    //  Hyperbolic sine integral or Shi functions.
    std::cout << "sinhint" << '\n';
    basename = "diff_sinhint";
    rundiff(sinhint, gsl::sinhint, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+5}),
			       std::make_pair(false, true), 101));

    //  Hyperbolic cosine integral or Chi functions.
    std::cout << "coshint" << '\n';
    basename = "diff_coshint";
    rundiff(coshint, gsl::coshint, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+5}),
			       std::make_pair(false, true), 101));

    // Dawson integral.
    std::cout << "dawson" << '\n';
    basename = "diff_dawson";
    rundiff(dawson, gsl::dawson, basename,
	     "x", fill_argument(std::make_pair(Real{0}, Real{+20}),
				std::make_pair(false, true), 201));

    // Jacobian elliptic integrals.
    std::cout << "jacobi_sn" << '\n';
    basename = "diff_jacobi_sn";
    rundiff(jacobi_sn, gsl::jacobi_sn, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101));

    // Jacobian elliptic integrals.
    std::cout << "jacobi_cn" << '\n';
    basename = "diff_jacobi_cn";
    rundiff(jacobi_cn, gsl::jacobi_cn, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101));

    // Jacobian elliptic integrals.
    std::cout << "jacobi_dn" << '\n';
    basename = "diff_jacobi_dn";
    rundiff(jacobi_dn, gsl::jacobi_dn, basename,
	     "k", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21),
	     "u", fill_argument(std::make_pair(Real{-5}, Real{+5}),
				std::make_pair(true, true), 101));

    //  Exponential integral E_n.
    std::cout << "expint" << '\n';
    basename = "diff_expint_en";
    rundiff(expint, gsl::expint, basename,
	    "n", {0, 1, 2, 3, 5},
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+5}),
			       std::make_pair(false, true), 101));


    //  Fresnel cosine integral.
    std::cout << "fresnel_c" << '\n';
    basename = "diff_fresnel_c";
    rundiff(fresnel_c, gsl::fresnel_c, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{-20}, _TpGSL{+20}),
			       std::make_pair(false, true), 401));

    //  Fresnel sine integral.
    std::cout << "fresnel_s" << '\n';
    basename = "diff_fresnel_s";
    rundiff(fresnel_s, gsl::fresnel_s, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{-20}, _TpGSL{+20}),
			       std::make_pair(false, true), 401));

    //  Dawson integral.
    std::cout << "dawson" << '\n';
    basename = "diff_dawson";
    rundiff(dawson, gsl::dawson, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+5}),
			       std::make_pair(false, true), 101));

    // Normalized sine cardinal function.
    std::cout << "sinc" << '\n';
    basename = "diff_sinc";
    rundiff(sinc, gsl::sinc, basename,
	    "x", fill_argument(std::make_pair(Real{-20}, Real{+20}),
			       std::make_pair(true, true), 401));

    //  Sine cardinal function.
    std::cout << "sinc_pi" << '\n';
    basename = "diff_sinc_pi";
    rundiff(sinc_pi, gsl::sinc_pi, basename,
	    "x", fill_argument(std::make_pair(_TpGSL{-20}, _TpGSL{+20}),
			       std::make_pair(false, true), 401));

    // Log rising factorial symbol.
    std::cout << "lrising_factorial" << '\n';
    basename = "diff_lrising_factorial";
    rundiff(lrising_factorial, beast::lrising_factorial, basename,
	    "a", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 21),
	    "x", {Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}});

    // Log falling factorials.
    std::cout << "lfalling_factorial" << '\n';
    basename = "diff_lfalling_factorial";
    rundiff(lfalling_factorial, beast::lfalling_factorial, basename,
	    "a", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 21),
	    "x", {Real{0}, Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}});

    // Rising factorials.
    std::cout << "rising_factorial" << '\n';
    basename = "diff_rising_factorial";
    rundiff(rising_factorial, beast::rising_factorial, basename,
	    "a", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 21),
	    "x", dvorder);

    // Falling factorials.
    std::cout << "falling_factorial" << '\n';
    basename = "diff_falling_factorial";
    rundiff(falling_factorial, beast::falling_factorial, basename,
	    "a", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 21),
	    "x", {Real{0}, Real{1}, Real{2}, Real{5}, Real{10}, Real{20}, Real{50}, Real{100}});

    // Regular modified spherical bessel functions.
    std::cout << "sph_bessel_i" << '\n';
    basename = "diff_sph_bessel_i";
    rundiff(sph_bessel_i, gsl::sph_bessel_i, basename,
	    "n", sph_order,
	    "x", fill_argument(
		       fill_argument(std::make_pair(Real{0}, Real{5}),
				     std::make_pair(true, true), 21),
		       fill_argument(std::make_pair(Real{0}, Real{100}),
				     std::make_pair(true, true), 21)));

    // Irregular modified spherical bessel functions.
    std::cout << "sph_bessel_k" << '\n';
    basename = "diff_sph_bessel_k";
    rundiff(sph_bessel_k, gsl::sph_bessel_k, basename,
	    "n", sph_order,
	    "x", fill_argument(
			fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(false, true), 21),
			fill_argument(std::make_pair(Real{0}, Real{100}),
			       std::make_pair(false, true), 21)));

    // Legendre functions of the second kind.
    std::cout << "legendre_q" << '\n';
    basename = "diff_legendre_q";
    rundiff(legendre_q, gsl::legendre_q, basename,
	    "l", vorder,
	    "x", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 21));

    // Factorial.
    std::cout << "factorial" << '\n';
    basename = "diff_factorial";
    rundiff(factorial<Real>, gsl::factorial, basename,
	    "n", fill_argument(std::make_pair(0U, 50U),
			       std::make_pair(true, true), 51));

    // Log factorial.
    std::cout << "lfactorial" << '\n';
    basename = "diff_lfactorial";
    rundiff(lfactorial<Real>, gsl::lfactorial, basename,
	    "n", fill_argument(std::make_pair(0U, 500U),
			       std::make_pair(true, true), 501));

    // Double factorial.
    std::cout << "double_factorial" << '\n';
    basename = "diff_double_factorial";
    rundiff(double_factorial<Real>, gsl::double_factorial, basename,
	    "n", fill_argument(std::make_pair(0, 50),
			       std::make_pair(true, true), 51));

    // Log double factorial.
    std::cout << "ldouble_factorial" << '\n';
    basename = "diff_ldouble_factorial";
    rundiff(ldouble_factorial<Real>, gsl::ldouble_factorial, basename,
	    "n", fill_argument(std::make_pair(0, 500),
			       std::make_pair(true, true), 501));

    // Binomial coefficient.
    std::cout << "binomial" << '\n';
    basename = "diff_binomial";
    rundiff(binomial<Real>, gsl::binomial, basename,
	    "n", fill_argument(std::make_pair(0U, 50U),
			       std::make_pair(true, true), 51),
	    "k", fill_argument(std::make_pair(0U, 50U),
			       std::make_pair(true, true), 51));

    // Log binomial coefficient.
    std::cout << "lbinomial" << '\n';
    basename = "diff_lbinomial";
    rundiff(lbinomial<Real>, gsl::lbinomial, basename,
	    "n", fill_argument(std::make_pair(0U, 200U),
			       std::make_pair(true, true), 201),
	    "k", fill_argument(std::make_pair(0U, 200U),
			       std::make_pair(true, true), 201));

    // Gegenbauer polynomials.
    std::cout << "gegenbauer" << '\n';
    basename = "diff_gegenbauer";
    rundiff(gegenbauer, gsl::gegenbauer, basename,
	    "n", vorder,
	    "alpha", fill_argument(std::make_pair(Real{0}, Real{5}),
				   std::make_pair(true, true), 11),
            "x", fill_argument(std::make_pair(Real{0}, Real{20}),
			       std::make_pair(true, true), 41));

    // Chebyshev polynomials of the first kind.
    std::cout << "chebyshev_t - UNTESTED" << '\n';
/*
    basename = "diff_chebyshev_t";
    rundiff(chebyshev_t, gsl::chebyshev_t, basename,
	     "n", {0U, 1U, 5U, 8U, 10U, 20U, 40U, 100U},
	     "x", fill_argument(std::make_pair(Real{-1}, Real{+1}),
				std::make_pair(true, true), 21));
*/
    // Chebyshev polynomials of the second kind.
    std::cout << "chebyshev_u - UNTESTED" << '\n';
//    basename = "diff_chebyshev_u";

    // Chebyshev polynomials of the third kind.
    std::cout << "chebyshev_v - UNTESTED" << '\n';
//    basename = "diff_chebyshev_v";

    // Chebyshev polynomials of the fourth kind.
    std::cout << "chebyshev_w - UNTESTED" << '\n';
//    basename = "diff_chebyshev_w";

    // Jacobi polynomials.
    std::cout << "jacobi" << '\n';
    basename = "diff_jacobi";
    rundiff(jacobi, gsl::jacobi, basename,
	    "n", vorder,
	    "alpha", fill_argument(std::make_pair(Real{0}, Real{5}),
				   std::make_pair(true, true), 11),
            "beta", fill_argument(std::make_pair(Real{0}, Real{5}),
				  std::make_pair(true, true), 21),
            "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
				  std::make_pair(true, true), 41));

    // Radial polynomials.
    std::cout << "radpoly" << '\n';
    basename = "diff_radpoly";
    rundiff(radpoly, gsl::radpoly, basename,
	    "n", vorder, "m", vorder,
            "rho", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, true), 21));

    // Zernike polynomials.
    std::cout << "zernike" << '\n';
    basename = "diff_zernike";
    rundiff(zernike, gsl::zernike, basename,
	    "n", vorder, "m", iorder,
            "rho", fill_argument(std::make_pair(Real{0}, Real{1}),
				 std::make_pair(true, true), 21),
            "phi", vphid);

    // Confluent hypergeometric limit functions.
    // Skip the singularity at c = 0.
    std::cout << "conf_hyperg_lim" << '\n';
    basename = "diff_conf_hyperg_lim";
    rundiff(conf_hyperg_lim, gsl::conf_hyperg_lim, basename,
	    "c", fill_argument(std::make_pair(Real{0}, Real{10}),
			       std::make_pair(false, true), 11),
	    "x", fill_argument(std::make_pair(Real{-10}, Real{10}),
			       std::make_pair(true, true), 21));

    // Heuman lambda functions.
    // Avoid poles at |x| = 1.
    std::cout << "heuman_lambda" << '\n';
    basename = "diff_heuman_lambda";
    rundiff(heuman_lambda, beast::heuman_lambda, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 21),
	    "phi", vphid);

    // Elliptic D integrals.
    // Avoid poles at |x| = 1.
    std::cout << "ellint_d" << '\n';
    basename = "diff_ellint_d";
    rundiff(ellint_d, beast::ellint_d, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 21),
	    "phi", vphid);

    // Complementary elliptic D integrals.
    // Avoid poles at |x| = 1.
    std::cout << "comp_ellint_d" << '\n';
    basename = "diff_comp_ellint_d";
    rundiff(comp_ellint_d, beast::comp_ellint_d, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 21));

    // Jacobi zeta functions.
    // Avoid poles at |x| = 1.
    std::cout << "jacobi_zeta" << '\n';
    basename = "diff_jacobi_zeta";
    rundiff(jacobi_zeta, beast::jacobi_zeta, basename,
	    "k", fill_argument(std::make_pair(Real{-1}, Real{1}),
			       std::make_pair(false, false), 21),
	    "phi", vphid);

    // Cylindrical Hankel functions of the first kind.
    std::cout << "cyl_hankel_1" << '\n';
    basename = "diff_cyl_hankel_1";
    rundiff(cyl_hankel_1, beast::cyl_hankel_1, basename,
	    "nu", cyl_neg_order,
	    "x", fill_argument(
		       fill_argument(std::make_pair(Real{0}, Real{5}),
				     std::make_pair(true, true), 21),
		       fill_argument(std::make_pair(Real{0}, Real{100}),
				     std::make_pair(true, true), 21)));

    // Cylindrical Hankel functions of the second kind.
    std::cout << "cyl_hankel_2" << '\n';
    basename = "diff_cyl_hankel_2";
    rundiff(cyl_hankel_2, beast::cyl_hankel_2, basename,
	    "nu", cyl_neg_order,
	    "x", fill_argument(
		       fill_argument(std::make_pair(Real{0}, Real{5}),
				     std::make_pair(true, true), 21),
		       fill_argument(std::make_pair(Real{0}, Real{100}),
				     std::make_pair(true, true), 21)));

    // Spherical Hankel functions of the first kind.
    std::cout << "sph_hankel_1" << '\n';
    basename = "diff_sph_hankel_1";
    rundiff(sph_hankel_1, beast::sph_hankel_1, basename,
	     "n", sph_order,
	     "x", fill_argument(
			fill_argument(std::make_pair(Real{0}, Real{5}),
				      std::make_pair(true, true), 21),
			fill_argument(std::make_pair(Real{0}, Real{100}),
				      std::make_pair(true, true), 21)));

    // Spherical Hankel functions of the second kind.
    std::cout << "sph_hankel_2" << '\n';
    basename = "diff_sph_hankel_2";
    rundiff(sph_hankel_2, beast::sph_hankel_2, basename,
	    "n", sph_order,
	    "x", fill_argument(
		       fill_argument(std::make_pair(Real{0}, Real{5}),
				     std::make_pair(true, true), 21),
		       fill_argument(std::make_pair(Real{0}, Real{100}),
				     std::make_pair(true, true), 21)));

    // Spherical harmonic functions.
    std::cout << "sph_harmonic" << '\n';
    basename = "diff_sph_harmonic";
    rundiff(sph_harmonic, beast::sph_harmonic, basename,
	    "l", vorder, "m", iorder,
	    "theta", fill_argument(std::make_pair(Real{0}, s_pi),
				   std::make_pair(true, true), 21),
	    "phi", vphid);

    // Dirichlet eta function.
    std::cout << "dirichlet_eta" << '\n';
    // Skip the pole at 1.
    basename = "diff_dirichlet_eta";
    rundiff(dirichlet_eta, gsl::dirichlet_eta, basename,
	    "s", fill_argument(
		       fill_argument(std::make_pair(Real{-10}, Real{1}),
				     std::make_pair(true, false), 56),
		       fill_argument(std::make_pair(Real{1}, Real{30}),
				     std::make_pair(false, true), 146)));

    // Owens T functions.
    std::cout << "owens_t" << '\n';
    basename = "diff_owens_t";
    rundiff(owens_t, beast::owens_t, basename,
	    "h", fill_argument(std::make_pair(Real{-5}, Real{5}),
			       std::make_pair(true, true), 41),
	    "a", fill_argument(std::make_pair(Real{0}, Real{5}),
			       std::make_pair(true, true), 21));

    // clausen_cl function.
    std::cout << "clausen_cl" << '\n';
    basename = "diff_clausen_cl";
    rundiff(clausen_cl, gsl::clausen_cl, basename,
	    "m", fill_argument(std::make_pair(2U, 2U),
			       std::make_pair(true, true), 1),
	    "w", fill_argument(std::make_pair(Real{-10}, Real{+10}),
			       std::make_pair(true, true), 41));

    // Bernoulli numbers.
    std::cout << "bernoulli" << '\n';
    basename = "diff_bernoulli";
    rundiff(bernoulli<Real>, beast::bernoulli, basename,
	    "n", fill_argument(std::make_pair(0U, 100U),
			       std::make_pair(true, true), 101));

    // Reperiodized sine function.
    std::cout << "sin_pi" << '\n';
    basename = "diff_sin_pi";
    rundiff(sin_pi, beast::sin_pi, basename,
	    "x", fill_argument(std::make_pair(Real{-20}, Real{+50}),
			       std::make_pair(false, true), 701));

    // Reperiodized cosine function.
    std::cout << "cos_pi" << '\n';
    basename = "diff_cos_pi";
    rundiff(cos_pi, beast::cos_pi, basename,
	    "x", fill_argument(std::make_pair(Real{-20}, Real{+50}),
			       std::make_pair(false, true), 701));

    // Chebyshev polynomials of the first kind.
    std::cout << "debye" << '\n';
    basename = "diff_debye";
    rundiff(debye, gsl::debye, basename,
	     "n", {1U, 2U, 3U, 4U, 5U, 6U},
	     "x", fill_argument(std::make_pair(Real{0}, Real{10}),
				std::make_pair(true, true), 41));

  return 0;
}

