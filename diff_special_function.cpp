
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
#include <bits/specfun.h>
#include "gsl_wrap.h"

#include "test_func.tcc"


///
///
///
int
main()
{
  using _TpGSL = double;

  //  Unsigned integer orders for various polynomials, harmonics, and spherical bessels.
  std::vector<unsigned int> vorder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

  //  ... corresponding "_TpGSL" integer orders for GSL.
  std::vector<_TpGSL> dvorder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

  //  Orders for cylindrical Bessel functions.
  std::vector<_TpGSL> cborderd{0, _TpGSL{1}/_TpGSL{3}, _TpGSL{0.5}, _TpGSL{2}/_TpGSL{3},
			       1, 2, 3, 5, 10, 20, 50, 100};

  //  Orders for spherical bessel functions.
  std::vector<unsigned int> sborder{0, 1, 2, 3, 4, 5, 10, 20, 50, 100};

  const unsigned int num_phi = 19; // 0 - 180 degrees.
  _TpGSL phi[num_phi];
  for (unsigned int i = 0; i < num_phi; ++i)
    phi[i] = _TpGSL{10} * i * M_PI / _TpGSL{180};
  std::vector<_TpGSL> vphid(phi, phi + num_phi);

  std::vector<_TpGSL> vab{0, 0.5, 1, 2, 5, 10, 20};

  std::string basename;

#if STD
  using __gnu_cxx::airy_ai;
  using __gnu_cxx::airy_bi;
  using       std::assoc_laguerre;
  using       std::assoc_legendre;
  using       std::beta;
  using       std::comp_ellint_1;
  using       std::comp_ellint_2;
  using       std::comp_ellint_3;
  using __gnu_cxx::conf_hyperg;
  using __gnu_cxx::coshint;
  using __gnu_cxx::cosint;
  using       std::cyl_bessel_i; 
  using       std::cyl_bessel_j; 
  using       std::cyl_bessel_k; 
  using       std::cyl_neumann;  
  using __gnu_cxx::dawson;
  using __gnu_cxx::dilog;
  using __gnu_cxx::double_factorial;
  using       std::ellint_1;	 
  using       std::ellint_2;	 
  using       std::ellint_3;	 
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
  using       std::hermite;	 
  using __gnu_cxx::hurwitz_zeta;
  using __gnu_cxx::hyperg;
  using __gnu_cxx::ibeta;
  using __gnu_cxx::jacobi_sn;
  using __gnu_cxx::jacobi_cn;
  using __gnu_cxx::jacobi_dn;
  using       std::laguerre;
  using __gnu_cxx::ldouble_factorial;
  using       std::legendre;
  using __gnu_cxx::legendre_q;
  using __gnu_cxx::lfactorial;
  using __gnu_cxx::lpochhammer_l;
  using __gnu_cxx::lpochhammer_u;
  using __gnu_cxx::pochhammer_l;
  using __gnu_cxx::pochhammer_u;
  using __gnu_cxx::psi;
  using       std::riemann_zeta;
  using __gnu_cxx::sinc;
  using __gnu_cxx::sinhint;
  using __gnu_cxx::sinint;
  using       std::sph_bessel;
  using __gnu_cxx::sph_bessel_i;
  using __gnu_cxx::sph_bessel_k;
  using       std::sph_legendre;
  using       std::sph_neumann;
#else
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
#endif // STD

#if STD
  //  Airy Ai functions.
  std::cout << "airy_a" << std::endl;
  basename = "diff_airy_ai";
  rundiff(airy_ai, gsl::airy_ai, basename,
	  "x", fill_argument(std::make_pair(-10.0, +10.0),
			     std::make_pair(true, true), 41));

  //  Airy Bi functions.
  std::cout << "airy_bi" << std::endl;
  basename = "diff_airy_bi";
  rundiff(airy_bi, gsl::airy_bi, basename,
	  "x", fill_argument(std::make_pair(-10.0, +10.0),
			     std::make_pair(true, true), 41));
#endif // STD

  //  Associated Laguerre polynomials.
  std::cout << "assoc_laguerre" << std::endl;
  basename = "diff_assoc_laguerre";
  rundiff(assoc_laguerre, gsl::laguerre_nm, basename,
	  "n", vorder, "m", vorder,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
	    		     std::make_pair(true, true)));


  //  Associated Legendre functions.
  std::cout << "assoc_legendre" << std::endl;
  basename = "diff_assoc_legendre";
  rundiff(assoc_legendre, gsl::legendre_Plm, basename,
	  "l", vorder, "m", vorder,
	  "x", fill_argument(std::make_pair(-1.0, 1.0),
	    		     std::make_pair(true, true), 1001));


  //  Beta function.
  std::cout << "beta" << std::endl;
  basename = "diff_beta";
  rundiff(beta, gsl::beta, basename,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
			     std::make_pair(false, true)),
	  "y", fill_argument(std::make_pair(0.0, 100.0),
			     std::make_pair(false, true)));


  //  Complete elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "comp_ellint_1" << std::endl;
  basename = "diff_comp_ellint_1";
  rundiff(comp_ellint_1, gsl::ellint_Kcomp, basename,
	  "k", fill_argument(std::make_pair(-1.0, 1.0),
			     std::make_pair(false, false)));


  //  Complete elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "comp_ellint_2" << std::endl;
  basename = "diff_comp_ellint_2";
  rundiff(comp_ellint_2, gsl::ellint_Ecomp, basename,
	  "k", fill_argument(std::make_pair(-1.0, 1.0),
			     std::make_pair(false, false)));


  //  Complete elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "comp_ellint_3" << std::endl;
  basename = "diff_comp_ellint_3";
  rundiff(comp_ellint_3, gsl::ellint_Pcomp, basename,
	  "k", fill_argument(std::make_pair(-1.0, 1.0),
			     std::make_pair(false, false)),
	  "nu", fill_argument(std::make_pair(0.0, 1.0),
			      std::make_pair(true, false), 11));


  //  Confluent hypergeometric functions.
  //  Skip the singularity at c = 0.
  std::cout << "conf_hyperg" << std::endl;
  basename = "diff_conf_hyperg";
  rundiff(conf_hyperg, gsl::hyperg_1F1, basename,
	  "a", vab,
	  "c", fill_argument(std::make_pair(0.0, 10.0),
			     std::make_pair(false, true), 11),
	  "x", fill_argument(std::make_pair(-10.0, 10.0),
			     std::make_pair(true, true), 201));


  //  Regular modified cylindrical Bessel functions.
  std::cout << "cyl_bessel_i" << std::endl;
  basename = "diff_cyl_bessel_i";
  rundiff(cyl_bessel_i, gsl::bessel_Inu, basename,
	  "nu", cborderd,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
			     std::make_pair(true, true), 1001));


  //  Cylindrical Bessel functions (of the first kind).
  std::cout << "cyl_bessel_j" << std::endl;
  basename = "diff_cyl_bessel_j";
  rundiff(cyl_bessel_j, gsl::bessel_Jnu, basename,
	  "nu", cborderd,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
			     std::make_pair(true, true), 1001));


  //  Irregular modified cylindrical Bessel functions.
  // Skip the pole at the origin.
  std::cout << "cyl_bessel_k" << std::endl;
  basename = "diff_cyl_bessel_k";
  rundiff(cyl_bessel_k, gsl::bessel_Knu, basename,
	  "nu", cborderd,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
			     std::make_pair(false, true), 1001));


  //  Cylindrical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "cyl_neumann" << std::endl;
  basename = "diff_cyl_neumann";
  rundiff(cyl_neumann, gsl::bessel_Ynu, basename,
	  "nu", cborderd,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
			     std::make_pair(false, true), 1001));


  //  Elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "ellint_1" << std::endl;
  basename = "diff_ellint_1";
  rundiff(ellint_1, gsl::ellint_F, basename,
	  "k", fill_argument(std::make_pair(-1.0, 1.0),
			     std::make_pair(false, false)),
	  "phi", vphid);


  //  Elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "ellint_2" << std::endl;
  basename = "diff_ellint_2";
  rundiff(ellint_2, gsl::ellint_E, basename,
	  "k", fill_argument(std::make_pair(-1.0, 1.0),
			     std::make_pair(false, false)),
	  "phi", vphid);


  //  Elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "ellint_3" << std::endl;
  basename = "diff_ellint_3";
  rundiff(ellint_3, gsl::ellint_P, basename,
	  "k", fill_argument(std::make_pair(-1.0, 1.0),
			     std::make_pair(false, false)),
	  "nu", fill_argument(std::make_pair(0.0, 1.0),
			      std::make_pair(true, false), 11),
	  "phi", vphid);


  //  Exponential integral.
  //  Skip the pole at 0.
  std::cout << "expint" << std::endl;
  basename = "diff_expint_neg";
  rundiff(expint, gsl::expint_Ei, basename,
	  "x", fill_argument(std::make_pair(-50.0, 0.0),
			     std::make_pair(true, false), 51));
  basename = "diff_expint_pos";
  rundiff(expint, gsl::expint_Ei, basename,
	  "x", fill_argument(std::make_pair(0.0, 50.0),
			     std::make_pair(false, true), 51));


  //  Hermite polynomials
  std::cout << "hermite" << std::endl;
  basename = "diff_hermite";
  rundiff(hermite, gsl::hermite, basename,
	  "n", vorder,
	  "x", fill_argument(std::make_pair(-10.0, 10.0),
	  		     std::make_pair(true, true)));


  //  Hypergeometric functions.
  //  Skip the singularity at c = 0.
  //  Skip the singularities at x = -1.
  std::cout << "hyperg" << std::endl;
  basename = "diff_hyperg";
  rundiff(hyperg, gsl::hyperg_2F1, basename,
	  "a", vab, "b", vab,
	  "c", fill_argument(std::make_pair(0.0, 10.0),
			     std::make_pair(false, true), 11),
	  "x", fill_argument(std::make_pair(-1.0, 1.0),
			     std::make_pair(true, false), 21));


  //  Laguerre polynomials.
  std::cout << "laguerre" << std::endl;
  basename = "diff_laguerre";
  rundiff(laguerre, gsl::laguerre_n, basename,
	  "n", vorder,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
	  		     std::make_pair(true, true), 1001));


  //  Legendre polynomials.
  std::cout << "legendre" << std::endl;
  basename = "diff_legendre";
  rundiff(legendre, gsl::legendre_Pl, basename,
	  "l", vorder,
	  "x", fill_argument(std::make_pair(-1.0, 1.0),
	  		     std::make_pair(true, true), 1001));


  //  Riemann zeta function.
  //  Skip the pole at 1.
  std::cout << "riemann_zeta" << std::endl;
  basename = "diff_riemann_zeta_neg";
  rundiff(riemann_zeta, gsl::zeta, basename,
	  "x", fill_argument(std::make_pair(-10.0, 1.0),
			     std::make_pair(true, false), 56));
  basename = "diff_riemann_zeta_pos";
  rundiff(riemann_zeta, gsl::zeta, basename,
	  "x", fill_argument(std::make_pair(1.0, 30.0),
			     std::make_pair(false, true), 146));

#if STD
  //  Hurwitz zeta function.
  std::cout << "hurwitz_zeta" << std::endl;
  //  Skip the pole at 1.
  basename = "diff_hurwitz_zeta";
  rundiff(hurwitz_zeta, gsl::hzeta, basename,
	  "s", fill_argument(std::make_pair(1.0, 30.0),
			     std::make_pair(false, true), 146),
	  "a", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(false, true), 26));
#endif // STD


  //  Spherical Bessel functions.
  std::cout << "sph_bessel" << std::endl;
  basename = "diff_sph_bessel";
  rundiff(sph_bessel, gsl::bessel_jl, basename,
	  "n", sborder,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
	  		     std::make_pair(true, true), 1001));


  //  Spherical Legendre functions.
  std::cout << "sph_legendre" << std::endl;
  basename = "diff_sph_legendre";
  rundiff(sph_legendre, gsl::legendre_sphPlm, basename,
	  "l", vorder, "m", vorder,
	  "theta", fill_argument(std::make_pair(0.0, static_cast<_TpGSL>(M_PI)),
	    			 std::make_pair(true, true), 1001));


  //  Spherical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "sph_neumann" << std::endl;
  basename = "diff_sph_neumann";
  rundiff(sph_neumann, gsl::bessel_yl, basename,
	  "n", sborder,
	  "x", fill_argument(std::make_pair(0.0, 100.0),
	  		     std::make_pair(false, true), 1001));

#if STD
  //  Carlson elliptic functions R_C.
  std::cout << "ellint_rc" << std::endl;
  basename = "diff_ellint_rc";
  rundiff(ellint_rc, gsl::ellint_RC, basename,
	  "x", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(false, true), 11),
	  "y", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(false, true), 11));

  //  Carlson elliptic functions R_D.
  std::cout << "ellint_rd" << std::endl;
  basename = "diff_ellint_rd";
  rundiff(ellint_rd, gsl::ellint_RD, basename,
	  "x", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(false, true), 11),
	  "y", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(true, true), 11),
	  "z", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(false, true), 11));

  //  Carlson elliptic functions R_F.
  std::cout << "ellint_rf" << std::endl;
  basename = "diff_ellint_rf";
  rundiff(ellint_rf, gsl::ellint_RF, basename,
	  "x", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(false, true), 11),
	  "y", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(true, true), 11),
	  "z", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(false, true), 11));

  //  Carlson elliptic functions R_J.
  std::cout << "ellint_rj" << std::endl;
  basename = "diff_ellint_rj";
  rundiff(ellint_rj, gsl::ellint_RJ, basename,
	  "x", fill_argument(std::make_pair(0.0, 5.0),
			     std::make_pair(false, true), 11),
	  "y", fill_argument(std::make_pair(0.0, 5.0),
		 	     std::make_pair(true, true), 11),
	  "z", fill_argument(std::make_pair(0.0, 5.0),
		 	     std::make_pair(false, true), 11),
	  "p", fill_argument(std::make_pair(0.0, 5.0),
		 	     std::make_pair(false, true), 11));

  //  Dilogarithm functions.
  std::cout << "dilog" << std::endl;
  basename = "diff_dilog";
  rundiff(dilog, gsl::dilog, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{-10}, _TpGSL{1}),
			     std::make_pair(true, true), 23));

  //  Upper incomplete Gamma functions.
  std::cout << "gamma_u" << std::endl;
  basename = "diff_gamma_u";
  rundiff(gamma_u, gsl::gamma_inc, basename,
	  "a", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{5}),
			     std::make_pair(false, true), 11),
	  "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{5}),
			     std::make_pair(true, true), 11));

  //  Incomplete Beta functions.
  std::cout << "ibeta" << std::endl;
  basename = "diff_ibeta";
  rundiff(ibeta, gsl::beta_inc, basename,
	  "a", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{5}),
			     std::make_pair(false, true), 11),
	  "b", fill_argument(std::make_pair(_TpGSL{5}, _TpGSL{0}),
			     std::make_pair(true, false), 11),
	  "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{1}),
			     std::make_pair(false, false), 21));

  //  Digamma or psi functions.
  std::cout << "psi" << std::endl;
  basename = "diff_psi";
  rundiff(psi, gsl::psi, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{-9.9375}, _TpGSL{10.0625}),
			     std::make_pair(true, true), 801));

  //  Sine integral or Si functions.
  std::cout << "sinint" << std::endl;
  basename = "diff_sinint";
  rundiff(sinint, gsl::Si, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+10}),
			     std::make_pair(false, true), 101));

  //  Cosine integral or Ci functions.
  std::cout << "cosint" << std::endl;
  basename = "diff_cosint";
  rundiff(cosint, gsl::Ci, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+10}),
			     std::make_pair(false, true), 101));

  //  Hyperbolic sine integral or Shi functions.
  std::cout << "sinhint" << std::endl;
  basename = "diff_sinhint";
  rundiff(sinhint, gsl::Shi, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+5}),
			     std::make_pair(false, true), 101));

  //  Hyperbolic cosine integral or Chi functions.
  std::cout << "coshint" << std::endl;
  basename = "diff_coshint";
  rundiff(coshint, gsl::Chi, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+5}),
			     std::make_pair(false, true), 101));

  //  Fresnel cosine integral.
  std::cout << "fresnel_c" << std::endl;
  basename = "diff_fresnel_c";
  rundiff(fresnel_c, gsl::fresnel_c, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{-20}, _TpGSL{+20}),
			     std::make_pair(false, true), 401));

  //  Fresnel sine integral.
  std::cout << "fresnel_s" << std::endl;
  basename = "diff_fresnel_s";
  rundiff(fresnel_s, gsl::fresnel_s, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{-20}, _TpGSL{+20}),
			     std::make_pair(false, true), 401));

  //  Dawson integral.
  std::cout << "dawson" << std::endl;
  basename = "diff_dawson";
  rundiff(dawson, gsl::dawson, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+5}),
			     std::make_pair(false, true), 101));

  //  Exponential integral E1.
  std::cout << "expint_e1" << std::endl;
  basename = "diff_expint_e1";
  rundiff(expint_e1, gsl::expint_E1, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{0}, _TpGSL{+5}),
			     std::make_pair(false, true), 101));

  //  Sine cardinal function.
  std::cout << "sinc" << std::endl;
  basename = "diff_sinc";
  rundiff(sinc, gsl::sinc, basename,
	  "x", fill_argument(std::make_pair(_TpGSL{-20}, _TpGSL{+20}),
			     std::make_pair(false, true), 401));

#endif // STD

  return 0;
}

