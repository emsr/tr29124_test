
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <stdexcept>

#define STD 1

#if STD
#  if LOCAL
#    include "cmath_local"
#  else
#    include <cmath>
#  endif
#  include <array>
#else
#  include <cmath>
#  include <tr1/cmath>
#  include <tr1/array>
#endif

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
#else
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
#endif // STD

#if STD
  //  Airy Ai functions.
  std::cout << "airy_a" << std::endl;
  basename = "diff_airy_ai";
  rundiff<_TpGSL, _TpGSL>(airy_ai, wrap_gsl_sf_airy_ai, basename,
			  "x", fill_argument(std::make_pair(-10.0, +10.0),
					     std::make_pair(true, true), 41));

  //  Airy Bi functions.
  std::cout << "airy_bi" << std::endl;
  basename = "diff_airy_bi";
  rundiff<_TpGSL, _TpGSL>(airy_bi, wrap_gsl_sf_airy_bi, basename,
			  "x", fill_argument(std::make_pair(-10.0, +10.0),
					     std::make_pair(true, true), 41));
#endif // STD

  //  Associated Laguerre polynomials.
  std::cout << "assoc_laguerre" << std::endl;
  basename = "diff_assoc_laguerre";
  rundiff<_TpGSL, unsigned int, unsigned int, _TpGSL>(assoc_laguerre, wrap_gsl_sf_laguerre_nm, basename,
						      "n", vorder, "m", vorder,
						      "x", fill_argument(std::make_pair(0.0, 100.0),
									 std::make_pair(true, true)));


  //  Associated Legendre functions.
  std::cout << "assoc_legendre" << std::endl;
  basename = "diff_assoc_legendre";
  rundiff<_TpGSL, unsigned int, unsigned int, _TpGSL>(assoc_legendre, wrap_gsl_sf_legendre_Plm, basename,
						      "l", vorder, "m", vorder,
						      "x", fill_argument(std::make_pair(-1.0, 1.0),
									 std::make_pair(true, true), 1001));


  //  Beta function.
  std::cout << "beta" << std::endl;
  basename = "diff_beta";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(beta, wrap_gsl_sf_beta, basename,
				  "x", fill_argument(std::make_pair(0.0, 100.0),
						     std::make_pair(false, true)),
				  "y", fill_argument(std::make_pair(0.0, 100.0),
						     std::make_pair(false, true)));


  //  Complete elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "comp_ellint_1" << std::endl;
  basename = "diff_comp_ellint_1";
  rundiff<_TpGSL, _TpGSL>(comp_ellint_1, wrap_gsl_sf_ellint_Kcomp, basename,
			  "k", fill_argument(std::make_pair(-1.0, 1.0),
					     std::make_pair(false, false)));


  //  Complete elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "comp_ellint_2" << std::endl;
  basename = "diff_comp_ellint_2";
  rundiff<_TpGSL, _TpGSL>(comp_ellint_2, wrap_gsl_sf_ellint_Ecomp, basename,
			  "k", fill_argument(std::make_pair(-1.0, 1.0),
					     std::make_pair(false, false)));


  //  Complete elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "comp_ellint_3" << std::endl;
  basename = "diff_comp_ellint_3";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(comp_ellint_3, wrap_gsl_sf_ellint_Pcomp, basename,
				  "k", fill_argument(std::make_pair(-1.0, 1.0),
						     std::make_pair(false, false)),
				  "nu", fill_argument(std::make_pair(0.0, 1.0),
						      std::make_pair(true, false), 11));


  //  Confluent hypergeometric functions.
  //  Skip the singularity at c = 0.
  std::cout << "conf_hyperg" << std::endl;
  basename = "diff_conf_hyperg";
  rundiff<_TpGSL, _TpGSL, _TpGSL, _TpGSL>(conf_hyperg, wrap_gsl_sf_hyperg_1F1, basename,
					  "a", vab,
					  "c", fill_argument(std::make_pair(0.0, 10.0),
							     std::make_pair(false, true), 11),
					  "x", fill_argument(std::make_pair(-10.0, 10.0),
							     std::make_pair(true, true), 201));


  //  Regular modified cylindrical Bessel functions.
  std::cout << "cyl_bessel_i" << std::endl;
  basename = "diff_cyl_bessel_i";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(cyl_bessel_i, wrap_gsl_sf_bessel_Inu, basename,
				  "nu", cborderd,
				  "x", fill_argument(std::make_pair(0.0, 100.0),
						     std::make_pair(true, true), 1001));


  //  Cylindrical Bessel functions (of the first kind).
  std::cout << "cyl_bessel_j" << std::endl;
  basename = "diff_cyl_bessel_j";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(cyl_bessel_j, wrap_gsl_sf_bessel_Jnu, basename,
				  "nu", cborderd,
				  "x", fill_argument(std::make_pair(0.0, 100.0),
						     std::make_pair(true, true), 1001));


  //  Irregular modified cylindrical Bessel functions.
  // Skip the pole at the origin.
  std::cout << "cyl_bessel_k" << std::endl;
  basename = "diff_cyl_bessel_k";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(cyl_bessel_k, wrap_gsl_sf_bessel_Knu, basename,
				  "nu", cborderd,
				  "x", fill_argument(std::make_pair(0.0, 100.0),
						     std::make_pair(false, true), 1001));


  //  Cylindrical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "cyl_neumann" << std::endl;
  basename = "diff_cyl_neumann";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(cyl_neumann, wrap_gsl_sf_bessel_Ynu, basename,
				  "nu", cborderd,
				  "x", fill_argument(std::make_pair(0.0, 100.0),
						     std::make_pair(false, true), 1001));


  //  Elliptic integrals of the first kind.
  //  Avoid poles at |x| = 1.
  std::cout << "ellint_1" << std::endl;
  basename = "diff_ellint_1";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(ellint_1, wrap_gsl_sf_ellint_F, basename,
				  "k", fill_argument(std::make_pair(-1.0, 1.0),
						     std::make_pair(false, false)),
				  "phi", vphid);


  //  Elliptic integrals of the second kind.
  //  Avoid poles at |x| = 1.
  std::cout << "ellint_2" << std::endl;
  basename = "diff_ellint_2";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(ellint_2, wrap_gsl_sf_ellint_E, basename,
				  "k", fill_argument(std::make_pair(-1.0, 1.0),
						     std::make_pair(false, false)),
				  "phi", vphid);


  //  Elliptic integrals of the third kind.
  //  Avoid poles at |x| = 1 and at nu = 1.
  std::cout << "ellint_3" << std::endl;
  basename = "diff_ellint_3";
  rundiff<_TpGSL, _TpGSL, _TpGSL, _TpGSL>(ellint_3, wrap_gsl_sf_ellint_P, basename,
					  "k", fill_argument(std::make_pair(-1.0, 1.0),
							     std::make_pair(false, false)),
					  "nu", fill_argument(std::make_pair(0.0, 1.0),
							      std::make_pair(true, false), 11),
					  "phi", vphid);


  //  Exponential integral.
  //  Skip the pole at 0.
  std::cout << "expint" << std::endl;
  basename = "diff_expint_neg";
  rundiff<_TpGSL, _TpGSL>(expint, wrap_gsl_sf_expint_Ei, basename,
			  "x", fill_argument(std::make_pair(-50.0, 0.0),
					     std::make_pair(true, false), 51));
  basename = "diff_expint_pos";
  rundiff<_TpGSL, _TpGSL>(expint, wrap_gsl_sf_expint_Ei, basename,
			  "x", fill_argument(std::make_pair(0.0, 50.0),
					     std::make_pair(false, true), 51));


  //  Hermite polynomials
  std::cout << "hermite  UNTESTED" << std::endl;
//  basename = "diff_hermite";
//  rundiff<_TpGSL, unsigned int, _TpGSL>(hermite, gsl_xxx, basename, vorder,
//					fill_argument(std::make_pair(-10.0, 10.0),
//						      std::make_pair(true, true)));


  //  Hypergeometric functions.
  //  Skip the singularity at c = 0.
  //  Skip the singularities at x = -1.
  std::cout << "hyperg" << std::endl;
  basename = "diff_hyperg";
  rundiff<_TpGSL, _TpGSL, _TpGSL, _TpGSL, _TpGSL>(hyperg, wrap_gsl_sf_hyperg_2F1, basename,
						  "a", vab, "b", vab,
						  "c", fill_argument(std::make_pair(0.0, 10.0),
								     std::make_pair(false, true), 11),
						  "x", fill_argument(std::make_pair(-1.0, 1.0),
								     std::make_pair(true, false), 21));


  //  Laguerre polynomials.
  std::cout << "laguerre" << std::endl;
  basename = "diff_laguerre";
  rundiff<_TpGSL, unsigned int, _TpGSL>(laguerre, wrap_gsl_sf_laguerre_n, basename,
					"n", vorder,
					"x", fill_argument(std::make_pair(0.0, 100.0),
							   std::make_pair(true, true), 1001));


  //  Legendre polynomials.
  std::cout << "legendre" << std::endl;
  basename = "diff_legendre";
  rundiff<_TpGSL, unsigned int, _TpGSL>(legendre, wrap_gsl_sf_legendre_Pl, basename,
					"l", vorder,
					"x", fill_argument(std::make_pair(-1.0, 1.0),
							   std::make_pair(true, true), 1001));


  //  Riemann zeta function.
  //  Skip the pole at 1.
  std::cout << "riemann_zeta" << std::endl;
  basename = "diff_riemann_zeta_neg";
  rundiff<_TpGSL, _TpGSL>(riemann_zeta, wrap_gsl_sf_zeta, basename,
			  "x", fill_argument(std::make_pair(-10.0, 1.0),
					     std::make_pair(true, false), 56));
  basename = "diff_riemann_zeta_pos";
  rundiff<_TpGSL, _TpGSL>(riemann_zeta, wrap_gsl_sf_zeta, basename,
			  "x", fill_argument(std::make_pair(1.0, 30.0),
					     std::make_pair(false, true), 146));

#if STD
  //  Hurwitz zeta function.
  std::cout << "hurwitz_zeta" << std::endl;
  //  Skip the pole at 1.
  basename = "diff_hurwitz_zeta";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(hurwitz_zeta, wrap_gsl_sf_hzeta, basename,
				  "s", fill_argument(std::make_pair(1.0, 30.0),
						     std::make_pair(false, true), 146),
				  "a", fill_argument(std::make_pair(0.0, 5.0),
						     std::make_pair(false, true), 26));
#endif // STD


  //  Spherical Bessel functions.
  std::cout << "sph_bessel" << std::endl;
  basename = "diff_sph_bessel";
  rundiff<_TpGSL, unsigned int, _TpGSL>(sph_bessel, wrap_gsl_sf_bessel_jl, basename,
					"n", sborder,
					"x", fill_argument(std::make_pair(0.0, 100.0),
							   std::make_pair(true, true), 1001));


  //  Spherical Legendre functions.
  std::cout << "sph_legendre" << std::endl;
  basename = "diff_sph_legendre";
  rundiff<_TpGSL, unsigned int, unsigned int, _TpGSL>(sph_legendre, wrap_gsl_sf_legendre_sphPlm, basename,
						      "l", vorder, "m", vorder,
						      "theta", fill_argument(std::make_pair(0.0, static_cast<_TpGSL>(M_PI)),
									     std::make_pair(true, true), 1001));


  //  Spherical Neumann functions.
  // Skip the pole at the origin.
  std::cout << "sph_neumann" << std::endl;
  basename = "diff_sph_neumann";
  rundiff<_TpGSL, unsigned int, _TpGSL>(sph_neumann, wrap_gsl_sf_bessel_yl, basename,
					"n", sborder,
					"x", fill_argument(std::make_pair(0.0, 100.0),
							   std::make_pair(false, true), 1001));

#if STD
  //  Carlson elliptic functions R_C.
  std::cout << "ellint_rc" << std::endl;
  basename = "diff_ellint_rc";
  rundiff<_TpGSL, _TpGSL, _TpGSL>(ellint_rc, wrap_gsl_sf_ellint_RC, basename,
				  "x", fill_argument(std::make_pair(0.0, 5.0),
						     std::make_pair(false, true), 11),
				  "y", fill_argument(std::make_pair(0.0, 5.0),
						     std::make_pair(false, true), 11));

  //  Carlson elliptic functions R_D.
  std::cout << "ellint_rd" << std::endl;
  basename = "diff_ellint_rd";
  rundiff<_TpGSL, _TpGSL, _TpGSL, _TpGSL>(ellint_rd, wrap_gsl_sf_ellint_RD, basename,
					  "x", fill_argument(std::make_pair(0.0, 5.0),
							     std::make_pair(false, true), 11),
					  "y", fill_argument(std::make_pair(0.0, 5.0),
							     std::make_pair(true, true), 11),
					  "z", fill_argument(std::make_pair(0.0, 5.0),
							     std::make_pair(false, true), 11));

  //  Carlson elliptic functions R_F.
  std::cout << "ellint_rf" << std::endl;
  basename = "diff_ellint_rf";
  rundiff<_TpGSL, _TpGSL, _TpGSL, _TpGSL>(ellint_rf, wrap_gsl_sf_ellint_RF, basename,
					  "x", fill_argument(std::make_pair(0.0, 5.0),
							     std::make_pair(false, true), 11),
					  "y", fill_argument(std::make_pair(0.0, 5.0),
							     std::make_pair(true, true), 11),
					  "z", fill_argument(std::make_pair(0.0, 5.0),
							     std::make_pair(false, true), 11));

  //  Carlson elliptic functions R_J.
  std::cout << "ellint_rj" << std::endl;
  basename = "diff_ellint_rj";
  rundiff<_TpGSL, _TpGSL, _TpGSL, _TpGSL, _TpGSL>(ellint_rj, wrap_gsl_sf_ellint_RJ, basename,
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
    basename = "dilog";
    rundiff<_TpGSL, _TpGSL>(dilog, wrap_gsl_sf_dilog, basename,
			    "x", fill_argument(std::make_pair(-_TpGSL{10}, _TpGSL{1}),
					       std::make_pair(true, true), 23));

    //  Upper incomplete Gamma functions.
    std::cout << "gamma_u" << std::endl;
    basename = "gamma_u";
    rundiff<_TpGSL, _TpGSL, _TpGSL>(gamma_u, wrap_gsl_sf_gamma_inc, basename,
				    "a", fill_argument(std::make_pair(_TpGSL{0}, +_TpGSL{5}),
						       std::make_pair(false, true), 11),
				    "x", fill_argument(std::make_pair(_TpGSL{0}, +_TpGSL{5}),
						       std::make_pair(true, true), 11));

    //  Incomplete Beta functions.
    std::cout << "ibeta" << std::endl;
    basename = "ibeta";
    rundiff<_TpGSL, _TpGSL, _TpGSL, _TpGSL>(ibeta, wrap_gsl_sf_beta_inc, basename,
				    "a", fill_argument(std::make_pair(_TpGSL{0}, +_TpGSL{5}),
						       std::make_pair(false, true), 11),
				    "b", fill_argument(std::make_pair(_TpGSL{5}, +_TpGSL{0}),
						       std::make_pair(false, true), 11),
				    "x", fill_argument(std::make_pair(_TpGSL{0}, +_TpGSL{1}),
						       std::make_pair(false, false), 21));

    //  Digamma or psi functions.
    std::cout << "psi" << std::endl;
    basename = "psi";
    rundiff<_TpGSL, _TpGSL>(psi, wrap_gsl_sf_psi, basename,
			    "x", fill_argument(std::make_pair(-_TpGSL{10}, _TpGSL{10}),
					       std::make_pair(true, true), 41));
#endif // STD

  return 0;
}

