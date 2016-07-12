#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <experimental/string_view>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <functional>
#include <utility>
#include <tuple>

#define STD 1

#include "specfun_testcase.h"
#include "gsl_wrap.h"
#include "boost_wrap.h"

#include "testcase.tcc"

std::string
get_filename(const std::string & path,
	     const std::string & prefix,
	     const std::string & basename,
	     const std::string & extra,
	     const std::string & suffix)
{
  auto filename = path + "/" + prefix;
  filename += basename + extra + suffix;

  return filename;
}


template<typename Real>
  void
  harness()
  {
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
    using __gnu_cxx::expint;
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

    // Unsigned integer orders for various polynomials, harmonics.
    std::vector<unsigned int> vorder{0, 1, 2, 5, 10, 20, 50, 100};

    // Integer orders for various polynomials, harmonics, and spherical bessels.
    std::vector<int> iorder{0, 1, 2, 5, 10, 20, 50, 100};

    // ... corresponding "Real" integer orders for GSL.
    std::vector<Real> dvorder{0, 1, 2, 5, 10, 20, 50, 100};

    // Orders for cylindrical Bessel functions.
    std::vector<Real> cyl_neg_order{-5, -2, -1, -Real{2.0L/3.0L},
				    -Real{0.5L}, -Real{1.0L/3.0L}};

    std::vector<Real> cyl_order{0, Real{1.0L/3.0L},
				Real{0.5L}, Real{2.0L/3.0L},
				1, 2, 5, 10, 20, 50, 100};

    // Orders for spherical bessel functions.
    std::vector<unsigned int> sph_order{0, 1, 2, 5, 10, 20, 50, 100};

    const unsigned int num_phi = 10;
    Real phi[num_phi];
    for (unsigned int i = 0; i < num_phi; ++i)
      phi[i] = Real{10} * i * __gnu_cxx::__math_constants<Real>::__pi / Real{180};
    std::vector<Real> vphid(phi, phi + num_phi);

    std::vector<Real> vab{0, Real{0.5L}, 1, 2, 5, 10, 20};

    const std::string path = ".";
    const std::string prefix = "/check_";

#if STD
    // Airy functions of the first kind.
    std::ofstream fairy_ai(path + prefix + "airy_ai" + ".cc");
    testcase<Real, Real>
      xairy_ai("testcase_airy",
	      test_function<Real, Real>("airy_ai", airy_ai),
	      baseline_function<Real, Real>("GSL", "gsl::airy_ai", gsl::airy_ai),
	      mask_function<Real>([](Real){ return true; }),
	      argument<Real>("x", fill_argument(std::make_pair(Real{-10}, Real{10}),
						std::make_pair(true, true), 41)));
    xairy_ai(fairy_ai);

#endif // STD

  }


int
main()
{
  harness<double>();

  return 0;
}


