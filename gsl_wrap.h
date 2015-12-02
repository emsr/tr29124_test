#include <gsl/gsl_sf.h>
//int gsl_sf_bessel_Jnu_asympx_e(const double nu, const double x, gsl_sf_result * result);
//int gsl_sf_bessel_Ynu_asympx_e(const double nu, const double x, gsl_sf_result * result);

///  Airy functions.
double wrap_gsl_sf_airy_ai(double x);
double wrap_gsl_sf_airy_bi(double x);

///  Associated Laguerre polynomials.
double wrap_gsl_sf_laguerre_nm(unsigned int n, unsigned int m, double x);


///  Associated Legendre functions.
double wrap_gsl_sf_legendre_Plm(unsigned int l, unsigned int m, double x);


///  Beta function.
double wrap_gsl_sf_beta(double x, double y);


///  Complete elliptic integrals of the first kind.
double wrap_gsl_sf_ellint_Kcomp(double k);


///  Complete elliptic integrals of the second kind.
double wrap_gsl_sf_ellint_Ecomp(double k);


///  Complete elliptic integrals of the third kind.
double wrap_gsl_sf_ellint_Pcomp(double k, double nu);


///  Confluent hypergeometric functions.
double wrap_gsl_sf_hyperg_1F1(double a, double c, double x);


///  Regular modified cylindrical Bessel functions.
double wrap_gsl_sf_bessel_Inu(double nu, double x);


///  Cylindrical Bessel functions (of the first kind).
double wrap_gsl_sf_bessel_Jnu(double nu, double x);

//double wrap_gsl_sf_bessel_Jnu_asymp(double nu, double x);


///  Irregular modified cylindrical Bessel functions.
double wrap_gsl_sf_bessel_Knu(double nu, double x);


///  Cylindrical Neumann functions.
double wrap_gsl_sf_bessel_Ynu(double nu, double x);

//double wrap_gsl_sf_bessel_Ynu_asymp(double nu, double x);


///  Elliptic integrals of the first kind.
double wrap_gsl_sf_ellint_F(double k, double phi);


///  Elliptic integrals of the second kind.
double wrap_gsl_sf_ellint_E(double k, double phi);


///  Elliptic integrals of the third kind.
double wrap_gsl_sf_ellint_P(double k, double nu, double phi);

///  Carlson elliptic integrals.
double wrap_gsl_sf_ellint_RC(double x, double y);
double wrap_gsl_sf_ellint_RD(double x, double y, double z);
double wrap_gsl_sf_ellint_RF(double x, double y, double z);
double wrap_gsl_sf_ellint_RJ(double x, double y, double z, double p);

///  Exponential integral.
double wrap_gsl_sf_expint_Ei(double x);


///  Hypergeometric functions.
double wrap_gsl_sf_hyperg_2F1(double a, double b, double c, double x);


///  Laguerre polynomials.
double wrap_gsl_sf_laguerre_n(unsigned int n, double x);


///  Legendre polynomials.
double wrap_gsl_sf_legendre_Pl(unsigned int l, double x);


///  Riemann zeta function.
double wrap_gsl_sf_zeta(double x);

///  Hurwitz zeta function.
double wrap_gsl_sf_hzeta(double s, double q);

///  Spherical Bessel functions.
double wrap_gsl_sf_bessel_jl(unsigned int n, double x);


///  Spherical Legendre functions.
double wrap_gsl_sf_legendre_sphPlm(unsigned int l, unsigned int m, double theta);


///  Spherical Neumann functions.
double wrap_gsl_sf_bessel_yl(unsigned int n, double x);

//  Normalized incomplete gamma functions.
double wrap_gsl_sf_gamma_inc_Q(double a, double x);

//  Complementary normalized incomplete gamma functions.
double wrap_gsl_sf_gamma_inc_P(double a, double x);

//  Non-normalized incomplete gamma functions.
double wrap_gsl_sf_gamma_inc(double a, double x);

//  Incomplete beta functions.
double wrap_gsl_sf_beta_inc(double a, double b, double x);

//  Dilogarithm functions.
double wrap_gsl_sf_dilog(double x);

//  Digamma or psi functions.
double wrap_gsl_sf_psi(double x);

