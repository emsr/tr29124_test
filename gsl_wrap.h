#include <gsl/gsl_sf.h>

#include "gslextras/Fresnel/fresnel.h"

namespace gsl
{

/// Airy functions.
double airy_ai(double x);
double airy_bi(double x);

/// Associated Laguerre polynomials.
double laguerre_nm(unsigned int n, unsigned int m, double x);


/// Associated Legendre functions.
double legendre_Plm(unsigned int l, unsigned int m, double x);


/// Beta function.
double beta(double x, double y);


/// Complete elliptic integrals of the first kind.
double ellint_Kcomp(double k);


/// Complete elliptic integrals of the second kind.
double ellint_Ecomp(double k);


/// Complete elliptic integrals of the third kind.
double ellint_Pcomp(double k, double nu);


/// Confluent hypergeometric functions.
double hyperg_1F1(double a, double c, double x);

/// Confluent hypergeometric limit functions.
double hyperg_0F1(double c, double x);


/// Regular modified cylindrical Bessel functions.
double bessel_Inu(double nu, double x);


/// Cylindrical Bessel functions (of the first kind).
double bessel_Jnu(double nu, double x);

//double bessel_Jnu_asymp(double nu, double x);


/// Irregular modified cylindrical Bessel functions.
double bessel_Knu(double nu, double x);


/// Cylindrical Neumann functions.
double bessel_Ynu(double nu, double x);

//double bessel_Ynu_asymp(double nu, double x);


/// Elliptic integrals of the first kind.
double ellint_F(double k, double phi);


/// Elliptic integrals of the second kind.
double ellint_E(double k, double phi);


/// Elliptic integrals of the third kind.
double ellint_P(double k, double nu, double phi);

/// Carlson elliptic integrals.
double ellint_RC(double x, double y);
double ellint_RD(double x, double y, double z);
double ellint_RF(double x, double y, double z);
double ellint_RJ(double x, double y, double z, double p);

/// Exponential integrals.
double expint_Ei(double x);

double expint_E1(double x);

double expint_En(unsigned int n, double x);

/// Hermite polynomials.
double hermite(unsigned int n, double x);

/// Hypergeometric functions.
double hyperg_2F1(double a, double b, double c, double x);


/// Laguerre polynomials.
double laguerre_n(unsigned int n, double x);


/// Legendre polynomials.
double legendre_Pl(unsigned int l, double x);


/// Riemann zeta function.
double zeta(double x);

/// Hurwitz zeta function.
double hzeta(double s, double q);

/// Spherical Bessel functions.
double bessel_jl(unsigned int n, double x);


/// Spherical Legendre functions.
double legendre_sphPlm(unsigned int l, unsigned int m, double theta);


/// Spherical Neumann functions.
double bessel_yl(unsigned int n, double x);

/// Normalized incomplete gamma functions.
double gamma_inc_Q(double a, double x);

/// Complementary normalized incomplete gamma functions.
double gamma_inc_P(double a, double x);

/// Non-normalized incomplete gamma functions.
double gamma_inc(double a, double x);

/// Incomplete beta functions.
double beta_inc(double a, double b, double x);

/// Dilogarithm functions.
double dilog(double x);

/// Digamma or psi functions.
double psi(double x);

/// Sine integral.
double Si(double x);

/// Cosine integral.
double Ci(double x);

/// Hyperbolic sine integral.
double Shi(double x);

/// Hyperbolic cosine integral.
double Chi(double x);

/// Gegenbauer polynomials.
double gegenpoly_n(int n, double lambda, double x);

/// Hydrogen wave functions.
double hydrogen(int n, double l, double Z, double r);

/// Dawson integral.
double dawson(double x);

/// Jacobian elliptic integrals.
double elljac_sn(double u, double m);

/// Jacobian elliptic integrals.
double elljac_cn(double u, double m);

/// Jacobian elliptic integrals.
double elljac_dn(double u, double m);

/// Fresnel cosine integral.
double fresnel_c(double x);

/// Fresnel sine integral.
double fresnel_s(double x);

/// Sinus cardinal function.
double sinc(double x);

} // namespace gsl

