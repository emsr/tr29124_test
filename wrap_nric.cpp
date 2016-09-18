#include "wrap_nric.h"


namespace nric
{

/// Airy Ai function.
double
airy_ai(double x)
{
  double ai, bi, aip, bip;
  airy(x, &ai, &bi, &aip, &bip );
  return ai;
}

/// Airy Bi function.
double
airy_bi(double x)
{
  double ai, bi, aip, bip;
  airy(x, &ai, &bi, &aip, &bip );
  return bi;
}

/// Associated Laguerre polynomials.
double
laguerre_nm(unsigned int n, unsigned int m, double x)
{ return laguerre_poly(n, double(m), x); }

/// Associated Legendre functions.
double
assoc_legendre(unsigned int l, unsigned int m, double x)
{ return legendre_poly(l, m, x); }

/// Beta functions.
double
beta(double x, double y)
{ return beta(x, y); }

/// Complete elliptic integrals of the first kind.
double
ellint_Kcomp(double k)
{ return 0.0; }

/// Complete elliptic integrals of the second kind.
double
ellint_Ecomp(double k)
{ return 0.0; }

/// Complete elliptic integrals of the third kind.
double
ellint_Pcomp(double k, double nu)
{ return 0.0; }

/// Confluent hypergeometric functions.
double
hyperg_1F1(double a, double c, double x)
{ return 0.0; }

/// Confluent hypergeometric limit functions.
double
hyperg_0F1(double c, double x)
{ return 0.0; }

/// Regular modified cylindrical Bessel functions.
double
bessel_Inu(double nu, double x)
{ return 0.0; }

/// Cylindrical Bessel functions (of the first kind).
double
bessel_Jnu(double nu, double x)
{ return 0.0; }

/// Irregular modified cylindrical Bessel functions.
double
bessel_Knu(double nu, double x)
{ return 0.0; }

/// Cylindrical Neumann functions.
double
bessel_Ynu(double nu, double x)
{ return 0.0; }

/// Elliptic integrals of the first kind.
double
ellint_F(double k, double phi)
{ return legendre_f(phi, k); }

/// Elliptic integrals of the second kind.
double
ellint_E(double k, double phi)
{ return legendre_e(phi, k); }

/// Elliptic integrals of the third kind.
double
ellint_P(double k, double nu, double phi)
{ return legendre_pi(phi, nu, k); }

/// Carlson elliptic integrals R_C.
double
ellint_RC(double x, double y)
{ return carlson_rc(x, y); }

/// Carlson elliptic integrals R_D.
double
ellint_RD(double x, double y, double z)
{ return carlson_rd(x, y, z); }

/// Carlson elliptic integrals R_F.
double
ellint_RF(double x, double y, double z)
{ return carlson_rf(x, y, z); }

/// Carlson elliptic integrals R_J.
double
ellint_RJ(double x, double y, double z, double p)
{ return carlson_rj(x, y, z, p); }

/// Exponential integral Ei.
double
expint_Ei(double x)
{ return ei(x); }

/// Exponential integral E_1.
double
expint_E1(double x)
{ return exp_int(1, x); }

/// Exponential integrals E_n.
double
expint_En(unsigned int n, double x)
{ return exp_int(n, x); }

/// Hermite polynomials.
double
hermite(unsigned int n, double x)
{ return hermite_h(n, x); }

/// Hypergeometric functions.
double
hyperg_2F1(double a, double b, double c, double x)
{ return 0.0; }

/// Laguerre polynomials.
double
laguerre_n(unsigned int n, double x)
{ return 0.0; }

/// Legendre polynomials.
double
legendre_Pl(unsigned int l, double x)
{ return legendre_p(l, x); }

/// Legendre polynomials of the second kind.
double
legendre_Ql(unsigned int l, double x)
{ return legendre_q(l, x); }

/// Riemann zeta function.
double
zeta(double x)
{ return 0.0; }

/// Hurwitz zeta functions.
double
hzeta(double s, double q)
{ return 0.0; }

/// Spherical Bessel functions.
double
sph_bessel(unsigned int n, double x)
{ return sph_bessel_j(n, x); }

/// Spherical Legendre functions.
double
legendre_sphPlm(unsigned int l, unsigned int m, double theta)
{ return spherical_harmonic(l, m, theta, phi); }

/// Spherical Neumann functions.
double
sph_neumann(unsigned int n, double x)
{ return sph_bessel_y(n, x); }

/// Normalized incomplete gamma functions.
double
qgamma(double a, double x)
{ return gamma_q(a, x); }

/// Complementary normalized incomplete gamma functions.
double
pgamma(double a, double x)
{ return gamma_p(a, x); }

/// Non-normalized incomplete gamma functions.
double
gamma_inc(double a, double x)
{ return 0.0; }

/// Incomplete beta functions.
double
beta_inc(double a, double b, double x)
{ return 0.0; }

/// Dilogarithm function.
double
dilog(double x)
{ return 0.0; }

/// Digamma or psi function.
double
psi(double x)
{ return 0.0; }

/// Sine integral.
double
Si(double x)
{
  double ci, si;
  cisi(x, &ci, &si);
  return si;
}

/// Cosine integral.
double
Ci(double x)
{
  double ci, si;
  cisi(x, &ci, &si);
  return ci;
}

/// Hyperbolic sine integral.
double
Shi(double x)
{ return 0.0; }

/// Hyperbolic cosine integral.
double
Chi(double x)
{ return 0.0; }

/// Gegenbauer polynomials.
double
gegenpoly_n(int n, double lambda, double x)
{ return gegenbauer_poly(n, lambda, x); }

/// Hydrogen wave functions.
double
hydrogen(int n, double l, double Z, double r)
{ return 0.0; }

/// Dawson integral.
double
dawson(double x)
{ return ::dawson(x); }

/// Jacobian elliptic integrals sn.
double
elljac_sn(double u, double k)
{
  double mc = k * k;
  double sn, cn, dn;
  jacobian_sncndn(u, mc, &sn, &cn, &dn);
  return sn;
}

/// Jacobian elliptic integrals cn.
double
elljac_cn(double u, double k)
{
  double mc = k * k;
  double sn, cn, dn;
  jacobian_sncndn(u, mc, &sn, &cn, &dn);
  return cn;
}

/// Jacobian elliptic integrals dn.
double
elljac_dn(double u, double k)
{
  double mc = k * k;
  double sn, cn, dn;
  jacobian_sncndn(u, mc, &sn, &cn, &dn);
  return dn;
}

/// Fresnel cosine integral.
double
fresnel_c(double x)
{ return ::fresnel_c(x); }

/// Fresnel sine integral.
double
fresnel_s(double x)
{ return ::fresnel_s(x); }

/// Sinus cardinal function.
double
sinc(double x)
{ return 0.0; }

/// Log upper Pochhammer symbol.
double
lnpoch(double a, double x)
{ return 0.0; }

/// Upper Pochhammer symbol.
double
poch(double a, double x)
{ return 0.0; }

/// Log factorial.
double
lnfact(unsigned int n)
{ return 0.0; }

/// Factorial.
double
fact(unsigned int n)
{ return 0.0; }

/// Log double factorial.
double
lndoublefact(unsigned int n)
{ return 0.0; }

/// Double factorial.
double
doublefact(unsigned int n)
{ return 0.0; }

/// Regular modified spherical bessel functions.
double
sph_bessel_i(unsigned int n, double x)
{ return sph_bessel_i(n, x); }

/// Irregular modified spherical bessel functions.
double
sph_bessel_k(unsigned int n, double x)
{ return sph_bessel_k(n, x); }

/// Chebyshev polynomials of the first kind.
double
chebyshev_t(unsigned int n, double x)
{ return ::chebyshev_t(n, x); }

/// Jacobi polynomials.
double
jacobi(unsigned int n, double alpha, double beta, double x)
{ double jacobi_poly(n, alpha, beta, x); }

} // namespace nric
