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
{ std::numeric_limits<double>::quiet_NaN(); }

/// Complete elliptic integrals of the second kind.
double
ellint_Ecomp(double k)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Complete elliptic integrals of the third kind.
double
ellint_Pcomp(double k, double nu)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Confluent hypergeometric functions.
double
hyperg_1F1(double a, double c, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Confluent hypergeometric limit functions.
double
conf_hyperg_lim(double c, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Regular modified cylindrical Bessel functions.
double
bessel_Inu(double nu, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Cylindrical Bessel functions (of the first kind).
double
bessel_Jnu(double nu, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Irregular modified cylindrical Bessel functions.
double
bessel_Knu(double nu, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Cylindrical Neumann functions.
double
bessel_Ynu(double nu, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

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
{ std::numeric_limits<double>::quiet_NaN(); }

/// Laguerre polynomials.
double
laguerre_n(unsigned int n, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

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
{ std::numeric_limits<double>::quiet_NaN(); }

/// Hurwitz zeta functions.
double
hzeta(double s, double q)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Spherical Bessel functions.
double
sph_bessel(unsigned int n, double x)
{ return sph_bessel_j(n, x); }

/// Spherical Legendre functions.
double
sph_legendre(unsigned int l, unsigned int m, double theta)
{ return spherical_harmonic(l, m, theta, phi); }

/// Spherical Neumann functions.
double
sph_neumann(unsigned int n, double x)
{ return sph_bessel_y(n, x); }

/// Normalized incomplete gamma functions.
double
gamma_q(double a, double x)
{ return gamma_q(a, x); }

/// Complementary normalized incomplete gamma functions.
double
gamma_p(double a, double x)
{ return gamma_p(a, x); }

/// Non-normalized incomplete gamma functions.
double
gamma_inc(double a, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Incomplete beta functions.
double
beta_inc(double a, double b, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Dilogarithm function.
double
dilog(double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Digamma or psi function.
double
digamma(double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Sine integral.
double
sinint(double x)
{
  double ci, si;
  cisi(x, &ci, &si);
  return si;
}

/// Cosine integral.
double
cosint(double x)
{
  double ci, si;
  cisi(x, &ci, &si);
  return ci;
}

/// Hyperbolic sine integral.
double
sinhint(double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Hyperbolic cosine integral.
double
coshint(double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Gegenbauer polynomials.
double
gegenbauer(unsigned int n, double lambda, double x)
{ return gegenbauer_poly(n, lambda, x); }

/// Hydrogen wave functions.
double
hydrogen(int n, double l, double Z, double r)
{ std::numeric_limits<double>::quiet_NaN(); }

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
{ std::numeric_limits<double>::quiet_NaN(); }

/// Reperiodized sinus cardinal function.
double
sinc_pi(double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Hyperbolic sinus cardinal function.
double
sinhc(double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Reperiodized hyperbolic sinus cardinal function.
double
sinhc_pi(double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Log upper Pochhammer symbol.
double
lnpoch(double a, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Upper Pochhammer symbol.
double
poch(double a, double x)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Log factorial.
double
lnfact(unsigned int n)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Factorial.
double
fact(unsigned int n)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Log double factorial.
double
lndoublefact(unsigned int n)
{ std::numeric_limits<double>::quiet_NaN(); }

/// Double factorial.
double
doublefact(unsigned int n)
{ std::numeric_limits<double>::quiet_NaN(); }

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