
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "wrap_pheces.h"

#include "cephes/cmath/cephes_cmath.h"
#include "cephes/cprob/cephes_cprob.h"
#include "cephes/misc/cephes_misc.h"
#include "cephes/ellf/cephes_ellf.h"
#include "cephes/polyn/cephes_polyn.h"
#include "cephes/bessel/cephes_bessel.h"

namespace pheces
{

/// Airy Ai function.
double
airy_ai(double x)
{
  double ai, aip, bi, bip;
  ::airy(x, &ai, &aip, &bi, &bip);
  return ai;
}

/// Airy Bi function.
double
airy_bi(double x)
{
  double ai, aip, bi, bip;
  ::airy(x, &ai, &aip, &bi, &bip);
  return bi;
}

/// Associated Laguerre polynomials.
double
assoc_laguerre(unsigned int /*n*/, unsigned int /*m*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Associated Legendre functions.
double
assoc_legendre(unsigned int /*l*/, unsigned int /*m*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Associated Legendre functions of the second kind.
double
assoc_legendre_q(unsigned int /*l*/, unsigned int /*m*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Beta functions.
double
beta(double x, double y)
{
  return ::beta(x, y);
}

/// Complete elliptic integrals of the first kind.
double
comp_ellint_1(double /*k*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complete elliptic integrals of the second kind.
double
comp_ellint_2(double /*k*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complete elliptic integrals of the third kind.
double
comp_ellint_3(double /*k*/, double /*nu*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complete Legendre elliptic D integrals.
double
comp_ellint_d(double /*k*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Confluent hypergeometric functions.
double
conf_hyperg(double a, double c, double x)
{
  return ::hyperg(a, c, x);
}

/// Confluent hypergeometric limit functions.
double
hyperg_0F1(double /*c*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hypergeometric functions.
double
hyperg_3F0(double a, double b, double c, double x)
{
  double err;
  return threef0(a, b, c, x, &err);
}

/// Regular modified cylindrical Bessel functions.
double
cyl_bessel_i(double /*nu*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cylindrical Bessel functions (of the first kind).
double
cyl_bessel_j(double nu, double x)
{
  return jv(nu, x);
}

/// Irregular modified cylindrical Bessel functions.
double
cyl_bessel_k(double nu, double x)
{
  return iv(nu, x);
}

/// Cylindrical Neumann functions.
double
cyl_neumann(double nu, double x)
{
  return yv(nu, x);
}

/// Elliptic integrals of the first kind.
double
ellint_1(double /*k*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Elliptic integrals of the second kind.
double
ellint_2(double /*k*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Elliptic integrals of the third kind.
double
ellint_3(double /*k*/, double /*nu*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Legendre elliptic D integrals.
double
ellint_d(double /*k*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_C.
double
ellint_rc(double /*x*/, double /*y*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_D.
double
ellint_rd(double /*x*/, double /*y*/, double /*z*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_F.
double
ellint_rf(double /*x*/, double /*y*/, double /*z*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_G.
double
ellint_rg(double /*x*/, double /*y*/, double /*z*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_J.
double
ellint_rj(double /*x*/, double /*y*/, double /*z*/, double /*p*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Exponential integral Ei.
double
expint(double x)
{
  return ei(x);
}

/// Exponential integrals E_n.
double
expint(unsigned int n, double x)
{
  return expn(n, x);
}

/// Hermite polynomials.
double
hermite(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hypergeometric functions.
double
hyperg(double a, double b, double c, double x)
{
  return ::hyp2f1(a, b, c, x);
}

/// Laguerre polynomials.
double
laguerre(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Legendre polynomials.
double
legendre_p(unsigned int /*l*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Legendre polynomials of the second kind.
double
legendre_q(unsigned int /*l*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Riemann zeta function.
double
riemann_zeta(double s)
{
  return ::zetac(s);
}

/// Hurwitz zeta functions.
double
hurwitz_zeta(double s, double q)
{
  return ::zeta(s, q);
}

/// Dirichlet eta function.
double
dirichlet_eta(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Bessel functions.
double
sph_bessel(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Legendre functions.
double
legendre_sphPlm(unsigned int /*l*/, unsigned int /*m*/, double /*theta*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Neumann functions.
double
sph_neumann(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Non-normalized lower incomplete gamma functions.
double
tgamma_lower(double a, double x)
{
  return ::gamma(x) * ::igam(a, x);
}

/// Non-normalized lower incomplete gamma functions.
double
tgamma(double a, double x)
{
  return ::gamma(x) * ::igamc(a, x);
}

/// Normalized incomlete gamma functions.
double
qgamma(double a, double x)
{
  return igamc(a, x);
}

/// Normalized incomlete gamma functions.
double
pgamma(double a, double x)
{
  return ::igam(a, x);
}

/// Incomlete beta functions.
double
ibeta(double a, double b, double x)
{
  return ::incbi(a, b, x);
}

/// Complementary incomlete beta functions.
double
ibetac(double a, double b, double x)
{
  return 1.0 - ::incbi(b, a, 1.0 - x);
}

/// Dilogarithm function.
double
dilog(double x)
{
  return ::polylog(0, x);
}

/// Digamma or psi function.
double
psi(double x)
{
  return ::spence(x);
}

/// Sine integral.
double
sinint(double x)
{
  double si, ci;
  ::sici(x, &si, &ci);
  return si;
}

/// Cosine integral.
double
cosint(double x)
{
  double si, ci;
  ::sici(x, &si, &ci);
  return ci;
}

/// Hyperbolic sine integral.
double
sinhint(double x)
{
  double shi, chi;
  ::shichi(x, &shi, &chi);
  return shi;
}

/// Hyperbolic cosine integral.
double
coshint(double x)
{
  double shi, chi;
  ::shichi(x, &shi, &chi);
  return chi;
}

/// Gegenbauer polynomials.
double
gegenpoly_n(unsigned int /*n*/, double /*lambda*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hydrogen wave functions.
double
hydrogen(unsigned int /*n*/, double /*l*/, double /*Z*/, double /*r*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Dawson integral.
double
dawson(double x)
{
  return ::dawsn(x);
}

/// Jacobian elliptic integrals sn.
double
jacobi_sn(double k, double u)
{
  double sn, cn, dn, ph;
  ::ellpj(u, k, &sn, &cn, &dn, &ph);
  return sn;
}

/// Jacobian elliptic integrals cn.
double
jacobi_cn(double k, double u)
{
  double sn, cn, dn, ph;
  ::ellpj(u, k, &sn, &cn, &dn, &ph);
  return cn;
}

/// Jacobian elliptic integrals dn.
double
jacobi_dn(double k, double u)
{
  double sn, cn, dn, ph;
  ::ellpj(u, k, &sn, &cn, &dn, &ph);
  return dn;
}

/// Fresnel cosine integral.
double
fresnel_c(double x)
{
  double S, C;
  ::fresnl(x, &S, &C);
  return C;
}

/// Fresnel sine integral.
double
fresnel_s(double x)
{
  double S, C;
  ::fresnl(x, &S, &C);
  return S;
}

/// Sinus cardinal function.
double
sinc(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Reperiodized sinus cardinal function.
double
sinc_pi(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hyperbolic sinus cardinal function.
double
sinhc(double x)
{
  return std::sinh(x) / x;
}

/// Reperiodized hyperbolic sinus cardinal function.
double
sinhc_pi(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log upper Pochhammer symbol.
double
lpochhammer(double /*a*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log lower Pochhammer symbol.
double
lpochhammer_lower(double /*a*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Upper Pochhammer symbol.
double
pochhammer(double /*a*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Lower Pochhammer symbol.
double
pochhammer_lower(double /*a*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log factorial.
double
lfactorial(unsigned int n)
{
  return std::log(fac(n));
}

/// Factorial.
double
factorial(unsigned int n)
{
  return ::fac(n);
}

/// Log double factorial.
double
ldouble_factorial(int /*n*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Double factorial.
double
double_factorial(int /*n*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Regular modified spherical bessel functions.
double
sph_bessel_i(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Irregular modified spherical bessel functions.
double
sph_bessel_k(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Chebyshev polynomials of the first kind.
double
chebyshev_t(unsigned int /*n*/, double /*x*/)
{ return std::numeric_limits<double>::quiet_NaN(); }

/// Binomial coefficients.
double
choose(unsigned int /*n*/, unsigned int /*k*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log binomial coefficients.
double
lnchoose(unsigned int /*n*/, unsigned int /*k*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Taylor coefficients.
double
taylorcoeff(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobi polynomials.
double
jacobi(unsigned int /*n*/, double /*alpha*/, double /*beta*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Radial polynomials
double
radpoly(unsigned int /*n*/, unsigned int /*m*/, double /*rho*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Zernike polynomials
double
zernike(unsigned int /*n*/, int /*m*/, double /*rho*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cylindrical Hankel functions of the first kind.
std::complex<double>
cyl_hankel_1(double /*nu*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cylindrical Hankel functions of the second kind.
std::complex<double>
cyl_hankel_2(double /*nu*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Hankel functions of the first kind.
std::complex<double>
sph_hankel_1(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Hankel functions of the second kind.
std::complex<double>
sph_hankel_2(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Heuman lambda functions.
double
heuman_lambda(double /*k*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobi zeta functions.
double
jacobi_zeta(double /*k*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse error function.
double
erf_inv(double /*p*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse complementary error function.
double
erfc_inv(double /*p*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical harmonic functions.
std::complex<double>
sph_harmonic(unsigned int /*l*/, int /*m*/, double /*theta*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Owen's T function.
double
owens_t(double /*h*/, double /*a*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Clausen function of order 2.
double
clausen_c(unsigned int /*m*/, double /*w*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Struve H function.
double
struve_h(double /*nu*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Struve L function.
double
struve_l(double /*nu*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Fermi-Dirac integrals
double
fermi_dirac(double /*s*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Bose-Einstein integrals
double
bose_einstein(double /*s*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

} // namespace pheces
