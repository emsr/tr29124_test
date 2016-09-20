#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/bessel.hpp> 
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/ellint_d.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/jacobi_zeta.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/owens_t.hpp>

#include "wrap_boost.h"

namespace beast
{

/// Airy Ai function.
double
airy_ai(double x)
{
  return boost::math::airy_ai(x);
}

/// Airy Bi function.
double
airy_bi(double x)
{
  return boost::math::airy_bi(x);
}

/// Associated Laguerre polynomials.
double
assoc_laguerre(unsigned int n, unsigned int m, double x)
{
  return boost::math::laguerre(n, m, x);
}

/// Associated Legendre functions.
double
assoc_legendre(unsigned int l, unsigned int m, double x)
{
  return boost::math::legendre_p(l, m, x);
}

/// Associated Legendre functions of the second kind.
double
assoc_legendre_q(unsigned int l, unsigned int m, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Beta functions.
double
beta(double x, double y)
{
  return boost::math::beta(x, y);
}

/// Complete elliptic integrals of the first kind.
double
comp_ellint_1(double k)
{
  return boost::math::ellint_1(k);
}

/// Complete elliptic integrals of the second kind.
double
comp_ellint_2(double k)
{
  return boost::math::ellint_2(k);
}

/// Complete elliptic integrals of the third kind.
double
comp_ellint_3(double k, double nu)
{
  return boost::math::ellint_3(k, nu);
}

/// Complete Legendre elliptic D integrals.
double
comp_ellint_d(double k)
{
  return boost::math::ellint_d(k);
}

/// Confluent hypergeometric functions.
double
conf_hyperg(double a, double c, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Confluent hypergeometric limit functions.
double
hyperg_0F1(double c, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Regular modified cylindrical Bessel functions.
double
cyl_bessel_i(double nu, double x)
{
  return boost::math::cyl_bessel_i(nu, x);
}

/// Cylindrical Bessel functions (of the first kind).
double
cyl_bessel_j(double nu, double x)
{
  return boost::math::cyl_bessel_j(nu, x);
}

/// Irregular modified cylindrical Bessel functions.
double
cyl_bessel_k(double nu, double x)
{
  return boost::math::cyl_bessel_k(nu, x);
}

/// Cylindrical Neumann functions.
double
cyl_neumann(double nu, double x)
{
  return boost::math::cyl_neumann(nu, x);
}

/// Elliptic integrals of the first kind.
double
ellint_1(double k, double phi)
{
  return boost::math::ellint_1(k, phi);
}

/// Elliptic integrals of the second kind.
double
ellint_2(double k, double phi)
{
  return boost::math::ellint_2(k, phi);
}

/// Elliptic integrals of the third kind.
double
ellint_3(double k, double nu, double phi)
{
  return boost::math::ellint_3(k, nu, phi);
}

/// Legendre elliptic D integrals.
double
ellint_d(double k, double phi)
{
  return boost::math::ellint_d(k, phi);
}

/// Carlson elliptic integrals R_C.
double
ellint_rc(double x, double y)
{
  return boost::math::ellint_rc(x, y);
}

/// Carlson elliptic integrals R_D.
double
ellint_rd(double x, double y, double z)
{
  return boost::math::ellint_rd(x, y, z);
}

/// Carlson elliptic integrals R_F.
double
ellint_rf(double x, double y, double z)
{
  return boost::math::ellint_rf(x, y, z);
}

/// Carlson elliptic integrals R_G.
double
ellint_rg(double x, double y, double z)
{
  return boost::math::ellint_rg(x, y, z);
}

/// Carlson elliptic integrals R_J.
double
ellint_rj(double x, double y, double z, double p)
{
  return boost::math::ellint_rj(x, y, z, p);
}

/// Exponential integral Ei.
double
expint(double x)
{
  return boost::math::expint(x);
}

/// Exponential integrals E_n.
double
expint(unsigned int n, double x)
{
  return boost::math::expint(n, x);
}

/// Hermite polynomials.
double
hermite(unsigned int n, double x)
{
  return boost::math::hermite(n, x);
}

/// Hypergeometric functions.
double
hyperg(double a, double b, double c, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Laguerre polynomials.
double
laguerre(unsigned int n, double x)
{
  return boost::math::laguerre(n, x);
}

/// Legendre polynomials.
double
legendre_p(unsigned int l, double x)
{
  return boost::math::legendre_p(l, x);
}

/// Legendre functions of the second kind.
double
legendre_q(unsigned int l, double x)
{
  return boost::math::legendre_q(l, x);
}

/// Riemann zeta function.
double
riemann_zeta(double x)
{
  return boost::math::zeta(x);
}

/// Hurwitz zeta functions.
double
hurwitz_zeta(double s, double q)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Dirichlet eta function.
double dirichlet_eta(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Bessel functions.
double
sph_bessel(unsigned int n, double x)
{
  return boost::math::sph_bessel(n, x);
}

/// Spherical Legendre functions.
double
legendre_sphPlm(unsigned int l, unsigned int m, double theta)
{
  return boost::math::spherical_harmonic_r(l, m, theta, 0.0);
}

/// Spherical Neumann functions.
double
sph_neumann(unsigned int n, double x)
{
  return boost::math::sph_neumann(n, x);
}

/// Non-normalized upper incomplete gamma functions.
double
tgamma(double a, double x)
{
  return boost::math::tgamma(a, x);
}

/// Non-normalized lower incomplete gamma functions.
double
tgamma_lower(double a, double x)
{
  return boost::math::tgamma_lower(a, x);
}

/// Normalized upper incomplete gamma functions.
double
qgamma(double a, double x)
{
  return boost::math::gamma_q(a, x);
}

/// Inverse normalized upper incomplete gamma functions.
double
qgamma_inv(double a, double q)
{
  return boost::math::gamma_q_inv(a, q);
}

/// Inverse parameter normalized upper incomplete gamma functions.
double
qgamma_inva(double a, double q)
{
  return boost::math::gamma_q_inva(a, q);
}

/// Normalized lower incomplete gamma functions.
double
pgamma(double a, double x)
{
  return boost::math::gamma_p(a, x);
}

/// Inverse normalized lower incomplete gamma functions.
double
pgamma_inv(double a, double q)
{
  return boost::math::gamma_p_inv(a, q);
}

/// Inverse parameter normalized lower incomplete gamma functions.
double
pgamma_inva(double a, double q)
{
  return boost::math::gamma_p_inva(a, q);
}

/// Incomplete beta functions.
double
ibeta(double a, double b, double x)
{
  return boost::math::ibeta(a, b, x);
}

double
ibeta_inv(double a, double b, double p)
{
  return boost::math::ibeta_inv(a, b, p);
}

double
ibetac_inv(double a, double b, double p)
{
  return boost::math::ibetac_inv(a, b, p);
}

double
ibeta_inva(double b, double x, double p)
{
  return boost::math::ibeta_inva(b, x, p);
}

double
ibetac_inva(double b, double x, double p)
{
  return boost::math::ibetac_inva(b, x, p);
}

double
ibeta_invb(double a, double x, double p)
{
  return boost::math::ibeta_invb(a, x, p);
}

double
ibetac_invb(double a, double x, double p)
{
  return boost::math::ibetac_invb(a, x, p);
}

/// Dilogarithm function.
double
dilog(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Digamma or psi function.
double
psi(double x)
{
  return boost::math::digamma(x);
}

/// Polygamma functions.
double
polygamma(int n, double x)
{
  return boost::math::polygamma(n, x);
}

/// Sine integral.
double
Si(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cosine integral.
double
Ci(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hyperbolic sine integral.
double
Shi(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hyperbolic cosine integral.
double
Chi(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Gegenbauer polynomials.
double
gegenpoly_n(unsigned int n, double lambda, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hydrogen wave functions.
double
hydrogen(int n, double l, double Z, double r)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Dawson integral.
double
dawson(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobian elliptic integrals sn.
double
jacobi_sn(double k, double u)
{
  return boost::math::jacobi_sn(k, u);
}

/// Jacobian elliptic integrals cn.
double
jacobi_cn(double k, double u)
{
  return boost::math::jacobi_cn(k, u);
}

/// Jacobian elliptic integrals dn.
double
jacobi_dn(double k, double u)
{
  return boost::math::jacobi_dn(k, u);
}

/// Fresnel cosine integral.
double
fresnel_c(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Fresnel sine integral.
double
fresnel_s(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Sinus cardinal function.
double
sinc_pi(double x)
{
  return boost::math::sinc_pi(x);
}

/// Normalized sinus cardinal function.
double
sinc(double x)
{
  return boost::math::sinc_pi(x / M_PI);
}

/// Hyperbolic sinus cardinal function.
double
sinhc_pi(double x)
{
  return std::sinh(x) / x;
}

/// Normalized hyperbolic sinus cardinal function.
double
sinhc(double x)
{
  return std::sinh(x / M_PI) / (x / M_PI);
}

/// Log upper Pochhammer symbol.
double
lpochhammer(double a, double x)
{
  return std::log(boost::math::rising_factorial(a, x));
}

/// Log lower Pochhammer symbol.
double
lpochhammer_lower(double a, double x)
{
  return std::log(boost::math::falling_factorial(a, x));
}

/// Upper Pochhammer symbol.
double
pochhammer(double a, double x)
{
  return boost::math::rising_factorial(a, x);
}

/// Lower Pochhammer symbol.
double
pochhammer_lower(double a, double x)
{
  return boost::math::falling_factorial(a, x);
}

/// Log factorial.
double
lfactorial(unsigned int n)
{
  return std::log(boost::math::factorial<double>(n));
}

/// Factorial.
double
factorial(unsigned int n)
{
  return boost::math::factorial<double>(n);
}

/// Log double factorial.
double
ldouble_factorial(int n)
{
  return std::log(boost::math::double_factorial<double>(n));
}

/// Double factorial.
double
double_factorial(unsigned int n)
{
  return boost::math::double_factorial<double>(n);
}

/// Regular modified spherical bessel functions.
double
sph_bessel_i(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Irregular modified spherical bessel functions.
double
sph_bessel_k(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Chebyshev polynomials of the first kind.
double
chebyshev_t(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Chebyshev polynomials of the second kind.
double
chebyshev_u(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Chebyshev polynomials of the third kind.
double
chebyshev_v(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Chebyshev polynomials of the fourth kind.
double
chebyshev_w(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobi polynomials.
double
jacobi(unsigned int n, double alpha, double beta, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Binomial coefficients.
double
choose(unsigned int n, unsigned int k)
{
  return std::numeric_limits<double>::quiet_NaN();//boost::math::binomial_coefficient(n, k);
}

/// Log binomial coefficients.
double
lnchoose(unsigned int n, unsigned int k)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Taylor coefficients.
double
taylorcoeff(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Radial polynomials
double
radpoly(unsigned int n, unsigned int m, double rho)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Zernike polynomials
double
zernike(unsigned int n, unsigned int m, double rho, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cylindrical Hankel functions of the first kind.
std::complex<double>
cyl_hankel_1(double nu, double x)
{
  return boost::math::cyl_hankel_1(nu, x);
}

/// Cylindrical Hankel functions of the second kind.
std::complex<double>
cyl_hankel_2(double nu, double x)
{
  return boost::math::cyl_hankel_2(nu, x);
}

/// Spherical Hankel functions of the first kind.
std::complex<double>
sph_hankel_1(unsigned int n, double x)
{
  return boost::math::sph_hankel_1(n, x);
}

/// Spherical Hankel functions of the second kind.
std::complex<double>
sph_hankel_2(unsigned int n, double x)
{
  return boost::math::sph_hankel_2(n, x);
}

/// Heuman lambda functions.
double
heuman_lambda(double k, double phi)
{
  return boost::math::heuman_lambda(k, phi);
}

/// Jacobi zeta functions.
double
jacobi_zeta(double k, double phi)
{
  return boost::math::jacobi_zeta(k, phi);
}

/// Inverse error function.
double
erf_inv(double p)
{
  return boost::math::erf_inv(p);
}

/// Inverse complementary error function.
double
erfc_inv(double p)
{
  return boost::math::erfc_inv(p);
}

/// Spherical harmonic functions.
std::complex<double>
sph_harmonic(unsigned int l, int m, double theta, double phi)
{
  return boost::math::spherical_harmonic(l, m, theta, phi);
}

/// Owen's T function.
double
owens_t(double h, double a)
{
  return boost::math::owens_t(h, a);
}

/// Clausen function of order 2.
double
clausen_c(unsigned int m, double w)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Struve H function.
double
struve_h(double nu, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Struve L function.
double
struve_l(double nu, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

} // namespace beast

