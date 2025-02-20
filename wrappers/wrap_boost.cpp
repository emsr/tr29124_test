
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/ellint_d.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/jacobi_zeta.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/zeta.hpp>

#include <wrap_boost.h>

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
/// Boost returns Condon-Shortley phase.
double
assoc_legendre(unsigned int l, unsigned int m, double x)
{
  return ((m & 1) ? -1 : +1) * boost::math::legendre_p(l, m, x);
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
conf_hyperg(double /*a*/, double /*c*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Tricomi confluent hypergeometric functions.
double
tricomi_u(double /*a*/, double /*c*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Confluent hypergeometric limit functions.
double
conf_hyperg_lim(double /*c*/, double /*x*/)
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
hyperg(double /*a*/, double /*b*/, double /*c*/, double /*x*/)
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
riemann_zeta(double s)
{
  return boost::math::zeta(s);
}

/// Hurwitz zeta functions.
double
hurwitz_zeta(double /*s*/, double /*q*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Dirichlet eta function.
double dirichlet_eta(double /*x*/)
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
sph_legendre(unsigned int l, unsigned int m, double theta)
{
  return boost::math::spherical_harmonic_r(l, m, theta, 0.0);
}

/// Spherical Neumann functions.
double
sph_neumann(unsigned int n, double x)
{
  return boost::math::sph_neumann(n, x);
}

/// Log gamma function.
double
lgamma(double a)
{
  if (a <= 0.0 && a == std::nearbyint(a))
    return std::numeric_limits<double>::infinity();
  else
    return boost::math::lgamma(a);
}

/// Gamma function.
double
tgamma(double a)
{
  return boost::math::tgamma(a);
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
gamma_q(double a, double x)
{
  return boost::math::gamma_q(a, x);
}

/// Inverse normalized upper incomplete gamma functions.
double
gamma_q_inv(double a, double q)
{
  return boost::math::gamma_q_inv(a, q);
}

/// Inverse parameter normalized upper incomplete gamma functions.
double
gamma_q_inva(double a, double q)
{
  return boost::math::gamma_q_inva(a, q);
}

/// Normalized lower incomplete gamma functions.
double
gamma_p(double a, double x)
{
  return boost::math::gamma_p(a, x);
}

/// Inverse normalized lower incomplete gamma functions.
double
gamma_p_inv(double a, double q)
{
  return boost::math::gamma_p_inv(a, q);
}

/// Inverse parameter normalized lower incomplete gamma functions.
double
gamma_p_inva(double a, double q)
{
  return boost::math::gamma_p_inva(a, q);
}

/// Incomplete beta functions.
double
ibeta(double a, double b, double x)
{
  return boost::math::ibeta(a, b, x);
}

/// Incomplete beta functions.
double
ibetac(double a, double b, double x)
{
  return boost::math::ibetac(a, b, x);
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
dilog(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Digamma or psi function.
double
digamma(double x)
{
  return boost::math::digamma(x);
}

/// Polygamma functions.
double
polygamma(unsigned int m, double x)
{
  return boost::math::polygamma(m, x);
}

/// Sine integral.
double
sinint(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cosine integral.
double
cosint(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hyperbolic sine integral.
double
sinhint(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hyperbolic cosine integral.
double
coshint(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Gegenbauer polynomials.
double
gegenbauer(unsigned int /*n*/, double /*lambda*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hydrogen wave functions.
double
hydrogen(int /*n*/, double /*l*/, double /*Z*/, double /*r*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Dawson integral.
double
dawson(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobian elliptic integrals sn.
double
jacobi_sn(double k, double u)
{
  return boost::math::jacobi_sn(k, u);
}

double
jacobi_sc(double k, double u)
{
  return boost::math::jacobi_sc(k, u);
}

double
jacobi_sd(double k, double u)
{
  return boost::math::jacobi_sd(k, u);
}

/// Jacobian elliptic integrals cn.
double
jacobi_cn(double k, double u)
{
  return boost::math::jacobi_cn(k, u);
}

double
jacobi_cd(double k, double u)
{
  return boost::math::jacobi_cd(k, u);
}

double
jacobi_cs(double k, double u)
{
  return boost::math::jacobi_cs(k, u);
}

/// Jacobian elliptic integrals dn.
double
jacobi_dn(double k, double u)
{
  return boost::math::jacobi_dn(k, u);
}

double
jacobi_ds(double k, double u)
{
  return boost::math::jacobi_ds(k, u);
}

double
jacobi_dc(double k, double u)
{
  return boost::math::jacobi_dc(k, u);
}

double
jacobi_ns(double k, double u)
{
  return boost::math::jacobi_ns(k, u);
}

double
jacobi_nc(double k, double u)
{
  return boost::math::jacobi_nc(k, u);
}

double
jacobi_nd(double k, double u)
{
  return boost::math::jacobi_nd(k, u);
}

/// Fresnel cosine integral.
double
fresnel_c(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Fresnel sine integral.
double
fresnel_s(double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Sinus cardinal function.
double
sinc(double x)
{
  return boost::math::sinc_pi(x);
}

/// Reperiodized sinus cardinal function.
double
sinc_pi(double x)
{
  using namespace boost::math::constants;
  return boost::math::sinc_pi(pi<double>() * x);
}

/// Hyperbolic sinus cardinal function.
double
sinhc(double x)
{
  return std::sinh(x) / x;
}

/// Reperiodized hyperbolic sinus cardinal function.
double
sinhc_pi(double x)
{
  using namespace boost::math::constants;
  auto arg = pi<double>() * x;
  return std::sinh(arg) / arg;
}

/// Log rising factorials.
double
lrising_factorial(double a, double x)
{
  return std::log(boost::math::rising_factorial(a, x));
}

/// Log falling factorials.
double
lfalling_factorial(double a, double x)
{
  auto ff = boost::math::falling_factorial(a, x);
  if (ff == 0)
    return std::numeric_limits<double>::infinity();
  else
    return std::log(std::abs(ff));
}

/// Rising factorials.
double
rising_factorial(double a, double x)
{
  return boost::math::rising_factorial(a, x);
}

/// Falling factorials.
double
falling_factorial(double a, double x)
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
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Chebyshev polynomials of the second kind.
double
chebyshev_u(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Chebyshev polynomials of the third kind.
double
chebyshev_v(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Chebyshev polynomials of the fourth kind.
double
chebyshev_w(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobi polynomials.
double
jacobi(unsigned int /*n*/, double /*alpha*/, double /*beta*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Binomial coefficients.
double
binomial(unsigned int /*n*/, unsigned int /*k*/)
{
  return std::numeric_limits<double>::quiet_NaN();//boost::math::binomial_coefficient(n, k);
}

/// Log binomial coefficients.
double
lbinomial(unsigned int /*n*/, unsigned int /*k*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Taylor coefficients.
double
taylorcoeff(unsigned int /*n*/, double /*x*/)
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
zernike(unsigned int /*n*/, unsigned int /*m*/, double /*rho*/, double /*phi*/)
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

/// Clausen Cl function of order 2.
double
clausen_cl(unsigned int /*m*/, double /*w*/)
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

/// Bernoulli numbers.
double
bernoulli(unsigned int n)
{
  if (n == 1)
    return -0.5;
  else if ((n & 1) == 1)
    return 0.0;
  else
    return boost::math::bernoulli_b2n<double>(n / 2);
}

/// Reperiodized sine.
double
sin_pi(double x)
{
  return boost::math::sin_pi(x);
}

/// Reperiodized cosine.
double
cos_pi(double x)
{
  return boost::math::cos_pi(x);
}

/// Fermi-Dirac integrals.
double
fermi_dirac(double /*s*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Bose-Einstein integrals.
double
bose_einstein(double /*s*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Debye functions.
double
debye(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Polylogarithms.
double
polylog(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complex polylogarithms.
std::complex<double>
polylog(unsigned int /*n*/, std::complex<double> /*x*/)
{
  const auto NaN = std::numeric_limits<double>::quiet_NaN();
  return std::complex<double>(NaN, NaN);
}

/// Reciprocal gamma.
double
gamma_reciprocal(unsigned int /*a*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

} // namespace beast

