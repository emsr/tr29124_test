
#include <limits>

#include "burkhardt/special_functions.h"
#include "wrap_burkhardt.h"


namespace burkhardt
{

/// Airy Ai function.
double
airy_ai(double x)
{
  double ai{}, bi{}, ad{}, bd{};
  airya_(&x, &ai, &bi, &ad, &bd);
  return ai;
}

/// Airy Bi function.
double
airy_bi(double x)
{
  double ai{}, bi{}, ad{}, bd{};
  airya_(&x, &ai, &bi, &ad, &bd);
  return bi;
}

/// Associated Laguerre polynomials.
double
assoc_laguerre(unsigned int n, unsigned int m, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Associated Legendre functions.
double
assoc_legendre(unsigned int l, unsigned int m, double x)
{
  //int mm = l * m;
  //double cpm[mm], cpd[mm];
  //double x{}, y{};
  //clpmn_(&mm, &l, &m, &x, &y, &cpm, &cpd);
  return std::numeric_limits<double>::quiet_NaN();
}

/// Associated Legendre functions of the second kind.
double
assoc_legendre_q(unsigned int l, unsigned int m, double x)
{
  //int mm = l * m;
  //double cpm[mm], cpd[mm];
  //double x{}, y{};
  //clqmn_(&mm, &m, &n, &x, &y, &cqm, &cqd);
  return std::numeric_limits<double>::quiet_NaN();
}

/// Beta functions.
double
beta(double x, double y)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complementary beta functions.
double
betac(double x, double y)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complete elliptic integrals of the first kind.
double
comp_ellint_1(double k)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complete elliptic integrals of the second kind.
double
comp_ellint_2(double k)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complete elliptic integrals of the third kind.
double
comp_ellint_3(double k, double nu)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complete Legendre elliptic D integrals.
double
comp_ellint_d(double k)
{
  return std::numeric_limits<double>::quiet_NaN();
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
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cylindrical Bessel functions (of the first kind).
double
cyl_bessel_j(double nu, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Irregular modified cylindrical Bessel functions.
double
cyl_bessel_k(double nu, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cylindrical Neumann functions.
double
cyl_neumann(double nu, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Elliptic integrals of the first kind.
double
ellint_1(double k, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Elliptic integrals of the second kind.
double
ellint_2(double k, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Elliptic integrals of the third kind.
double
ellint_3(double k, double nu, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Legendre elliptic D integrals.
double
ellint_d(double k, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_C.
double
ellint_rc(double x, double y)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_D.
double
ellint_rd(double x, double y, double z)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_F.
double
ellint_rf(double x, double y, double z)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_G.
double
ellint_rg(double x, double y, double z)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_J.
double
ellint_rj(double x, double y, double z, double p)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Exponential integral Ei.
double
expint(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Exponential integral E_n.
double
expint(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hermite polynomials.
double
hermite(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
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
  return std::numeric_limits<double>::quiet_NaN();
}

/// Legendre polynomials.
double
legendre_p(unsigned int l, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Legendre functions of the second kind.
double
legendre_q(unsigned int l, double x)
{
  //clqmn_(mm, m, n, x, y, cqm, cqd);
  return std::numeric_limits<double>::quiet_NaN();
}

/// Riemann zeta function.
double
riemann_zeta(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hurwitz zeta functions.
double
hurwitz_zeta(double s, double q)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Dirichlet eta function.
double
dirichlet_eta(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Bessel functions.
double
bessel_jl(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Legendre functions.
double
legendre_sphPlm(unsigned int l, unsigned int m, double theta)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Neumann functions.
double
bessel_yl(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Non-normalized lower incomplete gamma functions. (See Boost tgamma_lower(a, x)).
double
gamma_l(double a, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Normalized upper incomplete gamma functions.
double
qgamma(double a, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse normalized upper incomplete gamma functions.
double
qgamma_inv(double a, double q)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse parameter normalized upper incomplete gamma functions.
double
qgamma_inva(double x, double q)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Normalized lower incomplete gamma functions.
double
pgamma(double a, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse normalized lower incomplete gamma functions.
double
pgamma_inv(double a, double p)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse parameter normalized lower incomplete gamma functions.
double
pgamma_inva(double x, double p)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Non-normalized (upper) incomplete gamma functions. (See Boost tgamma(a, x)).
double
gamma_u(double a, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Incomplete beta functions.
double
ibeta(double a, double b, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Complementary incomplete beta functions.
double
ibetac(double a, double b, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
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
  return std::numeric_limits<double>::quiet_NaN();
}

/// Polygamma functions.
double
polygamma(int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
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
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobian elliptic integrals cn.
double
jacobi_cn(double k, double u)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobian elliptic integrals dn.
double
jacobi_dn(double k, double u)
{
  return std::numeric_limits<double>::quiet_NaN();
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
  return std::numeric_limits<double>::quiet_NaN();
}

/// Normalized sinus cardinal function.
double
sinc(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Hyperbolic sinus cardinal function.
double
sinhc_pi(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Normalized hyperbolic sinus cardinal function.
double
sinhc(double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log upper Pochhammer symbol.
double
lpochhammer_u(double a, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log lower Pochhammer symbol.
double
lpochhammer_l(double a, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Upper Pochhammer symbol.
double
pochhammer_u(double a, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Lower Pochhammer symbol.
double
pochhammer_l(double a, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log factorial.
double
lfactorial(unsigned int n)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Factorial.
double
factorial(unsigned int n)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log double factorial.
double
ldouble_factorial(int n)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Double factorial.
double
double_factorial(int n)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Regular modified spherical bessel functions.
double
bessel_il(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Irregular modified spherical bessel functions.
double
bessel_kl(unsigned int n, double x)
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

/// Binomial coefficients.
double
choose(unsigned int n, unsigned int k)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Log binomial coefficients.
double
lnchoose(unsigned int n, unsigned int k)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobi polynomials.
double
jacobi(unsigned int n, double alpha, double beta, double x)
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
zernike(unsigned int n, int m, double rho, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}
/*
/// Cylindrical Hankel functions of the first kind.
std::complex<double> cyl_hankel_1(double nu, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cylindrical Hankel functions of the second kind.
std::complex<double> cyl_hankel_2(double nu, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Hankel functions of the first kind.
std::complex<double> sph_hankel_1(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Hankel functions of the second kind.
std::complex<double> sph_hankel_2(unsigned int n, double x)
{
  return std::numeric_limits<double>::quiet_NaN();
}
*/
/// Heuman lambda functions.
double
heuman_lambda(double k, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobi zeta functions.
double
jacobi_zeta(double k, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse error function.
double
erf_inv(double p)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse complementary error function.
double
erfc_inv(double q)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse incomplete beta function.
double
ibeta_inv(double a, double b, double p)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse complementary incomplete beta function.
double
ibetac_inv(double a, double b, double q)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse parameter incomplete beta function.
double
ibeta_inva(double b, double x, double p)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse parameter complementary incomplete beta function.
double
ibetac_inva(double b, double x, double q)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse parameter incomplete beta function.
double
ibeta_invb(double a, double x, double p)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse parameter complementary incomplete beta function.
double
ibetac_invb(double a, double x, double q)
{
  return std::numeric_limits<double>::quiet_NaN();
}
/*
/// Spherical harmonic functions.
std::complex<double> sph_harmonic(unsigned int l, int m, double theta, double phi)
{
  return std::numeric_limits<double>::quiet_NaN();
}
*/
/// Owen's T function.
double
owens_t(double h, double a)
{
  return std::numeric_limits<double>::quiet_NaN();
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
  double sh{};
  stvhv_(&nu, &x, &sh);
  return sh;
}

/// Struve L function.
double
struve_l(double nu, double x)
{
  double sl{};
  stvlv_(&nu, &x, &sl);
  return sl;
}

} // namespace burkhardt

