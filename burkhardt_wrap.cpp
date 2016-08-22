
namespace burkhardt
{

/// Airy Ai function.
double
airy_ai(double x)
{
  double ai, bi, ad, bd;
  airya_( x, ai, bi, ad, bd )
  return ai;
}

/// Airy Bi function.
double
airy_bi(double x)
{
  double ai, bi, ad, bd;
  airya_( x, ai, bi, ad, bd )
  return bi;
}

/// Associated Laguerre polynomials.
double
assoc_laguerre(unsigned int n, unsigned int m, double x)
{
  return 0.0;
}

/// Associated Legendre functions.
double
assoc_legendre(unsigned int l, unsigned int m, double x)
{
  //int mm = l * m;
  //double cpm[mm], cpd[mm];
  //double x, y;
  //clpmn_(mm, l, m, x, y, cpm, cpd);
  return 0.0;
}

/// Associated Legendre functions of the second kind.
double
assoc_legendre_q(unsigned int l, unsigned int m, double x)
{
  //int mm = l * m;
  //double cpm[mm], cpd[mm];
  //double x, y;
  //clqmn_(mm, m, n, x, y, cqm, cqd);
  return 0.0;
}

/// Beta functions.
double
beta(double x, double y)
{
  return 0.0;
}

/// Complementary beta functions.
double
betac(double x, double y)
{
  return 0.0;
}

/// Complete elliptic integrals of the first kind.
double
comp_ellint_1(double k)
{
  return 0.0;
}

/// Complete elliptic integrals of the second kind.
double
comp_ellint_2(double k)
{
  return 0.0;
}

/// Complete elliptic integrals of the third kind.
double
comp_ellint_3(double k, double nu)
{
  return 0.0;
}

/// Complete Legendre elliptic D integrals.
double
comp_ellint_d(double k)
{
  return 0.0;
}

/// Confluent hypergeometric functions.
double
conf_hyperg(double a, double c, double x)
{
  return 0.0;
}

/// Confluent hypergeometric limit functions.
double
hyperg_0F1(double c, double x)
{
  return 0.0;
}

/// Regular modified cylindrical Bessel functions.
double
cyl_bessel_i(double nu, double x)
{
  return 0.0;
}

/// Cylindrical Bessel functions (of the first kind).
double
cyl_bessel_j(double nu, double x)
{
  return 0.0;
}

/// Irregular modified cylindrical Bessel functions.
double
cyl_bessel_k(double nu, double x)
{
  return 0.0;
}

/// Cylindrical Neumann functions.
double
cyl_neumann(double nu, double x)
{
  return 0.0;
}

/// Elliptic integrals of the first kind.
double
ellint_1(double k, double phi)
{
  return 0.0;
}

/// Elliptic integrals of the second kind.
double
ellint_2(double k, double phi)
{
  return 0.0;
}

/// Elliptic integrals of the third kind.
double
ellint_3(double k, double nu, double phi)
{
  return 0.0;
}

/// Legendre elliptic D integrals.
double
ellint_d(double k, double phi)
{
  return 0.0;
}

/// Carlson elliptic integrals R_C.
double
ellint_rc(double x, double y)
{
  return 0.0;
}

/// Carlson elliptic integrals R_D.
double
ellint_rd(double x, double y, double z)
{
  return 0.0;
}

/// Carlson elliptic integrals R_F.
double
ellint_rf(double x, double y, double z)
{
  return 0.0;
}

/// Carlson elliptic integrals R_G.
double
ellint_rg(double x, double y, double z)
{
  return 0.0;
}

/// Carlson elliptic integrals R_J.
double
ellint_rj(double x, double y, double z, double p)
{
  return 0.0;
}

/// Exponential integral Ei.
double
expint(double x)
{
  return 0.0;
}

/// Exponential integral E_n.
double
expint(unsigned int n, double x)
{
  return 0.0;
}

/// Hermite polynomials.
double
hermite(unsigned int n, double x)
{
  return 0.0;
}

/// Hypergeometric functions.
double
hyperg(double a, double b, double c, double x)
{
  return 0.0;
}

/// Laguerre polynomials.
double
laguerre(unsigned int n, double x)
{
  return 0.0;
}

/// Legendre polynomials.
double
legendre_p(unsigned int l, double x)
{
  return 0.0;
}

/// Legendre functions of the second kind.
double
legendre_q(unsigned int l, double x)
{
  //clqmn_(mm, m, n, x, y, cqm, cqd);
  return 0.0;
}

/// Riemann zeta function.
double
riemann_zeta(double x)
{
  return 0.0;
}

/// Hurwitz zeta functions.
double
hurwitz_zeta(double s, double q)
{
  return 0.0;
}

/// Dirichlet eta function.
double
dirichlet_eta(double x)
{
  return 0.0;
}

/// Spherical Bessel functions.
double
bessel_jl(unsigned int n, double x)
{
  return 0.0;
}

/// Spherical Legendre functions.
double
legendre_sphPlm(unsigned int l, unsigned int m, double theta)
{
  return 0.0;
}

/// Spherical Neumann functions.
double
bessel_yl(unsigned int n, double x)
{
  return 0.0;
}

/// Non-normalized lower incomplete gamma functions. (See Boost tgamma_lower(a, x)).
double
gamma_l(double a, double x)
{
  return 0.0;
}

/// Normalized upper incomplete gamma functions.
double
qgamma(double a, double x)
{
  return 0.0;
}

/// Inverse normalized upper incomplete gamma functions.
double
qgamma_inv(double a, double q)
{
  return 0.0;
}

/// Inverse parameter normalized upper incomplete gamma functions.
double
qgamma_inva(double x, double q)
{
  return 0.0;
}

/// Normalized lower incomplete gamma functions.
double
pgamma(double a, double x)
{
  return 0.0;
}

/// Inverse normalized lower incomplete gamma functions.
double
pgamma_inv(double a, double p)
{
  return 0.0;
}

/// Inverse parameter normalized lower incomplete gamma functions.
double
pgamma_inva(double x, double p)
{
  return 0.0;
}

/// Non-normalized (upper) incomplete gamma functions. (See Boost tgamma(a, x)).
double
gamma_u(double a, double x)
{
  return 0.0;
}

/// Incomplete beta functions.
double
ibeta(double a, double b, double x)
{
  return 0.0;
}

/// Complementary incomplete beta functions.
double
ibetac(double a, double b, double x)
{
  return 0.0;
}

/// Dilogarithm function.
double
dilog(double x)
{
  return 0.0;
}

/// Digamma or psi function.
double
psi(double x)
{
  return 0.0;
}

/// Polygamma functions.
double
polygamma(int n, double x)
{
  return 0.0;
}

/// Sine integral.
double
Si(double x)
{
  return 0.0;
}

/// Cosine integral.
double
Ci(double x)
{
  return 0.0;
}

/// Hyperbolic sine integral.
double
Shi(double x)
{
  return 0.0;
}

/// Hyperbolic cosine integral.
double
Chi(double x)
{
  return 0.0;
}

/// Gegenbauer polynomials.
double
gegenpoly_n(unsigned int n, double lambda, double x)
{
  return 0.0;
}

/// Hydrogen wave functions.
double
hydrogen(int n, double l, double Z, double r)
{
  return 0.0;
}

/// Dawson integral.
double
dawson(double x)
{
  return 0.0;
}

/// Jacobian elliptic integrals sn.
double
jacobi_sn(double k, double u)
{
  return 0.0;
}

/// Jacobian elliptic integrals cn.
double
jacobi_cn(double k, double u)
{
  return 0.0;
}

/// Jacobian elliptic integrals dn.
double
jacobi_dn(double k, double u)
{
  return 0.0;
}

/// Fresnel cosine integral.
double
fresnel_c(double x)
{
  return 0.0;
}

/// Fresnel sine integral.
double
fresnel_s(double x)
{
  return 0.0;
}

/// Sinus cardinal function.
double
sinc_pi(double x)
{
  return 0.0;
}

/// Normalized sinus cardinal function.
double
sinc(double x)
{
  return 0.0;
}

/// Hyperbolic sinus cardinal function.
double
sinhc_pi(double x)
{
  return 0.0;
}

/// Normalized hyperbolic sinus cardinal function.
double
sinhc(double x)
{
  return 0.0;
}

/// Log upper Pochhammer symbol.
double
lpochhammer_u(double a, double x)
{
  return 0.0;
}

/// Log lower Pochhammer symbol.
double
lpochhammer_l(double a, double x)
{
  return 0.0;
}

/// Upper Pochhammer symbol.
double
pochhammer_u(double a, double x)
{
  return 0.0;
}

/// Lower Pochhammer symbol.
double
pochhammer_l(double a, double x)
{
  return 0.0;
}

/// Log factorial.
double
lfactorial(unsigned int n)
{
  return 0.0;
}

/// Factorial.
double
factorial(unsigned int n)
{
  return 0.0;
}

/// Log double factorial.
double
ldouble_factorial(int n)
{
  return 0.0;
}

/// Double factorial.
double
double_factorial(int n)
{
  return 0.0;
}

/// Regular modified spherical bessel functions.
double
bessel_il(unsigned int n, double x)
{
  return 0.0;
}

/// Irregular modified spherical bessel functions.
double
bessel_kl(unsigned int n, double x)
{
  return 0.0;
}

/// Chebyshev polynomials of the first kind.
double
chebyshev_t(unsigned int n, double x)
{
  return 0.0;
}

/// Chebyshev polynomials of the second kind.
double
chebyshev_u(unsigned int n, double x)
{
  return 0.0;
}

/// Chebyshev polynomials of the third kind.
double
chebyshev_v(unsigned int n, double x)
{
  return 0.0;
}

/// Chebyshev polynomials of the fourth kind.
double
chebyshev_w(unsigned int n, double x)
{
  return 0.0;
}

/// Binomial coefficients.
double
choose(unsigned int n, unsigned int k)
{
  return 0.0;
}

/// Log binomial coefficients.
double
lnchoose(unsigned int n, unsigned int k)
{
  return 0.0;
}

/// Jacobi polynomials.
double
jacobi(unsigned int n, double alpha, double beta, double x)
{
  return 0.0;
}

/// Taylor coefficients.
double
taylorcoeff(unsigned int n, double x)
{
  return 0.0;
}

/// Radial polynomials
double
radpoly(unsigned int n, unsigned int m, double rho)
{
  return 0.0;
}

/// Zernike polynomials
double
zernike(unsigned int n, int m, double rho, double phi)
{
  return 0.0;
}

/// Cylindrical Hankel functions of the first kind.
std::complex<double> cyl_hankel_1(double nu, double x)
{
  return 0.0;
}

/// Cylindrical Hankel functions of the second kind.
std::complex<double> cyl_hankel_2(double nu, double x)
{
  return 0.0;
}

/// Spherical Hankel functions of the first kind.
std::complex<double> sph_hankel_1(unsigned int n, double x)
{
  return 0.0;
}

/// Spherical Hankel functions of the second kind.
std::complex<double> sph_hankel_2(unsigned int n, double x)
{
  return 0.0;
}

/// Heuman lambda functions.
double
heuman_lambda(double k, double phi)
{
  return 0.0;
}

/// Jacobi zeta functions.
double
jacobi_zeta(double k, double phi)
{
  return 0.0;
}

/// Inverse error function.
double
erf_inv(double p)
{
  return 0.0;
}

/// Inverse complementary error function.
double
erfc_inv(double q)
{
  return 0.0;
}

/// Inverse incomplete beta function.
double
ibeta_inv(double a, double b, double p)
{
  return 0.0;
}

/// Inverse complementary incomplete beta function.
double
ibetac_inv(double a, double b, double q)
{
  return 0.0;
}

/// Inverse parameter incomplete beta function.
double
ibeta_inva(double b, double x, double p)
{
  return 0.0;
}

/// Inverse parameter complementary incomplete beta function.
double
ibetac_inva(double b, double x, double q)
{
  return 0.0;
}

/// Inverse parameter incomplete beta function.
double
ibeta_invb(double a, double x, double p)
{
  return 0.0;
}

/// Inverse parameter complementary incomplete beta function.
double
ibetac_invb(double a, double x, double q)
{
  return 0.0;
}

/// Spherical harmonic functions.
std::complex<double> sph_harmonic(unsigned int l, int m, double theta, double phi)
{
  return 0.0;
}

/// Owen's T function.
double
owens_t(double h, double a)
{
  return 0.0;
}

/// Clausen function of order 2.
double
clausen_c(unsigned int m, double w)
{
  return 0.0;
}

/// Struve H function.
double
struve_h(double nu, double x)
{
  double sh{};
  stvhv_(nu, x, sh)
  return sh;
}

/// Struve L function.
double
struve_l(double nu, double x)
{
  double sl{};
  stvlv_(nu, x, sl)
  return sl;
}

} // namespace burkhardt

