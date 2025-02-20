#ifndef WRAP_BURKHARDT_H
#define WRAP_BURKHARDT_H 1

#include <vector>
#include <complex>

#include <emsr/quadrature_point.h>

namespace burkhardt
{

/// Airy Ai function.
double airy_ai(double x);

/// Airy Bi function.
double airy_bi(double x);

/// Associated Laguerre polynomials.
double assoc_laguerre(unsigned int n, unsigned int m, double x);

/// Associated Legendre functions.
double assoc_legendre(unsigned int l, unsigned int m, double x);

/// Associated Legendre functions of the second kind.
double assoc_legendre_q(unsigned int l, unsigned int m, double x);

/// Beta functions.
double beta(double x, double y);

/// Complementary beta functions.
double betac(double x, double y);

/// Complete elliptic integrals of the first kind.
double comp_ellint_1(double k);

/// Complete elliptic integrals of the second kind.
double comp_ellint_2(double k);

/// Complete elliptic integrals of the third kind.
double comp_ellint_3(double k, double nu);

/// Complete Legendre elliptic D integrals.
double comp_ellint_d(double k);

/// Confluent hypergeometric functions.
double conf_hyperg(double a, double c, double x);

/// Tricomi confluent hypergeometric functions.
double tricomi_u(double a, double c, double x);

/// Confluent hypergeometric limit functions.
double conf_hyperg_lim(double c, double x);

/// Regular modified cylindrical Bessel functions.
double cyl_bessel_i(double nu, double x);

/// Cylindrical Bessel functions (of the first kind).
double cyl_bessel_j(double nu, double x);

/// Irregular modified cylindrical Bessel functions.
double cyl_bessel_k(double nu, double x);

/// Cylindrical Neumann functions.
double cyl_neumann(double nu, double x);

/// Elliptic integrals of the first kind.
double ellint_1(double k, double phi);

/// Elliptic integrals of the second kind.
double ellint_2(double k, double phi);

/// Elliptic integrals of the third kind.
double ellint_3(double k, double nu, double phi);

/// Legendre elliptic D integrals.
double ellint_d(double k, double phi);

/// Carlson elliptic integrals R_C.
double ellint_rc(double x, double y);

/// Carlson elliptic integrals R_D.
double ellint_rd(double x, double y, double z);

/// Carlson elliptic integrals R_F.
double ellint_rf(double x, double y, double z);

/// Carlson elliptic integrals R_G.
double ellint_rg(double x, double y, double z);

/// Carlson elliptic integrals R_J.
double ellint_rj(double x, double y, double z, double p);

/// Exponential integral Ei.
double expint(double x);

/// Exponential integrals E_n.
double expint(unsigned int n, double x);

/// Hermite polynomials.
double hermite(unsigned int n, double x);

/// Hypergeometric functions.
double hyperg(double a, double b, double c, double x);

/// Laguerre polynomials.
double laguerre(unsigned int n, double x);

/// Legendre polynomials.
double legendre_p(unsigned int l, double x);

/// Legendre functions of the second kind.
double legendre_q(unsigned int l, double x);

/// Riemann zeta function.
double riemann_zeta(double s);

/// Hurwitz zeta functions.
double hurwitz_zeta(double s, double q);

/// Dirichlet eta function.
double dirichlet_eta(double x);

/// Spherical Bessel functions.
double sph_bessel(unsigned int n, double x);

/// Spherical Legendre functions.
double sph_legendre(unsigned int l, unsigned int m, double theta);

/// Spherical Neumann functions.
double sph_neumann(unsigned int n, double x);

/// Normalized upper incomplete gamma functions.
double gamma_q(double a, double x);

/// Inverse normalized upper incomplete gamma functions.
double gamma_q_inv(double a, double q);

/// Inverse parameter normalized upper incomplete gamma functions.
double gamma_q_inva(double x, double q);

/// Normalized lower incomplete gamma functions.
double gamma_p(double a, double x);

/// Inverse normalized lower incomplete gamma functions.
double gamma_p_inv(double a, double p);

/// Inverse parameter normalized lower incomplete gamma functions.
double gamma_p_inva(double x, double p);

/// Non-normalized (upper) incomplete gamma functions.
double tgamma(double a, double x);

/// Non-normalized lower incomplete gamma functions.
double tgamma_lower(double a, double x);

/// Incomplete beta functions.
double ibeta(double a, double b, double x);

/// Complementary incomplete beta functions.
double ibetac(double a, double b, double x);

/// Dilogarithm function.
double dilog(double x);

/// Digamma or psi function.
double digamma(double x);

/// Polygamma functions.
double polygamma(unsigned int m, double x);

/// Sine integral.
double sinint(double x);

/// Cosine integral.
double cosint(double x);

/// Hyperbolic sine integral.
double sinhint(double x);

/// Hyperbolic cosine integral.
double coshint(double x);

/// Gegenbauer polynomials.
double gegenbauer(unsigned int n, double lambda, double x);

/// Hydrogen wave functions.
double hydrogen(int n, double l, double Z, double r);

/// Dawson integral.
double dawson(double x);

/// Jacobian elliptic integrals sn.
double jacobi_sn(double k, double u);

/// Jacobian elliptic integrals cn.
double jacobi_cn(double k, double u);

/// Jacobian elliptic integrals dn.
double jacobi_dn(double k, double u);

/// Fresnel cosine integral.
double fresnel_c(double x);

/// Fresnel sine integral.
double fresnel_s(double x);

/// Sinus cardinal function.
double sinc(double x);

/// Reperiodized sinus cardinal function.
double sinc_pi(double x);

/// Hyperbolic sinus cardinal function.
double sinhc(double x);

/// Reperiodized hyperbolic sinus cardinal function.
double sinhc_pi(double x);

/// Log rising factorials.
double lrising_factorial(double a, double x);

/// Log falling factorials.
double lfalling_factorial(double a, double x);

/// Rising factorials.
double rising_factorial(double a, double x);

/// Falling factorials.
double falling_factorial(double a, double x);

/// Log factorial.
double lfactorial(unsigned int n);

/// Factorial.
double factorial(unsigned int n);

/// Log double factorial.
double ldouble_factorial(int n);

/// Double factorial.
double double_factorial(int n);

/// Regular modified spherical bessel functions.
double sph_bessel_i(unsigned int n, double x);

/// Irregular modified spherical bessel functions.
double sph_bessel_k(unsigned int n, double x);

/// Chebyshev polynomials of the first kind.
double chebyshev_t(unsigned int n, double x);

/// Chebyshev polynomials of the second kind.
double chebyshev_u(unsigned int n, double x);

/// Chebyshev polynomials of the third kind.
double chebyshev_v(unsigned int n, double x);

/// Chebyshev polynomials of the fourth kind.
double chebyshev_w(unsigned int n, double x);

/// Binomial coefficients.
double binomial(unsigned int n, unsigned int k);

/// Log binomial coefficients.
double lbinomial(unsigned int n, unsigned int k);

/// Jacobi polynomials.
double jacobi(unsigned int n, double alpha, double beta, double x);

/// Taylor coefficients.
double taylorcoeff(unsigned int n, double x);

/// Radial polynomials
double radpoly(unsigned int n, unsigned int m, double rho);

/// Zernike polynomials
double zernike(unsigned int n, int m, double rho, double phi);

/// Cylindrical Hankel functions of the first kind.
//std::complex<double> cyl_hankel_1(double nu, double x);

/// Cylindrical Hankel functions of the second kind.
//std::complex<double> cyl_hankel_2(double nu, double x);

/// Spherical Hankel functions of the first kind.
//std::complex<double> sph_hankel_1(unsigned int n, double x);

/// Spherical Hankel functions of the second kind.
//std::complex<double> sph_hankel_2(unsigned int n, double x);

/// Heuman lambda functions.
double heuman_lambda(double k, double phi);

/// Jacobi zeta functions.
double jacobi_zeta(double k, double phi);

/// Inverse error function.
double erf_inv(double p);

/// Inverse complementary error function.
double erfc_inv(double q);

/// Inverse incomplete beta function.
double ibeta_inv(double a, double b, double p);

/// Inverse complementary incomplete beta function.
double ibetac_inv(double a, double b, double q);

/// Inverse parameter incomplete beta function.
double ibeta_inva(double b, double x, double p);

/// Inverse parameter complementary incomplete beta function.
double ibetac_inva(double b, double x, double q);

/// Inverse parameter incomplete beta function.
double ibeta_invb(double a, double x, double p);

/// Inverse parameter complementary incomplete beta function.
double ibetac_invb(double a, double x, double q);

/// Spherical harmonic functions.
//std::complex<double> sph_harmonic(unsigned int l, int m, double theta, double phi);

/// Owen's T function.
double owens_t(double h, double a);

/// Clausen Cl function of order 2.
double clausen_cl(unsigned int m, double w);

/// Struve H function.
double struve_h(double nu, double x);

/// Struve L function.
double struve_l(double nu, double x);

/// Bernoulli numbers.
double bernoulli(unsigned int n);

/// Bernoulli polynomials.
double bernoulli(unsigned int n, double x);

/// Reperiodized sine.
double sin_pi(double x);

/// Reperiodized cosine.
double cos_pi(double x);

/// Fermi-Dirac integrals.
double fermi_dirac(double s, double x);

/// Bose-Einstein integrals.
double bose_einstein(double s, double x);

/// Debye functions.
double debye(unsigned int n, double x);

/// Polylogarithms.
double polylog(unsigned int n, double x);

/// Complex polylogarithms.
//std::complex<double> polylog(unsigned int n, std::complex<double> x);

/// Reciprocal gamma.
double gamma_reciprocal(unsigned int a);

/// Euler numbers.
double euler(unsigned int n);

/// Euler polynomials.
double euler(unsigned int n, double x);

/// Eulerian numbers of the first kind.
double eulerian_1(unsigned int n, unsigned int m);

/// Eulerian numbers of the second kind.
double eulerian_2(unsigned int n, unsigned int m);

/// Stirling numbers of the first kind.
double stirling_1(unsigned int n, unsigned int m);

/// Stirling numbers of the second kind.
double stirling_2(unsigned int n, unsigned int m);

/// Wright omega function.
std::complex<double> wright_omega(std::complex<double> z);

/// Lambert W function.
double lambert_w(double x);

/// Gauss-Chebyshev T rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_chebyshev_t_rule(std::size_t n);

/// Gauss-Chebyshev U rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_chebyshev_u_rule(std::size_t n);

/// Gauss-Chebyshev V rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_chebyshev_v_rule(std::size_t n);

/// Gauss-Radau rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_legendre_radau_rule(std::size_t n);

/// Gauss-Lobatto rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_legendre_lobatto_rule(std::size_t n);

/// Gauss-Laguerre rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_laguerre_rule(std::size_t n);

/// Gauss-Associative Laguerre rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_laguerre_rule(std::size_t n, double alpha);

/// Gauss-Jacobi rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_jacobi_rule(std::size_t n, double alpha, double beta);

/// Gauss-Gegenbauer rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_gegenbauer_rule(std::size_t n, double lambda);

/// Gauss-Hermite rule.
std::vector<emsr::QuadraturePoint<double>>
gauss_hermite_rule(std::size_t n);

/// Clenshaw-Curtis rule.
std::vector<emsr::QuadraturePoint<double>>
clenshaw_curtis_rule(std::size_t n);

/// Fejer 1 rule.
std::vector<emsr::QuadraturePoint<double>>
fejer_1_rule(std::size_t n);

/// Fejer 2 rule.
std::vector<emsr::QuadraturePoint<double>>
fejer_2_rule(std::size_t n);

// Bell numbers.
std::vector<unsigned int>
bell(unsigned int n);

/// Meixner polynomials.
double
meixner(int n, double beta, double c, double x);

/// Krawtchouk polynomials.
double
krawtchouk(int n, double p, double x, int m);

/// Charlier polynomials.
double
charlier(int n, double a, double x);

/// Gudermannian function.
double gd(double x);

/// Inverse Gudermannian function.
double agd(double g);
/*
/// Gamma probability density function.
gamma_pdf(double alpha, double beta, double x);

/// Gamma distribution.
gamma_p(double alpha, double beta, double x);

/// Normal probability density function.
normal_pdf(double mu, double sigma, double x);
normal_pdouble mu, double sigma, double x);

/// Lognormal probability density function.
lognormal_pdf(double mu, double sigma, double x);
lognormal_p(double mu, double sigma, double x);

/// Exponential probability density function.
exponential_pdf(double lambda, double x);
exponential_p(double lambda, double x);

/// Weibull probability density function.
weibull_pdf(double a, double b, double x);
weibull_p(double a, double b, double x);

/// Student-T probability density function.
student_t_pdf(double t, unsigned int nu);
student_t_p(double t, unsigned int nu);

/// Fisher-F probability density function.
fisher_f_pdf(double F, unsigned int nu1, unsigned int nu2);
fisher_f_p(double F, unsigned int nu1, unsigned int nu2);

/// Binomial probability density function.
binomial_pdf(double p, unsigned int n, unsigned int k);
binomial_p(double p, unsigned int n, unsigned int k);

/// Logistic probability density function.
logistic_pdf(double a, double b, double x);
logistic_p(double a, double b, double x);
*/
} // namespace burkhardt

#endif // WRAP_BURKHARDT_H

