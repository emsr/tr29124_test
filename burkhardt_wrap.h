
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

/// Confluent hypergeometric limit functions.
double hyperg_0F1(double c, double x);

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

/// Exponential integral E_n.
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
double riemann_zeta(double x);

/// Hurwitz zeta functions.
double hurwitz_zeta(double s, double q);

/// Dirichlet eta function.
double dirichlet_eta(double x);

/// Spherical Bessel functions.
double bessel_jl(unsigned int n, double x);

/// Spherical Legendre functions.
double legendre_sphPlm(unsigned int l, unsigned int m, double theta);

/// Spherical Neumann functions.
double bessel_yl(unsigned int n, double x);

/// Non-normalized lower incomplete gamma functions. (See Boost tgamma_lower(a, x)).
double gamma_l(double a, double x);

/// Normalized upper incomplete gamma functions.
double qgamma(double a, double x);

/// Inverse normalized upper incomplete gamma functions.
double qgamma_inv(double a, double q);

/// Inverse parameter normalized upper incomplete gamma functions.
double qgamma_inva(double x, double q);

/// Normalized lower incomplete gamma functions.
double pgamma(double a, double x);

/// Inverse normalized lower incomplete gamma functions.
double pgamma_inv(double a, double p);

/// Inverse parameter normalized lower incomplete gamma functions.
double pgamma_inva(double x, double p);

/// Non-normalized (upper) incomplete gamma functions. (See Boost tgamma(a, x)).
double gamma_u(double a, double x);

/// Incomplete beta functions.
double ibeta(double a, double b, double x);

/// Complementary incomplete beta functions.
double ibetac(double a, double b, double x);

/// Dilogarithm function.
double dilog(double x);

/// Digamma or psi function.
double psi(double x);

/// Polygamma functions.
double polygamma(int n, double x);

/// Sine integral.
double Si(double x);

/// Cosine integral.
double Ci(double x);

/// Hyperbolic sine integral.
double Shi(double x);

/// Hyperbolic cosine integral.
double Chi(double x);

/// Gegenbauer polynomials.
double gegenpoly_n(unsigned int n, double lambda, double x);

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
double sinc_pi(double x);

/// Normalized sinus cardinal function.
double sinc(double x);

/// Hyperbolic sinus cardinal function.
double sinhc_pi(double x);

/// Normalized hyperbolic sinus cardinal function.
double sinhc(double x);

/// Log upper Pochhammer symbol.
double lpochhammer_u(double a, double x);

/// Log lower Pochhammer symbol.
double lpochhammer_l(double a, double x);

/// Upper Pochhammer symbol.
double pochhammer_u(double a, double x);

/// Lower Pochhammer symbol.
double pochhammer_l(double a, double x);

/// Log factorial.
double lfactorial(unsigned int n);

/// Factorial.
double factorial(unsigned int n);

/// Log double factorial.
double ldouble_factorial(int n);

/// Double factorial.
double double_factorial(int n);

/// Regular modified spherical bessel functions.
double bessel_il(unsigned int n, double x);

/// Irregular modified spherical bessel functions.
double bessel_kl(unsigned int n, double x);

/// Chebyshev polynomials of the first kind.
double chebyshev_t(unsigned int n, double x);

/// Chebyshev polynomials of the second kind.
double chebyshev_u(unsigned int n, double x);

/// Chebyshev polynomials of the third kind.
double chebyshev_v(unsigned int n, double x);

/// Chebyshev polynomials of the fourth kind.
double chebyshev_w(unsigned int n, double x);

/// Binomial coefficients.
double choose(unsigned int n, unsigned int k);

/// Log binomial coefficients.
double lnchoose(unsigned int n, unsigned int k);

/// Jacobi polynomials.
double jacobi(unsigned int n, double alpha, double beta, double x);

/// Taylor coefficients.
double taylorcoeff(unsigned int n, double x);

/// Radial polynomials
double radpoly(unsigned int n, unsigned int m, double rho);

/// Zernike polynomials
double zernike(unsigned int n, int m, double rho, double phi);
/*
/// Cylindrical Hankel functions of the first kind.
std::complex<double> cyl_hankel_1(double nu, double x);

/// Cylindrical Hankel functions of the second kind.
std::complex<double> cyl_hankel_2(double nu, double x);

/// Spherical Hankel functions of the first kind.
std::complex<double> sph_hankel_1(unsigned int n, double x);

/// Spherical Hankel functions of the second kind.
std::complex<double> sph_hankel_2(unsigned int n, double x);
*/
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
/*
/// Spherical harmonic functions.
std::complex<double> sph_harmonic(unsigned int l, int m, double theta, double phi);
*/
/// Owen's T function.
double owens_t(double h, double a);

/// Clausen function of order 2.
double clausen_c(unsigned int m, double w);

} // namespace burkhardt

