
namespace gsl
{

/// Airy Ai function.
double airy_ai(double x);

/// Airy Bi function.
double airy_bi(double x);

/// Associated Laguerre polynomials.
double laguerre_nm(unsigned int n, unsigned int m, double x);

/// Associated Legendre functions.
double legendre_Plm(unsigned int l, unsigned int m, double x);

/// Associated Legendre functions of the second kind.
double legendre_Qlm(unsigned int l, unsigned int m, double x);

/// Beta functions.
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

/// Irregular modified cylindrical Bessel functions.
double bessel_Knu(double nu, double x);

/// Cylindrical Neumann functions.
double bessel_Ynu(double nu, double x);

/// Elliptic integrals of the first kind.
double ellint_F(double k, double phi);

/// Elliptic integrals of the second kind.
double ellint_E(double k, double phi);

/// Elliptic integrals of the third kind.
double ellint_P(double k, double nu, double phi);

/// Carlson elliptic integrals R_C.
double ellint_RC(double x, double y);

/// Carlson elliptic integrals R_D.
double ellint_RD(double x, double y, double z);

/// Carlson elliptic integrals R_F.
double ellint_RF(double x, double y, double z);

/// Carlson elliptic integrals R_J.
double ellint_RJ(double x, double y, double z, double p);

/// Exponential integral Ei.
double expint_Ei(double x);

/// Exponential integral E_1.
double expint_E1(double x);

/// Exponential integrals E_n.
double expint_En(unsigned int n, double x);

/// Hermite polynomials.
double hermite(unsigned int n, double x);

/// Hypergeometric functions.
double hyperg_2F1(double a, double b, double c, double x);

/// Laguerre polynomials.
double laguerre_n(unsigned int n, double x);

/// Legendre polynomials.
double legendre_p(unsigned int l, double x);

/// Legendre polynomials of the second kind.
double legendre_q(unsigned int l, double x);

/// Riemann zeta function.
double zeta(double x);

/// Hurwitz zeta functions.
double hzeta(double s, double q);

/// Spherical Bessel functions.
double bessel_jl(unsigned int n, double x);

/// Spherical Legendre functions.
double legendre_sphPlm(unsigned int l, unsigned int m, double theta);

/// Spherical Neumann functions.
double bessel_yl(unsigned int n, double x);

/// Non-normalized lower incomplete gamma functions.
double gamma_l(double a, double x);

/// Normalized incomplete gamma functions.
double gamma_q(double a, double x);

/// Complementary normalized incomplete gamma functions.
double gamma_p(double a, double x);

/// Non-normalized (upper) incomplete gamma functions.
double gamma_u(double a, double x);

/// Incomplete beta functions.
double ibeta(double a, double b, double x);

/// Dilogarithm function.
double dilog(double x);

/// Digamma or psi function.
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

/// Sinus cardinal function.
double sinc(double x);

/// Hyperbolic sinus cardinal function.
double sinhc_pi(double x);

/// Hyperbolic sinus cardinal function.
double sinhc(double x);

/// Log upper Pochhammer symbol.
double lnpoch(double a, double x);

/// Upper Pochhammer symbol.
double poch(double a, double x);

/// Log factorial.
double lnfact(unsigned int n);

/// Factorial.
double fact(unsigned int n);

/// Log double factorial.
double lndoublefact(int n);

/// Double factorial.
double doublefact(int n);

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

/// Zernicke polynomials
double zernicke(unsigned int n, int m, double rho, double phi);

/// Cylindrical Hankel functions of the first kind.
std::complex<double> cyl_hankel_1(double nu, double x);

/// Cylindrical Hankel functions of the second kind.
std::complex<double> cyl_hankel_2(double nu, double x);

/// Spherical Hankel functions of the first kind.
std::complex<double> sph_hankel_1(unsigned int n, double x);

/// Spherical Hankel functions of the second kind.
std::complex<double> sph_hankel_2(unsigned int n, double x);

/// Heuman lambda functions.
double heuman_lambda(double k, double phi);

/// Jacobi zeta functions.
double jacobi_zeta(double k, double phi);

} // namespace gsl

