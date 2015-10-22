

///
///  5.2.1.1  Associated Laguerre polynomials.
///
double wrap_gsl_sf_laguerre_nm(unsigned int n, unsigned int m, double x);


///
///  5.2.1.2  Associated Legendre functions.
///
double wrap_gsl_sf_legendre_Plm(unsigned int l, unsigned int m, double x);


///
///  5.2.1.3  Beta function.
///
double wrap_gsl_sf_beta(double x, double y);


///
///  5.2.1.4  Complete elliptic integrals of the first kind.
///
double wrap_gsl_sf_ellint_Kcomp(double k);


///
///  5.2.1.5  Complete elliptic integrals of the second kind.
///
double wrap_gsl_sf_ellint_Ecomp(double k);


///
///  5.2.1.6  Complete elliptic integrals of the third kind.
///
double wrap_gsl_sf_ellint_Pcomp(double k, double nu);


///
///  5.2.1.7  Confluent hypergeometric functions.
///
double wrap_gsl_sf_hyperg_1F1(double a, double c, double x);


///
///  5.2.1.8  Regular modified cylindrical Bessel functions.
///
double wrap_gsl_sf_bessel_Inu(double nu, double x);


///
///  5.2.1.9  Cylindrical Bessel functions (of the first kind).
///
double wrap_gsl_sf_bessel_Jnu(double nu, double x);

double wrap_gsl_sf_bessel_Jnu_asymp(double nu, double x);


///
///  5.2.1.10  Irregular modified cylindrical Bessel functions.
///
double wrap_gsl_sf_bessel_Knu(double nu, double x);


///
///  5.2.1.11  Cylindrical Neumann functions.
///
double wrap_gsl_sf_bessel_Ynu(double nu, double x);

double wrap_gsl_sf_bessel_Ynu_asymp(double nu, double x);


///
///  5.2.1.12  Elliptic integrals of the first kind.
///
double wrap_gsl_sf_ellint_F(double k, double phi);


///
///  5.2.1.13  Elliptic integrals of the second kind.
///
double wrap_gsl_sf_ellint_E(double k, double phi);


///
///  5.2.1.14  Elliptic integrals of the third kind.
///
double wrap_gsl_sf_ellint_P(double k, double nu, double phi);


///
///  5.2.1.15  Exponential integral.
///
double wrap_gsl_sf_expint_Ei(double x);


///
///  5.2.1.17  Hypergeometric functions.
///
double wrap_gsl_sf_hyperg_2F1(double a, double b, double c, double x);


///
///  5.2.1.18  Laguerre polynomials.
///
double wrap_gsl_sf_laguerre_n(unsigned int n, double x);


///
///  5.2.1.19  Legendre polynomials.
///
double wrap_gsl_sf_legendre_Pl(unsigned int l, double x);


///
///  5.2.1.20  Riemann zeta function.
///
double wrap_gsl_sf_zeta(double x);


///
///  5.2.1.21  Spherical Bessel functions.
///
double wrap_gsl_sf_bessel_jl(unsigned int n, double x);


///
///  5.2.1.22  Spherical Legendre functions.
///
double wrap_gsl_sf_legendre_sphPlm(unsigned int l, unsigned int m, double theta);


///
///  5.2.1.23  Spherical Neumann functions.
///
double wrap_gsl_sf_bessel_yl(unsigned int n, double x);

