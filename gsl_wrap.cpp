
#include <cmath>
#include <stdexcept>

//#include <gsl/gsl_sf.h>

#include "gsl_wrap.h"


///  Airy Ai functions.
double
wrap_gsl_sf_airy_ai(double x)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_airy_Ai_e(x, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_airy_ai");
  else
    return result.val;
}


///  Airy Bi functions.
double
wrap_gsl_sf_airy_bi(double x)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_airy_Bi_e(x, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_airy_bi");
  else
    return result.val;
}


///  5.2.1.1  Associated Laguerre polynomials.
double
wrap_gsl_sf_laguerre_nm(unsigned int n, unsigned int m, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_laguerre_n_e(static_cast<int>(n), static_cast<int>(m), x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_laguerre_nm");
  else
    return result.val;
}


///  5.2.1.2  Associated Legendre functions.
double
wrap_gsl_sf_legendre_Plm(unsigned int l, unsigned int m, double x)
{
  if (l < m)
   return 0.0;
  else
    {
      gsl_sf_result result;
      int stat = gsl_sf_legendre_Plm_e(static_cast<int>(l), static_cast<int>(m), x, &result);
      if (stat != GSL_SUCCESS)
        throw std::runtime_error("Error in wrap_gsl_sf_legendre_Plm");
      else
        return result.val;
    }
}


///  5.2.1.3  Beta function.
double
wrap_gsl_sf_beta(double x, double y)
{
  gsl_sf_result result;
  int stat = gsl_sf_beta_e(x, y, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_beta");
  else
    return result.val;
}


///  5.2.1.4  Complete elliptic integrals of the first kind.
double
wrap_gsl_sf_ellint_Kcomp(double k)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_Kcomp_e(k, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_Kcomp");
  else
    return result.val;
}


///  5.2.1.5  Complete elliptic integrals of the second kind.
double
wrap_gsl_sf_ellint_Ecomp(double k)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_Ecomp_e(k, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_Ecomp");
  else
    return result.val;
}


///  5.2.1.6  Complete elliptic integrals of the third kind.
double
wrap_gsl_sf_ellint_Pcomp(double k, double nu)
{
  //double phi = M_PI / 2.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  //int stat = gsl_sf_ellint_P_e(phi, k, nu, mode, &result);
  int stat = gsl_sf_ellint_Pcomp_e(k, nu, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_Pcomp");
  else
    return result.val;
}


///  5.2.1.7  Confluent hypergeometric functions.
double
wrap_gsl_sf_hyperg_1F1(double a, double c, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hyperg_1F1_e(a, c, x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_hyperg_1F1");
  else
    return result.val;
}


///  5.2.1.8  Regular modified cylindrical Bessel functions.
double
wrap_gsl_sf_bessel_Inu(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Inu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_bessel_Inu");
  else
    return result.val;
}


///  5.2.1.9  Cylindrical Bessel functions (of the first kind).
double
wrap_gsl_sf_bessel_Jnu(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Jnu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_bessel_Jnu");
  else
    return result.val;
}


///  5.2.1.9  Cylindrical Bessel functions (of the first kind).
//double
//wrap_gsl_sf_bessel_Jnu_asymp(double nu, double x)
//{
//  gsl_sf_result result;
//  int stat = gsl_sf_bessel_Jnu_asympx_e(nu, x, &result);
//  if (stat != GSL_SUCCESS)
//    throw std::runtime_error("Error in wrap_gsl_sf_bessel_Jnu_asymp");
//  else
//    return result.val;
//}


///  5.2.1.10  Irregular modified cylindrical Bessel functions.
double
wrap_gsl_sf_bessel_Knu(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Knu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_bessel_Knu");
  else
    return result.val;
}


///  5.2.1.11  Cylindrical Neumann functions.
double
wrap_gsl_sf_bessel_Ynu(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Ynu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_bessel_Ynu");
  else
    return result.val;
}


///  5.2.1.11  Cylindrical Neumann functions.
//double
//wrap_gsl_sf_bessel_Ynu_asymp(double nu, double x)
//{
//  gsl_sf_result result;
//  int stat = gsl_sf_bessel_Ynu_asympx_e(nu, x, &result);
//  if (stat != GSL_SUCCESS)
//    throw std::runtime_error("Error in wrap_gsl_sf_bessel_Ynu_asymp");
//  else
//    return result.val;
//}


///  5.2.1.12  Elliptic integrals of the first kind.
double
wrap_gsl_sf_ellint_F(double k, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_F_e(phi, k, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_F");
  else
    return result.val;
}


///  5.2.1.13  Elliptic integrals of the second kind.
double
wrap_gsl_sf_ellint_E(double k, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_E_e(phi, k, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_E");
  else
    return result.val;
}


///  5.2.1.14  Elliptic integrals of the third kind.
double
wrap_gsl_sf_ellint_P(double k, double nu, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_P_e(phi, k, nu, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_P");
  else
    return result.val;
}


///  Carlson elliptic integrals.
double
wrap_gsl_sf_ellint_RC(double x, double y)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RC_e(x, y, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_RC");
  else
    return result.val;
}

double
wrap_gsl_sf_ellint_RD(double x, double y, double z)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RD_e(x, y, z, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_RD");
  else
    return result.val;
}

double
wrap_gsl_sf_ellint_RF(double x, double y, double z)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RF_e(x, y, z, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_RF");
  else
    return result.val;
}

double
wrap_gsl_sf_ellint_RJ(double x, double y, double z, double p)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RJ_e(x, y, z, p, mode, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_ellint_RJ");
  else
    return result.val;
}


///  5.2.1.15  Exponential integral.
double
wrap_gsl_sf_expint_Ei(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_expint_Ei_e(x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_expint_Ei");
  else
    return result.val;
}


///  5.2.1.17  Hypergeometric functions.
double
wrap_gsl_sf_hyperg_2F1(double a, double b, double c, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hyperg_2F1_e(a, b, c, x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_hyperg_2F1");
  else
    return result.val;
}


///  5.2.1.18  Laguerre polynomials.
double
wrap_gsl_sf_laguerre_n(unsigned int n, double x)
{
  int m = 0;
  gsl_sf_result result;
  int stat = gsl_sf_laguerre_n_e(static_cast<int>(n), m, x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_laguerre_n");
  else
    return result.val;
}


///  5.2.1.19  Legendre polynomials.
double
wrap_gsl_sf_legendre_Pl(unsigned int l, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_legendre_Pl_e(static_cast<int>(l), x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_legendre_Pl");
  else
    return result.val;
}


///  5.2.1.20  Riemann zeta function.
double
wrap_gsl_sf_zeta(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_zeta_e(x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_zeta");
  else
    return result.val;
}


///  Hurwitz zeta function.
double
wrap_gsl_sf_hzeta(double s, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hzeta_e(s, x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_hzeta");
  else
    return result.val;
}


///  5.2.1.21  Spherical Bessel functions.
double
wrap_gsl_sf_bessel_jl(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_jl_e(static_cast<int>(n), x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_bessel_jl");
  else
    return result.val;
}


///  5.2.1.22  Spherical Legendre functions.
double
wrap_gsl_sf_legendre_sphPlm(unsigned int l, unsigned int m, double theta)
{
  double x = std::cos(theta);
  if (l < m)
    return 0.0;
  else
    {
      gsl_sf_result result;
      int stat = gsl_sf_legendre_sphPlm_e(static_cast<int>(l), static_cast<int>(m), x, &result);
      if (stat != GSL_SUCCESS)
        throw std::runtime_error("Error in wrap_gsl_sf_legendre_sphPlm");
      else
        return result.val;
    }
}


///  5.2.1.23  Spherical Neumann functions.
double
wrap_gsl_sf_bessel_yl(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_yl_e(static_cast<int>(n), x, &result);
  if (stat != GSL_SUCCESS)
    throw std::runtime_error("Error in wrap_gsl_sf_bessel_yl");
  else
    return result.val;
}
