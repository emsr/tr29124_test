

#if LOCAL
#  include "cmath_local"
#else
#  include <cmath>
#endif

#include <sstream>
#include <stdexcept>

#include "gsl_wrap.h"


///  Airy Ai functions.
double
wrap_gsl_sf_airy_ai(double x)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_airy_Ai_e(x, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_airy_ai:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_airy_bi:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_laguerre_nm:");
      msg << " n=" << n << " m=" << m << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
        {
          std::ostringstream msg("Error in wrap_gsl_sf_legendre_Plm:");
          msg << " l=" << l << " m=" << m << " x=" << x;
          throw std::runtime_error(msg.str());
        }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_beta:");
      msg << " x=" << x << " y=" << y;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_Kcomp:");
      msg << " k=" << k;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_Ecomp:");
      msg << " k=" << k;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_Pcomp:");
      msg << " k=" << k << " nu=" << nu;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_hyperg_1F1:");
      msg << " a=" << a << " c=" << c << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_bessel_Inu:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_bessel_Jnu:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_bessel_Knu:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_bessel_Ynu:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_F:");
      msg << " k=" << k << " phi=" << phi;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_E:");
      msg << " phi=" << phi << " k=" << k;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_P:");
      msg << " k=" << k << " nu=" << nu << " phi=" << phi;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


///  Carlson elliptic integrals.
double
wrap_gsl_sf_ellint_RC(double x, double y)
{
  if (x == 0.0 && y == 0.0)
    return 0.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RC_e(x, y, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_RC:");
      msg << " x=" << x << " y=" << y;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_ellint_RD(double x, double y, double z)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return 0.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RD_e(x, y, z, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_RD:");
      msg << " x=" << x << " y=" << y << " z=" << z;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_ellint_RF(double x, double y, double z)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return 0.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RF_e(x, y, z, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_RF:");
      msg << " x=" << x << " y=" << y << " z=" << z;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_ellint_RJ(double x, double y, double z, double p)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return 0.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RJ_e(x, y, z, p, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_ellint_RJ:");
      msg << " x=" << x << " y=" << y << " z=" << z << " p=" << p;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_expint_Ei:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_hyperg_2F1:");
      msg << " a=" << a << " b=" << b << " c=" << c << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_laguerre_n:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_legendre_Pl:");
      msg << " l=" << l << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_zeta:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_hzeta:");
      msg << " s=" << s << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_bessel_jl:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
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
        {
          std::ostringstream msg("Error in wrap_gsl_sf_legendre_sphPlm");
          msg << " l=" << l << " m=" << m << " theta=" << theta;
          throw std::runtime_error();
        }
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
    {
      std::ostringstream msg("Error in wrap_gsl_sf_bessel_yl:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_gamma_inc_Q(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_Q_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_gamma_inc_Q:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_gamma_inc_P(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_P_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_gamma_inc_P:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_gamma_inc(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_gamma_inc:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_beta_inc(double a, double b, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_beta_inc_e(a, b, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_beta_inc:");
      msg << " a=" << a << " b=" << b << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_dilog(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_dilog_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_dilog:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
wrap_gsl_sf_psi(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_psi_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in wrap_gsl_sf_psi:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}
