
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "wrap_gsl.h"

#include <gsl/gsl_sf.h>
#include "gslextras/Fresnel/fresnel.h"
#include "gslextras/Jacobi/jacobi-0.9.2/src/jacobi.h"
#include "gslextras/Hermite/gsl_sf_hermite.h"

namespace gsl
{

/// Airy Ai function.
double
airy_ai(double x)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_airy_Ai_e(x, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in airy_ai:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Airy Bi function.
double
airy_bi(double x)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_airy_Bi_e(x, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in airy_bi:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Associated Laguerre polynomials.
double
assoc_laguerre(unsigned int n, unsigned int m, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_laguerre_n_e(static_cast<int>(n), static_cast<int>(m), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in assoc_laguerre:");
      msg << " n=" << n << " m=" << m << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Associated Legendre functions.
double
assoc_legendre(unsigned int l, unsigned int m, double x)
{
  if (l < m)
   return std::numeric_limits<double>::quiet_NaN();
  else
    {
      gsl_sf_result result;
      int stat = gsl_sf_legendre_Plm_e(static_cast<int>(l), static_cast<int>(m), x, &result);
      if (stat != GSL_SUCCESS)
        {
          std::ostringstream msg("Error in assoc_legendre:");
          msg << " l=" << l << " m=" << m << " x=" << x;
          throw std::runtime_error(msg.str());
        }
      else
        return result.val;
    }
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
  gsl_sf_result result;
  int stat = gsl_sf_beta_e(x, y, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in beta:");
      msg << " x=" << x << " y=" << y;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Complete elliptic integrals of the first kind.
double
comp_ellint_1(double k)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_Kcomp_e(k, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_1:");
      msg << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Complete elliptic integrals of the second kind.
double
comp_ellint_2(double k)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_Ecomp_e(k, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_2:");
      msg << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Complete elliptic integrals of the third kind.
double
comp_ellint_3(double k, double nu)
{
  //double phi = M_PI / 2.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  //int stat = gsl_sf_ellint_P_e(phi, k, nu, mode, &result);
  int stat = gsl_sf_ellint_Pcomp_e(k, nu, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in comp_ellint_3:");
      msg << " k=" << k << " nu=" << nu;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Complete Legendre elliptic D integrals.
double
comp_ellint_d(double /*k*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Confluent hypergeometric functions.
double
conf_hyperg(double a, double c, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hyperg_1F1_e(a, c, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in conf_hyperg:");
      msg << " a=" << a << " c=" << c << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Confluent hypergeometric limit functions.
double
hyperg_0F1(double c, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hyperg_0F1_e(c, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in hyperg_0F1:");
      msg << " c=" << c << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Regular modified cylindrical Bessel functions.
double
cyl_bessel_i(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Inu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in cyl_bessel_i:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Cylindrical Bessel functions (of the first kind).
double
cyl_bessel_j(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Jnu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in cyl_bessel_j:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Irregular modified cylindrical Bessel functions.
double
cyl_bessel_k(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Knu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in cyl_bessel_k:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Cylindrical Neumann functions.
double
cyl_neumann(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Ynu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in cyl_neumann:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Elliptic integrals of the first kind.
double
ellint_1(double k, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_F_e(phi, k, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_1:");
      msg << " k=" << k << " phi=" << phi;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Elliptic integrals of the second kind.
double
ellint_2(double k, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_E_e(phi, k, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_2:");
      msg << " phi=" << phi << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Elliptic integrals of the third kind.
double
ellint_3(double k, double nu, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_P_e(phi, k, nu, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_3:");
      msg << " k=" << k << " nu=" << nu << " phi=" << phi;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Legendre elliptic D integrals.
double
ellint_d(double /*k*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_C.
double
ellint_rc(double x, double y)
{
  if (x == 0.0 && y == 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RC_e(x, y, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_rc:");
      msg << " x=" << x << " y=" << y;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Carlson elliptic integrals R_D.
double
ellint_rd(double x, double y, double z)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RD_e(x, y, z, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_rd:");
      msg << " x=" << x << " y=" << y << " z=" << z;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Carlson elliptic integrals R_F.
double
ellint_rf(double x, double y, double z)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RF_e(x, y, z, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_rf:");
      msg << " x=" << x << " y=" << y << " z=" << z;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Carlson elliptic integrals R_G.
double
ellint_rg(double /*x*/, double /*y*/, double /*z*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Carlson elliptic integrals R_J.
double
ellint_rj(double x, double y, double z, double p)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RJ_e(x, y, z, p, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_rj:");
      msg << " x=" << x << " y=" << y << " z=" << z << " p=" << p;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Exponential integral Ei.
double
expint(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_expint_Ei_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in expint:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Exponential integrals E_n.
double
expint(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_expint_En_e(n, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in expint:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Hermite polynomials.
double
hermite(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hermite_phys_e(n, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in hermite:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Hypergeometric functions.
double
hyperg(double a, double b, double c, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hyperg_2F1_e(a, b, c, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in hyperg:");
      msg << " a=" << a << " b=" << b << " c=" << c << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Laguerre polynomials.
double
laguerre(unsigned int n, double x)
{
  int m = 0;
  gsl_sf_result result;
  int stat = gsl_sf_laguerre_n_e(static_cast<int>(n), m, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in laguerre:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Legendre polynomials.
double
legendre_p(unsigned int l, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_legendre_Pl_e(static_cast<int>(l), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in legendre_p:");
      msg << " l=" << l << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Legendre polynomials of the second kind.
double
legendre_q(unsigned int l, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_legendre_Ql_e(static_cast<int>(l), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in legendre_q:");
      msg << " l=" << l << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Riemann zeta function.
double
riemann_zeta(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_zeta_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in riemann_zeta:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Hurwitz zeta functions.
double
hurwitz_zeta(double s, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hzeta_e(s, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in hurwitz_zeta:");
      msg << " s=" << s << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Dirichlet eta function.
double
dirichlet_eta(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_eta_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in dirichlet_eta:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Spherical Bessel functions.
double
sph_bessel(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_jl_e(static_cast<int>(n), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in sph_bessel:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Spherical Legendre functions.
double
sph_legendre(unsigned int l, unsigned int m, double theta)
{
  double x = std::cos(theta);
  if (l < m)
    return std::numeric_limits<double>::quiet_NaN();
  else
    {
      gsl_sf_result result;
      int stat = gsl_sf_legendre_sphPlm_e(static_cast<int>(l), static_cast<int>(m), x, &result);
      if (stat != GSL_SUCCESS)
        {
          std::ostringstream msg("Error in sph_legendre");
          msg << " l=" << l << " m=" << m << " theta=" << theta;
          throw std::runtime_error(msg.str());
        }
      else
        return result.val;
    }
}

/// Spherical Neumann functions.
double
sph_neumann(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_yl_e(static_cast<int>(n), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in sph_neumann:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Non-normalized lower incomplete gamma functions.
double
tgamma_lower(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_e(a, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in tgamma_lower:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val - tgamma(a, x);
}

/// Non-normalized lower incomplete gamma functions.
double
tgamma(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in tgamma:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Normalized incomlete gamma functions.
double
qgamma(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_Q_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in qgamma:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Normalized incomlete gamma functions.
double
pgamma(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_P_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in pgamma:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Incomlete beta functions.
double
ibeta(double a, double b, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_beta_inc_e(a, b, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ibeta:");
      msg << " a=" << a << " b=" << b << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Complementary incomlete beta functions.
double
ibetac(double a, double b, double x)
{
  return 1.0 - ibeta(b, a, 1.0 - x);
}

/// Dilogarithm function.
double
dilog(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_dilog_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in dilog:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Digamma or psi function.
double
psi(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_psi_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in psi:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Sine integral.
double
sinint(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_Si_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in sinint:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Cosine integral.
double
cosint(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_Ci_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in cosint:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Hyperbolic sine integral.
double
sinhint(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_Shi_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in sinhint:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Hyperbolic cosine integral.
double
coshint(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_Chi_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in coshint:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Gegenbauer polynomials.
double
gegenpoly_n(unsigned int n, double lambda, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gegenpoly_n_e(n, lambda, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in gegenpoly_n:");
      msg << " n=" << n << " lambda=" << lambda << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Hydrogen wave functions.
double
hydrogen(unsigned int n, double l, double Z, double r)
{
  gsl_sf_result result;
  int stat = gsl_sf_hydrogenicR_e(n, l, Z, r, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in hydrogen:");
      msg << " n=" << n << " l=" << l << " Z=" << Z << " r=" << r;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Dawson integral.
double
dawson(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_dawson_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in dawson:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Jacobian elliptic integrals sn.
double
jacobi_sn(double k, double u)
{
  double m = k * k;
  double sn, cn, dn;
  int stat = gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in jacobi_sn:");
      msg << " u=" << u << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return sn;
}

/// Jacobian elliptic integrals cn.
double
jacobi_cn(double k, double u)
{
  double m = k * k;
  double sn, cn, dn;
  int stat = gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in jacobi_cn:");
      msg << " u=" << u << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return cn;
}

/// Jacobian elliptic integrals dn.
double
jacobi_dn(double k, double u)
{
  double m = k * k;
  double sn, cn, dn;
  int stat = gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in jacobi_dn:");
      msg << " u=" << u << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return dn;
}

/// Fresnel cosine integral.
double
fresnel_c(double x)
{ return ::fresnel_c(x); }

/// Fresnel sine integral.
double
fresnel_s(double x)
{ return ::fresnel_s(x); }

/// Sinus cardinal function.
double
sinc(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_sinc_e(x / M_PI, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in sinc:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Reperiodized sinus cardinal function.
double
sinc_pi(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_sinc_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in sinc_pi:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
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
  return std::sinh(M_PI * x) / (M_PI * x);
}

/// Log upper Pochhammer symbol.
double
lpochhammer(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_lnpoch_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in lpochhammer:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Log lower Pochhammer symbol.
double
lpochhammer_lower(double a, double x)
{
  if (a == x)
    return std::numeric_limits<double>::infinity();
  gsl_sf_result result_num;
  int stat_num = gsl_sf_lngamma_e(std::abs(a - x), &result_num);
  gsl_sf_result result_den;
  int stat_den = gsl_sf_lngamma_e(std::abs(a), &result_den);
  if (stat_num != GSL_SUCCESS && stat_den != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in lpochhammer_lower:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result_num.val - result_den.val;
}

/// Upper Pochhammer symbol.
double
pochhammer(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_poch_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in pochhammer:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Lower Pochhammer symbol.
double
pochhammer_lower(double a, double x)
{
  if (a == x)
    return std::numeric_limits<double>::infinity();
  if (std::fmod(std::abs(a - x), M_PI) < 1.0e-12)
    return std::numeric_limits<double>::infinity();
  if (std::fmod(std::abs(a), M_PI) < 1.0e-12)
    return std::numeric_limits<double>::infinity();
  gsl_sf_result result_num;
  int stat_num = gsl_sf_gamma_e(std::abs(a - x), &result_num);
  gsl_sf_result result_den;
  int stat_den = gsl_sf_gamma_e(std::abs(a), &result_den);
  if (stat_num != GSL_SUCCESS && stat_den != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in lpochhammer_lower:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result_num.val / result_den.val;
}

/// Log factorial.
double
lfactorial(unsigned int n)
{
  gsl_sf_result result;
  int stat = gsl_sf_lnfact_e(n, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in lfactorial:");
      msg << " n=" << n;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Factorial.
double
factorial(unsigned int n)
{
  gsl_sf_result result;
  int stat = gsl_sf_fact_e(n, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in factorial:");
      msg << " n=" << n;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Log double factorial.
double
ldouble_factorial(int n)
{
  gsl_sf_result result;
  int stat = gsl_sf_lndoublefact_e(n, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ldouble_factorial:");
      msg << " n=" << n;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Double factorial.
double
double_factorial(int n)
{
  gsl_sf_result result;
  int stat = gsl_sf_doublefact_e(n, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in double_factorial:");
      msg << " n=" << n;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Regular modified spherical bessel functions.
double
sph_bessel_i(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_il_scaled_e(n, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in sph_bessel_i:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return std::exp(std::abs(x)) * result.val;
}

/// Irregular modified spherical bessel functions.
double
sph_bessel_k(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_kl_scaled_e(n, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in sph_bessel_k:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return std::exp(-x) * result.val;
}

/// Chebyshev polynomials of the first kind.
double
chebyshev_t(unsigned int /*n*/, double /*x*/)
{ return std::numeric_limits<double>::quiet_NaN(); }

/// Binomial coefficients.
double
choose(unsigned int n, unsigned int k)
{
  if (k > n)
    return std::numeric_limits<double>::quiet_NaN(); // GSL barfs on this for no reason.
  gsl_sf_result result;
  int stat = gsl_sf_choose_e(n, k, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in choose:");
      msg << " n=" << n << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Log binomial coefficients.
double
lnchoose(unsigned int n, unsigned int k)
{
  if (k > n)
    return -std::numeric_limits<double>::infinity(); // GSL barfs on this for no reason.
  gsl_sf_result result;
  int stat = gsl_sf_lnchoose_e(n, k, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in lnchoose:");
      msg << " n=" << n << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Taylor coefficients.
double
taylorcoeff(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_taylorcoeff_e(n, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in taylorcoeff:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Jacobi polynomials.
double
jacobi(unsigned int n, double alpha, double beta, double x)
{
  gsl_sf_result result;
  int stat = jac_jacobi_e(x, n, alpha, beta, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in jacobi:");
      msg << " n=" << n << " alpha=" << alpha << " beta=" << beta << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Radial polynomials
double
radpoly(unsigned int n, unsigned int m, double rho)
{
  if (m > n)
    throw std::range_error("radpoly: m > n");
  else if ((n - m) % 2 == 1)
    return std::numeric_limits<double>::quiet_NaN();
  else
    {
      auto k = (n - m) / 2;
      return (k % 2 == 0 ? +1 : -1) * std::pow(rho, double(m))
	   * jacobi(k, double(m), 0.0, 1.0 - 2.0 * rho * rho);
    }
}

/// Zernike polynomials
double
zernike(unsigned int n, int m, double rho, double phi)
{
  return radpoly(n, std::abs(m), rho)
       * (m >= 0 ? std::cos(m * phi) : std::sin(m * phi));
}

/// Cylindrical Hankel functions of the first kind.
std::complex<double>
cyl_hankel_1(double /*nu*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Cylindrical Hankel functions of the second kind.
std::complex<double>
cyl_hankel_2(double /*nu*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Hankel functions of the first kind.
std::complex<double>
sph_hankel_1(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical Hankel functions of the second kind.
std::complex<double>
sph_hankel_2(unsigned int /*n*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Heuman lambda functions.
double
heuman_lambda(double /*k*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Jacobi zeta functions.
double
jacobi_zeta(double /*k*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse error function.
double
erf_inv(double /*p*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Inverse complementary error function.
double
erfc_inv(double /*p*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Spherical harmonic functions.
std::complex<double>
sph_harmonic(unsigned int /*l*/, int /*m*/, double /*theta*/, double /*phi*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Owen's T function.
double
owens_t(double /*h*/, double /*a*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

/// Clausen function of order 2.
double
clausen_c(unsigned int m, double w)
{
  if (m != 2)
    return std::numeric_limits<double>::quiet_NaN();
  gsl_sf_result result;
  int stat = gsl_sf_clausen_e(w, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in clausen_c:");
      msg << " w=" << w;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
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

/// Fermi-Dirac integrals
double
fermi_dirac(double s, double x)
{
  int stat;
  gsl_sf_result result;
  if (s == -0.5)
    stat = gsl_sf_fermi_dirac_mhalf_e(x, &result) / std::tgamma(s + 1.0);
  else if (s == 0.5)
    stat = gsl_sf_fermi_dirac_half_e(x, &result) / std::tgamma(s + 1.0);
  else if (s == 1.5)
    stat = gsl_sf_fermi_dirac_3half_e(x, &result) / std::tgamma(s + 1.0);
  else if (int(s) == s)
    stat = gsl_sf_fermi_dirac_int_e(int(s), x, &result);
  else
    return std::numeric_limits<double>::quiet_NaN();

  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in fermi_dirac:");
      msg << " s=" << s << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Bose-Einstein integrals
double
bose_einstein(double /*s*/, double /*x*/)
{
  return std::numeric_limits<double>::quiet_NaN();
}

} // namespace gsl
