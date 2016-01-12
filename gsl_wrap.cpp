

#if LOCAL
#  include "cmath_local"
#else
#  include <cmath>
#endif

#include <sstream>
#include <stdexcept>

#include "gsl_wrap.h"

#include <jacobi.h>

namespace gsl
{

/// Airy Ai functions.
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


/// Airy Bi functions.
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
laguerre_nm(unsigned int n, unsigned int m, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_laguerre_n_e(static_cast<int>(n), static_cast<int>(m), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in laguerre_nm:");
      msg << " n=" << n << " m=" << m << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Associated Legendre functions.
double
legendre_Plm(unsigned int l, unsigned int m, double x)
{
  if (l < m)
   return 0.0;
  else
    {
      gsl_sf_result result;
      int stat = gsl_sf_legendre_Plm_e(static_cast<int>(l), static_cast<int>(m), x, &result);
      if (stat != GSL_SUCCESS)
        {
          std::ostringstream msg("Error in legendre_Plm:");
          msg << " l=" << l << " m=" << m << " x=" << x;
          throw std::runtime_error(msg.str());
        }
      else
        return result.val;
    }
}


/// Beta function.
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
ellint_Kcomp(double k)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_Kcomp_e(k, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_Kcomp:");
      msg << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Complete elliptic integrals of the second kind.
double
ellint_Ecomp(double k)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_Ecomp_e(k, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_Ecomp:");
      msg << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Complete elliptic integrals of the third kind.
double
ellint_Pcomp(double k, double nu)
{
  //double phi = M_PI / 2.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  //int stat = gsl_sf_ellint_P_e(phi, k, nu, mode, &result);
  int stat = gsl_sf_ellint_Pcomp_e(k, nu, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_Pcomp:");
      msg << " k=" << k << " nu=" << nu;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Confluent hypergeometric functions.
double
hyperg_1F1(double a, double c, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hyperg_1F1_e(a, c, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in hyperg_1F1:");
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
bessel_Inu(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Inu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in bessel_Inu:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Cylindrical Bessel functions (of the first kind).
double
bessel_Jnu(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Jnu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in bessel_Jnu:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Irregular modified cylindrical Bessel functions.
double
bessel_Knu(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Knu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in bessel_Knu:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Cylindrical Neumann functions.
double
bessel_Ynu(double nu, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_Ynu_e(nu, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in bessel_Ynu:");
      msg << " nu=" << nu << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Elliptic integrals of the first kind.
double
ellint_F(double k, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_F_e(phi, k, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_F:");
      msg << " k=" << k << " phi=" << phi;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Elliptic integrals of the second kind.
double
ellint_E(double k, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_E_e(phi, k, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_E:");
      msg << " phi=" << phi << " k=" << k;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Elliptic integrals of the third kind.
double
ellint_P(double k, double nu, double phi)
{
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_P_e(phi, k, nu, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_P:");
      msg << " k=" << k << " nu=" << nu << " phi=" << phi;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Carlson elliptic integrals.
double
ellint_RC(double x, double y)
{
  if (x == 0.0 && y == 0.0)
    return 0.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RC_e(x, y, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_RC:");
      msg << " x=" << x << " y=" << y;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}



double
ellint_RD(double x, double y, double z)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return 0.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RD_e(x, y, z, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_RD:");
      msg << " x=" << x << " y=" << y << " z=" << z;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}



double
ellint_RF(double x, double y, double z)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return 0.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RF_e(x, y, z, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_RF:");
      msg << " x=" << x << " y=" << y << " z=" << z;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


double
ellint_RJ(double x, double y, double z, double p)
{
  if (x == 0.0 && y == 0.0 && z == 0.0)
    return 0.0;
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  gsl_sf_result result;
  int stat = gsl_sf_ellint_RJ_e(x, y, z, p, mode, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in ellint_RJ:");
      msg << " x=" << x << " y=" << y << " z=" << z << " p=" << p;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Exponential integral.
double
expint_Ei(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_expint_Ei_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in expint_Ei:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
expint_E1(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_expint_E1_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in expint_E1:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

double
expint_En(int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_expint_En_e(n, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in expint_En:");
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
hyperg_2F1(double a, double b, double c, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hyperg_2F1_e(a, b, c, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in hyperg_2F1:");
      msg << " a=" << a << " b=" << b << " c=" << c << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Laguerre polynomials.
double
laguerre_n(unsigned int n, double x)
{
  int m = 0;
  gsl_sf_result result;
  int stat = gsl_sf_laguerre_n_e(static_cast<int>(n), m, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in laguerre_n:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Legendre polynomials.
double
legendre_Pl(unsigned int l, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_legendre_Pl_e(static_cast<int>(l), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in legendre_Pl:");
      msg << " l=" << l << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Legendre polynomials.
double
legendre_Ql(unsigned int l, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_legendre_Ql_e(static_cast<int>(l), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in legendre_Ql:");
      msg << " l=" << l << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Riemann zeta function.
double
zeta(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_zeta_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in zeta:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Hurwitz zeta function.
double
hzeta(double s, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_hzeta_e(s, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in hzeta:");
      msg << " s=" << s << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Spherical Bessel functions.
double
bessel_jl(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_jl_e(static_cast<int>(n), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in bessel_jl:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Spherical Legendre functions.
double
legendre_sphPlm(unsigned int l, unsigned int m, double theta)
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
          std::ostringstream msg("Error in legendre_sphPlm");
          msg << " l=" << l << " m=" << m << " theta=" << theta;
          throw std::runtime_error(msg.str());
        }
      else
        return result.val;
    }
}


/// Spherical Neumann functions.
double
bessel_yl(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_yl_e(static_cast<int>(n), x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in bessel_yl:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}



double
gamma_inc_Q(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_Q_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in gamma_inc_Q:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}



double
gamma_inc_P(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_P_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in gamma_inc_P:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}



double
gamma_inc(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_gamma_inc_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in gamma_inc:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}



double
beta_inc(double a, double b, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_beta_inc_e(a, b, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in beta_inc:");
      msg << " a=" << a << " b=" << b << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
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


/// Digamma function.
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
Si(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_Si_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in Si:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Cosine integral.
double
Ci(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_Ci_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in Ci:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Hyperbolic sine integral.
double
Shi(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_Shi_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in Shi:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Hyperbolic cosine integral.
double
Chi(double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_Chi_e(x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in Chi:");
      msg << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}


/// Gegenbauer polynomials.
double
gegenpoly_n(int n, double lambda, double x)
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
hydrogen(int n, double l, double Z, double r)
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


/// Jacobian elliptic integrals.
double
elljac_sn(double u, double k)
{
  double m = k * k;
  double sn, cn, dn;
  int stat = gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in elljac_sn:");
      msg << " u=" << u << " m=" << m;
      throw std::runtime_error(msg.str());
    }
  else
    return sn;
}

double
elljac_cn(double u, double k)
{
  double m = k * k;
  double sn, cn, dn;
  int stat = gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in elljac_cn:");
      msg << " u=" << u << " m=" << m;
      throw std::runtime_error(msg.str());
    }
  else
    return cn;
}

double
elljac_dn(double u, double k)
{
  double m = k * k;
  double sn, cn, dn;
  int stat = gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in elljac_dn:");
      msg << " u=" << u << " m=" << m;
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
  // Scale x by pi to match the deinition for the C++ proposal:
  // sinc(x) = sin(x)/x
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

/// Log upper Pochhammer symbol.
double
lnpoch(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_lnpoch_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in lnpoch:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Upper Pochhammer symbol.
double
poch(double a, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_poch_e(a, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in poch:");
      msg << " a=" << a << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Log factorial.
double
lnfact(unsigned int n)
{
  gsl_sf_result result;
  int stat = gsl_sf_lnfact_e(n, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in lnfact:");
      msg << " n=" << n;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Factorial.
double
fact(unsigned int n)
{
  gsl_sf_result result;
  int stat = gsl_sf_fact_e(n, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in fact:");
      msg << " n=" << n;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Log double factorial.
double
lndoublefact(unsigned int n)
{
  gsl_sf_result result;
  int stat = gsl_sf_lndoublefact_e(n, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in lndoublefact:");
      msg << " n=" << n;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Double factorial.
double
doublefact(unsigned int n)
{
  gsl_sf_result result;
  int stat = gsl_sf_doublefact_e(n, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in doublefact:");
      msg << " n=" << n;
      throw std::runtime_error(msg.str());
    }
  else
    return result.val;
}

/// Regular modified spherical bessel functions.
double
bessel_il(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_il_scaled_e(n, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in bessel_il:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return std::exp(std::abs(x)) * result.val;
}

/// Irregular modified spherical bessel functions.
double
bessel_kl(unsigned int n, double x)
{
  gsl_sf_result result;
  int stat = gsl_sf_bessel_kl_scaled_e(n, x, &result);
  if (stat != GSL_SUCCESS)
    {
      std::ostringstream msg("Error in bessel_kl:");
      msg << " n=" << n << " x=" << x;
      throw std::runtime_error(msg.str());
    }
  else
    return std::exp(-x) * result.val;
}

/// Chebyshev polynomials of the first kind.
double
chebyshev_t(unsigned int n, double x)
{ 0.0; }

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

} // namespace gsl
