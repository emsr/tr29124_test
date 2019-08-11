
#include <wrap_gst.h>

#include <scorer.h>
#include <parabolic_cylinder.h>
#include <incomplete_gamma.h>
#include <spheroidal_harmonic.h>
#include <toroidal_harmonic.h>

#include <vector>

namespace gst
{

/// Scorer functions.

double
scorer_gi(double x)
{
  int ifacg = 1;
  double y = 0.0, re_g, im_g, re_gp, im_gp;
  int ierrog = 0;
  scorergi(ifacg, x, &y, &re_g, &im_g, &re_gp, &im_gp, &ierrog);
  return re_g;
}

double
scorer_hi(double x)
{
  int ifach = 1;
  double y = 0.0, re_h, im_h, re_hp, im_hp;
  int ierroh = 0;
  scorergi(ifach, x, &y, &re_h, &im_h, &re_hp, &im_hp, &ierroh);
  return re_h;
}

/// Parabolic cylinder functions.

double
parab_cylinder_u(double a, double x)
{
  int mode = 1;
  double uaxx, vaxx;
  int ierr = 0;
  parab_cyl(a, x, mode, &uaxx, &vaxx, &ierr);
  return uaxx;
}

double
parab_cylinder_v(double a, double x)
{
  int mode = 1;
  double uaxx, vaxx;
  int ierr = 0;
  parab_cyl(a, x, mode, &uaxx, &vaxx, &ierr);
  return vaxx;
}

/// Incomplete gamma functions.

double
pgamma(double a, double x)
{
  double p, q;
  int ierr = 0;
  inc_gamma(a, x, &p, &q, &ierr);
  return p;
}

double
qgamma(double a, double x)
{
  double p, q;
  int ierr = 0;
  inc_gamma(a, x, &p, &q, &ierr);
  return q;
}

/// Toroidal harmonic functions

double
tor_harmonic_p(unsigned int l, unsigned int m, double x)
{
  std::vector<double> pl(l + 1), ql(l + 1);
  int newn = 0;
  tor_harmonic(x, m, l, pl.data(), ql.data(), &newn);
  return pl[l];
}

double
tor_harmonic_q(unsigned int l, unsigned int m, double x)
{
  std::vector<double> pl(l + 1), ql(l + 1);
  int newn = 0;
  tor_harmonic(x, m, l, pl.data(), ql.data(), &newn);
  return ql[l];
}

/// Spheroidal harmonic functions

double
pro_sph_harmonic_p(unsigned int l, unsigned int m, double x)
{
  std::vector<double> pl(l + 1), ql(l + 1);
  int mode = 1, nuevo = 0;
  pro_sph_harm(x, m, l, mode, pl.data(), ql.data(), nuevo);
  return pl[l];
}

double
pro_sph_harmonic_q(unsigned int l, unsigned int m, double x)
{
  std::vector<double> pl(l + 1), ql(l + 1);
  int mode = 1, nuevo = 0;
  pro_sph_harm(x, m, l, mode, pl.data(), ql.data(), nuevo);
  return ql[l];
}

double
obl_sph_harmonic_r(unsigned int l, unsigned int m, double x)
{
  std::vector<double> rl(l + 1), tl(l + 1);
  int mode = 1, nuevo = 0;
  obl_sph_harm(x, m, l, mode, rl.data(), tl.data(), nuevo);
  return rl[l];
}

double
obl_sph_harmonic_t(unsigned int l, unsigned int m, double x)
{
  std::vector<double> rl(l + 1), tl(l + 1);
  int mode = 1, nuevo = 0;
  obl_sph_harm(x, m, l, mode, rl.data(), tl.data(), nuevo);
  return tl[l];
}

} // namespace gst


