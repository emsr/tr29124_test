#ifndef WRAP_GST_H
#define WRAP_GST_H 1

#include <utility> // For pair.

namespace gst
{

/// Scorer functions.

double scorer_gi(double x);

double scorer_hi(double x);

/// Parabolic cylinder functions.

double parab_cylinder_u(double a, double x);

double parab_cylinder_v(double a, double x);

/// Incomplete gamma functions.

double pgamma(double a, double x);

double qgamma(double a, double x);

std::pair<double, double> gamma_cdf(double a, double x);

/// Inverse Incomplete gamma functions.

double inv_gamma_cdf(double a, std::pair<double, double> pq);

/// Toroidal harmonic functions

double tor_harmonic_p(unsigned int l, unsigned int m, double x);

double tor_harmonic_q(unsigned int l, unsigned int m, double x);

/// Spheroidal harmonic functions

double pro_sph_harmonic_p(unsigned int l, unsigned int m, double x);

double pro_sph_harmonic_q(unsigned int l, unsigned int m, double x);

double obl_sph_harmonic_r(unsigned int l, unsigned int m, double x);

double obl_sph_harmonic_t(unsigned int l, unsigned int m, double x);

} // namespace gst

#endif // WRAP_GST_H
