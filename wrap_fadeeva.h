#ifndef WRAP_FADEEVA_H
#define WRAP_FADEEVA_H 1

#include <complex>

namespace faddeeva
{

  /// Compute w(z) = exp(-z^2) erfc(-iz)
  /// the Faddeeva or scaled complex error function.
  std::complex<double> faddeeva(std::complex<double> z);

  /// Compute w(z) = exp(-z^2) erfc(-iz)
  /// the Faddeeva or scaled complex error function.
  double faddeeva(double x);

  /// Compute erfcx(z) = exp(z^2) erfc(z).
  std::complex<double> erfc_scaled(std::complex<double> z);

  /// Compute erfcx(z) = exp(z^2) erfc(z).
  double erfc_scaled(double x);

  /// Compute erf(z), the error function of complex arguments.
  std::complex<double> erf(std::complex<double> z);

  /// Compute erf(z), the error function of real arguments.
  double erf(double x);

  /// Compute erfi(z) = -i erf(iz), the imaginary error function.
  std::complex<double> erfi(std::complex<double> z);

  /// Compute erfi(z) = -i erf(iz), the imaginary error function.
  double erfi(double x);

  /// Compute erfc(z) = 1 - erf(z), the complementary error function.
  std::complex<double> erfc(std::complex<double> z);

  /// Compute erfc(z) = 1 - erf(z), the complementary error function.
  double erfc(double x);

  /// Compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z).
  std::complex<double> dawson(std::complex<double> z);

  /// Compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z).
  double dawson(double x);

} // namespace faddeeva

#endif // WRAP_FADEEVA_H

