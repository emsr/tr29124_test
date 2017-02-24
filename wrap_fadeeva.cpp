
#include "wrap_faddeeva.h"

#include "Faddeeva/Faddeeva.h"

namespace faddeeva
{

  /// Compute w(z) = exp(-z^2) erfc(-iz),
  /// the Faddeeva or scaled complex error function.
  std::complex<double>
  faddeeva(std::complex<double> z)
  { return Faddeeva::w(z); }

  /// Compute w(z) = exp(-z^2) erfc(-iz),
  /// the Faddeeva or scaled complex error function.
  double
  faddeeva(double x)
  { return Faddeeva::w_im(x); }

  /// Compute erfcx(z) = exp(z^2) erfc(z).
  std::complex<double>
  erfc_scaled(std::complex<double> z)
  { return Faddeeva::erfcx(z); }

  /// Compute erfcx(z) = exp(z^2) erfc(z).
  double
  erfc_scaled(double x)
  { return Faddeeva::erfcx(double x); }

  /// Compute erf(z), the error function of complex arguments.
  std::complex<double>
  erf(std::complex<double> z)
  { return Faddeeva::erf(z); }

  /// Compute erf(z), the error function of complex arguments.
  double
  erf(double x)
  { return Faddeeva::erf(double x); }

  /// Compute erfi(z) = -i erf(iz), the imaginary error function.
  std::complex<double>
  erfi(std::complex<double> z)
  { return Faddeeva::erfi(z); }

  /// Compute erfi(z) = -i erf(iz), the imaginary error function.
  double
  erfi(double x)
  { return Faddeeva::erfi(x); }

  /// Compute erfc(z) = 1 - erf(z), the complementary error function.
  std::complex<double>
  erfc(std::complex<double> z)
  { return Faddeeva::erfc(z); }

  /// Compute erfc(z) = 1 - erf(z), the complementary error function.
  double
  erfc(double x)
  { return Faddeeva::erfc(x); }

  /// Compute Dawson(z) = sqrt(pi)/2 exp(-z^2) erfi(z).
  std::complex<double>
  dawson(std::complex<double> z)
  { return Faddeeva::Dawson(z); }

  /// Compute Dawson(z) = sqrt(pi)/2 exp(-z^2) erfi(z).
  double
  dawson(double x)
  { return Faddeeva::Dawson(x); }

} // namespace faddeeva
