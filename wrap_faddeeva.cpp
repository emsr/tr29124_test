
#include "wrap_faddeeva.h"

#include <ext/cmath>

#include "Faddeeva/Faddeeva.h"

namespace faddeeva
{

  /**
   * Compute the Faddeeva or scaled complex complementary error function
   * for complex argument:
   * @f[
   *    w(z) = exp(-z^2) erfc(-iz)
   * @f]
   */
  std::complex<double>
  faddeeva(std::complex<double> z)
  { return Faddeeva::w(z); }

  /**
   * Compute the Faddeeva or scaled complex complementary error function
   * for real argument:
   * @f[
   *    w(z) = exp(-z^2) erfc(-iz)
   * @f]
   */
  double
  faddeeva(double x)
  { return Faddeeva::w_im(x); }

  /**
   * Compute the scaled complementary error function of complex argument:
   * @f[
   *    erfcx(z) = exp^{z^2} erfc(z).
   * @f]
   */
  std::complex<double>
  erfc_scaled(std::complex<double> z)
  { return Faddeeva::erfcx(z); }

  /**
   * Compute the scaled complementary error function of real argument:
   * @f[
   *    erfcx(z) = exp^{z^2} erfc(z).
   * @f]
   */
  double
  erfc_scaled(double x)
  { return Faddeeva::erfcx(x); }

  /**
   * Compute the error function of complex argument:
   * @f[
   *    erf(z) = \frac{2}{\sqrt{\pi}}\int_{0}^{z} e^{-t^2} dt
   * @f]
   */
  std::complex<double>
  erf(std::complex<double> z)
  { return Faddeeva::erf(z); }

  /**
   * Compute the error function of real argument:
   * @f[
   *    erf(x) = \frac{2}{\sqrt{\pi}}\int_{0}^{x} e^{-t^2} dt
   * @f]
   */
  double
  erf(double x)
  { return Faddeeva::erf(x); }

  /**
   * Compute the imaginary error function for complex argument:
   * @f[
   *   erfi(z) = -i erf(iz)
   * @f]
   */
  std::complex<double>
  erfi(std::complex<double> z)
  { return Faddeeva::erfi(z); }

  /**
   * Compute the imaginary error function of real argument:
   * @f[
   *    erfi(x) = -i erf(ix)
   * @f]
   */
  double
  erfi(double x)
  { return Faddeeva::erfi(x); }

  /**
   * Compute the complementary error function for complex argument:
   * @f[
   *    erfc(z) = 1 - erf(z) = \frac{2}{\sqrt{\pi}}\int_{z}^{\infty} e^{-t^2} dt
   * @f]
   */
  std::complex<double>
  erfc(std::complex<double> z)
  { return Faddeeva::erfc(z); }

  /**
   * Compute the complementary error function for real argument:
   * @f[
   *   erfc(z) = 1 - erf(z)
   * @f]
   */
  double
  erfc(double x)
  { return Faddeeva::erfc(x); }

  /**
   * Compute the Dawson integral for complex argument:
   * @f[
   *    D(z) = e^{-z^2} \int_{0}^{z} e^{t^2}dt
   *         = \frac{\sqrt{\pi}}{2} exp^{-z^2} erfi(z).
   * @f]
   */
  std::complex<double>
  dawson(std::complex<double> z)
  { return Faddeeva::Dawson(z); }

  /**
   * Compute the Dawson integral for real argument:
   * @f[
   *    D(x) = e^{-x^2} \int_{0}^{x} e^{t^2}dt
   *         = \frac{\sqrt{\pi}}{2} exp^{-x^2} erfi(x).
   * @f]
   */
  double
  dawson(double x)
  { return Faddeeva::Dawson(x); }

  /**
   * Compute Voigt function for real shape parameter t and argument x:
   * @f[
   *   W(x,t) = U(x,t) + iV(x,t)
   * @f]
   * where
   * @f[
   *   U(x,t) = \frac{1}{\sqrt{4t}} \int_{-\infty}^{+\infty}
   *            \frac{e^{(x - y)^2 / (4t)}}{1 + y^2} dy
   * @f]
   * and
   * @f[
   *   V(x,t) = \frac{1}{\sqrt{4t}} \int_{-\infty}^{+\infty}
   *            \frac{y e^{(x - y)^2 / (4t)}}{1 + y^2} dy
   * @f]
   */
  std::complex<double>
  voigt(double x, double t)
  {
    const auto _S_sqrt_pi = 2 * __gnu_cxx::__math_constants<double>::__root_pi_div_2;
    auto s = 1.0 / (2.0 * std::sqrt(t));
    std::complex z(s, -s * x);
    return _S_sqrt_pi * s * faddeeva(z);
  }

  /**
   * Compute the Voigt U function for real shape parameter t and argument x:
   * @f[
   *   U(x,t) = \frac{1}{\sqrt{4t}} \int_{-\infty}^{+\infty}
   *            \frac{e^{(x - y)^2 / (4t)}}{1 + y^2} dy
   * @f]
   */
  std::complex<double>
  voigt_u(double x, double t)
  { return std::real(voigt(x, t)); }

  /**
   * Compute the Voigt V function for real shape parameter t and argument x:
   * @f[
   *   V(x,t) = \frac{1}{\sqrt{4t}} \int_{-\infty}^{+\infty}
   *            \frac{y e^{(x - y)^2 / (4t)}}{1 + y^2} dy
   * @f]
   */
  std::complex<double>
  voigt_v(double x, double t)
  { return std::imag(voigt(x, t)); }

} // namespace faddeeva
