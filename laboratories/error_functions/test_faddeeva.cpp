/**
 *
 */

#include <iostream>
#include <limits>

#include <wrap_faddeeva.h>

namespace emsr
{

  /**
   * Compute the Faddeeva or scaled complex complementary error function
   * for complex argument:
   * @f[
   *    w(z) = exp(-z^2) erfc(-iz)
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    faddeeva(std::complex<Tp> z);

  /**
   * Compute the Faddeeva or scaled complex complementary error function
   * for real argument:
   * @f[
   *    w(x) = exp(-x^2) erfc(-ix)
   * @f]
   */
  template<typename Tp>
    Tp
    faddeeva(Tp x);

  /**
   * Compute the scaled complementary error function of complex argument:
   * @f[
   *    erfcx(z) = exp^{z^2} erfc(z).
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    erfc_scaled(std::complex<Tp> z);

  /**
   * Compute the scaled complementary error function of real argument:
   * @f[
   *    erfcx(x) = exp^{x^2} erfc(x).
   * @f]
   */
  template<typename Tp>
    Tp
    erfc_scaled(Tp x);

  /**
   * Compute the error function of complex argument:
   * @f[
   *    erf(z) = \frac{2}{\sqrt{\pi}}\int_{0}^{z} e^{-t^2} dt
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    erf(std::complex<Tp> z);

  /**
   * Compute the error function of real argument:
   * @f[
   *    erf(x) = \frac{2}{\sqrt{\pi}}\int_{0}^{x} e^{-t^2} dt
   * @f]
   */
  template<typename Tp>
    Tp
    erf(Tp x);

  /**
   * Compute the imaginary error function for complex argument:
   * @f[
   *   erfi(z) = -i erf(iz)
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    erfi(std::complex<Tp> z);

  /**
   * Compute the imaginary error function of real argument:
   * @f[
   *    erfi(x) = -i erf(ix)
   * @f]
   */
  template<typename Tp>
    Tp
    erfi(Tp x);

  /**
   * Compute the complementary error function for complex argument:
   * @f[
   *    erfc(z) = 1 - erf(z) = \frac{2}{\sqrt{\pi}}\int_{z}^{\infty} e^{-t^2} dt
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    erfc(std::complex<Tp> z);

  /**
   * Compute the complementary error function for real argument:
   * @f[
   *   erfc(x) = 1 - erf(x) = \frac{2}{\sqrt{\pi}}\int_{x}^{\infty} e^{-t^2} dt
   * @f]
   */
  template<typename Tp>
    Tp
    erfc(Tp x);

  /**
   * Compute the Dawson integral for complex argument:
   * @f[
   *    D(z) = e^{-z^2} \int_{0}^{z} e^{t^2}dt
   *         = \frac{\sqrt{\pi}}{2} exp^{-z^2} erfi(z).
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    dawson(std::complex<Tp> z);

  /**
   * Compute the Dawson integral for real argument:
   * @f[
   *    D(x) = e^{-x^2} \int_{0}^{x} e^{t^2}dt
   *         = \frac{\sqrt{\pi}}{2} exp^{-x^2} erfi(x).
   * @f]
   */
  template<typename Tp>
    Tp
    dawson(Tp x);

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
  template<typename Tp>
    std::complex<Tp>
    voigt(Tp x, Tp t);

  /**
   * Compute the Voigt U function for real shape parameter t and argument x:
   * @f[
   *   U(x,t) = \frac{1}{\sqrt{4t}} \int_{-\infty}^{+\infty}
   *            \frac{e^{(x - y)^2 / (4t)}}{1 + y^2} dy
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    voigt_u(Tp x, Tp t);

  /**
   * Compute the Voigt V function for real shape parameter t and argument x:
   * @f[
   *   V(x,t) = \frac{1}{\sqrt{4t}} \int_{-\infty}^{+\infty}
   *            \frac{y e^{(x - y)^2 / (4t)}}{1 + y^2} dy
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    voigt_v(Tp x, Tp t);

} // namespace emsr

bool
not_zero(double x)
{
  return std::abs(x) > std::numeric_limits<double>::epsilon();
}

bool
is_zero(double x)
{
  return !not_zero(x);
}

/**
 * Return the Faddeeva function - a complex scaled error function.
 * Given a complex number z, this function returns
 * @f[
 *   w(z) = exp(-z^2)erfc(-iz)
 * @f]
 */
std::complex<double>
faddeeva_old(const std::complex<double>& z)
{
  constexpr double dconst = 1.12837916709551;

  double x = std::real(z), y = std::imag(z);

  int capn = 0;
  int nu = 0;
  double h = 0.0, h2 = 0.0, lambda = 0.0;
  if (y < 4.29 && x < 5.33)
    {
      auto s = 2.0 * (1.0 - y / 4.29) * std::sqrt(1.0 - x * x / 28.41);
      h = 1.6 * s;
      h2 = 2.0 * h;
      capn = int(12 + 23 * s);
      lambda = std::pow(h2, capn);
      nu = int(17 + 21 * s);
    }
  else
    nu = 16;
  bool b = is_zero(h) || is_zero(lambda);
  auto r1 = 0.0;
  auto r2 = 0.0;
  auto s1 = 0.0;
  auto s2 = 0.0;
  for (auto n = nu; n >= 0; --n)
    {
      const auto np1 = n + 1;
      auto t1 = y + h + np1 * r1;
      auto t2 = x - np1 * r2;
      const auto c = 0.5 / (t1 * t1 + t2 * t2);
      r1 = c * t1;
      r2 = c * t2;
      if (h > 0.0 && n <= capn)
	{
	  t1 = lambda + s1;
	  s1 = r1 * t1 - r2 * s2;
	  s2 = r2 * t1 + r1 * s2;
	  lambda /= h2;
	}
    }

  double re, im;
  if (is_zero(y))
    re = std::exp(-x * x);
  else
    {
      if (b)
	re = dconst * r1;
      else
	re = dconst * s1;
    }
  if (b)
    im = dconst * r2;
  else
    im = dconst * s2;

  return std::complex<double>(re, im);
}

int
main()
{
  for (int i = -100; i <= +100; ++i)
    {
      double x = 0.05 * i;
      for (int j = -100; j <= +100; ++j)
	{
	  double y = 0.05 * j;
	  std::complex<double> z(x, y);
	  std::complex<double> w = faddeeva_old(z);
	  std::cout << ' ' << x
		    << ' ' << y
		    << ' ' << std::real(w)
		    << ' ' << std::imag(w)
		    << '\n';
	}
      std::cout << '\n';
    }
}

