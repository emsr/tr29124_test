
#include <wrap_faddeeva.h>

namespace __gnu_cxx
{

  /**
   * Compute the Faddeeva or scaled complex complementary error function
   * for complex argument:
   * @f[
   *    w(z) = exp(-z^2) erfc(-iz)
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    faddeeva(std::complex<_Tp> __z);

  /**
   * Compute the Faddeeva or scaled complex complementary error function
   * for real argument:
   * @f[
   *    w(z) = exp(-z^2) erfc(-iz)
   * @f]
   */
  template<typename _Tp>
    _Tp
    faddeeva(_Tp __x);

  /**
   * Compute the scaled complementary error function of complex argument:
   * @f[
   *    erfcx(z) = exp^{z^2} erfc(z).
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    erfc_scaled(std::complex<_Tp> __z);

  /**
   * Compute the scaled complementary error function of real argument:
   * @f[
   *    erfcx(z) = exp^{z^2} erfc(z).
   * @f]
   */
  template<typename _Tp>
    _Tp
    erfc_scaled(_Tp __x);

  /**
   * Compute the error function of complex argument:
   * @f[
   *    erf(z) = \frac{2}{\sqrt{\pi}}\int_{0}^{z} e^{-t^2} dt
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    erf(std::complex<_Tp> __z);

  /**
   * Compute the error function of real argument:
   * @f[
   *    erf(x) = \frac{2}{\sqrt{\pi}}\int_{0}^{x} e^{-t^2} dt
   * @f]
   */
  template<typename _Tp>
    _Tp
    erf(_Tp __x);

  /**
   * Compute the imaginary error function for complex argument:
   * @f[
   *   erfi(z) = -i erf(iz)
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    erfi(std::complex<_Tp> __z);

  /**
   * Compute the imaginary error function of real argument:
   * @f[
   *    erfi(x) = -i erf(ix)
   * @f]
   */
  template<typename _Tp>
    _Tp
    erfi(_Tp __x);

  /**
   * Compute the complementary error function for complex argument:
   * @f[
   *    erfc(z) = 1 - erf(z) = \frac{2}{\sqrt{\pi}}\int_{z}^{\infty} e^{-t^2} dt
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    erfc(std::complex<_Tp> __z);

  /**
   * Compute the complementary error function for real argument:
   * @f[
   *   erfc(x) = 1 - erf(x) = \frac{2}{\sqrt{\pi}}\int_{x}^{\infty} e^{-t^2} dt
   * @f]
   */
  template<typename _Tp>
    _Tp
    erfc(_Tp __x);

  /**
   * Compute the Dawson integral for complex argument:
   * @f[
   *    D(z) = e^{-z^2} \int_{0}^{z} e^{t^2}dt
   *         = \frac{\sqrt{\pi}}{2} exp^{-z^2} erfi(z).
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    dawson(std::complex<_Tp> __z);

  /**
   * Compute the Dawson integral for real argument:
   * @f[
   *    D(x) = e^{-x^2} \int_{0}^{x} e^{t^2}dt
   *         = \frac{\sqrt{\pi}}{2} exp^{-x^2} erfi(x).
   * @f]
   */
  template<typename _Tp>
    _Tp
    dawson(_Tp __x);

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
  template<typename _Tp>
    std::complex<_Tp>
    voigt(_Tp __x, _Tp __t);

  /**
   * Compute the Voigt U function for real shape parameter t and argument x:
   * @f[
   *   U(x,t) = \frac{1}{\sqrt{4t}} \int_{-\infty}^{+\infty}
   *            \frac{e^{(x - y)^2 / (4t)}}{1 + y^2} dy
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    voigt_u(_Tp __x, _Tp __t);

  /**
   * Compute the Voigt V function for real shape parameter t and argument x:
   * @f[
   *   V(x,t) = \frac{1}{\sqrt{4t}} \int_{-\infty}^{+\infty}
   *            \frac{y e^{(x - y)^2 / (4t)}}{1 + y^2} dy
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    voigt_v(_Tp __x, _Tp __t);

} // namespace __gnu_cxx

int
main()
{
}

