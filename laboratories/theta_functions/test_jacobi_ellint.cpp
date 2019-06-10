/**
 *
 */

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>

#include <wrap_boost.h>
#include <wrap_gsl.h>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return the Jacobi amplitude function of modulus k and argument u.
   * @f[
   *    am(k, u) = arcsin(sn(k, u))
   * @f]
   * @see Consult jacobi_sn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_am(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__am();
    }

  /**
   * Return the Jacobi elliptic sc function of modulus k and argument u.
   * @f[
   *    sc(k, u) = \frac{sn(k,u)}{cn(k,u)}
   * @f]
   * @see Consult jacobi_sn and jacobi_cn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_sc(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__sc();
    }

  /**
   * Return the Jacobi elliptic  function of modulus k and argument u.
   * @f[
   *    sd(k, u) = \frac{sn(k,u)}{dn(k,u)}
   * @f]
   * @see Consult jacobi_sn and jacobi_dn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_sd(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__sd();
    }

  /**
   * Return the Jacobi elliptic cd function of modulus k and argument u.
   * @f[
   *    cd(k, u) = \frac{cn(k,u)}{dn(k,u)}
   * @f]
   * @see Consult jacobi_cn and jacobi_dn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_cd(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__cd();
    }

  /**
   * Return the Jacobi elliptic cs function of modulus k and argument u.
   * @f[
   *    cs(k, u) = \frac{cn(k,u)}{sn(k,u)}
   * @f]
   * @see Consult jacobi_cn and jacobi_sn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_cs(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__cs();
    }

  /**
   * Return the Jacobi elliptic dc function of modulus k and argument u.
   * @f[
   *    dc(k, u) = \frac{dn(k,u)}{cn(k,u)}
   * @f]
   * @see Consult jacobi_dn and jacobi_cn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_dc(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__dc();
    }

  /**
   * Return the Jacobi elliptic ds function of modulus k and argument u.
   * @f[
   *    ds(k, u) = \frac{dn(k,u)}{sn(k,u)}
   * @f]
   * @see Consult jacobi_dn and jacobi_sn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_ds(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__ds();
    }

  /**
   * Return the Jacobi elliptic nc function of modulus k and argument u.
   * @f[
   *    nc(k, u) = \frac{1}{cn(k,u)}
   * @f]
   * @see Consult jacobi_cn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_nc(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__nc();
    }

  /**
   * Return the Jacobi elliptic nd function of modulus k and argument u.
   * @f[
   *    nd(k, u) = \frac{1}{dn(k,u)}
   * @f]
   * @see Consult jacobi_dn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_nd(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__nd();
    }

  /**
   * Return the Jacobi elliptic ns function of modulus k and argument u.
   * @f[
   *    ns(k, u) = \frac{1}{sn(k,u)}
   * @f]
   * @see Consult jacobi_sn for more details.
   */
  template<typename _Kp, typename _Up>
    inline __gnu_cxx::fp_promote_t<_Kp, _Up>
    jacobi_ns(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::fp_promote_t<_Kp, _Up>;
      return std::__detail::__jacobi_ellint<__type>(__k, __u).__ns();
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

void
test_gsl()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  std::vector<double> kvals{0.0, 0.5, 0.75, 0.95, 1.0};

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic sine amplitude function sn(k,u)\n";
  for (auto k : kvals)
    for (int i = -50; i <= 50; ++i)
      {
	auto u = 0.2 * i;
	auto fgnu = __gnu_cxx::jacobi_sn(k, u);
	auto fgsl = gsl::jacobi_sn(k, u);
	std::cout << ' ' << std::setw(width) << u;
	  std::cout << ' ' << std::setw(width) << fgnu
		    << ' ' << std::setw(width) << fgsl
		    << ' ' << fgnu - fgsl;
	std::cout << '\n';
      }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic cosine amplitude function cn(k,u)\n";
  for (auto k : kvals)
    for (int i = -50; i <= 50; ++i)
      {
	auto u = 0.2 * i;
	auto fgnu = __gnu_cxx::jacobi_cn(k, u);
	auto fgsl = gsl::jacobi_cn(k, u);
	std::cout << ' ' << std::setw(width) << u;
	  std::cout << ' ' << std::setw(width) << fgnu
		    << ' ' << std::setw(width) << fgsl
		    << ' ' << fgnu - fgsl;
	std::cout << '\n';
      }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic delta amplitude function dn(k,u)\n";
  for (auto k : kvals)
    for (int i = -50; i <= 50; ++i)
      {
	auto u = 0.2 * i;
	auto fgnu = __gnu_cxx::jacobi_dn(k, u);
	auto fgsl = gsl::jacobi_dn(k, u);
	std::cout << ' ' << std::setw(width) << u;
	  std::cout << ' ' << std::setw(width) << fgnu
		    << ' ' << std::setw(width) << fgsl
		    << ' ' << fgnu - fgsl;
	std::cout << '\n';
      }
}

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  std::vector<double> kvals{0.0, 0.5, 0.75, 0.95, 1.0};

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic sine amplitude function sn(k,u)\n";
  for (int i = -200; i <= 200; ++i)
    {
      auto u = 0.05 * i;
      std::cout << ' ' << std::setw(width) << u;
      for (auto k : kvals)
	std::cout << ' ' << std::setw(width) << __gnu_cxx::jacobi_sn(k, u);
      std::cout << '\n';
    }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic cosine amplitude function cn(k,u)\n";
  for (int i = -200; i <= 200; ++i)
    {
      auto u = 0.05 * i;
      std::cout << ' ' << std::setw(width) << u;
      for (auto k : kvals)
	std::cout << ' ' << std::setw(width) << __gnu_cxx::jacobi_cn(k, u);
      std::cout << '\n';
    }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic delta amplitude function dn(k,u)\n";
  for (int i = -200; i <= 200; ++i)
    {
      auto u = 0.05 * i;
      std::cout << ' ' << std::setw(width) << u;
      for (auto k : kvals)
	std::cout << ' ' << std::setw(width) << __gnu_cxx::jacobi_dn(k, u);
      std::cout << '\n';
    }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic amplitude function am(k,u)\n";
  for (int i = -200; i <= 200; ++i)
    {
      auto u = 0.05 * i;
      std::cout << ' ' << std::setw(width) << u;
      for (auto k : kvals)
	std::cout << ' ' << std::setw(width) << __gnu_cxx::jacobi_am(k, u);
      std::cout << '\n';
    }

  test_gsl();
}
