/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

  // Jacobi elliptic cosine amplitude inverse functions.

  /**
   * Return the inverse of the Jacobi elliptic cosine amplitude function
   * @f$ cn(k,u) @f$ (@see jacobi_cn) for @c float arguments.
   *
   * @see jacobi_acn for more information.
   */
  float
  jacobi_acnf(float k, float v)
  { return emsr::ellint_1(k, std::acos(v)); }

  /**
   * Return the inverse of the Jacobi elliptic cosine amplitude function
   * @f$ cn(k,u) @f$ (@see jacobi_cn) for <tt>long double</tt> arguments.
   *
   * @see jacobi_acn for more information.
   */
  long double
  jacobi_acnl(long double k, long double v)
  { return emsr::ellint_1(k, std::acos(v)); }

  /**
   * Return the inverse of the Jacobi elliptic cosine amplitude function
   * @f$ cn(k,u) @f$ (@see jacobi_cn) for real arguments.
   * The inverse is given by
   * @f[
   *    arccn(k,v) = F(k, arccos(v))
   * @f]
   * where @f$ F(k, \phi) @f$ is the Legendre elliptic function
   * of the first kind (@see ellint_1).
   */
  template<typename _Tk, typename _Tv>
    emsr::fp_promote_t<_Tk, _Tv>
    jacobi_acn(_Tk k, _Tv v)
    {
      using type = emsr::fp_promote_t<_Tk, _Tv>;
      return emsr::ellint_1<type>(k, std::acos(v));
    }

  // Jacobi elliptic sine amplitude inverse functions.

  /**
   * Return the inverse of the Jacobi elliptic sine amplitude function
   * @f$ sn(k,u) @f$ (@see jacobi_sn) for @c float arguments.
   *
   * @see jacobi_asn for more information.
   */
  float
  jacobi_asnf(float k, float v)
  { return emsr::ellint_1(k, std::asin(v)); }

  /**
   * Return the inverse of the Jacobi elliptic sine amplitude function
   * @f$ sn(k,u) @f$ (@see jacobi_sn) for <tt>long double</tt> arguments.
   *
   * @see jacobi_asn for more information.
   */
  long double
  jacobi_asnl(long double k, long double v)
  { return emsr::ellint_1(k, std::asin(v)); }

  /**
   * Return the inverse of the Jacobi elliptic sine amplitude function
   * @f$ sn(k,u) @f$ (@see jacobi_sn) for real arguments.
   * The inverse is given by
   * @f[
   *    arcsn(k,v) = F(k, arcsin(v))
   * @f]
   * where @f$ F(k, \phi) @f$ is the Legendre elliptic function
   * of the first kind (@see ellint_1).
   */
  template<typename _Tk, typename _Tv>
    emsr::fp_promote_t<_Tk, _Tv>
    jacobi_asn(_Tk k, _Tv v)
    {
      using type = emsr::fp_promote_t<_Tk, _Tv>;
      return emsr::ellint_1<type>(k, std::asin(v));
    }

  // Jacobi elliptic delta amplitude inverse functions.

  /**
   * Return the inverse of the Jacobi elliptic delta amplitude function
   * @f$ dn(k,u) @f$ (@see jacobi_dn) for @c float arguments.
   *
   * @see jacobi_adn for more information.
   */
  float
  jacobi_adnf(float k, float v)
  { return emsr::ellint_1(k, std::asin(std::sqrt(1.0F - v * v) / k)); }

  /**
   * Return the inverse of the Jacobi elliptic delta amplitude function
   * @f$ dn(k,u) @f$ (@see jacobi_dn) for <tt>long double</tt> arguments.
   *
   * @see jacobi_adn for more information.
   */
  long double
  jacobi_adnl(long double k, long double v)
  { return emsr::ellint_1(k, std::asin(std::sqrt(1.0L - v * v) / k)); }

  /**
   * Return the inverse of the Jacobi elliptic delta amplitude function
   * @f$ dn(k,u) @f$ (@see jacobi_dn) for real arguments.
   * The inverse is given by
   * @f[
   *    arcdn(k,v) = F(k, arcsin(\sqrt{1-v^2}/2))
   * @f]
   * where @f$ F(k, \phi) @f$ is the Legendre elliptic function
   * of the first kind (@see ellint_1).
   */
  template<typename _Tk, typename _Tv>
    emsr::fp_promote_t<_Tk, _Tv>
    jacobi_adn(_Tk k, _Tv v)
    {
      using type = emsr::fp_promote_t<_Tk, _Tv>;
      auto root = std::sqrt(type{1} - v * v);
      return emsr::ellint_1<type>(k, std::asin(root / k));
    }


int
main()
{
  using Tp = double;
  std::cout.precision(std::numeric_limits<Tp>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  auto k = 0.5;

  std::cout << "\n\nJacobi cn inverse\n";
  std::cout << "++++++++++++++++++++++++++++++\n";
  std::cout << std::setw(width) << "u"
	    << std::setw(width) << "v = cn(k,u)"
	    << std::setw(width) << "|acn(k,v) - u|"
	      << '\n';
  std::cout << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << '\n';
  for (auto iu = 0; iu < 150; ++iu)
    {
      auto u = Tp(iu * 0.01L);
      auto v = emsr::jacobi_cn(k, u);
      auto w = jacobi_acn(k, v);
      std::cout << std::setw(width) << u
		<< std::setw(width) << v
		<< std::setw(width) << std::abs(w - u)
		<< '\n';
    }

  std::cout << "\n\nJacobi sn inverse\n";
  std::cout << "++++++++++++++++++++++++++++++\n";
  std::cout << std::setw(width) << "u"
	    << std::setw(width) << "v = sn(k,u)"
	    << std::setw(width) << "|asn(k,v) - u|"
	      << '\n';
  std::cout << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << '\n';
  for (auto iu = 0; iu < 150; ++iu)
    {
      auto u = Tp(iu * 0.01L);
      auto v = emsr::jacobi_sn(k, u);
      auto w = jacobi_asn(k, v);
      std::cout << std::setw(width) << u
		<< std::setw(width) << v
		<< std::setw(width) << std::abs(w - u)
		<< '\n';
    }

  std::cout << "\n\nJacobi dn inverse\n";
  std::cout << "++++++++++++++++++++++++++++++\n";
  std::cout << std::setw(width) << "u"
	    << std::setw(width) << "v = dn(k,u)"
	    << std::setw(width) << "|adn(k,v) - u|"
	      << '\n';
  std::cout << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << std::setw(width) << "=============="
	    << '\n';
  for (auto iu = 0; iu < 150; ++iu)
    {
      auto u = Tp(iu * 0.01L);
      auto v = emsr::jacobi_dn(k, u);
      auto w = jacobi_adn(k, v);
      std::cout << std::setw(width) << u
		<< std::setw(width) << v
		<< std::setw(width) << std::abs(w - u)
		<< '\n';
    }
}
