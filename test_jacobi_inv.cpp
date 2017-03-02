/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_inv test_jacobi_inv.cpp -L$HOME/bin/lib64 -lquadmath
./test_jacobi_inv > test_jacobi_inv.txt

g++ -std=gnu++17 -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_inv test_jacobi_inv.cpp -lquadmath
./test_jacobi_inv > test_jacobi_inv.txt
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
  jacobi_acnf(float __k, float __v)
  { return std::ellint_1(__k, std::acos(__v)); }

  /**
   * Return the inverse of the Jacobi elliptic cosine amplitude function
   * @f$ cn(k,u) @f$ (@see jacobi_cn) for <tt>long double</tt> arguments.
   *
   * @see jacobi_acn for more information.
   */
  long double
  jacobi_acnl(long double __k, long double __v)
  { return std::ellint_1(__k, std::acos(__v)); }

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
    __gnu_cxx::__promote_fp_t<_Tk, _Tv>
    jacobi_acn(_Tk __k, _Tv __v)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tk, _Tv>;
      return std::ellint_1<__type>(__k, std::acos(__v));
    }

  // Jacobi elliptic sine amplitude inverse functions.

  /**
   * Return the inverse of the Jacobi elliptic sine amplitude function
   * @f$ sn(k,u) @f$ (@see jacobi_sn) for @c float arguments.
   *
   * @see jacobi_asn for more information.
   */
  float
  jacobi_asnf(float __k, float __v)
  { return std::ellint_1(__k, std::asin(__v)); }

  /**
   * Return the inverse of the Jacobi elliptic sine amplitude function
   * @f$ sn(k,u) @f$ (@see jacobi_sn) for <tt>long double</tt> arguments.
   *
   * @see jacobi_asn for more information.
   */
  long double
  jacobi_asnl(long double __k, long double __v)
  { return std::ellint_1(__k, std::asin(__v)); }

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
    __gnu_cxx::__promote_fp_t<_Tk, _Tv>
    jacobi_asn(_Tk __k, _Tv __v)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tk, _Tv>;
      return std::ellint_1<__type>(__k, std::asin(__v));
    }

  // Jacobi elliptic delta amplitude inverse functions.

  /**
   * Return the inverse of the Jacobi elliptic delta amplitude function
   * @f$ dn(k,u) @f$ (@see jacobi_dn) for @c float arguments.
   *
   * @see jacobi_adn for more information.
   */
  float
  jacobi_adnf(float __k, float __v)
  { return std::ellint_1(__k, std::asin(std::sqrt(1.0F - __v * __v) / __k)); }

  /**
   * Return the inverse of the Jacobi elliptic delta amplitude function
   * @f$ dn(k,u) @f$ (@see jacobi_dn) for <tt>long double</tt> arguments.
   *
   * @see jacobi_adn for more information.
   */
  long double
  jacobi_adnl(long double __k, long double __v)
  { return std::ellint_1(__k, std::asin(std::sqrt(1.0L - __v * __v) / __k)); }

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
    __gnu_cxx::__promote_fp_t<_Tk, _Tv>
    jacobi_adn(_Tk __k, _Tv __v)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Tk, _Tv>;
      auto __root = std::sqrt(__type{1} - __v * __v);
      return std::ellint_1<__type>(__k, std::asin(__root / __k));
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
      auto v = __gnu_cxx::jacobi_cn(k, u);
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
      auto v = __gnu_cxx::jacobi_sn(k, u);
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
      auto v = __gnu_cxx::jacobi_dn(k, u);
      auto w = jacobi_adn(k, v);
      std::cout << std::setw(width) << u
		<< std::setw(width) << v
		<< std::setw(width) << std::abs(w - u)
		<< '\n';
    }
}
