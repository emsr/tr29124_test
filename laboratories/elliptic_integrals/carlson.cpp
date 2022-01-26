/**
 *
 */

#include <cmath>
#include <iostream>
#include <complex>
#include <limits>

// From bits/specfun_util.h

namespace emsr
{
namespace detail
{

  /// A class to abstract the scalar data type in a generic way.
  template<typename _Tp>
    struct num_traits
    {
      using value_type = _Tp;
    };

  template<>
    template<typename _Tp>
      struct num_traits<std::complex<_Tp>>
      {
	using value_type = typename std::complex<_Tp>::value_type;
      };

  template<typename _Tp>
    using num_traits_t = typename num_traits<_Tp>::value_type;

  /// A class to encapsulate type dependent floating point
  /// constants.  Not everything will be able to be expressed as
  /// type logic.
  template<typename _Tp>
    struct floating_point_constant
    {
      static const _Tp value;
    };

  /// A structure for numeric constants.
  template<typename _Tp>
    struct numeric_constants
    {
      ///  Constant @f$ \pi @f$.
      static constexpr _Tp pi() noexcept
      { return static_cast<_Tp>(3.1415926535897932384626433832795029L); }
      ///  Constant @f$ \pi / 2 @f$.
      static constexpr _Tp pi_2() noexcept
      { return static_cast<_Tp>(1.5707963267948966192313216916397514L); }
      ///  Constant @f$ \pi / 3 @f$.
      static constexpr _Tp pi_3() noexcept
      { return static_cast<_Tp>(1.0471975511965977461542144610931676L); }
      ///  Constant @f$ \pi / 4 @f$.
      static constexpr _Tp pi_4() noexcept
      { return static_cast<_Tp>(0.7853981633974483096156608458198757L); }
      ///  Constant @f$ 1 / \pi @f$.
      static constexpr _Tp __1_pi() noexcept
      { return static_cast<_Tp>(0.3183098861837906715377675267450287L); }
      ///  Constant @f$ 2 / \sqrt(\pi) @f$.
      static constexpr _Tp __2_sqrtpi() noexcept
      { return static_cast<_Tp>(1.1283791670955125738961589031215452L); }
      ///  Constant @f$ \sqrt(2) @f$.
      static constexpr _Tp sqrt2() noexcept
      { return static_cast<_Tp>(1.4142135623730950488016887242096981L); }
      ///  Constant @f$ \sqrt(3) @f$.
      static constexpr _Tp sqrt3() noexcept
      { return static_cast<_Tp>(1.7320508075688772935274463415058723L); }
      ///  Constant @f$ \sqrt(\pi/2) @f$.
      static constexpr _Tp sqrtpio2() noexcept
      { return static_cast<_Tp>(1.2533141373155002512078826424055226L); }
      ///  Constant @f$ 1 / sqrt(2) @f$.
      static constexpr _Tp sqrt1_2() noexcept
      { return static_cast<_Tp>(0.7071067811865475244008443621048490L); }
      ///  Constant @f$ \log(\pi) @f$.
      static constexpr _Tp lnpi() noexcept
      { return static_cast<_Tp>(1.1447298858494001741434273513530587L); }
      ///  Constant Euler's constant @f$ \gamma_E @f$.
      static constexpr _Tp gamma_e() noexcept
      { return static_cast<_Tp>(0.5772156649015328606065120900824024L); }
      ///  Constant Euler-Mascheroni @f$ e @f$
      static constexpr _Tp euler() noexcept
      { return static_cast<_Tp>(2.7182818284590452353602874713526625L); }
    };

  template<typename _Tp>
    inline bool
    isnan(const std::complex<_Tp>& x)
    { return std::isnan(std::real(x)) || std::isnan(std::imag(x)); }


  template<typename _Tp>
    inline _Tp
    quiet_NaN(_Tp)
    { return std::numeric_limits<_Tp>::quiet_NaN(); }

  template<typename _Tp>
    inline std::complex<_Tp>
    quiet_NaN(std::complex<_Tp>)
    {
      auto nan = std::numeric_limits<_Tp>::quiet_NaN();
      return std::complex<_Tp>(nan, nan);
    }


  template<typename _Tp>
    inline _Tp
    infinity(_Tp)
    { return std::numeric_limits<_Tp>::infinity(); }

  template<typename _Tp>
    inline std::complex<_Tp>
    infinity(std::complex<_Tp>)
    {
      auto inf = std::numeric_limits<_Tp>::infinity();
      return std::complex<_Tp>(inf, inf);
    }

} // namespace detail
} // namespace emsr

// End bits/specfun_util.h

namespace std
{
namespace detail
{

  /**
   *   @brief Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   *          of the first kind.
   *
   *   The Carlson elliptic function of the first kind is defined by:
   *   @f[
   *       R_F(x,y,z) = \frac{1}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}}
   *   @f]
   *
   *   @param  x  The first of three symmetric arguments.
   *   @param  y  The second of three symmetric arguments.
   *   @param  z  The third of three symmetric arguments.
   *   @return  The Carlson elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    ellint_rf(_Tp x, _Tp y, _Tp z)
    {
      using _Val = num_traits_t<_Tp>;
      const _Val r = std::numeric_limits<_Val>::epsilon();
      _Tp xt = x;
      _Tp yt = y;
      _Tp zt = z;
      _Tp a0 = (x + y + z) / _Val(3);
      _Val q = std::pow( _Val(3) * r, -_Val(1) / _Val(6) )
	       * std::max(std::abs(a0 - z),
			  std::max(std::abs(a0 - x),
				   std::abs(a0 - y)));
      _Tp a = a0;
      _Val f = _Val(1);

      while (true)
	{
	  _Tp lambda = std::sqrt(xt) * std::sqrt(yt)
		       + std::sqrt(yt) * std::sqrt(zt)
		       + std::sqrt(zt) * std::sqrt(xt);
	  a = (a + lambda) / _Val(4);
	  xt = (xt + lambda) / _Val(4);
	  yt = (yt + lambda) / _Val(4);
	  zt = (zt + lambda) / _Val(4);
	  f *= _Val(4);
	  if (q < f * std::abs(a))
	    {
	      _Tp xf = (a0 - x) / (f * a);
	      _Tp yf = (a0 - y) / (f * a);
	      _Tp zf = -(xf + yf);
	      _Tp e2 = xf * yf - zf * zf;
	      _Tp e3 = xf * yf * zf;
	      return (_Val(1)
		    - e2 / _Val(10)
		    + e3 / _Val(14)
		    + e2 * e2 / _Val(24)
		    - _Val(3) * e2 * e3 / _Val(44)) / std::sqrt(a);
	    }
	}

      return _Tp(0);
    }

  /**
   *   @brief  Return the Carlson elliptic function
   *           @f$ R_C(x,y) = R_F(x,y,y) @f$ where @f$ R_F(x,y,z) @f$
   *           is the Carlson elliptic function of the first kind.
   *
   *   The Carlson elliptic function is defined by:
   *   @f[
   *       R_C(x,y) = \frac{1}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)}
   *   @f]
   *
   *   Based on Carlson's algorithms:
   *   -  B. C. Carlson Numer. Math. 33, 1 (1979)
   *   -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   *   -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *      by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   *   @param  x  The first argument.
   *   @param  y  The second argument.
   *   @return  The Carlson elliptic function.
   */
  template<typename _Tp>
    _Tp
    ellint_rc(_Tp x, _Tp y)
    {
      using _Val = num_traits_t<_Tp>;
      if (std::imag(y) == _Val(0) && std::real(y) < _Val(0))
	return std::sqrt(x / (x - y)) * ellint_rc(x - y, -y);
      const _Val r = std::numeric_limits<_Val>::epsilon();
      _Tp xt = x;
      _Tp yt = y;
      _Tp a0 = (x + _Val(2) * y) / _Val(3);
      _Val q = std::pow( _Val(3) * r, -_Val(1) / _Val(8) )
	       * std::abs(a0 - x);
      _Tp a = a0;
      _Val f = _Val(1);

      while (true)
	{
	  _Tp lambda = _Val(2) * std::sqrt(xt) * std::sqrt(yt) + yt;
	  a = (a + lambda) / _Val(4);
	  xt = (xt + lambda) / _Val(4);
	  yt = (yt + lambda) / _Val(4);
	  f *= _Val(4);
	  if (q < f * std::abs(a))
	    {
	      _Tp s = (y - a0) / (f * a);
	      return (_Val(1) + s * s * (_Val(3) / _Val(10)
		    + s * (_Val(1) / _Val(7)
		    + s * (_Val(3) / _Val(8)
		    + s * (_Val(9) / _Val(22)
		    + s * (_Val(159) / _Val(208)
		    + s * (_Val(9) / _Val(8)))))))) / std::sqrt(a);
	    }
	}

      return _Tp(0);
    }

  /**
   *   @brief  Return the Carlson elliptic function @f$ R_J(x,y,z,p) @f$
   *           of the third kind.
   *
   *   The Carlson elliptic function of the third kind is defined by:
   *   @f[
   *       R_J(x,y,z,p) = \frac{3}{2} \int_0^\infty
   *       \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}(t + p)}
   *   @f]
   *
   *   Based on Carlson's algorithms:
   *   -  B. C. Carlson Numer. Math. 33, 1 (1979)
   *   -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   *   -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *      by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   *   @param  x  The first of three symmetric arguments.
   *   @param  y  The second of three symmetric arguments.
   *   @param  z  The third of three symmetric arguments.
   *   @param  p  The fourth argument.
   *   @return  The Carlson elliptic function of the fourth kind.
   */
  template<typename _Tp>
    _Tp
    ellint_rj(_Tp x, _Tp y, _Tp z, _Tp p)
    {
      using _Val = num_traits_t<_Tp>;
      const _Val r = std::numeric_limits<_Val>::epsilon();
      _Tp xt = x;
      _Tp yt = y;
      _Tp zt = z;
      _Tp pt = p;
      _Tp a0 = (x + y + z + _Val(2) * p) / _Val(5);
      _Tp delta = (p - x) * (p - y) * (p - z);
      _Val q = std::pow(r / _Val(4), -_Val(1) / _Val(6))
	       * std::max(std::abs(a0 - z),
			  std::max(std::abs(a0 - x),
				   std::max(std::abs(a0 - y), std::abs(a0 - p))));
      _Tp a = a0;
      _Val f = _Val(1);
      _Val fe = _Val(1);
      _Tp sum = _Tp();

      while (true)
	{
	  _Tp xroot = std::sqrt(xt);
	  _Tp yroot = std::sqrt(yt);
	  _Tp zroot = std::sqrt(zt);
	  _Tp proot = std::sqrt(pt);
	  _Tp lambda = xroot * yroot
		       + yroot * zroot
		       + zroot * xroot;
	  a = (a + lambda) / _Val(4);
	  xt = (xt + lambda) / _Val(4);
	  yt = (yt + lambda) / _Val(4);
	  zt = (zt + lambda) / _Val(4);
	  pt = (pt + lambda) / _Val(4);
	  _Tp d = (proot + xroot) * (proot + yroot) * (proot + zroot);
	  _Tp e = delta / (fe * d * d);
	  sum += ellint_rc(_Tp(1), _Tp(1) + e) / (f * d);
	  f *= _Val(4);
	  fe *= _Val(64);
	  if (q < f * std::abs(a))
	    {
	      _Tp xf = (a0 - x) / (f * a);
	      _Tp yf = (a0 - y) / (f * a);
	      _Tp zf = (a0 - z) / (f * a);
	      _Tp xyz = xf * yf * zf;
	      _Tp pf = -(xf + yf + zf) / _Val(2);
	      _Tp pp = pf * pf;
	      _Tp ppp = pp * pf;
	      _Tp e2 = xf * yf + yf * zf + zf * xf - _Val(3) * pp;
	      _Tp e3 = xyz + _Val(2) * e2 * pf + _Tp(4) * ppp;
	      _Tp e4 = (_Val(2) * xyz + e2 * pf + _Val(3) * ppp) * pf;
	      _Tp e5 = xyz * pp;
	      return (_Val(1) - _Val(3) * e2 / _Val(14)
			      + e3 / _Val(6)
			      + _Val(9) * e2 * e2 / _Val(88)
			      - _Val(3) * e4 / _Val(22)
			      - _Val(9) * e2 * e3 / _Val(52)
			      + _Val(3) * e5 / _Val(26)) / f / a / std::sqrt(a)
			      + _Val(6) * sum;
	    }
	}

      return _Tp(0);
    }

  /**
   *   @brief  Return the Carlson elliptic function of the second kind
   *           @f$ R_D(x,y,z) = R_J(x,y,z,z) @f$ where
   *           @f$ R_J(x,y,z,p) @f$ is the Carlson elliptic function
   *           of the third kind.
   *
   *   The Carlson elliptic function of the second kind is defined by:
   *   @f[
   *       R_D(x,y,z) = \frac{3}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{3/2}}
   *   @f]
   *
   *   Based on Carlson's algorithms:
   *   -  B. C. Carlson Numer. Math. 33, 1 (1979)
   *   -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   *   -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *      by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   *   @param  x  The first of two symmetric arguments.
   *   @param  y  The second of two symmetric arguments.
   *   @param  z  The third argument.
   *   @return  The Carlson elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    ellint_rd(_Tp x, _Tp y, _Tp z)
    {
      using _Val = num_traits_t<_Tp>;
      const _Val r = std::numeric_limits<_Val>::epsilon();
      _Tp xt = x;
      _Tp yt = y;
      _Tp zt = z;
      _Tp a0 = (x + y + _Val(3) * z) / _Val(5);
      _Val q = std::pow(r / _Val(4), -_Val(1) / _Val(6))
	       * std::max(std::abs(a0 - z),
			  std::max(std::abs(a0 - x),
			  std::abs(a0 - y)));
      _Tp a = a0;
      _Val f = _Val(1);
      _Tp sum = _Tp();

      while (true)
	{
	  _Tp lambda = std::sqrt(xt) * std::sqrt(yt)
		       + std::sqrt(yt) * std::sqrt(zt)
		       + std::sqrt(zt) * std::sqrt(xt);
	  sum += _Val(1) / f / std::sqrt(zt) / (zt + lambda);
	  a = (a + lambda) / _Val(4);
	  xt = (xt + lambda) / _Val(4);
	  yt = (yt + lambda) / _Val(4);
	  zt = (zt + lambda) / _Val(4);
	  f *= _Val(4);
	  if (q < f * std::abs(a))
	    {
	      _Tp xf = (a0 - x) / (f * a);
	      _Tp yf = (a0 - y) / (f * a);
	      _Tp zf = -(xf + yf) / _Val(3);
	      _Tp zz = zf * zf;
	      _Tp xy = xf * yf;
	      _Tp e2 = xy - _Val(6) * zz;
	      _Tp e3 = (_Val(3) * xy - _Val(8) * zz) * zf;
	      _Tp e4 = _Val(3) * (xy - zz) * zz;
	      _Tp e5 = xy * zf * zz;
	      return (_Val(1)
		    - _Val(3) * e2 / _Val(14)
		    + e3 / _Val(6)
		    + _Val(9) * e2 * e2 / _Val(88)
		    - _Val(3) * e4 / _Val(22)
		    - _Val(9) * e2 * e3 / _Val(52)
		    + _Val(3) * e5 / _Val(26)) / f / a / std::sqrt(a)
		    + _Val(3) * sum;
	    }
	}

      return _Tp(0);
    }

  template<typename _Tp>
    _Tp
    comp_ellint_rf(_Tp x, _Tp y)
    {
      using _Val = num_traits_t<_Tp>;
      const _Val r = std::numeric_limits<_Val>::epsilon();
      const _Val tolfact = _Val(2.7L) * std::sqrt(r);
      x = std::sqrt(x);
      y = std::sqrt(y);
      while (true)
	{
	  _Tp xt = x;
	  x = (x + y) / _Tp(2);
	  y = std::sqrt(xt * y);
	  if (std::abs(x - y) < tolfact * std::abs(x))
	    return _Val(numeric_constants<_Val>::pi()) / (x + y);
	}
    }


  template<typename _Tp>
    _Tp
    comp_ellint_rg(_Tp x, _Tp y);

  /**
   *   @brief  Return the symmetric Carlson elliptic function of the second kind
   *           @f$ R_G(x,y,z) @f$.
   *
   *   The Carlson symmetric elliptic function of the second kind is defined by:
   *   @f[
   *       R_G(x,y,z) = \frac{1}{4} \int_0^\infty
   *                 dt t [(t + x)(t + y)(t + z)]^{-1/2}
   *                 (\frac{x}{t + x} + \frac{y}{t + y} + \frac{z}{t + z})
   *   @f]
   *
   *   Based on Carlson's algorithms:
   *   -  B. C. Carlson Numer. Math. 33, 1 (1979)
   *   -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   *   -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *      by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   *   @param  x  The first of three symmetric arguments.
   *   @param  y  The second of three symmetric arguments.
   *   @param  z  The third of three symmetric arguments.
   *   @return  The Carlson symmetric elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    ellint_rg(_Tp x, _Tp y, _Tp z)
    {
      using _Val = num_traits_t<_Tp>;
      if (z == _Tp())
	{
	  if (x == _Tp())
	    return std::sqrt(y);
	  else if (y == _Tp())
	    return std::sqrt(x);
	  else
	    return comp_ellint_rg(x, y);
	}
      else if (x == _Tp())
	{
	  if (y == _Tp())
	    return std::sqrt(z);
	  else if (z == _Tp())
	    return std::sqrt(y);
	  else
	    return comp_ellint_rg(y, z);
	}
      else if (y == _Tp())
	{
	  if (z == _Tp())
	    return std::sqrt(x);
	  else if (x == _Tp())
	    return std::sqrt(z);
	  else
	    return comp_ellint_rg(z, x);
	}
      else
	return (z * ellint_rf(x, y, z)
	      - (x - z) * (y - z) * ellint_rd(x, y, z) / _Val(3)
	      + std::sqrt(x * y / z)) / _Val(2);
    }

  template<typename _Tp>
    _Tp
    comp_ellint_rg(_Tp x, _Tp y)
    {
      using _Val = num_traits_t<_Tp>;
      const _Val r = std::numeric_limits<_Val>::epsilon();
      const _Val tolfact = _Val(2.7L) * std::sqrt(r);
      _Tp xt = std::sqrt(x);
      _Tp yt = std::sqrt(y);
      const _Tp a = (xt + yt) / _Val(2);
      _Tp sum = _Tp();
      _Val sf = _Val(1) / _Val(2);
      while (true)
	{
	  _Tp xtt = xt;
	  xt = (xt + yt) / _Tp(2);
	  yt = std::sqrt(xtt * yt);
	  _Tp del = xt - yt;
	  if (std::abs(del) < tolfact * std::abs(xt))
	    return (a * a - sum) * _Val(numeric_constants<_Val>::pi()) / (xt + yt) / _Val(2);
	  sum += sf * del * del;
	  sf *= _Val(2);
	}
    }

  /**
   *   @brief  Return the complete elliptic integral of the first kind
   *           @f$ K(k) @f$ using the Carlson formulation.
   *
   *   The complete elliptic integral of the first kind is defined as
   *   @f[
   *     K(k) = F(k,\pi/2) = \int_0^{\pi/2}\frac{d\theta}
   *                                           {\sqrt{1 - k^2 sin^2\theta}}
   *   @f]
   *   where @f$ F(k,\phi) @f$ is the incomplete elliptic integral of the
   *   first kind.
   *
   *   @param  k  The argument of the complete elliptic function.
   *   @return  The complete elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    comp_ellint_1(const _Tp k)
    {
      using _Val = num_traits_t<_Tp>;
      if (std::isnan(k))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(k) >= _Val(1))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return ellint_rf(_Tp(0), _Tp(1) - k * k, _Tp(1));
    }

  /**
   *   @brief  Return the incomplete elliptic integral of the first kind
   *           @f$ F(k,\phi) @f$ using the Carlson formulation.
   *
   *   The incomplete elliptic integral of the first kind is defined as
   *   @f[
   *     F(k,\phi) = \int_0^{\phi}\frac{d\theta}
   *                                   {\sqrt{1 - k^2 sin^2\theta}}
   *   @f]
   *
   *   @param  k  The argument of the elliptic function.
   *   @param  phi  The integral limit argument of the elliptic function.
   *   @return  The elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    ellint_1(const _Tp k, const _Tp phi)
    {
      using _Val = num_traits_t<_Tp>;
      if (std::isnan(k) || std::isnan(phi))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(k) > _Val(1))
	throw std::domain_error("Bad argument in ellint_1.");
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int n = std::floor(phi / numeric_constants<_Tp>::pi()
				   + _Tp(0.5L));
	  const _Tp phi_red = phi
			      - n * numeric_constants<_Tp>::pi();

	  const _Tp s = std::sin(phi_red);
	  const _Tp c = std::cos(phi_red);

	  const _Tp F = s
			* ellint_rf(c * c,
				_Tp(1) - k * k * s * s, _Tp(1));

	  if (n == 0)
	    return F;
	  else
	    return F + _Tp(2) * n * comp_ellint_1(k);
	}
    }

  /**
   *   @brief  Return the complete elliptic integral of the second kind
   *           @f$ E(k) @f$ using the Carlson formulation.
   *
   *   The complete elliptic integral of the second kind is defined as
   *   @f[
   *     E(k,\pi/2) = \int_0^{\pi/2}\sqrt{1 - k^2 sin^2\theta}
   *   @f]
   *
   *   @param  k  The argument of the complete elliptic function.
   *   @return  The complete elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    comp_ellint_2(const _Tp k)
    {
      using _Val = num_traits_t<_Tp>;
      if (std::isnan(k))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(k) == _Val{1})
	return _Tp(1);
      else if (std::abs(k) > _Tp(1))
	throw std::domain_error("Bad argument in comp_ellint_2.");
      else
	{
	  const _Tp kk = k * k;

	  return ellint_rf(_Tp(0), _Tp(1) - kk, _Tp(1))
	       - kk * ellint_rd(_Tp(0), _Tp(1) - kk, _Tp(1)) / _Tp(3);
	}
    }

  /**
   *   @brief  Return the incomplete elliptic integral of the second kind
   *           @f$ E(k,\phi) @f$ using the Carlson formulation.
   *
   *   The incomplete elliptic integral of the second kind is defined as
   *   @f[
   *     E(k,\phi) = \int_0^{\phi} \sqrt{1 - k^2 sin^2\theta}
   *   @f]
   *
   *   @param  k  The argument of the elliptic function.
   *   @param  phi  The integral limit argument of the elliptic function.
   *   @return  The elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    ellint_2(const _Tp k, const _Tp phi)
    {
      using _Val = num_traits_t<_Tp>;
      if (std::isnan(k) || std::isnan(phi))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(k) > _Val(1))
	throw std::domain_error("Bad argument in ellint_2.");
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int n = std::floor(phi / numeric_constants<_Tp>::pi()
				   + _Tp(0.5L));
	  const _Tp phi_red = phi
			      - n * numeric_constants<_Tp>::pi();

	  const _Tp kk = k * k;
	  const _Tp s = std::sin(phi_red);
	  const _Tp ss = s * s;
	  const _Tp sss = ss * s;
	  const _Tp c = std::cos(phi_red);
	  const _Tp cc = c * c;

	  const _Tp E = s
			* ellint_rf(cc, _Tp(1) - kk * ss, _Tp(1))
			- kk * sss
			* ellint_rd(cc, _Tp(1) - kk * ss, _Tp(1))
			/ _Tp(3);

	  if (n == 0)
	    return E;
	  else
	    return E + _Tp(2) * n * comp_ellint_2(k);
	}
    }

  /**
   *   @brief Return the complete elliptic integral of the third kind
   *          @f$ \Pi(k,\nu) = \Pi(k,\nu,\pi/2) @f$ using the
   *          Carlson formulation.
   *
   *   The complete elliptic integral of the third kind is defined as
   *   @f[
   *     \Pi(k,\nu) = \int_0^{\pi/2}
   *                   \frac{d\theta}
   *                 {(1 - \nu \sin^2\theta)\sqrt{1 - k^2 \sin^2\theta}}
   *   @f]
   *
   *   @param  k  The argument of the elliptic function.
   *   @param  nu  The second argument of the elliptic function.
   *   @return  The complete elliptic function of the third kind.
   */
  template<typename _Tp>
    _Tp
    comp_ellint_3(const _Tp k, const _Tp nu)
    {
      using _Val = num_traits_t<_Tp>;
      if (std::isnan(k) || std::isnan(nu))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (nu == _Tp(1))
	return std::numeric_limits<_Tp>::infinity();
      else if (std::abs(k) > _Val(1))
	throw std::domain_error("Bad argument in comp_ellint_3.");
      else
	{
	  const _Tp kk = k * k;

	  return ellint_rf(_Tp(0), _Tp(1) - kk, _Tp(1))
	       - nu
	       * ellint_rj(_Tp(0), _Tp(1) - kk, _Tp(1), _Tp(1) + nu)
	       / _Tp(3);
	}
    }

  /**
   *   @brief Return the incomplete elliptic integral of the third kind
   *          @f$ \Pi(k,\nu,\phi) @f$ using the Carlson formulation.
   *
   *   The incomplete elliptic integral of the third kind is defined as
   *   @f[
   *     \Pi(k,\nu,\phi) = \int_0^{\phi}
   *                       \frac{d\theta}
   *                            {(1 - \nu \sin^2\theta)
   *                             \sqrt{1 - k^2 \sin^2\theta}}
   *   @f]
   *
   *   @param  k  The argument of the elliptic function.
   *   @param  nu  The second argument of the elliptic function.
   *   @param  phi  The integral limit argument of the elliptic function.
   *   @return  The elliptic function of the third kind.
   */
  template<typename _Tp>
    _Tp
    ellint_3(const _Tp k, const _Tp nu, const _Tp phi)
    {
      using _Val = num_traits_t<_Tp>;
      if (std::isnan(k) || std::isnan(nu) || std::isnan(phi))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(k) > _Val(1))
	throw std::domain_error("Bad argument in ellint_3.");
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int n = std::floor(phi / numeric_constants<_Tp>::pi()
				   + _Tp(0.5L));
	  const _Tp phi_red = phi
			      - n * numeric_constants<_Tp>::pi();

	  const _Tp kk = k * k;
	  const _Tp s = std::sin(phi_red);
	  const _Tp ss = s * s;
	  const _Tp sss = ss * s;
	  const _Tp c = std::cos(phi_red);
	  const _Tp cc = c * c;

	  const _Tp Pi = s
			 * ellint_rf(cc, _Tp(1) - kk * ss, _Tp(1))
			 - nu * sss
			 * ellint_rj(cc, _Tp(1) - kk * ss, _Tp(1),
				       _Tp(1) + nu * ss) / _Tp(3);

	  if (n == 0)
	    return Pi;
	  else
	    return Pi + _Tp(2) * n * comp_ellint_3(k, nu);
	}
    }

} // namespace detail
//_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

namespace emsr
{
//_GLIBCXX_BEGIN_NAMESPACE_VERSION

  inline float
  ellint_rff(float x, float y, float z)
  { return emsr::detail::ellint_rf<float>(x, y, z); }

  inline long double
  ellint_rfl(long double x, long double y, long double z)
  { return emsr::detail::ellint_rf<long double>(x, y, z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline typename emsr::fp_promote_t<_Tp, _Up, _Vp>
    ellint_rf(_Tp x, _Up y, _Vp z)
    {
      typedef typename emsr::fp_promote_t<_Tp, _Up, _Vp> type;
      return emsr::detail::ellint_rf<type>(x, y, z);
    }

  inline float
  ellint_rcf(float x, float y)
  { return emsr::detail::ellint_rc<float>(x, y); }

  inline long double
  ellint_rcl(long double x, long double y)
  { return emsr::detail::ellint_rc<long double>(x, y); }

  template<typename _Tp, typename _Up>
    inline typename emsr::fp_promote_t<_Tp, _Up>
    ellint_rc(_Tp x, _Up y)
    {
      typedef typename emsr::fp_promote_t<_Tp, _Up> type;
      return emsr::detail::ellint_rc<type>(x, y);
    }

  inline float
  ellint_rjf(float x, float y, float z, float p)
  { return emsr::detail::ellint_rj<float>(x, y, z, p); }

  inline long double
  ellint_rjl(long double x, long double y, long double z, long double p)
  { return emsr::detail::ellint_rj<long double>(x, y, z, p); }

  template<typename _Tp, typename _Up, typename _Vp, typename _Wp>
    inline typename emsr::fp_promote_t<_Tp, _Up, _Vp, _Wp>
    ellint_rj(_Tp x, _Up y, _Vp z, _Wp p)
    {
      typedef typename emsr::fp_promote_t<_Tp, _Up, _Vp, _Wp> type;
      return emsr::detail::ellint_rj<type>(x, y, z, p);
    }

  inline float
  ellint_rdf(float x, float y, float z)
  { return emsr::detail::ellint_rd<float>(x, y, z); }

  inline long double
  ellint_rdl(long double x, long double y, long double z)
  { return emsr::detail::ellint_rd<long double>(x, y, z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline typename emsr::fp_promote_t<_Tp, _Up, _Vp>
    ellint_rd(_Tp x, _Up y, _Vp z)
    {
      typedef typename emsr::fp_promote_t<_Tp, _Up, _Vp> type;
      return emsr::detail::ellint_rd<type>(x, y, z);
    }

  inline float
  ellint_rgf(float x, float y, float z)
  { return emsr::detail::ellint_rg<float>(x, y, z); }

  inline long double
  ellint_rgl(long double x, long double y, long double z)
  { return emsr::detail::ellint_rg<long double>(x, y, z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline typename emsr::fp_promote_t<_Tp, _Up, _Vp>
    ellint_rg(_Tp x, _Up y, _Vp z)
    {
      typedef typename emsr::fp_promote_t<_Tp, _Up, _Vp> type;
      return emsr::detail::ellint_rg<type>(x, y, z);
    }

//_GLIBCXX_END_NAMESPACE_VERSION
} // namespace emsr

template<typename _Tp>
  void
  test_carlson()
  {
    using _Cmplx = std::complex<_Tp>;

    std::cout.precision(std::numeric_limits<_Tp>::digits10);

    std::cout << '\n';
    std::cout << '\n';

    std::cout << "R_F(1, 2, 0) = " << emsr::ellint_rf(_Tp{1}, _Tp{2}, _Tp{0}) << '\n';  //  1.3110287771461
    std::cout << "R_F(i, -i, 0) = " << emsr::ellint_rf(_Cmplx(_Tp{0}, _Tp{1}),
						   _Cmplx(_Tp{0}, _Tp{-1}),
						   _Cmplx(_Tp{0}, _Tp{0})) << '\n';  //  1.8540746773014;
    std::cout << "R_F(i - 1, i, 0) = " << emsr::ellint_rf(_Cmplx(_Tp{-1}, _Tp{1}),
						      _Cmplx(_Tp{0}, _Tp{1}),
						      _Cmplx(_Tp{0}, _Tp{0})) << '\n';  //  0.79612586584234 - i 1.2138566698365
    std::cout << "R_F(0.5, 1, 0) = " << emsr::ellint_rf(0.5, _Tp{1}, _Tp{0}) << '\n';  //  1.8540746773014
    std::cout << "R_F(2, 3, 4) = " << emsr::ellint_rf(_Tp{2}, _Tp{3}, _Tp{4}) << '\n';  //  0.58408284167715
    std::cout << "R_F(i, -i, 2) = " << emsr::ellint_rf(_Cmplx(_Tp{0}, _Tp{1}),
						   _Cmplx(_Tp{0}, _Tp{-1}),
						   _Cmplx(_Tp{2}, _Tp{0})) << '\n';  //  1.0441445654064
    std::cout << "R_F(i - 1, 1 - i, i) = " << emsr::ellint_rf(_Cmplx(_Tp{-1}, _Tp{1}),
							  _Cmplx(_Tp{1}, _Tp{-1}),
							  _Cmplx(_Tp{0}, _Tp{1})) << '\n';  //  0.93912050218619 - i 0.53296252018635

    std::cout << '\n';

    std::cout << "R_C(0, 1/4) = " << emsr::ellint_rc(_Tp{0}, 0.25) << '\n';  //  pi = 3.1415926535898
    std::cout << "R_C(9/4, 2) = " << emsr::ellint_rc(2.25, _Tp{2}) << '\n';  //  ln 2 = 0.69314718055995
    std::cout << "R_C(0, i) = " << emsr::ellint_rc(_Cmplx(_Tp{0}, _Tp{0}),
							 _Cmplx(_Tp{0}, _Tp{1})) << '\n';  //  (1 - i) 1.1107207345396
    std::cout << "R_C(-i, i) = " << emsr::ellint_rc(_Cmplx(_Tp{0}, _Tp{-1}),
							  _Cmplx(_Tp{0}, _Tp{1})) << '\n';  //  1.2260849569072 - i 0.34471136988768
    std::cout << "R_C(1/4, -2) = " << emsr::ellint_rc(0.25, _Tp{-2}) << '\n';  //  ln 2/3 = 0.23104906018665
    std::cout << "R_C(i, -1) = " << emsr::ellint_rc(_Cmplx(_Tp{0}, _Tp{1}),
							  _Cmplx(_Tp{-1}, _Tp{0})) << '\n';  //  0.77778596920447 + i 0.19832484993429

    std::cout << '\n';

    std::cout << "R_J(0, 1, 2, 3) = " << emsr::ellint_rj(_Tp{0}, _Tp{1}, _Tp{2}, _Tp{3}) << '\n';  //  0.77688623778582
    std::cout << "R_J(2, 3, 4, 5) = " << emsr::ellint_rj(_Tp{2}, _Tp{3}, _Tp{4}, _Tp{5}) << '\n';  //  0.14297579667157
    std::cout << "R_J(2, 3, 4, -1 + i) = " << emsr::ellint_rj(_Cmplx(_Tp{2},_Tp{0}),
								    _Cmplx(_Tp{3},_Tp{0}),
								    _Cmplx(_Tp{4},_Tp{0}),
								    _Cmplx(_Tp{-1},_Tp{1})) << '\n';  //  0.13613945827771 - i 0.38207561624427
    std::cout << "R_J(i, -i, 0, 2) = " << emsr::ellint_rj(_Cmplx(_Tp{0},_Tp{1}),
								_Cmplx(_Tp{0},_Tp{-1}),
								_Cmplx(_Tp{0},_Tp{0}),
					     			_Cmplx(_Tp{2},_Tp{0})) << '\n';  //  1.6490011662711
    std::cout << "R_J(-1 + i, -1 - i, 1, 2) = " << emsr::ellint_rj(_Cmplx(_Tp{-1},_Tp{1}),
							       _Cmplx(_Tp{-1},_Tp{-1}),
							       _Cmplx(_Tp{1},_Tp{0}),
							       _Cmplx(_Tp{2},_Tp{0})) << '\n';  //  0.94148358841220
    std::cout << "R_J(i, -i, 0, 1 - i) = " << emsr::ellint_rj(_Cmplx(_Tp{0},_Tp{1}),
							  _Cmplx(_Tp{0},_Tp{-1}),
							  _Cmplx(_Tp{0},_Tp{0}),
							  _Cmplx(_Tp{1},_Tp{-1})) << '\n';  //  1.8260115229009 + i 1.2290661908643
    std::cout << "R_J(-1 + i, -1 - i, 1, -3 + i) = " << emsr::ellint_rj(_Cmplx(_Tp{-1},_Tp{1}),
								    _Cmplx(_Tp{-1},_Tp{-1}),
								    _Cmplx(_Tp{1},_Tp{0}),
								    _Cmplx(_Tp{-3},_Tp{1})) << '\n';  //  -0.61127970812028 - i 1.0684038390007
    std::cout << "R_J(-1 + i, -2 - i, -i, -1 + i) = " << emsr::ellint_rj(_Cmplx(_Tp{-1},_Tp{1}),
								     _Cmplx(_Tp{-2},_Tp{-1}),
								     _Cmplx(_Tp{0},_Tp{-1}),
								     _Cmplx(_Tp{-1},_Tp{1})) << '\n';  //  1.8249027393704 - i 1.2218475784827

    std::cout << '\n';

    std::cout << "R_D(0, 2, 1) = " << emsr::ellint_rd(_Tp{0}, _Tp{2}, _Tp{1}) << '\n';  //  1.7972103521034
    std::cout << "R_D(2, 3, 4) = " << emsr::ellint_rd(_Tp{2}, _Tp{3}, _Tp{4}) << '\n';  //  0.16510527294261
    std::cout << "R_D(i, -i, 2) = " << emsr::ellint_rd(_Cmplx(_Tp{0},_Tp{1}),
						   _Cmplx(_Tp{0},_Tp{-1}),
						   _Cmplx(_Tp{2},_Tp{0})) << '\n';  //  0.65933854154220
    std::cout << "R_D(0, i, -i) = " << emsr::ellint_rd(_Cmplx(_Tp{0},_Tp{0}),
						   _Cmplx(_Tp{0},_Tp{1}),
						   _Cmplx(_Tp{0},_Tp{-1})) << '\n';  //  1.2708196271910 + i 2.7811120159521
    std::cout << "R_D(0, i - 1, i) = " << emsr::ellint_rd(_Cmplx(_Tp{0},_Tp{0}),
						      _Cmplx(_Tp{-1},_Tp{1}),
						      _Cmplx(_Tp{0},_Tp{1})) << '\n';  //  -1.8577235439239 - i 0.96193450888839
    std::cout << "R_D(-2 - i, -i, -1 + i) = " << emsr::ellint_rd(_Cmplx(_Tp{-2},_Tp{-1}),
							     _Cmplx(_Tp{0},_Tp{-1}),
							     _Cmplx(_Tp{-1},_Tp{1})) << '\n';  //  1.8249027393704 - i 1.2218475784827

    std::cout << '\n';

    std::cout << "R_G(0, 16, 16) = 2E(0) = pi = " << emsr::ellint_rg(_Tp{0}, 16.0, 16.0) << '\n';  //  3.1415926535898
    std::cout << "R_G(2, 3, 4) = " << emsr::ellint_rg(2, 3, 4) << '\n';  //  1.7255030280692
    std::cout << "R_G(0, i, -i) = " << emsr::ellint_rg(_Cmplx(_Tp{0},_Tp{0}),
						   _Cmplx(_Tp{0},_Tp{1}),
						   _Cmplx(_Tp{0},_Tp{-1})) << '\n';  //  0.42360654239699
    std::cout << "R_G(i - 1, i, 0) = " << emsr::ellint_rg(_Cmplx(_Tp{-1},_Tp{1}),
						    _Cmplx(_Tp{0},_Tp{1}),
						    _Cmplx(_Tp{0},_Tp{0})) << '\n';  //  0.44660591677018 + i 0.70768352357515
    std::cout << "R_G(-i, i - 1, i) = " << emsr::ellint_rg(_Cmplx(_Tp{0},_Tp{-1}),
						    _Cmplx(_Tp{-1},_Tp{1}),
						    _Cmplx(_Tp{0},_Tp{1})) << '\n';  //  0.36023392184473 + i 0.40348623401722
    std::cout << "R_G(0, 0.0796, 4) = E(0.99) = " << emsr::ellint_rg(0,0.0796,4) << '\n';  //  1.0284758090288

    std::cout << '\n';

    std::cout << "R_F<int>(1, 2, 0) = " << emsr::ellint_rf<int>(_Tp{1}, _Tp{2}, _Tp{0}) << '\n';  //  1.3110287771461

    _Tp pi_2 = emsr::detail::numeric_constants<_Tp>::pi_2();
    std::cout << "K(pi/2) = " << emsr::detail::comp_ellint_1(_Cmplx(pi_2, _Tp{0})) << '\n';  //  1.5887715763658593607082818553065 - 1.3986463677643598308560440635658*i
  }

int
main()
{
  test_carlson<double>();

  return 0;
}

