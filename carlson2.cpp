
// ../bin/bin/g++ -g -std=c++11 -o carlson2 carlson2.cpp

#include <cmath>
#include <iostream>
#include <complex>
#include <limits>

//  This should be added to <complex>.  Or not?

namespace __gnu_cxx { //_GLIBCXX_BEGIN_NAMESPACE(__gnu_cxx)

  //  They recently poisoned __promote in ext/type_traits.h to prevent usage of math functions with weird types
  template<>
    struct __promote<std::complex<long double>>
    { typedef std::complex<long double> __type; };

  template<>
    struct __promote<std::complex<double>>
    { typedef std::complex<double> __type; };

  template<>
    struct __promote<std::complex<float>>
    { typedef std::complex<float> __type; };

} //_GLIBCXX_END_NAMESPACE


// From bits/special_function_util.h

namespace std
{

  namespace __detail
  {

    /// A class to encapsulate type dependent floating point
    /// constants.  Not everything will be able to be expressed as
    /// type logic.
    template<typename _Tp>
    struct __floating_point_constant
    {
      static const _Tp __value;
    };


    /// A structure for numeric constants.
    template<typename _Tp>
      struct __numeric_constants
      {
	///  Constant @f$ \pi @f$.
	static _Tp __pi() throw()
	{ return static_cast<_Tp>(3.1415926535897932384626433832795029L); }
	///  Constant @f$ \pi / 2 @f$.
	static _Tp __pi_2() throw()
	{ return static_cast<_Tp>(1.5707963267948966192313216916397514L); }
	///  Constant @f$ \pi / 3 @f$.
	static _Tp __pi_3() throw()
	{ return static_cast<_Tp>(1.0471975511965977461542144610931676L); }
	///  Constant @f$ \pi / 4 @f$.
	static _Tp __pi_4() throw()
	{ return static_cast<_Tp>(0.7853981633974483096156608458198757L); }
	///  Constant @f$ 1 / \pi @f$.
	static _Tp __1_pi() throw()
	{ return static_cast<_Tp>(0.3183098861837906715377675267450287L); }
	///  Constant @f$ 2 / \sqrt(\pi) @f$.
	static _Tp __2_sqrtpi() throw()
	{ return static_cast<_Tp>(1.1283791670955125738961589031215452L); }
	///  Constant @f$ \sqrt(2) @f$.
	static _Tp __sqrt2() throw()
	{ return static_cast<_Tp>(1.4142135623730950488016887242096981L); }
	///  Constant @f$ \sqrt(3) @f$.
	static _Tp __sqrt3() throw()
	{ return static_cast<_Tp>(1.7320508075688772935274463415058723L); }
	///  Constant @f$ \sqrt(\pi/2) @f$.
	static _Tp __sqrtpio2() throw()
	{ return static_cast<_Tp>(1.2533141373155002512078826424055226L); }
	///  Constant @f$ 1 / sqrt(2) @f$.
	static _Tp __sqrt1_2() throw()
	{ return static_cast<_Tp>(0.7071067811865475244008443621048490L); }
	///  Constant @f$ \log(\pi) @f$.
	static _Tp __lnpi() throw()
	{ return static_cast<_Tp>(1.1447298858494001741434273513530587L); }
	///  Constant Euler's constant @f$ \gamma_E @f$.
	static _Tp __gamma_e() throw()
	{ return static_cast<_Tp>(0.5772156649015328606065120900824024L); }
	///  Constant Euler-Mascheroni @f$ e @f$
	static _Tp __euler() throw()
	{ return static_cast<_Tp>(2.7182818284590452353602874713526625L); }
      };


#if _GLIBCXX_USE_C99_MATH && !_GLIBCXX_USE_C99_FP_MACROS_DYNAMIC

    /// This is a wrapper for the isnan function. Otherwise, for NaN,
    /// all comparisons result in false. If/when we build a std::isnan
    /// out of intrinsics, this will disappear completely in favor of
    /// std::isnan.
    template<typename _Tp>
    inline bool __isnan(const _Tp __x)
    {
      return std::isnan(__x);
    }

#else

    template<typename _Tp>
    inline bool __isnan(const _Tp __x)
    {
      return __builtin_isnan(__x);
    }

    template<>
    inline bool __isnan<float>(const float __x)
    {
      return __builtin_isnanf(__x);
    }

    template<>
    inline bool __isnan<long double>(const long double __x)
    {
      return __builtin_isnanl(__x);
    }

#endif

  } // namespace __detail

}

// End bits/special_function_util.h

namespace std
{

  namespace __detail
  {

    template<typename _Tp>
      struct __ellint_traits
      {
	typedef _Tp __value_type;
      };

    template<>
      template<typename _Tp>
	struct __ellint_traits<std::complex<_Tp> >
	{
	  typedef _Tp __value_type;
	};

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
     *   @param  __x  The first of three symmetric arguments.
     *   @param  __y  The second of three symmetric arguments.
     *   @param  __z  The third of three symmetric arguments.
     *   @return  The Carlson elliptic function of the first kind.
     */
    template<typename _Tp>
      _Tp
      __ellint_rf(_Tp __x, _Tp __y, _Tp __z)
      {
	typedef typename __ellint_traits<_Tp>::__value_type _Val;
	const _Val __r = std::numeric_limits<_Val>::epsilon();
	_Tp __xt = __x;
	_Tp __yt = __y;
	_Tp __zt = __z;
	_Tp __a0 = (__x + __y + __z) / _Val(3);
	_Val __q = std::pow( _Val(3) * __r, -_Val(1) / _Val(6) )
		 * std::max(std::abs(__a0 - __z),
			    std::max(std::abs(__a0 - __x),
				     std::abs(__a0 - __y)));
	_Tp __a = __a0;
	_Val __f = _Val(1);

	for (;;)
	  {
	    _Tp __lambda = std::sqrt(__xt) * std::sqrt(__yt)
			 + std::sqrt(__yt) * std::sqrt(__zt)
			 + std::sqrt(__zt) * std::sqrt(__xt);
	    __a = (__a + __lambda) / _Val(4);
	    __xt = (__xt + __lambda) / _Val(4);
	    __yt = (__yt + __lambda) / _Val(4);
	    __zt = (__zt + __lambda) / _Val(4);
	    __f *= _Val(4);
	    if (__q < __f * std::abs(__a))
	      {
		_Tp __xf = (__a0 - __x) / (__f * __a);
		_Tp __yf = (__a0 - __y) / (__f * __a);
		_Tp __zf = -(__xf + __yf);
		_Tp __e2 = __xf * __yf - __zf * __zf;
		_Tp __e3 = __xf * __yf * __zf;
		return (_Val(1)
		      - __e2 / _Val(10)
		      + __e3 / _Val(14)
		      + __e2 * __e2 / _Val(24)
		      - _Val(3) * __e2 * __e3 / _Val(44)) / std::sqrt(__a);
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
     *   @param  __x  The first argument.
     *   @param  __y  The second argument.
     *   @return  The Carlson elliptic function.
     */
    template<typename _Tp>
      _Tp
      __ellint_rc(_Tp __x, _Tp __y)
      {
	typedef typename __ellint_traits<_Tp>::__value_type _Val;
	if (std::imag(__y) == _Val(0) && std::real(__y) < _Val(0))
	  return std::sqrt(__x / (__x - __y)) * __ellint_rc(__x - __y, -__y);
	const _Val __r = std::numeric_limits<_Val>::epsilon();
	_Tp __xt = __x;
	_Tp __yt = __y;
	_Tp __a0 = (__x + _Val(2) * __y) / _Val(3);
	_Val __q = std::pow( _Val(3) * __r, -_Val(1) / _Val(8) )
		 * std::abs(__a0 - __x);
	_Tp __a = __a0;
	_Val __f = _Val(1);

	for (;;)
	  {
	    _Tp __lambda = _Val(2) * std::sqrt(__xt) * std::sqrt(__yt) + __yt;
	    __a = (__a + __lambda) / _Val(4);
	    __xt = (__xt + __lambda) / _Val(4);
	    __yt = (__yt + __lambda) / _Val(4);
	    __f *= _Val(4);
	    if (__q < __f * std::abs(__a))
	      {
		_Tp __s = (__y - __a0) / (__f * __a);
		return (_Val(1) + __s * __s * (_Val(3) / _Val(10)
		      + __s * (_Val(1) / _Val(7)
		      + __s * (_Val(3) / _Val(8)
		      + __s * (_Val(9) / _Val(22)
		      + __s * (_Val(159) / _Val(208)
		      + __s * (_Val(9) / _Val(8)))))))) / std::sqrt(__a);
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
     *   @param  __x  The first of three symmetric arguments.
     *   @param  __y  The second of three symmetric arguments.
     *   @param  __z  The third of three symmetric arguments.
     *   @param  __p  The fourth argument.
     *   @return  The Carlson elliptic function of the fourth kind.
     */
    template<typename _Tp>
      _Tp
      __ellint_rj(_Tp __x, _Tp __y, _Tp __z, _Tp __p)
      {
	typedef typename __ellint_traits<_Tp>::__value_type _Val;
	const _Val __r = std::numeric_limits<_Val>::epsilon();
	_Tp __xt = __x;
	_Tp __yt = __y;
	_Tp __zt = __z;
	_Tp __pt = __p;
	_Tp __a0 = (__x + __y + __z + _Val(2) * __p) / _Val(5);
	_Tp __delta = (__p - __x) * (__p - __y) * (__p - __z);
	_Val __q = std::pow(__r / _Val(4), -_Val(1) / _Val(6))
		 * std::max(std::abs(__a0 - __z),
			    std::max(std::abs(__a0 - __x),
				     std::max(std::abs(__a0 - __y), std::abs(__a0 - __p))));
	_Tp __a = __a0;
	_Val __f = _Val(1);
	_Val __fe = _Val(1);
	_Tp __sum = _Tp();

	for (;;)
	  {
	    _Tp __xroot = std::sqrt(__xt);
	    _Tp __yroot = std::sqrt(__yt);
	    _Tp __zroot = std::sqrt(__zt);
	    _Tp __proot = std::sqrt(__pt);
	    _Tp __lambda = __xroot * __yroot
			 + __yroot * __zroot
			 + __zroot * __xroot;
	    __a = (__a + __lambda) / _Val(4);
	    __xt = (__xt + __lambda) / _Val(4);
	    __yt = (__yt + __lambda) / _Val(4);
	    __zt = (__zt + __lambda) / _Val(4);
	    __pt = (__pt + __lambda) / _Val(4);
	    _Tp __d = (__proot + __xroot) * (__proot + __yroot) * (__proot + __zroot);
	    _Tp __e = __delta / (__fe * __d * __d);
	    __sum += __ellint_rc(_Tp(1), _Tp(1) + __e) / (__f * __d);
	    __f *= _Val(4);
	    __fe *= _Val(64);
	    if (__q < __f * std::abs(__a))
	      {
		_Tp __xf = (__a0 - __x) / (__f * __a);
		_Tp __yf = (__a0 - __y) / (__f * __a);
		_Tp __zf = (__a0 - __z) / (__f * __a);
		_Tp __xyz = __xf * __yf * __zf;
		_Tp __pf = -(__xf + __yf + __zf) / _Val(2);
		_Tp __pp = __pf * __pf;
		_Tp __ppp = __pp * __pf;
		_Tp __e2 = __xf * __yf + __yf * __zf + __zf * __xf - _Val(3) * __pp;
		_Tp __e3 = __xyz + _Val(2) * __e2 * __pf + _Tp(4) * __ppp;
		_Tp __e4 = (_Val(2) * __xyz + __e2 * __pf + _Val(3) * __ppp) * __pf;
		_Tp __e5 = __xyz * __pp;
		return (_Val(1) - _Val(3) * __e2 / _Val(14)
				+ __e3 / _Val(6)
				+ _Val(9) * __e2 * __e2 / _Val(88)
				- _Val(3) * __e4 / _Val(22)
				- _Val(9) * __e2 * __e3 / _Val(52)
				+ _Val(3) * __e5 / _Val(26)) / __f / __a / std::sqrt(__a)
				+ _Val(6) * __sum;
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
     *   @param  __x  The first of two symmetric arguments.
     *   @param  __y  The second of two symmetric arguments.
     *   @param  __z  The third argument.
     *   @return  The Carlson elliptic function of the second kind.
     */
    template<typename _Tp>
      _Tp
      __ellint_rd(_Tp __x, _Tp __y, _Tp __z)
      {
	typedef typename __ellint_traits<_Tp>::__value_type _Val;
	const _Val __r = std::numeric_limits<_Val>::epsilon();
	_Tp __xt = __x;
	_Tp __yt = __y;
	_Tp __zt = __z;
	_Tp __a0 = (__x + __y + _Val(3) * __z) / _Val(5);
	_Val __q = std::pow(__r / _Val(4), -_Val(1) / _Val(6))
		 * std::max(std::abs(__a0 - __z),
			    std::max(std::abs(__a0 - __x),
			    std::abs(__a0 - __y)));
	_Tp __a = __a0;
	_Val __f = _Val(1);
	_Tp __sum = _Tp();

	for (;;)
	  {
	    _Tp __lambda = std::sqrt(__xt) * std::sqrt(__yt)
			 + std::sqrt(__yt) * std::sqrt(__zt)
			 + std::sqrt(__zt) * std::sqrt(__xt);
	    __sum += _Val(1) / __f / std::sqrt(__zt) / (__zt + __lambda);
	    __a = (__a + __lambda) / _Val(4);
	    __xt = (__xt + __lambda) / _Val(4);
	    __yt = (__yt + __lambda) / _Val(4);
	    __zt = (__zt + __lambda) / _Val(4);
	    __f *= _Val(4);
	    if (__q < __f * std::abs(__a))
	      {
		_Tp __xf = (__a0 - __x) / (__f * __a);
		_Tp __yf = (__a0 - __y) / (__f * __a);
		_Tp __zf = -(__xf + __yf) / _Val(3);
		_Tp __zz = __zf * __zf;
		_Tp __xy = __xf * __yf;
		_Tp __e2 = __xy - _Val(6) * __zz;
		_Tp __e3 = (_Val(3) * __xy - _Val(8) * __zz) * __zf;
		_Tp __e4 = _Val(3) * (__xy - __zz) * __zz;
		_Tp __e5 = __xy * __zf * __zz;
		return (_Val(1)
		      - _Val(3) * __e2 / _Val(14)
		      + __e3 / _Val(6)
		      + _Val(9) * __e2 * __e2 / _Val(88)
		      - _Val(3) * __e4 / _Val(22)
		      - _Val(9) * __e2 * __e3 / _Val(52)
		      + _Val(3) * __e5 / _Val(26)) / __f / __a / std::sqrt(__a)
		      + _Val(3) * __sum;
	      }
	  }

	return _Tp(0);
      }

    template<typename _Tp>
      _Tp
      __comp_ellint_rf(_Tp __x, _Tp __y)
      {
	typedef typename __ellint_traits<_Tp>::__value_type _Val;
	const _Val __r = std::numeric_limits<_Val>::epsilon();
	const _Val __tolfact = _Val(2.7L) * std::sqrt(__r);
	__x = std::sqrt(__x);
	__y = std::sqrt(__y);
	for (;;)
	  {
	    _Tp __xt = __x;
	    __x = (__x + __y) / _Tp(2);
	    __y = std::sqrt(__xt * __y);
	    if (std::abs(__x - __y) < __tolfact * std::abs(__x))
	      return _Val(__numeric_constants<_Val>::__pi()) / (__x + __y);
	  }
      }

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
     *   @param  __x  The first of three symmetric arguments.
     *   @param  __y  The second of three symmetric arguments.
     *   @param  __z  The third of three symmetric arguments.
     *   @return  The Carlson symmetric elliptic function of the second kind.
     */

    template<typename _Tp>
      _Tp
      __comp_ellint_rg(_Tp __x, _Tp __y);

    template<typename _Tp>
      _Tp
      __ellint_rg(_Tp __x, _Tp __y, _Tp __z)
      {
	typedef typename __ellint_traits<_Tp>::__value_type _Val;
	if (__z == _Tp())
	  {
	    if (__x == _Tp())
	      return std::sqrt(__y);
	    else if (__y == _Tp())
	      return std::sqrt(__x);
	    else
	      return __comp_ellint_rg(__x, __y);
	  }
	else if (__x == _Tp())
	  {
	    if (__y == _Tp())
	      return std::sqrt(__z);
	    else if (__z == _Tp())
	      return std::sqrt(__y);
	    else
	      return __comp_ellint_rg(__y, __z);
	  }
	else if (__y == _Tp())
	  {
	    if (__z == _Tp())
	      return std::sqrt(__x);
	    else if (__x == _Tp())
	      return std::sqrt(__z);
	    else
	      return __comp_ellint_rg(__z, __x);
	  }
	else
	  return (__z * __ellint_rf(__x, __y, __z)
		- (__x - __z) * (__y - __z) * __ellint_rd(__x, __y, __z) / _Val(3)
		+ std::sqrt(__x * __y / __z)) / _Val(2);
      }

    template<typename _Tp>
      _Tp
      __comp_ellint_rg(_Tp __x, _Tp __y)
      {
	typedef typename __ellint_traits<_Tp>::__value_type _Val;
	const _Val __r = std::numeric_limits<_Val>::epsilon();
	const _Val __tolfact = _Val(2.7L) * std::sqrt(__r);
	_Tp __xt = std::sqrt(__x);
	_Tp __yt = std::sqrt(__y);
	const _Tp __a = (__xt + __yt) / _Val(2);
	_Tp __sum = _Tp();
	_Val __sf = _Val(1) / _Val(2);
	for (;;)
	  {
	    _Tp __xtt = __xt;
	    __xt = (__xt + __yt) / _Tp(2);
	    __yt = std::sqrt(__xtt * __yt);
	    _Tp __del = __xt - __yt;
	    if (std::abs(__del) < __tolfact * std::abs(__xt))
	      return (__a * __a - __sum) * _Val(__numeric_constants<_Val>::__pi()) / (__xt + __yt) / _Val(2);
	    __sum += __sf * __del * __del;
	    __sf *= _Val(2);
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
     *   @param  __k  The argument of the complete elliptic function.
     *   @return  The complete elliptic function of the first kind.
     */
    template<typename _Tp>
    _Tp
    __comp_ellint_1(const _Tp __k)
    {
      if (__isnan(std::real(__k)) || __isnan(std::imag(__k)))
	return std::numeric_limits<_Tp>::quiet_NaN();
      //else if (std::abs(__k) >= _Tp(1))
	//return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return __ellint_rf(_Tp(0), _Tp(1) - __k * __k, _Tp(1));
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
     *   @param  __k  The argument of the elliptic function.
     *   @param  __phi  The integral limit argument of the elliptic function.
     *   @return  The elliptic function of the first kind.
     */
    template<typename _Tp>
    _Tp
    __ellint_1(const _Tp __k, const _Tp __phi)
    {
      if (__isnan(__k) || __isnan(__phi))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(__k) > _Tp(1))
	std::__throw_domain_error(__N("Bad argument in __ellint_1."));
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(__phi / __numeric_constants<_Tp>::__pi()
				   + _Tp(0.5L));
	  const _Tp __phi_red = __phi
			      - __n * __numeric_constants<_Tp>::__pi();

	  const _Tp __s = std::sin(__phi_red);
	  const _Tp __c = std::cos(__phi_red);

	  const _Tp __F = __s
			* __ellint_rf(__c * __c,
				_Tp(1) - __k * __k * __s * __s, _Tp(1));

	  if (__n == 0)
	    return __F;
	  else
	    return __F + _Tp(2) * __n * __comp_ellint_1(__k);
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
     *   @param  __k  The argument of the complete elliptic function.
     *   @return  The complete elliptic function of the second kind.
     */
    template<typename _Tp>
    _Tp
    __comp_ellint_2(const _Tp __k)
    {
      if (__isnan(__k))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(__k) == 1)
	return _Tp(1);
      else if (std::abs(__k) > _Tp(1))
	std::__throw_domain_error(__N("Bad argument in __comp_ellint_2."));
      else
	{
	  const _Tp __kk = __k * __k;

	  return __ellint_rf(_Tp(0), _Tp(1) - __kk, _Tp(1))
	       - __kk * __ellint_rd(_Tp(0), _Tp(1) - __kk, _Tp(1)) / _Tp(3);
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
     *   @param  __k  The argument of the elliptic function.
     *   @param  __phi  The integral limit argument of the elliptic function.
     *   @return  The elliptic function of the second kind.
     */
    template<typename _Tp>
    _Tp
    __ellint_2(const _Tp __k, const _Tp __phi)
    {
      if (__isnan(__k) || __isnan(__phi))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(__k) > _Tp(1))
	std::__throw_domain_error(__N("Bad argument in __ellint_2."));
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(__phi / __numeric_constants<_Tp>::__pi()
				   + _Tp(0.5L));
	  const _Tp __phi_red = __phi
			      - __n * __numeric_constants<_Tp>::__pi();

	  const _Tp __kk = __k * __k;
	  const _Tp __s = std::sin(__phi_red);
	  const _Tp __ss = __s * __s;
	  const _Tp __sss = __ss * __s;
	  const _Tp __c = std::cos(__phi_red);
	  const _Tp __cc = __c * __c;

	  const _Tp __E = __s
			* __ellint_rf(__cc, _Tp(1) - __kk * __ss, _Tp(1))
			- __kk * __sss
			* __ellint_rd(__cc, _Tp(1) - __kk * __ss, _Tp(1))
			/ _Tp(3);

	  if (__n == 0)
	    return __E;
	  else
	    return __E + _Tp(2) * __n * __comp_ellint_2(__k);
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
     *   @param  __k  The argument of the elliptic function.
     *   @param  __nu  The second argument of the elliptic function.
     *   @return  The complete elliptic function of the third kind.
     */
    template<typename _Tp>
    _Tp
    __comp_ellint_3(const _Tp __k, const _Tp __nu)
    {
      if (__isnan(__k) || __isnan(__nu))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__nu == _Tp(1))
	return std::numeric_limits<_Tp>::infinity();
      else if (std::abs(__k) > _Tp(1))
	std::__throw_domain_error(__N("Bad argument in __comp_ellint_3."));
      else
	{
	  const _Tp __kk = __k * __k;

	  return __ellint_rf(_Tp(0), _Tp(1) - __kk, _Tp(1))
	       - __nu
	       * __ellint_rj(_Tp(0), _Tp(1) - __kk, _Tp(1), _Tp(1) + __nu)
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
     *   @param  __k  The argument of the elliptic function.
     *   @param  __nu  The second argument of the elliptic function.
     *   @param  __phi  The integral limit argument of the elliptic function.
     *   @return  The elliptic function of the third kind.
     */
    template<typename _Tp>
    _Tp
    __ellint_3(const _Tp __k, const _Tp __nu, const _Tp __phi)
    {
      if (__isnan(__k) || __isnan(__nu) || __isnan(__phi))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::abs(__k) > _Tp(1))
	std::__throw_domain_error(__N("Bad argument in __ellint_3."));
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(__phi / __numeric_constants<_Tp>::__pi()
				   + _Tp(0.5L));
	  const _Tp __phi_red = __phi
			      - __n * __numeric_constants<_Tp>::__pi();

	  const _Tp __kk = __k * __k;
	  const _Tp __s = std::sin(__phi_red);
	  const _Tp __ss = __s * __s;
	  const _Tp __sss = __ss * __s;
	  const _Tp __c = std::cos(__phi_red);
	  const _Tp __cc = __c * __c;

	  const _Tp __Pi = __s
			 * __ellint_rf(__cc, _Tp(1) - __kk * __ss, _Tp(1))
			 - __nu * __sss
			 * __ellint_rj(__cc, _Tp(1) - __kk * __ss, _Tp(1),
				       _Tp(1) + __nu * __ss) / _Tp(3);

	  if (__n == 0)
	    return __Pi;
	  else
	    return __Pi + _Tp(2) * __n * __comp_ellint_3(__k, __nu);
	}
    }

  } // __detail

}

namespace __gnu_cxx
{

  inline float
  ellint_rff(float __x, float __y, float __z)
  { return std::__detail::__ellint_rf<float>(__x, __y, __z); }

  inline long double
  ellint_rfl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rf<long double>(__x, __y, __z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline typename __gnu_cxx::__promote_3<_Tp, _Up, _Vp>::__type
    ellint_rf(_Tp __x, _Up __y, _Vp __z)
    {
      typedef typename __gnu_cxx::__promote_3<_Tp, _Up, _Vp>::__type __type;
      return std::__detail::__ellint_rf<__type>(__x, __y, __z);
    }

  inline float
  ellint_rcf(float __x, float __y)
  { return std::__detail::__ellint_rc<float>(__x, __y); }

  inline long double
  ellint_rcl(long double __x, long double __y)
  { return std::__detail::__ellint_rc<long double>(__x, __y); }

  template<typename _Tp, typename _Up>
    inline typename __gnu_cxx::__promote_2<_Tp, _Up>::__type
    ellint_rc(_Tp __x, _Up __y)
    {
      typedef typename __gnu_cxx::__promote_2<_Tp, _Up>::__type __type;
      return std::__detail::__ellint_rc<__type>(__x, __y);
    }

  inline float
  ellint_rjf(float __x, float __y, float __z, float __p)
  { return std::__detail::__ellint_rj<float>(__x, __y, __z, __p); }

  inline long double
  ellint_rjl(long double __x, long double __y, long double __z, long double __p)
  { return std::__detail::__ellint_rj<long double>(__x, __y, __z, __p); }

  template<typename _Tp, typename _Up, typename _Vp, typename _Wp>
    inline typename __gnu_cxx::__promote_4<_Tp, _Up, _Vp, _Wp>::__type
    ellint_rj(_Tp __x, _Up __y, _Vp __z, _Wp __p)
    {
      typedef typename __gnu_cxx::__promote_4<_Tp, _Up, _Vp, _Wp>::__type __type;
      return std::__detail::__ellint_rj<__type>(__x, __y, __z, __p);
    }

  inline float
  ellint_rdf(float __x, float __y, float __z)
  { return std::__detail::__ellint_rd<float>(__x, __y, __z); }

  inline long double
  ellint_rdl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rd<long double>(__x, __y, __z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline typename __gnu_cxx::__promote_3<_Tp, _Up, _Vp>::__type
    ellint_rd(_Tp __x, _Up __y, _Vp __z)
    {
      typedef typename __gnu_cxx::__promote_3<_Tp, _Up, _Vp>::__type __type;
      return std::__detail::__ellint_rd<__type>(__x, __y, __z);
    }

  inline float
  ellint_rgf(float __x, float __y, float __z)
  { return std::__detail::__ellint_rg<float>(__x, __y, __z); }

  inline long double
  ellint_rgl(long double __x, long double __y, long double __z)
  { return std::__detail::__ellint_rg<long double>(__x, __y, __z); }

  template<typename _Tp, typename _Up, typename _Vp>
    inline typename __gnu_cxx::__promote_3<_Tp, _Up, _Vp>::__type
    ellint_rg(_Tp __x, _Up __y, _Vp __z)
    {
      typedef typename __gnu_cxx::__promote_3<_Tp, _Up, _Vp>::__type __type;
      return std::__detail::__ellint_rg<__type>(__x, __y, __z);
    }
}


int
main()
{
  std::cout.precision(14);
  std::cout << "R_F(1, 2, 0) = " << __gnu_cxx::ellint_rf(1.0, 2.0, 0.0) << std::endl;  //  1.3110287771461
  std::cout << "R_F(i, -i, 0) = " << __gnu_cxx::ellint_rf(std::complex<double>(0.0, 1.0),
						 std::complex<double>(0.0, -1.0),
						 0.0) << std::endl;  //  1.8540746773014;
  std::cout << "R_F(i - 1, i, 0) = " << __gnu_cxx::ellint_rf(std::complex<double>(-1.0, 1.0),
						    std::complex<double>(0.0, 1.0),
						    0.0) << std::endl;  //  0.79612586584234 - i 1.2138566698365
  std::cout << "R_F(0.5, 1, 0) = " << __gnu_cxx::ellint_rf(0.5, 1.0, 0.0) << std::endl;  //  1.8540746773014
  std::cout << "R_F(2, 3, 4) = " << __gnu_cxx::ellint_rf(2.0, 3.0, 4.0) << std::endl;  //  0.58408284167715
  std::cout << "R_F(i, -i, 2) = " << __gnu_cxx::ellint_rf(std::complex<double>(0.0, 1.0),
						 std::complex<double>(0.0, -1.0),
						 2.0) << std::endl;  //  1.0441445654064
  std::cout << "R_F(i - 1, 1 - i, i) = " << __gnu_cxx::ellint_rf(std::complex<double>(-1.0, 1.0),
							std::complex<double>(1.0, -1.0),
							std::complex<double>(0.0, 1.0)) << std::endl;  //  0.93912050218619 - i 0.53296252018635

  std::cout << std::endl;

  std::cout << "R_C(0, 1/4) = " << __gnu_cxx::ellint_rc(0.0, 0.25) << std::endl;  //  pi = 3.1415926535898
  std::cout << "R_C(9/4, 2) = " << __gnu_cxx::ellint_rc(2.25, 2.0) << std::endl;  //  ln 2 = 0.69314718055995
  std::cout << "R_C(0, i) = " << __gnu_cxx::ellint_rc(0.0,
						       std::complex<double>(0.0, 1.0)) << std::endl;  //  (1 - i) 1.1107207345396
  std::cout << "R_C(-i, i) = " << __gnu_cxx::ellint_rc(std::complex<double>(0.0, -1.0),
							std::complex<double>(0.0, 1.0)) << std::endl;  //  1.2260849569072 - i 0.34471136988768
  std::cout << "R_C(1/4, -2) = " << __gnu_cxx::ellint_rc(0.25, -2.0) << std::endl;  //  ln 2/3 = 0.23104906018665
  std::cout << "R_C(i, -1) = " << __gnu_cxx::ellint_rc(std::complex<double>(0.0, 1.0),
							-1.0) << std::endl;  //  0.77778596920447 + i 0.19832484993429

  std::cout << std::endl;

  std::cout << "R_J(0, 1, 2, 3) = " << __gnu_cxx::ellint_rj(0.0, 1.0, 2.0, 3.0) << std::endl;  //  0.77688623778582
  std::cout << "R_J(2, 3, 4, 5) = " << __gnu_cxx::ellint_rj(2.0, 3.0, 4.0, 5.0) << std::endl;  //  0.14297579667157
  std::cout << "R_J(2, 3, 4, -1 + i) = " << __gnu_cxx::ellint_rj(2.0,
								  3.0,
								  4.0,
								  std::complex<double>(-1.0,1.0)) << std::endl;  //  0.13613945827771 - i 0.38207561624427
  std::cout << "R_J(i, -i, 0, 2) = " << __gnu_cxx::ellint_rj(std::complex<double>(0.0,1.0),
							      std::complex<double>(0.0,-1.0),
							      0.0,
					     		      2.0) << std::endl;  //  1.6490011662711
  std::cout << "R_J(-1 + i, -1 - i, 1, 2) = " << __gnu_cxx::ellint_rj(std::complex<double>(-1.0,1.0),
							     std::complex<double>(-1.0,-1.0),
							     1.0,
							     2.0) << std::endl;  //  0.94148358841220
  std::cout << "R_J(i, -i, 0, 1 - i) = " << __gnu_cxx::ellint_rj(std::complex<double>(0.0,1.0),
							std::complex<double>(0.0,-1.0),
							0.0,
							std::complex<double>(1.0,-1.0)) << std::endl;  //  1.8260115229009 + i 1.2290661908643
  std::cout << "R_J(-1 + i, -1 - i, 1, -3 + i) = " << __gnu_cxx::ellint_rj(std::complex<double>(-1.0,1.0),
								  std::complex<double>(-1.0,-1.0),
								  1.0,
								  std::complex<double>(-3.0,1.0)) << std::endl;  //  -0.61127970812028 - i 1.0684038390007
  std::cout << "R_J(-1 + i, -2 - i, -i, -1 + i) = " << __gnu_cxx::ellint_rj(std::complex<double>(-1.0,1.0),
								   std::complex<double>(-2.0,-1.0),
								   std::complex<double>(0.0,-1.0),
								   std::complex<double>(-1.0,1.0)) << std::endl;  //  1.8249027393704 - i 1.2218475784827

  std::cout << std::endl;

  std::cout << "R_D(0, 2, 1) = " << __gnu_cxx::ellint_rd(0.0, 2.0, 1.0) << std::endl;  //  1.7972103521034
  std::cout << "R_D(2, 3, 4) = " << __gnu_cxx::ellint_rd(2.0, 3.0, 4.0) << std::endl;  //  0.16510527294261
  std::cout << "R_D(i, -i, 2) = " << __gnu_cxx::ellint_rd(std::complex<double>(0.0,1.0),
						 std::complex<double>(0.0,-1.0),
						 2.0) << std::endl;  //  0.65933854154220
  std::cout << "R_D(0, i, -i) = " << __gnu_cxx::ellint_rd(0.0,
						 std::complex<double>(0.0,1.0),
						 std::complex<double>(0.0,-1.0)) << std::endl;  //  1.2708196271910 + i 2.7811120159521
  std::cout << "R_D(0, i - 1, i) = " << __gnu_cxx::ellint_rd(0.0,
						    std::complex<double>(-1.0,1.0),
						    std::complex<double>(0.0,1.0)) << std::endl;  //  -1.8577235439239 - i 0.96193450888839
  std::cout << "R_D(-2 - i, -i, -1 + i) = " << __gnu_cxx::ellint_rd(std::complex<double>(-2.0,-1.0),
							   std::complex<double>(0.0,-1.0),
							   std::complex<double>(-1.0,1.0)) << std::endl;  //  1.8249027393704 - i 1.2218475784827

  std::cout << std::endl;

  std::cout << "R_G(0, 16, 16) = 2E(0) = pi = " << __gnu_cxx::ellint_rg(0.0, 16.0, 16.0) << std::endl;  //  3.1415926535898
  std::cout << "R_G(2, 3, 4) = " << __gnu_cxx::ellint_rg(2, 3, 4) << std::endl;  //  1.725503020692
  std::cout << "R_G(0, i, -i) = " << __gnu_cxx::ellint_rg(0.0,
						 std::complex<double>(0.0,1.0),
						 std::complex<double>(0.0,-1.0)) << std::endl;  //  0.42360654239699
  std::cout << "R_G(i - 1, i, 0) = " << __gnu_cxx::ellint_rg(std::complex<double>(-1.0,1.0),
						  std::complex<double>(0.0,1.0),
						  0.0) << std::endl;  //  0.44660591677018 + i 0.70768352357515
  std::cout << "R_G(-i, i - 1, i) = " << __gnu_cxx::ellint_rg(std::complex<double>(0.0,-1.0),
						  std::complex<double>(-1.0,1.0),
						  std::complex<double>(0.0,1.0)) << std::endl;  //  0.36023392184473 + i 0.40348623401722
  std::cout << "R_G(0, 0.0796, 4) = E(0.99) = " << __gnu_cxx::ellint_rg(0,0.0796,4) << std::endl;  //  1.0284758090288

  std::cout << std::endl;

  std::cout << "R_F<int>(1, 2, 0) = " << __gnu_cxx::ellint_rf<int>(1.0, 2.0, 0.0) << std::endl;  //  1.3110287771461

  std::cout << std::endl;

  double __pi_2 = std::__detail::__numeric_constants<double>::__pi_2();
  std::cout << "K(pi/2) = " << std::__detail::__comp_ellint_1(std::complex<double>(__pi_2, 0.0)) << std::endl;  //  1.5887715763658593607082818553065 - 1.3986463677643598308560440635658*i

  return 0;
}

