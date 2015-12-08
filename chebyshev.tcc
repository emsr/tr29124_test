#ifndef _CHEBYSHEV_TCC
#define _CHEBYSHEV_TCC 1

#include <utility>
#include <iterator>
#include <algorithm> // For find_if.

  /**
   *  @brief  Construct a Chebyshev fit of a function.
   *  Given a function @c func, the lower limit @c a, the upper limit @c b,
   *  and the number of points @c n this routine computes the coefficients
   *  of a Chebyshev polynomial expansion such that
   *    func(x) =~ kSum_0^n-1 c_k T_k(x) - c_0/2.
   *  Use this routine with moderately large n of 30 or 50.
   *  Then truncate the series to a smaller number of terms to satisfy
   *  the accuracy requirements.
   */
  template<typename _Tp>
    _Chebyshev<_Tp>::_Chebyshev(_Tp __a, _Tp __b, unsigned __n, _Tp __func(_Tp))
    : _M_lower(__a),
      _M_upper(__b),
      _M_coef(__n)
    {
      constexpr _Tp _S_pi = _Tp{3.1415926535897932384626433832795029L};

      if (this->upper() == this->lower())
	std::__throw_domain_error("_Chebyshev: input domain has zero length");

      std::vector<_Tp> __f(__n);

      _Tp __bma = (_M_upper - _M_lower) / _Tp{2};
      _Tp __bpa = (_M_upper + _M_lower) / _Tp{2};
      for (unsigned int __k = 0; __k < __n; ++__k)
	{
	  _Tp __y = std::cos(_S_pi * (__k + _Tp(0.5L)) / __n);
	  __f[__k] = __func(__bpa + __y * __bma);
	}

      _Tp __fac = _Tp{2} / __n;
      for (unsigned int __j = 0; __j < __n; ++__j)
	{
	  _Tp __sum = _Tp{};
	  for (unsigned int __k = 0; __k < __n; ++__k)
	    __sum += __f[__k] * std::cos(_S_pi * __j * (__k + _Tp(0.5L)) / __n);
	  _M_coef[__j] = __fac * __sum;
	}
      return;
    }


  /**
   *  Constructs a Chebyshev expansion from a polynomial.
   *
   *  A polynomial valid in the range x=a to x=b is specified by coefficients d[0..m-1].
   *  An equivalent array of Chebyshev coefficients, c[0..m-1], is output.
   *  The index mfew will be set to the index of the first nonzero Chebyshev coefficient smaller than err.
   *  Then the polynomial is approximated by the first mfew cefficients c[0..mfew-1].
   */
  template<typename _Tp>
    template<typename _Up>
    _Chebyshev<_Tp>::_Chebyshev(_Up __a, _Up __b,
				const std::polynomial<_Up> & __poly)
    : _M_lower(__a),
      _M_upper(__b),
      _M_coef(__poly.order())
    {
      if (this->upper() == this->lower())
	std::__throw_domain_error("_Chebyshev: input domain has zero length");

      //  Map the interval of the polynomial from x:[a,b] to y:[-1,+1].
      poly_shift_coeff((-_Tp{2} - this->upper() - this->lower())
		     / (this->upper() - this->lower()),
		       (+_Tp{2} - this->upper() - this->lower())
		     / (this->upper() - this->lower()), __poly);

      _Tp __pow = _Tp{1};
      this->_M_coef[0] = _Tp{2} * __poly[0];
      for (int __k = 1; __k < __poly.size(); ++__k)
	{
	  this->_M_coef[__k] = _Tp{};
	  _Tp __fac = __poly[__k] / __pow;
	  int __jm = __k;
	  int __jp = 1;
	  for (int __j = __k; __j >= 0; __j -= 2, --__jm, ++__jp)
	    {
	      //  Increment this and lower orders of Chebyshev
	      //  with the combinatorial coefficient times d[k].
	      this->_M_coef[__j] += __fac;
	      __fac *= _Tp(__jm) / _Tp(__jp);
	    }
	  __pow += __pow;
	}

      //  Map the interval of the polynomial from y:[-1,+1] to x:[a,b].
      poly_shift_coeff(this->lower(), this->upper(), __poly);
    }


  /**
   *  Chebyshev evaluation.
   */
  template<typename _Tp>
    _Tp
    _Chebyshev<_Tp>::operator()(_Tp __x) const
    {
      if ((__x - this->lower()) * (__x - this->upper()) > _Tp{})
	std::__throw_domain_error("_Chebyshev::operator(): argument out of range");

      if (this->_M_coef.size() == 0)
	return _Tp{};

      auto __d = _Tp{};
      auto __dd = _Tp{};
      auto __y = (_Tp{2} * __x - this->lower() - this->upper())
	       / (this->upper() - this->lower());
      auto __y2 = _Tp{2} * __y;
      for (unsigned int __j = this->_M_coef.size() - 1; __j >= 1; --__j)
	__dd = std::exchange(__d, __y2 * __d - __dd + this->_M_coef[__j]);

      return __y * __d - __dd + this->_M_coef[0] / _Tp{2};
    }


  /**
   *  This routine returns a Chebyshev fit of the derivative of this Chebyshev fit.
   */
  template<typename _Tp>
    _Chebyshev<_Tp>
    _Chebyshev<_Tp>::derivative() const
    {
      if (this->_M_coef.size() < 3)
	return _Chebyshev<_Tp>{};

      unsigned int __n = this->_M_coef.size();
      std::vector<_Tp> __cder(__n);
      __cder[__n - 1] = _Tp{};
      __cder[__n - 2] = _Tp(2 * (__n - 1)) * this->_M_coef[__n - 1];
      for (unsigned int __j = __n - 3; __j >= 0; --__j)
	__cder[__j] = __cder[__j + 2] + _Tp(2 * (__j + 1)) * this->_M_coef[__j + 1];

      _Tp __con = _Tp{2} / (this->upper() - this->lower());
      for (unsigned int __j = 0; __j < __n; ++__j)
	__cder[__j] *= __con;

      return _Chebyshev<_Tp>(this->lower(), this->upper(), std::move(__cder));
    }


  /**
   *  This routine returns the Chebyshev fit of the integral of this Chebyshev fit.
   *  The constant of integration is set so that the integral vanishes at a.
   */
  template<typename _Tp>
    _Chebyshev<_Tp>
    _Chebyshev<_Tp>::integral() const
    {
      unsigned int __n = this->_M_coef.size();
      std::vector<_Tp> __cint(__n);
      _Tp __sum = _Tp{};
      _Tp __fact = _Tp{1};
      _Tp __con = (this->upper() - this->lower()) / _Tp{4};
      for (unsigned int __j = 1; __j < __n - 1; ++__j)
	{
	  //  Accumulate the constant of integration in sum.
	  //  fact will equal +/-1 to alternate the sign of the terms.
	  __cint[__j] = __con
		* (this->_M_coef[__j - 1] - this->_M_coef[__j + 1]) / __j;
	  __sum += __fact * __cint[__j];
	  __fact = -__fact;
	}
      __cint[__n - 1] = __con * this->_M_coef[__n - 2] / (__n - 1);
      __sum += __fact * __cint[__n - 1];

      //  Set the constant of integration.
      __cint[0] = _Tp{2} * __sum;

      return _Chebyshev<_Tp>(this->lower(), this->upper(), __cint);
    }


  /**
   *  This routine returns the array d[0..n-1], of coefficients of a polynomial
   *  expansion which is equivalent to the Chebyshev fit.
   */
  template<typename _Tp>
    std::polynomial<_Tp>
    _Chebyshev<_Tp>::to_polynomial() const
    {
      unsigned int __n = this->_M_coef.size();

      std::vector<_Tp> __d(__n);
      std::vector<_Tp> __dd(__n);

      __d[0] = this->_M_coef[__n - 1];
      for (unsigned int __j = __n - 2; __j >= 1; --__j)
	{
	  for (unsigned int __k = __n - __j; __k >= 1; --__k)
	    __dd[__k] = std::exchange(__d[__k], _Tp{2} * __d[__k - 1] - __dd[__k]);
	  __dd[0] = std::exchange(__d[0], -__dd[0] + this->_M_coef[__j]);
	}
      for (unsigned int __j = __n - 1; __j >= 1; --__j)
	__d[__j] = __d[__j - 1] - __dd[__j];
      __d[0] = -__dd[0] + this->_M_coef[0] / _Tp{2};

      //  Map the interval [-1,+1] to [a,b].
      poly_shift_coeff(this->lower(), this->upper(), __d);

      return std::polynomial<_Tp>(__d);
    }


  /**
   *  Truncate the Chebyshev series so that the first nonzero neglected term
   *  is less than eps.  Skipping zero coefficients takes care of series
   *  with alternating zero and nonzero terms.
   */
  template<typename _Tp>
    template<typename _Up>
      void
      _Chebyshev<_Tp>::truncate(_Up __eps)
      {
	// @todo Uniform container erasure would be better.
	auto __beg = std::find_if(std::begin(this->_M_coef),
				  std::end(this->_M_coef),
			[__eps](_Tp __c)
			{ return __c != _Tp{} && std::abs(__c) < __eps; });
	this->_M_coef.erase(__beg, std::end(this->_M_coef));
      }


  /**
   *  Write the Chebyshev expansion to an output stream.
   */
  template<typename _Tp,
	   typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& __os,
	       const _Chebyshev<_Tp>& __cheb)
    {
      auto __width = std::numeric_limits<_Tp>::max_digits10;
      auto __prec = __os.precision(__width);
      __os << '{';
      auto __begin = std::begin(__cheb._M_coef);
      auto __end = std::end(__cheb._M_coef);
      if (__begin != __end)
	{
	  __os << *__begin;
	  for (++__begin; __begin != __end; ++__begin)
	    __os << ',' << *__begin;
	}
      __os << '}';
      __os.precision(__prec);
      return __os;
    }

  /**
   *  Economize a polynomial poly over a range [a, b] by returning
   *  a new lower-order polynomial so that evaluations of the new polynomial
   *  will equal those of the old polynomial within the specified tolerance eps
   *  over the requested range.
   */
  template<typename _Tp>
    std::polynomial<_Tp>
    economize_polynomial(_Tp __a, _Tp __b, const std::polynomial<_Tp>& __poly,
			 _Tp __eps)
    {
      if (__poly.order() < 3)
	throw std::domain_error("poly_economize: input arrays too small");

      //  Convert the polynomial into a Chebyshev series.
      _Chebyshev<_Tp> cpoly(__a, __b, __poly);

      //  Truncate the Chebyshev series to the first term less than eps.
      cpoly.truncate(__eps);

      //  Convert the truncated Chebyshev series back into a polynomial.
      return cpoly.to_polynomial();
    }

#endif // _CHEBYSHEV_TCC
