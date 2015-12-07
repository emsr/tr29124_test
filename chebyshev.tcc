#ifndef _CHEBYSHEV_TCC
#define _CHEBYSHEV_TCC 1

#include <utility>
#include <iterator>

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
    {
      constexpr _Tp _S_pi = _Tp{3.1415926535897932384626433832795029L};
      _M_lower = __a;
      _M_upper = __b;
      _M_coef.resize(__n);

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
	  _Tp __sum = _Tp{0};
	  for (unsigned int __k = 0; __k < __n; ++__k)
	    __sum += __f[__k] * std::cos(_S_pi * __j * (__k + _Tp(0.5L)) / __n);
	  _M_coef[__j] = __fac * __sum;
	}
      return;
    }


  /**
   *  Chebyshev evaluation.  All arguments are input.  c[0..m-1] is an array of Chebyshev
   *  coefficients - the first m elements of c output from chebyshev_fit (previously called
   *  with the same a and b).  The Chebyshev polynomial is evaluated at a point y determined
   *  from x, a, b.  The result is returned as the function value.
   */
  template<typename _Tp>
    _Tp
    _Chebyshev<_Tp>::operator()(_Tp __x) const
    {
      if ((__x - this->_M_lower) * (__x - this->_M_upper) > _Tp{0})
	std::__throw_domain_error("_Chebyshev::operator(): Argument x out of range");

      if (this->_M_coef.size() == 0)
	return _Tp{0};

      auto __d = _Tp{0};
      auto __dd = _Tp{0};
      auto __y = (_Tp{2} * __x - this->_M_lower - this->_M_upper)
	       / (this->_M_upper - this->_M_lower);
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
      __cder[__n - 1] = _Tp{0};
      __cder[__n - 2] = _Tp(2 * (__n - 1)) * this->_M_coef[__n - 1];
      for (unsigned int __j = __n - 3; __j >= 0; --__j)
	__cder[__j] = __cder[__j + 2] + _Tp(2 * (__j + 1)) * this->_M_coef[__j + 1];

      _Tp __con = _Tp{2} / (this->_M_upper - this->_M_lower);
      for (unsigned int __j = 0; __j < __n; ++__j)
	__cder[__j] *= __con;

      return _Chebyshev<_Tp>(this->_M_lower, this->_M_upper, std::move(__cder));
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
      _Tp __sum = _Tp{0};
      _Tp __fact = _Tp{1};
      _Tp __con = (this->_M_upper - this->_M_lower) / _Tp{4};
      for (unsigned int __j = 1; __j < __n - 1; ++__j)
	{
	  //  Accumulate the constant of integration in sum.
	  //  fact will equal +/-1 to alternate the sign of the terms.
	  __cint[__j] = __con * (this->_M_coef[__j - 1] - this->_M_coef[__j + 1]) / __j;
	  __sum += __fact * __cint[__j];
	  __fact = -__fact;
	}
      __cint[__n - 1] = __con * this->_M_coef[__n - 2] / (__n - 1);
      __sum += __fact * __cint[__n - 1];

      //  Set the constant of integration.
      __cint[0] = _Tp{2} * __sum;

      return _Chebyshev<_Tp>(this->_M_lower, this->_M_upper, __cint);
    }


  /**
   *  This routine returns the array d[0..n-1], of coefficients of a polynomial
   *  expansion which is equivalent to the Chebyshev fit.
   */
  template<typename _Tp>
    std::polynomial<_Tp>
    _Chebyshev<_Tp>::to_polynomial()
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

      //    Map the interval [-1,+1] to [a,b].
      poly_shift_coeff(this->_M_lower, this->_M_upper, __d);

      return std::polynomial<_Tp>(__d);
    }


  /**
   *  Converts a polynomial to a Chebyshev expansion.
   *
   *  A polynomial valid in the range x=a to x=b is specified by coefficients d[0..m-1].
   *  An equivalent array of Chebyshev coefficients, c[0..m-1], is output.
   *  The index mfew will be set to the index of the first nonzero Chebyshev coefficient smaller than err.
   *  Then the polynomial is approximated by the first mfew cefficients c[0..mfew-1].
   */
  template<typename _Tp>
    template<typename _Up>
    _Chebyshev<_Tp>::_Chebyshev(_Up __a, _Up __b, const std::polynomial<_Up> & __poly)
    {
      this->_M_lower = _Tp(__a);
      this->_M_upper = _Tp(__b);
      this->_M_coef.resize(__poly.order());

      if (this->_M_upper == this->_M_lower)
	std::__throw_domain_error("poly_to_chebyshev: input range has zero length");

      //  Map the interval of the polynomial from x:[a,b] to y:[-1,+1].
      poly_shift_coeff((-_Tp{2} - this->_M_upper - this->_M_lower)
		     / (this->_M_upper - this->_M_lower),
		       (+_Tp{2} - this->_M_upper - this->_M_lower)
		     / (this->_M_upper - this->_M_lower), __poly);

      _Tp __pow = _Tp{1};
      this->_M_coef[0] = _Tp{2} * __poly[0];
      for (int __k = 1; __k < __poly.size(); ++__k)
	{
	  this->_M_coef[__k] = _Tp{0};
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
      poly_shift_coeff(this->_M_lower, this->_M_upper, __poly);
    }

  tempate<typename _Tp,
	  typename _CharT, typename _Traits = std::char_traits<_CharT>>
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& __os,
	       const _Chebyshev<_Tp>& __cheb)
    {
      __os << '{' << '\n';
      for (__cc : __cheb._M_coef)
	__os << ' ' << __cc << ',' << '\n';
      __os << '}' << '\n';
    }

  /**
   *  Economize the polynomial d[0..m-1] by a new polynomial d[0..m-1]
   *  so that the new polynomial evaluated with fewer terms d[0..mfew-1]
   *  will equal the old polynomial within the specified tolerance err.
   *  The index mfew will be set to the index of the first nonzero Chebyshev
   *  coefficient smaller than err.
   *
   *  The polynomial is approximated by the first mfew cefficients d[0..mfew-1].
   */
  template<typename _Tp>
    void
    poly_economize(_Tp __a, _Tp __b, std::vector<_Tp> & __d, std::vector<_Tp> & __c, _Tp __err)
    {
      if (__d.size() != __c.size())
	throw std::domain_error("Bad input arrays in poly_economize");

      //  Convert the shifted polynomial into a Chebyshev series.
      poly_to_chebyshev(__a, __b, __d);

      //  Truncate the Chebyshev series so that the first nonzero neglected term is less than err.
      //  Skipping zero coefficients takes care of series with alternating 0 and nonzero terms.
      //  Truncate the polynomial as well.
      unsigned int __mfew = 0;
      for (unsigned int __k = 0; __k < __c.size(); ++__k)
	{
	  if (__c[__k] == _Tp{0})
	    continue;
	  else if (__c[__k] < __err)
	    {
	      __mfew = __k;
	      break;
	    }
	  else
	    __mfew = __k + 1;
	}
      // @todo Uniform container erasure would be better.
      __d.erase(std::begin(__d) + __mfew, std::end(__d));

      //  Convert the truncated Chebyshev series back into a polynomial.
      chebyshev_to_poly(__a, __b, __c, __d, __mfew);

      return;
    }

#endif // _CHEBYSHEV_TCC
