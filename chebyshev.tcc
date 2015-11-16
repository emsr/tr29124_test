

#include <vector>
#include <iosfwd>

#include "polynomial"

  template<typename _Tp>
    class _Chebyshev
    {
      _Chebyshev(_Tp __a, _Tp __b, unsigned __n, _Tp __func(_Tp));

      _Chebyshev(_Tp __a, _Tp __b, std::initializer_list<_Tp> __il)
      : _M_lower{__a}, _M_upper{__b},
	_M_coef{__il}
      { }

      _Tp operator()(_Tp __x) const;

    private:
      _Tp _M_lower;
      _Tp _M_upper;
      std::vector<_Tp> & _M_coef;
    };


  /**
   *  @brief  Construct a Chebyshev fit.  Given a function, the lower limit a, the upper limit b,
   *    and the number of points n this routine computes the coefficients c[0..n-1] of a 
   *    Chebyshev polynomial expansion such that func(x) =~ kSum_0^n-1 c_k T_k(x) - c_0/2.
   *    Use this routine with moderately large n of 30 or 50.  Then you can truncate the
   *    series to a smaller number of terms to satisfy the accuracy requirements.
   */
  template<typename _Tp>
    void
    _Chebyshev<_Tp>::_Chebyshev(_Tp __a, _Tp __b, unsigned __n, _Tp __func(_Tp))
    {
      _M_lower = __a;
      _M_upper = __b;
      _M_coef.resize(__n);

      std::vector<_Tp> __f(__n);

      _Tp __bma = (_M_upper - _M_lower) / _Tp(2);
      _Tp __bpa = (_M_upper + _M_lower) / _Tp(2);
      for (unsigned int __k = 0; __k < __n; ++__k)
	{
	  _Tp __y = std::cos(PI * (__k + _Tp(0.5L)) / __n);
	  f[__k] = func(__bpa + __y * __bma);
	}

      _Tp __fac = _Tp(2) / __n;
      for (unsigned int __j = 0; __j < __n; ++__j)
	{
	  _Tp __sum = _Tp(0);
	  for (unsigned int __k = 0; __k < __n; ++__k)
	    sum += __f[__k] * std::cos(PI * __j * (__k + _Tp(0.5L)) / __n);
	  _M_coef[__j] = __fac * __sum;
	}
      return;
    }


  /**
   *    Chebyshev evaluation.  All arguments are input.  c[0..m-1] is an array of Chebyshev
   *    coefficients - the first m elements of c output from chebyshev_fit (previously called
   *    with the same a and b).  The Chebyshev polynomial is evaluated at a point y determined
   *    from x, a, b.  The result is returned as the function value.
   */
  template<typename _Tp>
    _Tp
    _Chebyshev<_Tp>::operator()(_Tp x)
    {
      if ((__x - _M_lower) * (__x - _M_upper) > _Tp(0))
	throw std::domain_error("Argument x out of range in _Chebyshev::operator().");

      if (_M_coef.size() == 0)
	return _Tp(0);

      _Tp __d = _Tp(0);
      _Tp __dd = _Tp(0);
      _Tp __y = (_Tp(2) * __x - _M_lower - _M_upper) / (_M_upper - _M_lower);
      _Tp __y2 = _Tp(2) * __y;
      for (unsigned int __j = _M_coef.size() - 1; __j >= 1; --__j)
	{
	  _Tp __sv = __d;
	  __d = __y2 * __d - __dd + coef[__j];
	  __dd = __sv;
	}

      return __y * __d - __dd + _M_coef[0] / _Tp(2);
    }


  /**
   *  This routine returns the array cder[0..n-1], the Chebyshev coefficients
   *  of the derivative of the function.
   */
  template<typename _Tp>
    void
    _Chebyshev<_Tp>::differentiate(std::vector<_Tp> & __cder)
    {
      __cder.resize(_M_coef.size());

      if (_M_coef.size() < 3)
	return;

      unsigned int __n = _M_coef.size();
      cder[__n - 1] = _Tp(0);
      cder[__n - 2] = _Tp(2 * (__n - 1)) * _M_coef[__n - 1];
      for (unsigned int __j = __n - 3; __j >= 0; --__j)
	__cder[j] = __cder[__j + 2] + _Tp(2 * (__j + 1)) * _M_coef[__j + 1];

      _Tp __con = _Tp(2) / (_M_upper - _M_lower);
      for (unsigned int __j = 0; __j < __n; ++__j)
	__cder[__j] *= __con;

      return;
    }


  /**
   *  This routine returns the array cint[0..m-1], the Chebyshev coefficients
   *  of the integral of the function.
   *  The constant of integration is set so that the integral vanishes at a.
   */
  template<typename _Tp>
    void
    _Chebyshev<_Tp>::integrate(std::vector<_Tp> & __cint)
    {
      __cint.resize(_M_coef.size());

      unsigned int __n = _M_coef.size();
      _Tp __sum = _Tp(0);
      _Tp __fact = _Tp(1);
      _Tp __con = (_M_upper - _M_lower) / _Tp(4);
      for (unsigned int __j = 1; __j < __n - 1; ++__j)
	{
	  //  Accumulate the constant of integration in sum.
	  //  fact will equal +/-1 to alternate the sign of the terms.
	  __cint[j] = __con * (_M_coef[__j - 1] - _M_coef[__j + 1]) / __j;
	  __sum += __fact * __cint[__j];
	  __fact = -__fact;
	}
      __cint[m - 1] = __con * _M_coef[__n - 2] / (__n - 1);
      __sum += __fact * __cint[__n - 1];

      //  Set the constant of integration.
      __cint[0] = _Tp(2) * __sum;

      return;
    }


  /**
   *  This routine returns the array d[0..n-1], of coefficients of a polynomial
   *  expansion which is equivalent to the Chebyshev fit.
   */
  template<typename _Tp>
    void
    _Chebyshev<_Tp>::to_polynomial(std::vector<_Tp> & __d)
    {
      unsigned int __n = _M_coef.size()

      __d.resize(__n);
      std::vector<_Tp> __dd(__n);

      __d[0] = _M_coef[__n - 1];
      for (unsigned int __j = __n - 2; __j >= 1; --__j)
	{
	  _Tp __sv;
	  for (unsigned int __k = __n - __j; __k >= 1; --__k)
	    {
	      __sv = __d[__k];
	      __d[k] = 2.0 * __d[__k-1] - __dd[__k];
	      __dd[k] = __sv;
	    }
	  __sv = __d[0];
	  __d[0] = -__dd[0] + _M_coef[__j];
	  __dd[0] = __sv;
	}
      for (unsigned int __j = __n - 1; __j >= 1; --__j)
	__d[j] = __d[__j - 1] - __dd[__j];
      __d[0] = -__dd[0] + _M_coef[0] / _Tp(2);

      //    Map the interval [-1,+1] to [a,b].
      poly_shift_coeff(a, b, d);

      return;
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
    void
    poly_to_chebyshev(_Tp a, _Tp b, const std::vector<_Tp> & d, std::vector<_Tp> & c)
    {
      if (b == a)
	nrerror("poly_to_chebyshev: input range has zero length");

      //  Map the interval of the polynomial from x:[a,b] to y:[-1,+1].
      poly_shift_coeff((-_Tp(2) - b - a) / (b - a), (+_Tp(2) - b - a) / (b - a), d, m);

      _Tp pow = _Tp(1);
      c[0] = _Tp(2) * d[0];
      for (int k = 1; k < d.size(); ++k)
	{
	  c[k] = _Tp(0);
	  _Tp fac = d[k] / pow;
	  int jm = k;
	  int jp = 1;
	  for (int j = k; j >= 0; j -= 2, --jm, ++jp)
	    {
	      //  Increment this and lower orders of Chebyshev with the combinatorial
	      //  coefficient times d[k].
	      c[j] += fac;
	      fac *= _Tp(jm) / _Tp(jp);
	    }
	  pow += pow;
	}

      //  Map the interval of the polynomial from y:[-1,+1] to x:[a,b].
      poly_shift_coeff(a, b, d);
    }


  /**
   *  Economize the polynomial d[0..m-1] by a new polynomial d[0..m-1] so that the new polynomial
   *  evaluated with fewer terms d[0..mfew-1] will equal the old polynomial within the specified tolerance err.
   *  The index mfew will be set to the index of the first nonzero Chebyshev coefficient smaller than err.
   *  Then the polynomial is approximated by the first mfew cefficients d[0..mfew-1].
   */
  template<typename _Tp>
    void
    poly_economize(_Tp a, _Tp b, std::vector<_Tp> & d, std::vector<_Tp> & c, int & mfew, _Tp err)
    {
      if (d.size() != c.size())
	throw std::domain_error("Bad input arrays in poly_economize");

      //    Convert the shifted polynomial into a Chebyshev series.
      poly_to_chebyshev(a, b, d, c);

      //  Truncate the Chebyshev series so that the first nonzero neglected term is less than err.
      //  Skipping zero coefficients takes care of series with alternating 0 and nonzero terms.
      //  Truncate the polynomial as well.
      for (mfew = 0, k = 0; k < c.size(); ++k)
	{
	  if (c[k] == _Tp(0))
	    continue;
	  else if (c[k] < err)
	    {
	      mfew = k;
	      break;
	    }
	  else
	   mfew = k + 1;
	}
      for (int k = mfew; k < c.size(); ++k)
	d[k] = c[k] = _Tp(0);

      //  Convert the truncated Chebyshev series back into a polynomial.
      chebyshev_to_poly(a, b, c, d, mfew);

      return;
    }


  /**
   *  Evaluate the definite integral of a function fitted by Chebyshev polynomials
   *  over the interval from a to b.  The output coefficients from chebyshev_fit
   *  are in c.  The maximum order used for the integral is m.
   */
  template<typename _Tp>
    _Tp
    clenshaw_curtis_quad(_Tp a, _Tp b, const std::vector<_Tp> & c, _Tp eps)
    {
      if (eps <= _Tp(0))
	nrerror("Error tolerance eps must be greater than 0 in clenshaw_curtis_quad.");

      _Tp sum = _Tp(0);
      for (int k = 0; (2*k + 1) < c.size(); ++k)
	sum -= c[2 * k + 1] / _Tp((2 * k + 1) * (2 * k - 1));
      sum -= c[1] / _Tp(2);

      return (b - a) * sum;
    }


