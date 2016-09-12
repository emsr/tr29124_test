// Math extensions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file ext/roots.tcc
 *  This file is a GNU extension to the Standard C++ Library.
 */

#ifndef _EXT_ROOTS_TCC
#define _EXT_ROOTS_TCC 1

#pragma GCC system_header

#include <stdexcept>
#include <cmath>

#include <ext/math_const.h>

namespace __gnu_cxx
{

  /**
   * Given a function @c __func and an initial range @c __x_lower
   * to @c __x_upper, the routine expands the range geometrically until a root
   * is bracketed by the returned values @c __x_lower and @c __x_upper
   * (in which case __root_bracket returns true)
   * or until the range becomes unacceptably large
   * (in which case __root_bracket returns false).
   * Success is guaranteed for a function which has opposite signs
   * for sufficiently large and small arguments.
   *
   * @param __func  A function
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp, typename _Func>
    bool
    __root_bracket(_Func __func, _Tp& __x_lower, _Tp& __x_upper,
		   std::size_t __max_iter)
    {
      const _Tp __golden = __gnu_cxx::__math_constants<_Tp>::__phi;

      if (__x_lower >= __x_upper)
	std::__throw_logic_error(__N("__root_bracket: bad initial range"));
      auto __f_lower = __func(__x_lower);
      auto __f_upper = __func(__x_upper);
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  if (__f_lower * __f_upper < _Tp{0})
	    return true;
	  if (std::abs(__f_lower) < std::abs(__f_upper))
	    __f_lower = __func(__x_lower += __golden * (__x_lower - __x_upper));
	  else
	    __f_upper = __func(__x_upper += __golden * (__x_upper - __x_lower));
	}
      return false;
    }


  /**
   * Given a function @c __func defined on an interval @c __x_lower
   * to @c __x_upper, the routine subdivides the interval
   * into @c n equally-spaced segments and searches for zero crossings
   * of the function.
   *
   * @param __func  A function
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __n  The number of subdivisions of the interval
   * @param __xb  The output vector of root bounds
   */
  template<typename _Tp, typename _Func>
    std::vector<std::pair<_Tp, _Tp>>
    __root_brackets(_Func __func,
		    _Tp __x_lower, _Tp __x_upper, std::size_t __n)
    {
      std::vector<std::pair<_Tp, _Tp>> __xb;
      auto __dx = (__x_upper - __x_lower) / _Tp(__n);
      auto __x = __x_lower;
      auto __f_lo = __func(__x);
      for (std::size_t __i = 0; __i < __n; ++__i)
	{
	  __x += __dx;
	  auto __f_hi = __func(__x);
	  if (__f_lo * __f_hi <= _Tp{0})
	    __xb.emplace_back(__x - __dx, __x);
	  __f_lo = __f_hi;
	}

      return __xb;
    }


  /**
   * Using bisection, find the root of a function @c __func
   * known to lie between @c __x_lower and @c __x_upper.
   * The returned root is refined until its accuracy is <tt>+/- __eps</tt>.
   *
   * @param __func  A function
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp, typename _Func>
    _Tp
    __root_bisect(_Func __func, _Tp __x_lower, _Tp __x_upper,
		  _Tp __eps, std::size_t __max_iter)
    {
      auto __f = __func(__x_lower);
      auto __f_mid = __func(__x_upper);
      if (__f * __f_mid > _Tp{0})
	std::__throw_logic_error(__N("__root_bisect: "
				 "Root must be bracketed for bisection"));
      //  Orient search so that f > _Tp{0} lies at x + dx.
      _Tp __dx;
      auto __x = __f < _Tp{0}
	       ? (__dx = __x_upper - __x_lower, __x_lower)
	       : (__dx = __x_lower - __x_upper, __x_upper);
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  __dx /= _Tp{2};
	  auto __x_mid = __x + __dx;
	  __f_mid = __func(__x_mid);
	  if (__f_mid < _Tp{0})
	    __x = __x_mid;
	  if (std::abs(__dx) < __eps || __f_mid == _Tp{0})
	    return __x;
	}

      std::__throw_logic_error(__N("__root_bisect: "
			       "Maximum number of bisections exceeded"));
    }


  /**
   * Using the secant method, find the root of a function @c __func
   * known to lie between @c __x_lower and @c __x_upper.
   * The returned root is refined until its accuracy is <tt>+/- __eps</tt>.
   *
   * @param __func  A function
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp, typename _Func>
    _Tp
    __root_secant(_Func __func, _Tp __x_lower, _Tp __x_upper,
		  _Tp __eps, std::size_t __max_iter)
    {
      _Tp __x_lo, __x;

      auto __f_lo = __func(__x_lower);
      auto __f = __func(__x_upper);
      if (std::abs(__f_lo) < std::abs(__f))
	{
	  __x = __x_lower;
	  __x_lo = __x_upper;
	  std::swap(__f_lo, __f);
	}
      else
	{
	  __x_lo = __x_lower;
	  __x = __x_upper;
	}

      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  auto __dx = (__x_lo - __x) * __f / (__f - __f_lo);
	  __x_lo = __x;
	  __f_lo = __f;
	  __x += __dx;
	  __f = __func(__x);
	  if (std::abs(__dx) < __eps || __f == _Tp{0})
	    return __x;
	}

      std::__throw_logic_error(__N("__root_secant: "
			       "Maximum number of iterations exceeded"));
    }


  /**
   * Using the false position method, find the root of a function @c __func
   * known to lie between @c __x_lower and @c __x_upper.
   * The returned root is refined until its accuracy is <tt>+/- __eps</tt>.
   *
   * @param __func  A function
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp, typename _Func>
    _Tp
    __root_false_position(_Func __func, _Tp __x_lower, _Tp __x_upper,
			  _Tp __eps, std::size_t __max_iter)
    {
      auto __f_lo = __func(__x_lower);
      auto __f_hi = __func(__x_upper);
      if (__f_lo * __f_hi > _Tp{0})
	std::__throw_logic_error(__N("__root_false_position: "
				    "Root must be bracketed"));

      _Tp __x_lo, __x_hi;
      if (__f_lo < _Tp{0})
	{
	  __x_lo = __x_lower;
	  __x_hi = __x_upper;
	}
      else
	{ 
	  __x_lo = __x_upper;
	  __x_hi = __x_lower;
	  std::swap(__f_lo, __f_hi);
	}

      auto __dx = __x_hi - __x_lo;
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  auto __x = __x_lo + __dx * __f_lo / (__f_lo - __f_hi);
	  auto __f = __func(__x);
	  _Tp __del;
	  if (__f < _Tp{0})
	    {
	      __del = __x_lo - __x;
	      __x_lo = __x;
	      __f_lo = __f;
	    }
	  else
	    {
	      __del = __x_hi - __x;
	      __x_hi = __x;
	      __f_hi = __f;
	    }
	  __dx = __x_hi - __x_lo;
	  if (std::abs(__del) < __eps || __f == _Tp{0})
	    return __x;
	}

      std::__throw_logic_error(__N("__root_false_position: "
			       "Maximum number of iterations exceeded"));
    }


  /**
   * Using Ridder's method, find the root of a function @c __func known
   * to lie between @c __x_lower and @c __x_upper.
   * The returned root is refined until its accuracy is <tt>+/- __eps</tt>.
   *
   * @param __func  A function
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp, typename _Func>
    _Tp
    __root_ridder(_Func __func, _Tp __x_lower, _Tp __x_upper,
		  _Tp __eps, std::size_t __max_iter)
    {
      auto __f_hi = __func(__x_upper);
      auto __f_lo = __func(__x_lower);
      auto __x_hi = __x_upper;
      auto __x_lo = __x_lower;

      if (__f_lo * __f_hi > _Tp{0})
	std::__throw_logic_error(__N("__root_ridder: Root must be bracketed"));

      if (__f_lo == _Tp{0})
	return __x_lower;

      if (__f_hi == _Tp{0})
	return __x_upper;

      const _Tp __UNUSED = -1.0e30; // an exceedingly unlikely answer.
      auto __ans = __UNUSED;
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  auto __xm = (__x_lo + __x_hi) / _Tp{2};
	  auto __fm = __func(__xm);
	  auto __s = std::sqrt(__fm * __fm - __f_lo * __f_hi);
	  if (__s == _Tp{0})
	    return __ans;
	  auto __xnew = __xm + (__xm - __x_lo)
		      * (__f_lo >= __f_hi ? 1.0 : -1.0) * __fm / __s;
	  if (std::abs(__xnew - __ans) < __eps)
	    return __ans;
	  __ans = __xnew;
	  auto __fnew = __func(__ans);
	  if (__fnew == _Tp{0})
	    return __ans;
	  if (std::copysign(__fm, __fnew) != __fm)
	    {
	      __x_lo = __xm;
	      __f_lo = __fm;
	      __x_hi = __xnew;
	      __f_hi = __fnew;
	    }
	  else if (std::copysign(__f_lo, __fnew) != __f_lo)
	    {
	      __x_hi = __ans;
	      __f_hi = __fnew;
	    }
	  else if (std::copysign(__f_hi, __fnew) != __f_hi)
	    {
	      __x_lo = __ans;
	      __f_lo = __fnew;
	    }
	  else
	    std::__throw_logic_error(__N("__root_ridder: "
				     "Some major malfunction"));

	  if (std::abs(__x_hi - __x_lo) < __eps)
	    return __ans;
	}

      std::__throw_logic_error(__N("__root_ridder: "
			       "Maximum number of iterations exceeded"));
    }


  /**
   * Using Brent's method, find the root of a function @c __func known
   * to lie between @c __x_lower and @c __x_upper.
   * The returned root is refined until it's accuracy is <tt>+/- __eps</tt>.
   *
   * @param __func  A function
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp, typename _Func>
    _Tp
    __root_brent(_Func __func, _Tp __x_lower, _Tp __x_upper,
		 _Tp __eps, std::size_t __max_iter)
    {
      auto __x_a = __x_lower;
      auto __x_b = __x_upper;
      auto __x_c = __x_upper;
      auto __f_a = __func(__x_a);
      auto __f_b = __func(__x_b);

      const _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();

      if (__f_a * __f_b > _Tp{0})
	std::__throw_logic_error(__N("__root_brent: Root must be bracketed"));

      auto __f_c = __f_b;
      for (std::size_t __iter = 0; __iter < __max_iter; ++__iter)
	{
	  _Tp __d, __e;
	  if (__f_b * __f_c > _Tp{0})
	    {
	      __x_c = __x_a;
	      __f_c = __f_a;
	      __e = __d = __x_b - __x_a;
	    }
	  if (std::abs(__f_c) < std::abs(__f_b))
	    {
	      __x_a = __x_b;
	      __x_b = __x_c;
	      __x_c = __x_a;
	      __f_a = __f_b;
	      __f_b = __f_c;
	      __f_c = __f_a;
	    }
	  auto __toler = _Tp{2} * _S_eps * std::abs(__x_b) + __eps / _Tp{2};
	  auto __xm = (__x_c - __x_b) / _Tp{2};
	  if (std::abs(__xm) <= __toler || __f_b == _Tp{0})
	    return __x_b;
	  if (std::abs(__e) >= __toler && std::abs(__f_a) > std::abs(__f_b))
	    {
	      _Tp __p, __q, __r;
	      auto __s = __f_b / __f_a;
	      if (__x_a == __x_c)
		{
		  __p = _Tp{2} * __xm * __s;
		  __q = _Tp{1} - __s;
		}
	      else
		{
		  __q = __f_a / __f_c;
		  __r = __f_b / __f_c;
		  __p = __s * (_Tp{2} * __xm * __q * (__q - __r)
		      - (__x_b - __x_a) * (__r - _Tp{1}));
		  __q = (__q - _Tp{1}) * (__r - _Tp{1}) * (__s - _Tp{1});
		}
	      if (__p > _Tp{0})
		__q = -__q;
	      __p = std::abs(__p);
	      auto __min1 = _Tp{3} * __xm * __q - std::abs(__toler * __q);
	      auto __min2 = std::abs(__e * __q);
	      if (_Tp{2} * __p < std::min(__min1, __min2))
		{
		  __e = __d;
		  __d = __p / __q;
		}
	      else
		{
		  __d = __xm;
		  __e = __d;
		}
	    }
	  else
	    {
	      __d = __xm;
	      __e = __d;
	    }
	  __x_a = __x_b;
	  __f_a = __f_b;
	  if (std::abs(__d) > __toler)
	    __x_b += __d;
	  else
	    __x_b += std::copysign(__toler, __xm);
	  __f_b = __func(__x_b);
	}

      std::__throw_logic_error(__N("__root_brent: "
			       "Maximum number of iterations exceeded"));
    }


  /**
   * Using the Newton-Raphson method, find the root of a function @c __func
   * known to lie in the interval @c __x_lower to @c __x_upper.
   * The returned root is refined until its accuracy is
   * within <tt>+/- __eps</tt>.
   *
   * @param __func A routine that provides both the function
   *               and the first derivative of the function at the point x.
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp, typename _StateFunc>
    _Tp
    __root_newton(_StateFunc __func, _Tp __x_lower, _Tp __x_upper,
		  _Tp __eps, std::size_t __max_iter)
    {
      auto __x = (__x_lower + __x_upper) / _Tp{2};
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  auto __sf = __func(__x);
	  auto __dx = __sf.__value / __sf.__deriv;
	  __x -= __dx;
	  if ((__x_lower - __x) * (__x - __x_upper) < _Tp{0})
	    std::__throw_logic_error(__N("__root_newton: "
				     "Jumped out of brackets"));
	  if (std::abs(__dx) < __eps)
	    return __x;
	}
      std::__throw_logic_error(__N("__root_newton: "
			       "Maximum number of iterations exceeded"));
    }


  /**
   * Using a combination of Newton-Raphson and bisection, find the root
   * of a function @c __func known to lie in the interval @c __x_lower
   * to @c __x_upper.
   * The returned root is refined until its accuracy is
   * within <tt>+/- __eps</tt>.
   *
   * @param __func  A routine that provides both the function
   *                and the first derivative of the function at the point x.
   * @param __x_lower  The lower end of the interval
   * @param __x_upper  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp, typename _StateFunc>
    _Tp
    __root_safe(_StateFunc __func, _Tp __x_lower, _Tp __x_upper,
		_Tp __eps, std::size_t __max_iter)
    {
      auto __sf_lo = __func(__x_lower);
      auto __sf_hi = __func(__x_upper);

      if (__sf_lo.__value * __sf_hi.__value > _Tp{0})
	std::__throw_logic_error(__N("__root_safe: Root must be bracketed"));

      if (__sf_lo.__value == _Tp{0})
	return __x_lower;
      if (__sf_hi.__value == _Tp{0})
	return __x_upper;

      _Tp __x_hi, __x_lo;
      if (__sf_lo.__value < _Tp{0})
	{
	  __x_lo = __x_lower;
	  __x_hi = __x_upper;
	}
      else
	{
	  __x_hi = __x_lower;
	  __x_lo = __x_upper;
	}

      auto __x = (__x_lower + __x_upper) / _Tp{2};
      auto __dxold = std::abs(__x_upper - __x_lower);
      auto __dx = __dxold;
      auto __sf = __func(__x);
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  if (((__x - __x_hi) * __sf.__deriv - __sf.__value)
	    * ((__x - __x_lo) * __sf.__deriv - __sf.__value) > _Tp{0}
	   || std::abs(_Tp{2} * __sf.__value)
	    > std::abs(__dxold * __sf.__deriv))
	    {
	      __dxold = __dx;
	      __dx = (__x_hi - __x_lo) / _Tp{2};
	      __x = __x_lo + __dx;
	      if (__x == __x_lo)
		return __x;
	    }
	  else
	    {
	      __dxold = __dx;
	      __dx = __sf.__value / __sf.__deriv;
	      auto __temp = __x;
	      __x -= __dx;
	      if (__temp == __x)
		return __x;
	    }
	  if (std::abs(__dx) < __eps)
	    return __x;

	  __sf = __func(__x);
	  if (__sf.__value < _Tp{0})
	    __x_lo = __x;
	  else
	    __x_hi = __x;
	}

      std::__throw_logic_error(__N("__root_safe: "
			       "Maximum number of iterations exceeded"));
    }

/*
  template<typename _Tp>
    void
    __root_laguerre(_Polymomial<std::complex<_Tp>>& __a, std::complex<_Tp>& __x, int& __its)
    {
      // Estimated fractional roundoff error.
      constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();

      // Number of fractional values.
      constexpr int MR = 8;
      // Fractions used to break a limit cycle.
      static const _Tp
      _S_frac[MR + 1]
      {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};

      // Number of steps taken before trying a new fraction.
      constexpr int MT = 10;

      constexpr int _S_max_iter = MT * MR;

      Complex d, f;
      int __m = __a.degree();
      for (int __iter = 1;__iter <= _S_max_iter; ++__iter)
	{ // Loop over iterations up to allowed maximum.
	  __its = __iter;
	  auto __b = __a[__m];
	  auto __err = std::abs(__b);
	  std::complex<_Tp> __d{}, __f{};
	  auto __abx = std::abs(__x);
	  for (int __j = __m - 1; __j >= 0; --__j)
	    {
	      // Efficient computation of the polynomial and its first two derivatives.
	      // f stores P''(x)/2.
	      __f = __x * __f + __d;
	      __d = __x * __d + __b;
	      __b = __x * __b + __a[__j];
	      __err = __abx * __err + std::abs(__b);
	    }
	  __err *= ___S_eps;
	  // Estimate of roundoff error in evaluating polynomial.
	  if (std::abs(__b) <= __err) // We have the root.
	    return;
	  // Use Laguerre's formula.
	  auto __g = __d / __b;
	  auto __g2 = __g * __g;
	  auto __h = __g2 - _Tp{2} * __f / __b;
	  auto __sq = std::sqrt(_Tp(__m - 1) * (_Tp(__m) * __h - __g2));
	  auto __gp = __g + __sq;
	  auto __gm = __g - __sq;
	  auto __abp = std::abs(__gp);
	  auto __abm = std::abs(__gm);
	  if (__abp < __abm)
	    __gp = __gm;
	  auto __dx = std::max(__abp, __abm) > _Tp{0}
		    ? _Tp(__m) / __gp
		    : std::polar(_Tp{1} + __abx, _Tp(__iter));
	  auto __x1 = __x - __dx;
	  if (__x == __x1)
	    return;
	  if (__iter % MT != 0)
	    __x = __x1;
	  else
	    __x -= ___S_frac[__iter / MT] * __dx;
	}
      std::__throwlogic_error("__root_laguerre: "
			      "Maximum number of iterations exceeded");
    }

  template<typename _Tp>
    void
    qroot(_Polynomial<std::complex<_Tp>>& __p, _Tp& __b, _Tp& __c, _Tp __eps)
    {
      using _Poly = _Polynomial<std::complex<_Tp>>;
      constexpr int _S_max_iter = 20;
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_tiny = _Tp{100} * _S_eps;
      auto __n = __p.order();
      _Poly __q, __qq, __rem;
      for (int __iter = 0; __iter < _S_max_iter; ++__iter)
	{
	  _Poly __d(__c, __b, _Tp{1});

	  // First division: r, s.
	  divmod(__p, __d, __q, __rem);
	  auto __s = __rem[0];
	  auto __r = __rem[1];
	  // Second division: partial r, s with respect to c.
	  divmod(__q, __d, __qq, __rem);
	  auto __sc = -__rem[0];
	  auto __rc = -__rem[1];
	  auto __sb = -__c * __rc;
	  auto __rb = -__b * __rc + __sc;
	  // Solve 2x2 equation.
	  auto __dv = _Tp{1} / (__sb * __rc - __sc * __rb);
	  auto __delb = ( __r * __sc - __s * __rc) * __dv;
	  auto __delc = (-__r * __sb + __s * __rb) * __dv;
	  __b += __delb;
	  __delc = (-__r * __sb + __s * __rb) * __dv;
	  __c += __delc;
	  if ((std::abs(__delb) <= eps * std::abs(__b) || std::abs(__b) < _S_tiny)
           && (std::abs(__delc) <= eps * std::abs(__c) || std::abs(__c) < _S_tiny))
	    return;
	}
      std::__throw_logic_error("qroot: "
			       "Maximum number of iterations exceeded");
    }
*/

} // namespace __gnu_cxx

#endif // _EXT_ROOTS_TCC
