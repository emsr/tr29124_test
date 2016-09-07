#ifndef ROOTS_TCC
#define ROOTS_TCC 1


#include <stdexcept>
#include <cmath>

#include "ext/math_const.h"

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
  template<typename _Tp>
    bool
    __root_bracket(_Tp (*__func)(_Tp), _Tp& __x_lower, _Tp& __x_upper,
		   std::size_t __max_iter)
    {
      const _Tp __golden = __gnu_cxx::__math_constants<_Tp>::__phi;

      if (__x_lower >= __x_upper)
	std::__throw_logic_error(__N("__root_bracket: bad initial range"));
      auto __f1 = __func(__x_lower);
      auto __f2 = __func(__x_upper);
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  if (__f1 * __f2 < _Tp{0})
	    return true;
	  if (std::abs(__f1) < std::abs(__f2))
	    __f1 = __func(__x_lower += __golden * (__x_lower - __x_upper));
	  else
	    __f2 = __func(__x_upper += __golden * (__x_upper - __x_lower));
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
  template<typename _Tp>
    std::vector<std::pair<_Tp, _Tp>>
    __root_brackets(_Tp (*__func)(_Tp),
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
  template<typename _Tp>
    _Tp
    __root_bisect(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper,
		  _Tp __eps, std::size_t __max_iter)
    {
      auto __f = __func(__x_lower);
      auto __f_mid = __func(__x_upper);
      if (__f * __f_mid >= _Tp{0})
	std::__throw_logic_error(__N("__root_bisect: "
				 "Root must be bracketed for bisection"));
      //  Orient search so that f > _Tp{0} lies at x + dx.
      _Tp __dx;
      auto __x = __f < _Tp{0}
	       ? (__dx = __x_upper - __x_lower, __x_lower)
	       : (__dx = __x_lower - __x_upper, __x_upper);
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  auto __x_mid = __x + (__dx *= 0.5);
	  __f_mid = __func(__x_mid);
	  if (__f_mid < _Tp{0})
	    __x = __x_mid;
	  if (std::abs(__dx) < __eps || __f_mid == _Tp{0})
	    return __x;
	}
      std::__throw_logic_error(__N("__root_bisect: "
			       "Maximum number of bisections exceeded"));

      return _Tp{0};
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
  template<typename _Tp>
    _Tp
    __root_secant(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper,
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
	  __x_lo = __x;
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

      return  _Tp{0};
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
  template<typename _Tp>
    _Tp
    __root_false_position(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper,
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

      return _Tp{0};
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
  template<typename _Tp>
    _Tp
    __root_ridder(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper,
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
	  auto __xm = 0.5 * (__x_lo + __x_hi);
	  auto __fm = __func(__xm);
	  auto __s = std::sqrt(__fm * __fm - __f_lo * __f_hi);
	  if (__s == _Tp{0})
	    return __ans;
	  auto __xnew = __xm + (__xm - __x_lo)
		      * (__f_lo >= __f_hi ? 1.0 : -1.0) * __fm / __s;
	  if (std::abs(__xnew - __ans) < __eps)
	    return __ans;
	  auto __fnew = __func(__ans = __xnew);
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
	      __x_hi = __xnew;
	      __f_hi = __fnew;
	    }
	  else if (std::copysign(__f_hi, __fnew) != __f_hi)
	    {
	      __x_lo = __xnew;
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
  template<typename _Tp>
    _Tp
    __root_brent(_Tp (*__func)(_Tp), _Tp __x_lower, _Tp __x_upper,
		 _Tp __eps, std::size_t __max_iter)
    {
      auto __x_a = __x_lower;
      auto __x_b = __x_upper;
      auto __x_c = __x_upper;
      auto __f_a = __func(__x_a);
      auto __f_b = __func(__x_b);

      const _Tp __EPS = 1.0e-12;

      if (__f_b * __f_a > _Tp{0})
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
	  auto __tol1 = _Tp{2} * __EPS * std::abs(__x_b) + __eps / _Tp{2};
	  auto __xm = (__x_c - __x_b) / _Tp{2};
	  if (std::abs(__xm) <= __tol1 || __f_b == _Tp{0})
	    return __x_b;
	  if (std::abs(__e) >= __tol1 && std::abs(__f_a) > std::abs(__f_b))
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
	      auto __min1 = _Tp{3} * __xm * __q - std::abs(__tol1 * __q);
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
	  if (std::abs(__d) > __tol1)
	    __x_b += __d;
	  else
	    __x_b += std::copysign(__tol1, __xm);
	  __f_b = __func(__x_b);
	}
      std::__throw_logic_error(__N("__root_brent: "
			       "Maximum number of iterations exceeded"));

      return  _Tp{0};
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
  template<typename _Tp>
    _Tp
    __root_newton(void (*__func)(_Tp, _Tp*, _Tp*), _Tp __x_lower, _Tp __x_upper,
		  _Tp __eps, std::size_t __max_iter)
    {
      auto __x = (__x_lower + __x_upper) / _Tp{2};
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  _Tp __df, __f;
	  __func(__x, &__f, &__df);
	  auto __dx = __f / __df;
	  __x -= __dx;
	  if ((__x_lower - __x) * (__x - __x_upper) < _Tp{0})
	    std::__throw_logic_error(__N("__root_newton: "
				     "Jumped out of brackets"));
	  if (std::abs(__dx) < __eps)
	    return __x;
	}
      std::__throw_logic_error(__N("__root_newton: "
			       "Maximum number of iterations exceeded"));

      return  _Tp{0};
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
  template<typename _Tp>
    _Tp
    __root_safe(void (*__func)(_Tp, _Tp*, _Tp*), _Tp __x_lower, _Tp __x_upper,
		_Tp __eps, std::size_t __max_iter)
    {
      _Tp __df, __f_hi, __f_lo;
      __func(__x_lower, &__f_lo, &__df);
      __func(__x_upper, &__f_hi, &__df);

      if (__f_lo * __f_hi > _Tp{0})
	std::__throw_logic_error(__N("__root_safe: Root must be bracketed"));

      if (__f_lo == _Tp{0})
	return __x_lower;
      if (__f_hi == _Tp{0})
	return __x_upper;

      _Tp __x_hi, __x_lo;
      if (__f_lo < _Tp{0})
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
      _Tp __f;
      __func(__x, &__f, &__df);
      for (std::size_t __i = 0; __i < __max_iter; ++__i)
	{
	  if (((__x - __x_hi) * __df - __f)
	    * ((__x - __x_lo) * __df - __f) >= _Tp{0}
	   || std::abs(_Tp{2} * __f) > std::abs(__dxold * __df))
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
	      __dx = __f / __df;
	      auto __temp = __x;
	      __x -= __dx;
	      if (__temp == __x)
		return __x;
	    }
	  if (std::abs(__dx) < __eps)
	    return __x;

	  __func(__x, &__f, &__df);
	  if (__f < _Tp{0})
	    __x_lo = __x;
	  else
	    __x_hi = __x;
	}

      std::__throw_logic_error(__N("__root_safe: "
			       "Maximum number of iterations exceeded"));

      return  _Tp{0};
    }

} // namespace __gnu_cxx

#endif // ROOTS_TCC
