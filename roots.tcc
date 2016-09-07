#ifndef ROOTS_TCC
#define ROOTS_TCC 1


#include <cmath>
#include <vector>
#include <utility>

#include "ext/math_const.h"


  /**
   * Given a function func and an initial guessd range x1 to x2, the routine
   * expands the range geometrically until a root is bracketed by the
   * returned values x1 and x2 (in which case bracket returns 1)
   * or until the range becomes unacceptably large (in which case zbrak
   * returns 0).  Success is guaranteed for a function which has opposite signs
   * for sufficiently large and small arguments.
   * @param __func  A function
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp>
    bool
    __root_bracket(_Tp __func(_Tp), _Tp& __x1, _Tp& __x2, int __max_iter = 50)
    {
      const _Tp __factor = __gnu_cxx::__math_constants<_Tp>::__phi;

      if (__x1 >= __x2)
	std::__throw_logic_error("__root_bracket: bad initial range");
      auto __f1 = __func(__x1);
      auto __f2 = __func(__x2);
      for (int __i = 0; __i < __max_iter; ++__i)
	{
	  if (__f1 * __f2 < _Tp{0})
	    return true;
	  if (std::abs(__f1) < std::abs(__f2))
	    __f1 = __func(__x1 += __factor * (__x1 - __x2));
	  else
	    __f2 = __func(__x2 += __factor * (__x2 - __x1));
	}
      return false;
    }


  /**
   * Given a function func defined on an interval x1 to x2, the routine
   * subdivides the interval into n equally spaced segments, and searches
   * for zero crossings of the function.
   * @param __func  A function
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __n  The number of subdivisions of the interval
   * @param __xb  The output vector of root bounds
   */
  template<typename _Tp>
    void
    root_brackets(_Tp __func(_Tp),
		  _Tp __x1, _Tp __x2, int __n,
		  std::vector<std::pair<_Tp, _Tp>>& __xb)
    {
      __xb.clear();
      auto __dx = (__x2 - __x1) / __n;
      auto __x = __x1;
      auto __fp = __func(__x);
      for (int __i = 0; __i < __n; ++__i)
	{
	  __x += __dx;
	  auto __fc = __func(__x);
	  if (__fc * __fp <= _Tp{0})
	    __xb.emplace_back(__x - __dx, __x);
	  __fp = __fc;
	}

      return;
    }


  /**
   * Using bisection, find the root of a function func known to lie between x1
   * and x2.
   * The root (which is returned) is refined until its accuracy is +/-eps.
   * @param __func  A function
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp>
    _Tp
    __root_bisect(_Tp __func(_Tp), _Tp __x1, _Tp __x2, _Tp __eps,
		  int __max_iter = 55)
    {
      auto __f = __func(__x1);
      auto __fmid = __func(__x2);
      if (__f * __fmid >= _Tp{0})
	std::__throw_logic_error("__root_bisect:"
				 " Root must be bracketed for bisection");
      //  Orient search so that f > _Tp{0} lies at x + dx.
      _Tp __dx;
      auto __x = __f < _Tp{0}
	       ? (__dx = __x2 - __x1, __x1)
	       : (__dx = __x1 - __x2, __x2);
      for (int __i = 0; __i < __max_iter; ++__i)
	{
	  auto __xmid = __x + (__dx *= 0.5);
	  __fmid = __func(__xmid);
	  if (__fmid < _Tp{0})
	    __x = __xmid;
	  if (std::abs(__dx) < __eps || __fmid == _Tp{0})
	    return __x;
	}
      std::__throw_logic_error("__root_bisect: Too many bisections");

      return _Tp{0};
    }


  /**
   * Using the secant method, find the root of a function func thought to
   * lie between x1 and x2.  The root, returned as secant, is refined
   * until its accuracy is +/- eps.
   * @param __func  A function
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp>
    _Tp
    __root_secant(_Tp __func(_Tp), _Tp __x1, _Tp __x2, _Tp __eps,
		  int __max_iter = 40)
    {
      _Tp __xl, __x;

      auto __fl = __func(__x1);
      auto __f = __func(__x2);
      if (std::abs(__fl) < std::abs(__f))
	{
	  __x = __x1;
	  __xl = __x2;
	  std::swap(__fl, __f);
	}
      else
	{
	  __xl = __x;
	  __x = __x2;
	}

      for (int __i = 0; __i < __max_iter; ++__i)
	{
	  auto __dx = (__xl - __x) * __f / (__f - __fl);
	  __xl = __x;
	  __fl = __f;
	  __x += __dx;
	  __f = __func(__x);
	  if (std::abs(__dx) < __eps || __f == _Tp{0})
	    return __x;
	}

      std::__throw_logic_error("__root_secant: "
			       "Maximum number of iterations exceeded");

      return  _Tp{0};
    }


  /**
   * Using the false position method, find the root of a function func known
   * to lie between x1 and x2.  The root, returned as root_false_pos,
   * is refined until its accuracy is +/- eps.
   * @param __func  A function
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp>
    _Tp
    __root_false_pos(_Tp __func(_Tp), _Tp __x1, _Tp __x2, _Tp __eps,
		     int __max_iter = 40)
    {
      auto __fl = __func(__x1);
      auto __fh = __func(__x2);
      if (__fl * __fh > _Tp{0})
	std::__throw_logic_error("__root_false_pos: Root must be bracketed");

      _Tp __xl, __xh;
      if (__fl < _Tp{0})
	{
	  __xl = __x1;
	  __xh = __x2;
	}
      else
	{ 
	  __xl = __x2;
	  __xh = __x1;
	  std::swap(__fl, __fh);
	}

      auto __dx = __xh - __xl;
      for (int __i = 0; __i < __max_iter; ++__i)
	{
	  auto __x = __xl + __dx * __fl / (__fl - __fh);
	  auto __f = func(__x);
	  _Tp __del;
	  if (__f < _Tp{0})
	    {
	      __del = __xl - __x;
	      __xl = __x;
	      __fl = __f;
	    }
	  else
	    {
	      __del = __xh - __x;
	      __xh = __x;
	      __fh = __f;
	    }
	  __dx = __xh - __xl;
	  if (std::abs(__del) < __eps || __f == _Tp{0})
	    return __x;
	}

      std::__throw_logic_error("__root_false_pos: "
			       "Maximum number of iterations exceeded");

      return _Tp{0};
    }


  /**
   * Using Ridder's method, find the root of a function func known
   * to lie between x1 and x2.  The root, returned as root_brent, is refined
   * until its accuracy is +/- eps.
   * @param __func  A function
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename _Tp>
    _Tp
    __root_ridder(_Tp __func(_Tp), _Tp __x1, _Tp __x2, _Tp __eps,
		  int __max_iter = 100)
    {
      auto __fh = __func(__x2);
      auto __fl = __func(__x1);
      auto __xh = __x2;
      auto __xl = __x1;

      const _Tp __UNUSED = -1.0e30; // an exceedingly unlikely answer.

      if (__fl * __fh < _Tp{0})
	{
	  auto __ans = __UNUSED;
	  for (int __i = 0; __i < __max_iter; ++__i)
	    {
	      auto __xm = 0.5 * (__xl + __xh);
	      auto __fm = __func(__xm);
	      auto __s = std::sqrt(__fm * __fm - __fl * __fh);
	      if (__s == _Tp{0})
		return __ans;
	      auto __xnew = __xm + (__xm - __xl)
			  * (__fl >= __fh ? 1.0 : -1.0) * __fm / __s;
	      if (std::abs(__xnew - __ans) < __eps)
		return __ans;
	      auto __fnew = __func(__ans = __xnew);
	      if (__fnew == _Tp{0})
		return __ans;
	      if (dsign(__fm, __fnew) != __fm)
		{
		  __xl = __xm;
		  __fl = __fm;
		  __xh = __xnew;
		  __fh = __fnew;
		}
	      else if (dsign(__fl, __fnew) != __fl)
		{
		  __xh = __xnew;
		  __fh = __fnew;
		}
	      else if (dsign(__fh, __fnew) != __fh)
		{
		  __xl = __xnew;
		  __fl = __fnew;
		}
	      else
		std::__throw_logic_error("__root_ridder: "
					 "Some major malfunction");

	      if (std::abs(__xh - __xl) < __eps)
		return __ans;
	    }

	  std::__throw_logic_error("__root_ridder: "
				   "Maximum number of iterations exceeded");
	}
      else
	{
	  if (__fl == _Tp{0})
	    return __x1;
	  if (__fh == _Tp{0})
	    return __x2;

	  std::__throw_logic_error("__root_ridder: Root must be bracketed");
	}

      return _Tp{0};
    }


  /**
   * Using Brent's method, find the root of a function func known
   * to lie between x1 and x2. The root, returned as brent, will be refined
   * until it's accuracy is eps.
   * @param __func  A function
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __eps  The tolerance
   */
  template<typename _Tp>
    _Tp
    __root_brent(_Tp __func(_Tp), _Tp __x1, _Tp __x2, _Tp __eps)
    {
      auto __a = __x1;
      auto __b = __x2;
      auto __c = __x2;
      auto __fa = __func(__a);
      auto __fb = __func(__b);

      const int __ITMAX = 100;
      const _Tp __EPS = 1.0e-12;

      if (__fb * __fa > _Tp{0})
	std::__throw_logic_error("__root_brent: Root must be bracketed");
      auto __fc = __fb;
      for (int __iter = 0; __iter < __ITMAX; ++__iter)
	{
	  _Tp  __d, __e;
	  if (__fb * __fc > _Tp{0})
	    {
	      __c = __a;
	      __fc = __fa;
	      __e = __d = __b - __a;
	    }
	  if (std::abs(__fc) < std::abs(__fb))
	    {
	      __a = __b;
	      __b = __c;
	      __c = __a;
	      __fa = __fb;
	      __fb = __fc;
	      __fc = __fa;
	    }
	  auto __tol1 = 2 * __EPS * std::abs(__b) + 0.5 * __eps;
	  auto __xm = 0.5 * (__c - __b);
	  if (std::abs(__xm) <= __tol1 || __fb == _Tp{0})
	    return __b;
	  if (std::abs(__e) >= __tol1 && std::abs(__fa) > std::abs(__fb))
	    {
	      _Tp __p, __q, __r;
	      auto __s = __fb / __fa;
	      if (__a == __c)
		{
		  __p = 2 * __xm * __s;
		  __q = 1 - __s;
		}
	      else
		{
		  __q = __fa / __fc;
		  __r = __fb / __fc;
		  __p = __s * (2 * __xm * __q * (__q - __r)
		      - (__b - __a) * (__r - 1));
		  __q = (__q - 1) * (__r - 1) * (__s - 1);
		}
	      if (__p > _Tp{0})
		__q = -__q;
	      __p = std::abs(__p);
	      auto __min1 = 3 * __xm * __q - std::abs(__tol1 * __q);
	      auto __min2 = std::abs(__e * __q);
	      if (2 * __p < std::min(__min1, __min2))
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
	  __a = __b;
	  __fa = __fb;
	  if (std::abs(__d) > __tol1)
	    __b += __d;
	  else
	    __b += dsign(__tol1, __xm);
	  __fb = func(__b);
	}
      std::__throw_logic_error("__root_brent: "
			       "Maximum number of iterations exceeded");

      return  _Tp{0};
    }


  /**
   * Using the Newton-Raphson method, find the root of a function known to lie
   * in the interval x1 to x2.  The root will be refined until its accuracy
   * is known within +/- eps.
   * @param __func A routine that provides both the function
   *               and the first derivative of the function at the point x.
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __eps  The tolerance
   */
  template<typename _Tp>
    _Tp
    __root_newton(void func(_Tp, _Tp*, _Tp*), _Tp __x1, _Tp __x2, _Tp __eps)
    {
      const int _S_max_iter = 40;

      auto __x = 0.5 * (__x1 + __x2);
      for (int __i = 0; __i < _S_max_iter; ++__i)
	{
	  _Tp __df, __f;
	  func(__x, &__f, &__df);
	  auto __dx = __f / __df;
	  __x -= __dx;
	  if ((__x1 - __x) * (__x - __x2) < _Tp{0})
	    std::__throw_logic_error("__root_newton: Jumped out of brackets");
	  if (std::abs(__dx) < __eps)
	    return __x;
	}
      std::__throw_logic_error("__root_newton: Maximum number of iterations");

      return  _Tp{0};
    }


  /**
   * Using a combination of Newton-Raphson and bisection,
   * find the root of a function known to lie in the interval
   * x1 to x2.  The root will be refined until its accuracy
   * is known within +/- eps.
   * @param __func  A routine that provides both the function
   *                and the first derivative of the function at the point x.
   * @param __x1  The lower end of the interval
   * @param __x2  The upper end of the interval
   * @param __eps  The tolerance
   */
  template<typename _Tp>
    _Tp
    __root_safe(void __func(_Tp, _Tp*, _Tp*), _Tp __x1, _Tp __x2, _Tp __eps)
    {
      const int _S_max_iter = 100;

      _Tp __df, __fh, __fl;
      __func(__x1, &__fl, &__df);
      __func(__x2, &__fh, &__df);

      if (__fl * __fh < _Tp{0})
	std::__throw_logic_error("__root_safe: Root must be bracketed");

      if (__fl == _Tp{0})
	return __x1;
      if (__fh == _Tp{0})
	return __x2;

      _Tp __xh, __xl;
      if (__fl < _Tp{0})
	{
	  __xl = __x1;
	  __xh = __x2;
	}
      else
	{
	  __xh = __x1;
	  __xl = __x2;
	}

      auto __x = 0.5 * (__x1 + __x2);
      auto __dxold = std::abs(__x2 - __x1);
      auto __dx = __dxold;
      _Tp __f;
      __func(__x, &__f, &__df);
      for (int __i = 0; __i < _S_max_iter; ++__i)
	{
	  if (((__x - __xh) * __df - __f)
	    * ((__x - __xl) * __df - __f) >= _Tp{0}
	   || std::abs(2 * __f) > std::abs(__dxold * __df))
	    {
	      __dxold = __dx;
	      __dx = 0.5 * (__xh - __xl);
	      __x = __xl + __dx;
	      if (__x == __xl)
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
	    __xl = __x;
	  else
	    __xh = __x;
	}

      std::__throw_logic_error("__root_safe: Maximum number of iterations");

      return  _Tp{0};
    }


#endif  //  ROOTS_TCC
