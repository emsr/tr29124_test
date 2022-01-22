// Math extensions -*- C++ -*-

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

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

/** @file ext/root_finding.tcc
 */

#ifndef ROOT_SEARCH_TCC
#define ROOT_SEARCH_TCC 1


#include <stdexcept>
#include <cmath>
#include <tuple> // For tie.

#include <emsr/math_constants.h>

namespace emsr
{

  /**
   * Given a function @c func and an initial range @c x_lower
   * to @c x_upper, the routine expands the range geometrically until a root
   * is bracketed by the returned values @c x_lower and @c x_upper
   * (in which case root_bracket returns true)
   * or until the range becomes unacceptably large
   * (in which case root_bracket returns false).
   * Success is guaranteed for a function which has opposite signs
   * for sufficiently large and small arguments.
   *
   * @param func  A simple one-argument function returning a value.
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename Func>
    bool
    root_bracket(Func func, Tp& x_lower, Tp& x_upper,
		 std::size_t max_iter)
    {
      const Tp golden = emsr::phi_v<Tp>;

      if (x_lower >= x_upper)
	throw std::domain_error("__root_bracket: Initial range must be ordered");
      auto f_lower = func(x_lower);
      auto f_upper = func(x_upper);
      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  if (f_lower * f_upper < Tp{0})
	    return true;
	  if (std::abs(f_lower) < std::abs(f_upper))
	    f_lower = func(x_lower += golden * (x_lower - x_upper));
	  else
	    f_upper = func(x_upper += golden * (x_upper - x_lower));
	}
      return false;
    }


  /**
   * Given a function @c func defined on an interval @c x_lower
   * to @c x_upper, the routine subdivides the interval
   * into @c n equally-spaced segments and searches for zero crossings
   * of the function.
   *
   * @param func  A simple one-argument function returning a value.
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param n  The number of subdivisions of the interval
   * @param xb  The output vector of root bounds
   */
  template<typename Tp, typename Func>
    std::vector<std::pair<Tp, Tp>>
    root_brackets(Func func,
		    Tp x_lower, Tp x_upper, std::size_t n)
    {
      std::vector<std::pair<Tp, Tp>> xb;
      auto dx = (x_upper - x_lower) / Tp(n);
      auto x = x_lower;
      auto f_lo = func(x);
      for (std::size_t i = 0; i < n; ++i)
	{
	  x += dx;
	  auto f_hi = func(x);
	  if (f_lo * f_hi <= Tp{0})
	    xb.emplace_back(x - dx, x);
	  f_lo = f_hi;
	}

      return xb;
    }


  /**
   * Using bisection, find the root of a function @c func
   * known to lie between @c x_lower and @c x_upper.
   * The returned root is refined until its accuracy is <tt>+/- eps</tt>.
   *
   * @param func  A simple one-argument function returning a value.
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param eps  The tolerance
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename Func>
    Tp
    root_bisect(Func func, Tp x_lower, Tp x_upper,
		  Tp eps, std::size_t max_iter)
    {
      auto f = func(x_lower);
      auto f_mid = func(x_upper);
      if (f * f_mid > Tp{0})
	throw std::domain_error("root_bisect: Root must be bracketed for bisection");
      //  Orient search so that f > 0 lies at x + dx.
      auto dx = f < Tp{0} ? x_upper - x_lower : x_lower - x_upper;
      auto x = f < Tp{0} ? x_lower : x_upper;
      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  dx /= Tp{2};
	  auto x_mid = x + dx;
	  f_mid = func(x_mid);
	  if (f_mid < Tp{0})
	    x = x_mid;
	  if (std::abs(dx) < eps || f_mid == Tp{0})
	    return x_mid;
	}

      throw std::runtime_error("root_bisect: Maximum number of bisections exceeded");
    }


  /**
   * Using the secant method, find the root of a function @c func
   * known to lie between @c x_lower and @c x_upper.
   * The returned root is refined until its accuracy is <tt>+/- eps</tt>.
   *
   * @param func  A simple one-argument function returning a value.
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param eps  The tolerance
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename Func>
    Tp
    root_secant(Func func, Tp x_lower, Tp x_upper,
		  Tp eps, std::size_t max_iter)
    {
      Tp x_lo, x;

      auto f_lo = func(x_lower);
      auto f = func(x_upper);
      if (std::abs(f_lo) < std::abs(f))
	{
	  x = x_lower;
	  x_lo = x_upper;
	  std::swap(f_lo, f);
	}
      else
	{
	  x_lo = x_lower;
	  x = x_upper;
	}

      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  auto dx = (x_lo - x) * f / (f - f_lo);
	  x_lo = x;
	  f_lo = f;
	  x += dx;
	  f = func(x);
	  if (std::abs(dx) < eps || f == Tp{0})
	    return x;
	}

      throw std::runtime_error("root_secant: Maximum number of iterations exceeded");
    }


  /**
   * Using the false position method, find the root of a function @c func
   * known to lie between @c x_lower and @c x_upper.
   * The returned root is refined until its accuracy is <tt>+/- eps</tt>.
   *
   * @param func  A simple one-argument function returning a value.
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param eps  The tolerance
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename Func>
    Tp
    root_false_position(Func func, Tp x_lower, Tp x_upper,
			  Tp eps, std::size_t max_iter)
    {
      auto f_lo = func(x_lower);
      auto f_hi = func(x_upper);
      if (f_lo * f_hi > Tp{0})
	throw std::domain_error("root_false_position: Root must be bracketed");

      Tp x_lo, x_hi;
      if (f_lo < Tp{0})
	{
	  x_lo = x_lower;
	  x_hi = x_upper;
	}
      else
	{ 
	  x_lo = x_upper;
	  x_hi = x_lower;
	  std::swap(f_lo, f_hi);
	}

      auto dx = x_hi - x_lo;
      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  auto x = x_lo + dx * f_lo / (f_lo - f_hi);
	  auto f = func(x);
	  Tp del;
	  if (f < Tp{0})
	    {
	      del = x_lo - x;
	      x_lo = x;
	      f_lo = f;
	    }
	  else
	    {
	      del = x_hi - x;
	      x_hi = x;
	      f_hi = f;
	    }
	  dx = x_hi - x_lo;
	  if (std::abs(del) < eps || f == Tp{0})
	    return x;
	}

      throw std::runtime_error("root_false_position: "
				     "Maximum number of iterations exceeded");
    }


  /**
   * Using Ridder's method, find the root of a function @c func known
   * to lie between @c x_lower and @c x_upper.
   * The returned root is refined until its accuracy is <tt>+/- eps</tt>.
   *
   * @param func  A simple one-argument function returning a value.
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param eps  The tolerance
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename Func>
    Tp
    root_ridder(Func func, Tp x_lower, Tp x_upper,
		  Tp eps, std::size_t max_iter)
    {
      auto f_hi = func(x_upper);
      auto f_lo = func(x_lower);
      auto x_hi = x_upper;
      auto x_lo = x_lower;

      if (f_lo * f_hi > Tp{0})
	throw std::domain_error("root_ridder: Root must be bracketed");

      if (f_lo == Tp{0})
	return x_lower;

      if (f_hi == Tp{0})
	return x_upper;

      const Tp UNUSED = -1.0e30; // an exceedingly unlikely answer.
      auto ans = UNUSED;
      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  auto xm = (x_lo + x_hi) / Tp{2};
	  auto fm = func(xm);
	  auto s = std::sqrt(fm * fm - f_lo * f_hi);
	  if (s == Tp{0})
	    return ans;
	  auto xnew = xm + (xm - x_lo)
		      * (f_lo >= f_hi ? 1.0 : -1.0) * fm / s;
	  if (std::abs(xnew - ans) < eps)
	    return ans;
	  ans = xnew;
	  auto fnew = func(ans);
	  if (fnew == Tp{0})
	    return ans;
	  if (std::copysign(fm, fnew) != fm)
	    {
	      x_lo = xm;
	      f_lo = fm;
	      x_hi = xnew;
	      f_hi = fnew;
	    }
	  else if (std::copysign(f_lo, fnew) != f_lo)
	    {
	      x_hi = ans;
	      f_hi = fnew;
	    }
	  else if (std::copysign(f_hi, fnew) != f_hi)
	    {
	      x_lo = ans;
	      f_lo = fnew;
	    }
	  else
	    throw std::runtime_error("root_ridder: Some major malfunction");

	  if (std::abs(x_hi - x_lo) < eps)
	    return ans;
	}

      throw std::runtime_error("root_ridder: Maximum number of iterations exceeded");
    }


  /**
   * Using Brent's method, find the root of a function @c func known
   * to lie between @c x_lower and @c x_upper.
   * The returned root is refined until it's accuracy is <tt>+/- eps</tt>.
   *
   * @param func  A simple one-argument function returning a value.
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param eps  The tolerance
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename Func>
    Tp
    root_brent(Func func, Tp x_lower, Tp x_upper,
	       Tp eps, std::size_t max_iter)
    {
      auto x_a = x_lower;
      auto x_b = x_upper;
      auto x_c = x_upper;
      auto f_a = func(x_a);
      auto f_b = func(x_b);

      const Tp s_eps = std::numeric_limits<Tp>::epsilon();

      if (f_a * f_b > Tp{0})
	throw std::domain_error("root_brent: Root must be bracketed");

      auto f_c = f_b;
      for (std::size_t iter = 0; iter < max_iter; ++iter)
	{
	  Tp d, e;
	  if (f_b * f_c > Tp{0})
	    {
	      x_c = x_a;
	      f_c = f_a;
	      e = d = x_b - x_a;
	    }
	  if (std::abs(f_c) < std::abs(f_b))
	    {
	      x_a = x_b;
	      x_b = x_c;
	      x_c = x_a;
	      f_a = f_b;
	      f_b = f_c;
	      f_c = f_a;
	    }
	  auto toler = Tp{2} * s_eps * std::abs(x_b) + eps / Tp{2};
	  auto xm = (x_c - x_b) / Tp{2};
	  if (std::abs(xm) <= toler || f_b == Tp{0})
	    return x_b;
	  if (std::abs(e) >= toler && std::abs(f_a) > std::abs(f_b))
	    {
	      Tp p, q, r;
	      auto s = f_b / f_a;
	      if (x_a == x_c)
		{
		  p = Tp{2} * xm * s;
		  q = Tp{1} - s;
		}
	      else
		{
		  q = f_a / f_c;
		  r = f_b / f_c;
		  p = s * (Tp{2} * xm * q * (q - r)
		      - (x_b - x_a) * (r - Tp{1}));
		  q = (q - Tp{1}) * (r - Tp{1}) * (s - Tp{1});
		}
	      if (p > Tp{0})
		q = -q;
	      p = std::abs(p);
	      auto min1 = Tp{3} * xm * q - std::abs(toler * q);
	      auto min2 = std::abs(e * q);
	      if (Tp{2} * p < std::min(min1, min2))
		{
		  e = d;
		  d = p / q;
		}
	      else
		{
		  d = xm;
		  e = d;
		}
	    }
	  else
	    {
	      d = xm;
	      e = d;
	    }
	  x_a = x_b;
	  f_a = f_b;
	  if (std::abs(d) > toler)
	    x_b += d;
	  else
	    x_b += std::copysign(toler, xm);
	  f_b = func(x_b);
	}

      throw std::runtime_error("root_brent: Maximum number of iterations exceeded");
    }


  /**
   * Using the Newton-Raphson method, find the root of a function @c func
   * known to lie in the interval @c x_lower to @c x_upper.
   * The returned root is refined until its accuracy is
   * within <tt>+/- eps</tt>.
   *
   * @f[
   *    x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
   * @f]
   *
   * @param func A routine that provides both the function value
   *               and the first derivative of the function at the point x.
   *               The return type must be decomposable into [value, deriv].
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param eps  The tolerance
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename StateFunc>
    Tp
    root_newton(StateFunc func, Tp x_lower, Tp x_upper,
		Tp eps, std::size_t max_iter)
    {
      auto x = (x_lower + x_upper) / Tp{2};
      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  const auto [value, deriv] = func(x);
	  const auto dx = value / deriv;
	  x -= dx;
	  if ((x_lower - x) * (x - x_upper) < Tp{0})
	    throw std::runtime_error("root_newton: Jumped out of brackets");
	  if (std::abs(dx) < eps)
	    return x;
	}
      throw std::runtime_error("root_newton: Maximum number of iterations exceeded");
    }


  /**
   * Using the Halley method, find the root of a function @c func
   * known to lie in the interval @c x_lower to @c x_upper.
   * The returned root is refined until its accuracy is
   * within <tt>+/- eps</tt>.
   *
   * @f[
   *    x_{n+1} - x_n = -\frac{2 f(x_n) f'(x_n)}
   *                    {2 [f'(x_n)]^2 - f(x_n) f''(x_n)}
   * @f]
   * This form indicates the close relationship to the Newton method:
   * @f[
   *    x_{n+1} - x_n = -\frac{f'(x_n)}
   *                    {f'(x_n) - [f(x_n) f''(x_n)]/[2f'(x_n)]}
   * @f]
   *
   * @param func A routine that provides both the function value, and the
   *               first and second derivative of the function at the point x.
   *               The return type must be decomposable into [val, der, der2].
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param eps  The tolerance
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename StateFunc>
    Tp
    root_halley(StateFunc func, Tp x_lower, Tp x_upper,
		Tp eps, std::size_t max_iter)
    {
      auto x = (x_lower + x_upper) / Tp{2};
      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  const auto [f, df, d2f] = func(x);
	  const auto dx = Tp{2} * f * df
			   / (Tp{2} * df * df - f * d2f);
	  x -= dx;
	  if ((x_lower - x) * (x - x_upper) < Tp{0})
	    throw std::runtime_error("root_halley: Jumped out of brackets");
	  if (std::abs(dx) < eps)
	    return x;
	}
      throw std::runtime_error("root_halley: Maximum number of iterations exceeded");
    }


  /**
   * Using the Steffensen method, find the root of a function @c __func
   * known to lie in the interval @c __x_lower to @c __x_upper.
   * The returned root is refined until its accuracy is
   * within <tt>+/- __eps</tt>.
   *
   * @f[
   *    x_{n+1} = x_n - \frac{f(x_n)}{\frac{f(x_n + f(x_n))}{f(x_n)} - 1}
   * @f]
   *
   * The Steffensen method amounts to estimating the derivative using the
   * previous function value as a stepsize and applying the Newton iteration.
   *
   * This method is supposed to have quadratic convergence like Newton's
   * method but it only required te function values.  On the other hand,
   * the initial estimate generally must be tighter than for Newton's method.
   *
   * @param __func A routine that provides the function value at the point x.
   * @param __x_lower  The lower end of the search interval.
   * @param __x_upper  The upper end of the search interval.
   * @param __eps  The tolerance
   * @param __max_iter  The maximum number of iterations
   */
  template<typename Tp, typename Func>
    Tp
    root_steffensen(Func func, Tp x_lower, Tp x_upper,
    		    Tp eps, std::size_t max_iter)
    {
      auto x = (x_lower + x_upper) / Tp{2};
      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  auto fval = func(x);
	  auto gval = func(x + fval) / fval - Tp{1};
	  auto dx = fval / gval;
	  x -= dx;
	  if ((x_lower - x) * (x - x_upper) < Tp{0})
	    throw std::runtime_error("root_steffensen: Jumped out of brackets");
	  if (std::abs(dx) < eps)
	    return x;
	}
      throw std::runtime_error("root_steffensen: Maximum number of iterations exceeded");
    }


  /**
   * Using a combination of Newton-Raphson and bisection, find the root
   * of a function @c func known to lie in the interval @c x_lower
   * to @c x_upper.
   * The returned root is refined until its accuracy is
   * within <tt>+/- eps</tt>.
   *
   * @param func  A routine that provides both the function
   *                and the first derivative of the function at the point x.
   * @param x_lower  The lower end of the search interval.
   * @param x_upper  The upper end of the search interval.
   * @param eps  The tolerance
   * @param max_iter  The maximum number of iterations
   */
  template<typename Tp, typename StateFunc>
    Tp
    root_safe(StateFunc func, Tp x_lower, Tp x_upper,
		Tp eps, std::size_t max_iter)
    {
      auto [value_lower, deriv_lower] = func(x_lower);
      auto [value_upper, deriv_upper] = func(x_upper);

      if (value_lower * value_upper > Tp{0})
	throw std::domain_error("root_safe: Root must be bracketed");

      if (value_lower == Tp{0})
	return x_lower;
      if (value_upper == Tp{0})
	return x_upper;

      Tp x_hi, x_lo;
      if (value_lower < Tp{0})
	{
	  x_lo = x_lower;
	  x_hi = x_upper;
	}
      else
	{
	  x_hi = x_lower;
	  x_lo = x_upper;
	}

      auto x = (x_lower + x_upper) / Tp{2};
      auto dxold = std::abs(x_upper - x_lower);
      auto dx = dxold;
      auto [value, deriv] = func(x);
      for (std::size_t i = 0; i < max_iter; ++i)
	{
	  if (((x - x_hi) * deriv - value)
	    * ((x - x_lo) * deriv - value) > Tp{0}
	   || std::abs(Tp{2} * value)
	    > std::abs(dxold * deriv))
	    {
	      dxold = dx;
	      dx = (x_hi - x_lo) / Tp{2};
	      x = x_lo + dx;
	      if (x == x_lo)
		return x;
	    }
	  else
	    {
	      dxold = dx;
	      dx = value / deriv;
	      auto temp = x;
	      x -= dx;
	      if (temp == x)
		return x;
	    }
	  if (std::abs(dx) < eps)
	    return x;

	  auto [val, der] = func(x);
	  value = val;
	  deriv = der;

	  if (value < Tp{0})
	    x_lo = x;
	  else
	    x_hi = x;
	}

      throw std::runtime_error("root_safe: Maximum number of iterations exceeded");
    }

} // namespace emsr

#endif // ROOT_SEARCH_TCC
