// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2017 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.
//
// Written by Jason Dick.
//
// Provides readable access to integration functions

#ifndef INTEGRATE_H
#define INTEGRATE_H 1

#include <limits>

namespace __gnu_test
{

/*
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp, _Tp, _Tp>
    integration_rule(_FuncTp __func, _Tp __a, _Tp __b,
		     const qk_intrule __qkintrule);
*/

  // Integrates a smooth function from a to b
  // absolute_error_lim and relative_error_lim are the absolute
  // and relative error limits.
  // max_iter is the maximum number of iterations allowed
  // qkintrule is the Gauss-Kronrod integration rule, and can be either:
  // QK_15, QK_21, QK_31, QK_41, QK_51, QK_61
  // Higher-order rules converge more rapidly for most functions,
  // but may slow convergence for less well-behaved ones.
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_smooth(const _FuncTp& __func, _Tp __a, _Tp __b,
		     _Tp __max_abs_error,
		     _Tp __max_rel_error,
		     const std::size_t __max_iter = 1024,
		     const qk_intrule __qkintrule = QK_61)
    {
      return qag_integrate(__func, __a, __b, __max_abs_error,
			   __max_rel_error, __max_iter, __qkintrule);
    }

  // Recursive Gauss-Kronrod integration optimized for
  // discontinuous or singular functions
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_singular(const _FuncTp& __func, _Tp __a, _Tp __b,
		       _Tp __max_abs_error,
		       _Tp __max_rel_error,
		       const std::size_t __max_iter = 1024)
    {
      return qags_integrate(__func, __a, __b, __max_abs_error,
			    __max_rel_error, __max_iter);
    }

  // Integrates function from -infinity to +infinity
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_infinite(const _FuncTp& __func,
		       _Tp __max_abs_error,
		       _Tp __max_rel_error,
		       const std::size_t __max_iter = 1024)
    {
      return qagi_integrate(__func, __max_abs_error,
			    __max_rel_error, __max_iter);
    }

  // Integrations function from -infinity to b
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_from_infinity(const _FuncTp& __func, _Tp __b,
			    _Tp __max_abs_error,
			    _Tp __max_rel_error,
			    const std::size_t __max_iter = 1024)
    {
      return qagil_integrate(__func, __b, __max_abs_error,
			     __max_rel_error, __max_iter);
    }

  // Integrations function from a to +infinity
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_to_infinity(const _FuncTp& __func, _Tp __a,
			  _Tp __max_abs_error,
			  _Tp __max_rel_error,
			  const std::size_t __max_iter = 1024)
    {
      return qagiu_integrate(__func, __a, __max_abs_error,
			     __max_rel_error, __max_iter);
    }

  // Integrates function, allows setting of limits as being +/- infinity
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
	      _Tp __max_abs_error,
	      _Tp __max_rel_error,
	      const std::size_t __max_iter = 1024)
    {
      using std::isnan; // To avoid ambiguous overload
      const _Tp __infty = std::numeric_limits<_Tp>::infinity();
      const _Tp __NaN = std::numeric_limits<_Tp>::quiet_NaN();

      if (isnan(__a) || isnan(__b))
	return std::make_tuple(__NaN, __NaN);
      else if (__a == -__infty)
	{
	  if (__b == __infty) // Integration from -inf to +inf
	    return integrate_infinite(__func, __max_abs_error,
				      __max_rel_error, __max_iter);
	  else if (__b == -__infty)
	    std::__throw_runtime_error("Attempted to integrate from -infinity"
				" to -infinity in integrate()");
	  else // Integration from -inf to finite value
	    return integrate_from_infinity(__func, __b, __max_abs_error,
					   __max_rel_error, __max_iter);
	}
      else if (__a == __infty)
	{
	  if (__b == __infty)
	    std::__throw_runtime_error("Attempted to integrate from +infinity"
				" to +infinity in integrate()");
	  else if (__b == -__infty) // Integration from +inf to -inf,
	    {		     // Call integrate_infinite() and flip sign
	      std::tuple<_Tp, _Tp> __res
		= integrate_infinite(__func, __max_abs_error,
				     __max_rel_error, __max_iter);
	      return std::make_tuple(-std::get<0>(__res), std::get<1>(__res));
	    }
	  else // Integration from +inf to finite value,
	    { // Call integrate_to_infinity and flip sign
	      std::tuple<_Tp, _Tp>
		__res = integrate_to_infinity(__func, __b, __max_abs_error,
					      __max_rel_error, __max_iter);
	      return std::make_tuple(-std::get<0>(__res), std::get<1>(__res));
	    }
	}
      else // a is finite
	{
	  if (__b == __infty)
	    return integrate_to_infinity(__func, __a, __max_abs_error,
					 __max_rel_error, __max_iter);
	  else if (__b == -__infty) // Integration from finite value to -inf,
	    { // Call integrate_from_infinity and flip sign
	      std::tuple<_Tp, _Tp>
		__res = integrate_from_infinity(__func, __a,
						__max_abs_error,
						__max_rel_error,
						__max_iter);
	      return std::make_tuple(-std::get<0>(__res), std::get<1>(__res));
	    }
	  else // Both a and b finite, call integrate_singular
	    return integrate_singular(__func, __a, __b, __max_abs_error,
				      __max_rel_error, __max_iter);
	}
    }

} // namespace __gnu_test

#endif // INTEGRATE_H
