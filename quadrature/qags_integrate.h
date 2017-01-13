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
// Ported from GSL by Jason Dick
// Originally written by Brian Gaugh
//
// This file implements Gauss-Kronrod integration with singularities
// Based upon gsl-2.3/integration/qags.c

#ifndef QAGS_INTEGRATE_H
#define QAGS_INTEGRATE_H 1

#include <utility>
#include <limits>
#include <tuple>
#include <stdexcept>

#include "integration_error.h"
#include "integration_transform.h"
#include "qk_integrate.h"
#include "integration_workspace.h"
#include "extrapolation_table.h"

namespace __gnu_test
{

  /**
   * Integrate potentially singular function from a to b using recursive
   * Gauss-Kronrod algorithm.
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qags_integrate(integration_workspace<_Tp>& __workspace,
		   const _FuncTp& func,
		   const _Tp a, const _Tp b,
		   const _Tp epsabs,
		   const _Tp epsrel);

  /**
   * Integrate potentially singular function from a to b using recursive
   * Gauss-Kronrod algorithm.
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qags_integrate(const _FuncTp& __func,
		   const _Tp __a, const _Tp __b,
		   const _Tp __epsabs,
		   const _Tp __epsrel,
		   const std::size_t __limit)
    {
      integration_workspace<_Tp> __workspace(__limit);
      return qags_integrate(__workspace, __func, __a, __b, __epsabs, __epsrel);
    }

  /**
   * Integrate a function defined over (-\infty, +\infty).
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qagi_integrate(integration_workspace<_Tp>& __workspace,
		   const _FuncTp& __func,
		   _Tp __epsabs, _Tp __epsrel)
    {
      return qags_integrate(__workspace, i_transform<_FuncTp, _Tp>(__func),
			    _Tp{0}, _Tp{1}, __epsabs, __epsrel);
    }

  /**
   * Integrate a function defined over (-\infty, +\infty).
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qagi_integrate(const _FuncTp& __func,
		   const _Tp __epsabs,
		   const _Tp __epsrel,
		   const std::size_t __limit)
    {
      integration_workspace<_Tp> __workspace(__limit);
      return qagi_integrate(__workspace, __func, __epsabs, __epsrel);
    }

  /**
   * Integrate a function defined over (-\infty, b].
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qagil_integrate(integration_workspace<_Tp>& __workspace,
		    const _FuncTp& __func, _Tp __b,
		    _Tp __epsabs, _Tp __epsrel)
    {
      return qags_integrate(__workspace, il_transform(__func, __b),
			    _Tp{0}, _Tp{1}, __epsabs, __epsrel);
    }

  /**
   * Integrate a function defined over (-\infty, b].
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qagil_integrate(const _FuncTp& __func, const _Tp __b,
		    const _Tp __epsabs,
		    const _Tp __epsrel,
		    const std::size_t __limit)
    {
      integration_workspace<_Tp> __workspace(__limit);
      return qagil_integrate(__workspace, __func, __b, __epsabs, __epsrel);
    }

  /**
   * Integrate a function defined over [a, +\infty).
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qagiu_integrate(integration_workspace<_Tp>& __workspace,
		    const _FuncTp& __func, _Tp __a,
		    _Tp __epsabs, _Tp __epsrel)
    {
      return qags_integrate(__workspace, iu_transform(__func, __a),
			    _Tp{0}, _Tp{1}, __epsabs, __epsrel);
    }

  /**
   * Integrate a function defined over [a, +\infty).
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qagiu_integrate(const _FuncTp& func, const _Tp a,
		    const _Tp epsabs,
		    const _Tp epsrel,
		    const std::size_t limit)
    {
      integration_workspace<_Tp> __workspace(limit);
      return qagiu_integrate(__workspace, func, a, epsabs, epsrel);
    }

  /**
   * Integrate potentially singular function from a to b using recursive
   * Gauss-Kronrod algorithm.
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qags_integrate(integration_workspace<_Tp>& __workspace,
		   const _FuncTp& __func,
		   const _Tp __a, const _Tp __b,
		   _Tp __epsabs,
		   _Tp __epsrel)
    {
      const qk_intrule __qk_rule = QK_21;

      _Tp __tolerance;

      _Tp __ertest = 0;
      _Tp __error_over_large_intervals = 0;
      _Tp __reseps = 0, __abseps = 0, __correc = 0;
      std::size_t __ktmin = 0;
      int __roundoff_type1 = 0, __roundoff_type2 = 0, __roundoff_type3 = 0;
      int __error_type = NO_ERROR, __error_type2 = NO_ERROR;

      const auto _S_max = std::numeric_limits<_Tp>::max();
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      std::size_t __iteration = 0;

      int __extrapolate = 0;
      int __disallow_extrapolation = 0;

      __workspace.set_initial_limits(__a, __b);

      const auto __limit = __workspace.capacity();

      if (__epsabs <= 0
	  && (__epsrel < 50 * _S_eps || __epsrel < 0.5e-28))
	std::__throw_logic_error("qags_integrate: "
				 "Tolerance cannot be achieved "
				 "with given absolute "
				 "and relative error limits.");

      // Perform the first integration.

      _Tp __result0, __abserr0, __resabs0, __resasc0;
      std::tie(__result0, __abserr0, __resabs0, __resasc0)
	  = qk_integrate(__func, __a, __b, __qk_rule);

      __workspace.set_initial_results(__result0, __abserr0);

      __tolerance = std::max(__epsabs, __epsrel * std::abs(__result0));

      if (__abserr0 <= 100 * _S_eps * __resabs0
	  && __abserr0 > __tolerance)
	__throw__IntegrationError("qags_integrate: "
				  "cannot reach tolerance because "
				  "of roundoff error on first attempt",
				  ROUNDOFF_ERROR, __result0, __abserr0);
      else if ((__abserr0 <= __tolerance && __abserr0 != __resasc0)
		|| __abserr0 == _Tp{0})
	return std::make_tuple(__result0, __abserr0);
      else if (__limit == 1)
	__throw__IntegrationError("qagp_integrate: "
				  "a maximum of one iteration was insufficient",
				  MAX_ITER_ERROR, __result0, __abserr0);

      extrapolation_table<_Tp> __table;
      __table.append(__result0);

      auto __area = __result0;
      auto __errsum = __abserr0;

      auto __res_ext = __result0;
      auto __err_ext = _S_max;

      auto __positive_integrand = __test_positivity(__result0, __resabs0);

      __iteration = 1;

      do
	{
	  // Bisect the subinterval with the largest error estimate.

	  _Tp __a_i, __b_i, __r_i, __e_i;
	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  auto __current_level = __workspace.current_level() + 1;

	  auto __a1 = __a_i;
	  auto __b1 = (__a_i + __b_i) / _Tp{2};
	  auto __a2 = __b1;
	  auto __b2 = __b_i;

	  ++__iteration;

	  _Tp __area1, __error1, __resabs1, __resasc1;
	  std::tie(__area1, __error1, __resabs1, __resasc1)
	    = qk_integrate(__func, __a1, __b1, __qk_rule);

	  _Tp __area2, __error2, __resabs2, __resasc2;
	  std::tie(__area2, __error2, __resabs2, __resasc2)
	    = qk_integrate(__func, __a2, __b2, __qk_rule);

	  auto __area12 = __area1 + __area2;
	  auto __error12 = __error1 + __error2;
	  auto __last_e_i = __e_i;

	  // Improve previous approximations to the integral and test for
	  // accuracy.

	  __errsum += __error12 - __e_i;
	  __area += __area12 - __r_i;

	  __tolerance = std::max(__epsabs, __epsrel * std::abs(__area));

	  if (__resasc1 != __error1 && __resasc2 != __error2)
	    {
	      auto __delta = __r_i - __area12;

	      if (std::abs(__delta) <= 1.0e-5 * std::abs(__area12)
		  && __error12 >= 0.99 * __e_i)
		{
		  if (!__extrapolate)
		    ++__roundoff_type1;
		  else
		    ++__roundoff_type2;
		}
	      if (__iteration > 10 && __error12 > __e_i)
		++__roundoff_type3;
	    }

	  // Test for roundoff and eventually set error flag.

	  if (__roundoff_type1 + __roundoff_type2 >= 10
	      || __roundoff_type3 >= 20)
	    __error_type = ROUNDOFF_ERROR;

	  if (__roundoff_type2 >= 5)
	    __error_type2 = MAX_ITER_ERROR;

	  // Set error flag in the case of bad integrand behaviour at
	  // a point of the integration range.

	  if (__workspace.subinterval_too_small(__a1, __a2, __b2))
	    __error_type = EXTRAP_ROUNDOFF_ERROR;

	  // Append the newly-created intervals to the list.
	  __workspace.update(__a1, __b1, __area1, __error1,
			     __a2, __b2, __area2, __error2);

	  if (__errsum <= __tolerance)
	    {
	      __check_error(__func__, __error_type,
			    __workspace.sum_results(), __errsum);
	      return std::make_tuple(__workspace.sum_results(), __errsum);
	    }

	  if (__error_type)
	    break;

	  if (__iteration >= __limit - 1)
	    {
	      __error_type = MAX_ITER_ERROR;
	      break;
	    }

	  if (__iteration == 2) // Set up variables on first iteration.
	    {
	      __error_over_large_intervals = __errsum;
	      __ertest = __tolerance;
	      __table.append(__area);
	      continue;
	    }

	  if (__disallow_extrapolation)
	    continue;

	  __error_over_large_intervals += -__last_e_i;

	  if (__current_level < __workspace.max_level())
	    __error_over_large_intervals += __error12;

	  if (!__extrapolate)
	    {
	      // Test whether the interval to be bisected next is the
	      // smallest interval.

	      if (__workspace.large_interval())
		continue;

	      __extrapolate = 1;
	      __workspace.set_nrmax(1);
	    }

	  if (!__error_type2 && __error_over_large_intervals > __ertest)
	    if (__workspace.increase_nrmax())
	      continue;

	  // Perform extrapolation.

	  __table.append(__area);

	  std::tie(__reseps, __abseps) = __table.qelg();

	  ++__ktmin;

	  if (__ktmin > 5 && __err_ext < 0.001 * __errsum)
	    __error_type = DIVERGENCE_ERROR;

	  if (__abseps < __err_ext)
	    {
	      __ktmin = 0;
	      __err_ext = __abseps;
	      __res_ext = __reseps;
	      __correc = __error_over_large_intervals;
	      __ertest = std::max(__epsabs, __epsrel * std::abs(__reseps));
	      if (__err_ext <= __ertest)
		break;
	    }

	  // Prepare bisection of the smallest interval.
	  if (__table.get_nn() == 1)
	    __disallow_extrapolation = 1;

	  if (__error_type == DIVERGENCE_ERROR)
	    break;

	  // Work on interval with largest error.
	  __workspace.reset_nrmax();
	  __extrapolate = 0;
	  __error_over_large_intervals = __errsum;

	}
      while (__iteration < __limit);

      auto __result = __res_ext;
      auto __abserr = __err_ext;

      if (__err_ext == _S_max)
	{
	  __check_error(__func__, __error_type,
			__workspace.sum_results(), __errsum);
	  return std::make_tuple(__workspace.sum_results(), __errsum);
	}

      if (__error_type || __error_type2)
	{
	  if (__error_type2)
	    __err_ext += __correc;

	  if (__error_type == 0)
	    __error_type = SINGULAR_ERROR;

	  if (__res_ext != _Tp{0} && __area != _Tp{0})
	    {
	      if (__err_ext / std::abs(__res_ext) > __errsum / std::abs(__area))
		{
		  __check_error(__func__, __error_type,
				__workspace.sum_results(), __errsum);
		  return std::make_tuple(__workspace.sum_results(), __errsum);
		}
	    }
	  else if (__err_ext > __errsum)
	    {
	      __check_error(__func__, __error_type,
			    __workspace.sum_results(), __errsum);
	      return std::make_tuple(__workspace.sum_results(), __errsum);
	    }
	  else if (__area == _Tp{0})
	    {
	      __check_error<_Tp>(__func__, __error_type);
	      std::__throw_runtime_error("qags_integrate: Unknown error.");
	    }
	}

      // Test on divergence.
      auto __max_area = std::max(std::abs(__res_ext), std::abs(__area));
      if (!__positive_integrand && __max_area < 0.01 * __resabs0)
	{
	  __check_error<_Tp>(__func__, __error_type);
	  std::__throw_runtime_error("qags_integrate: Unknown error.");
	}

      auto __ratio = __res_ext / __area;
      if (__ratio < 0.01 || __ratio > 100.0 || __errsum > std::abs(__area))
	__error_type = UNKNOWN_ERROR;

      if (__error_type == 0)
	return std::make_tuple(__result, __abserr);

      __check_error<_Tp>(__func__, __error_type);
      __throw__IntegrationError("qags_integrate: Unknown error.",
				UNKNOWN_ERROR, __result, __abserr);
    }

} // namespace __gnu_test

#endif // QAGS_INTEGRATE_H
