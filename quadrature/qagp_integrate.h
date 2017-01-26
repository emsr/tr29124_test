// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
// Copyright (C) 2016-2017 Free Software Foundation, Inc.
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
// Ported from GSL by Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements integration using a recursive Gauss-Kronrod algorithm
// Based on gsl/integration/qagp.c

#ifndef QAGP_INTEGRATE_H
#define QAGP_INTEGRATE_H 1

#include "qk_integrate.h"
#include "integration_workspace.h"
#include "extrapolation_table.h"

namespace __gnu_test
{

template<typename _Tp>
void
dump_ws(integration_workspace<_Tp>& workspace, const char* cmp, const char* msg)
{
  std::cerr.precision(8);
  std::cerr << '\n' << cmp << ": " << msg << '\n';
  std::cerr << std::setw(15) << workspace << '\n';
}

  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qagp_integrate(integration_workspace<_Tp>& __workspace,
		   const _FuncTp& __func,
		   std::vector<_Tp> __pts,
		   _Tp __epsabs, _Tp __epsrel)
    {
      return qagp_integrate(__workspace, __func, __pts,
			    __epsabs, __epsrel, QK_21);
    }

  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qagp_integrate(integration_workspace<_Tp>& __workspace,
		   const _FuncTp& __func,
		   std::vector<_Tp> __pts,
		   _Tp __epsabs, _Tp __epsrel,
		   const qk_intrule __qk_rule)
    {
      const auto _S_max = std::numeric_limits<_Tp>::max();
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto __limit = __workspace.capacity();
      const auto __n_ivals = __pts.size() - 1;

      bool __extrapolate = false;
      bool __disallow_extrapolation = false;

      if (__epsabs <= 0
	  && (__epsrel < 50 * _S_eps || __epsrel < 0.5e-28))
	std::__throw_runtime_error("qagp_integrate: "
				   "Tolerance cannot be achieved "
				   "with given absolute "
				   "and relative error limits.");

      if (__pts.size() > __workspace.capacity())
	std::__throw_runtime_error("qagp_integrate: "
				   "number of pts exceeds size of workspace");

      // Check that the integration range and break points are an
      // ascending sequence.

      for (std::size_t __i = 0; __i < __n_ivals; ++__i)
	if (__pts[__i + 1] < __pts[__i])
	  std::__throw_runtime_error("qagp_integrate: "
				     "points are not in an ascending sequence");

      __workspace.clear();

      // Perform the first integration.
      auto __result0 = _Tp{0};
      auto __abserr0 = _Tp{0};
      auto __resabs0 = _Tp{0};
      for (std::size_t __i = 0; __i < __n_ivals; ++__i)
	{
	  const auto __a = __pts[__i];
	  const auto __b = __pts[__i + 1];

	  _Tp __area1, __error1, __resabs1, __resasc1;
	  std::tie(__area1, __error1, __resabs1, __resasc1)
	    = qk_integrate(__func, __a, __b, __qk_rule);

	  __result0 += __area1;
	  __abserr0 += __error1;
	  __resabs0 += __resabs1;
	  std::size_t __level = (__error1 == __resasc1 && __error1 != _Tp{0})
				? 1 : 0;
	  __workspace.append(__a, __b, __area1, __error1, __level);
	}

      // Compute the initial error estimate.
      auto __errsum = _Tp{0};
      for (std::size_t __i = 0; __i < __n_ivals; ++__i)
	{
	  if (__workspace.depth(__i) == 1)
	    {
	      __workspace.set_abs_error(__i, __abserr0);
	      __workspace.set_depth(__i, 0);
	    }
	  __errsum += __workspace.abs_error(__i);
	}	

      // We must re-sort because the errors were reassigned.
      __workspace.sort_error();

      // Test on accuracy.
      auto __tolerance = std::max(__epsabs, __epsrel * std::abs(__result0));

      if (__abserr0 <= 100 * _S_eps * __resabs0
	   && __abserr0 > __tolerance)
	__throw__IntegrationError("qagp_integrate: "
				  "cannot reach tolerance because "
				  "of roundoff error on first attempt",
				  ROUNDOFF_ERROR, __result0, __abserr0);
      else if (__abserr0 <= __tolerance)
	return std::make_tuple(__result0, __abserr0);
      else if (__limit == 1)
	__throw__IntegrationError("qagp_integrate: "
				  "a maximum of one iteration was insufficient",
				  MAX_ITER_ERROR, __result0, __abserr0);

      extrapolation_table<_Tp> __table;
      __table.append(__result0);

      auto __res_ext = __result0;
      auto __err_ext = _S_max;

      auto __area = __result0;
      auto __iteration = __n_ivals - 1;
      auto __ktmin = 0u;
      auto __ertest = __tolerance;
      auto __error_over_large_intervals = __errsum;
      auto __reseps = _Tp{0}, __abseps = _Tp{0}, __correc = _Tp{0};
      int __error_type = NO_ERROR, __error_type2 = NO_ERROR;
      int __roundoff_type1 = 0, __roundoff_type2 = 0, __roundoff_type3 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.

	  _Tp __a_i, __b_i, __r_i, __e_i;
	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  const auto __current_depth = __workspace.current_depth() + 1;

	  const auto __a1 = __a_i;
	  const auto __b1 = (__a_i + __b_i) / _Tp{2};
	  const auto __a2 = __b1;
	  const auto __b2 = __b_i;

	  ++__iteration;

	  _Tp __area1, __error1, __resabs1, __resasc1;
	  std::tie(__area1, __error1, __resabs1, __resasc1)
	    = qk_integrate(__func, __a1, __b1, __qk_rule);

	  _Tp __area2, __error2, __resabs2, __resasc2;
	  std::tie(__area2, __error2, __resabs2, __resasc2)
	    = qk_integrate(__func, __a2, __b2, __qk_rule);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;
	  const auto __last_e_i = __e_i;

	  // Improve previous approximations to the integral and test for
	  // accuracy.

	  __errsum += __error12 - __e_i;
	  __area += __area12 - __r_i;

	  __tolerance = std::max(__epsabs, __epsrel * std::abs(__area));

	  if (__resasc1 != __error1 && __resasc2 != __error2)
	    {
	      const auto __delta = __r_i - __area12;

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
	  // a point of the integration range

	  if (__workspace.subinterval_too_small(__a1, __a2, __b2))
	    __error_type = EXTRAP_ROUNDOFF_ERROR;

	  // Split the current interval in two.
	  __workspace.split(__b1, __area1, __error1, __area2, __error2);

	  if (__errsum <= __tolerance)
	    {
	      const auto __result = __workspace.total_integral();
	      __check_error(__func__, __error_type, __result, __errsum);
	      return std::make_tuple(__result, __errsum);
	    }

	  if (__error_type != NO_ERROR)
	    break;

	  if (__iteration >= __limit - 1)
	    {
	      __error_type = MAX_ITER_ERROR;
	      break;
	    }

	  if (__disallow_extrapolation)
	    continue;

	  __error_over_large_intervals -= __last_e_i;

	  if (__current_depth < __workspace.max_depth())
	    __error_over_large_intervals += __error12;

	  if (!__extrapolate)
	    {
	      // Test whether the interval to be bisected next is the
	      // smallest interval.
	      if (__workspace.large_interval())
		continue;
	      else
		{
		  __extrapolate = true;
		  __workspace.set_start(1);
		}
	    }

	  // The smallest interval has the largest error.  Before
	  // bisecting decrease the sum of the errors over the larger
	  // intervals (error_over_large_intervals) and perform
	  // extrapolation.

	  if (__error_type2 == NO_ERROR
	   && __error_over_large_intervals > __ertest)
	    if (__workspace.increment_start())
	      continue;

	  // Perform extrapolation.
	  __table.append(__area);

	  if (__table.get_nn() < 3)
	    {
	      __workspace.reset_start();
	      __extrapolate = false;
	      __error_over_large_intervals = __errsum;
	      continue;
	    }

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
	    __disallow_extrapolation = true;

	  if (__error_type == DIVERGENCE_ERROR)
	    break;

	  // Work on interval with largest error.
	  __workspace.reset_start();
	  __extrapolate = false;
	  __error_over_large_intervals = __errsum;
	}
      while (__iteration < __limit);

      auto __result = __res_ext;
      auto __abserr = __err_ext;

      if (__err_ext == _S_max)
	{
	  const auto __result = __workspace.total_integral();
	  __check_error(__func__, __error_type, __result, __errsum);
	  return std::make_tuple(__result, __errsum);
	}

      if (__error_type != NO_ERROR || __error_type2 != NO_ERROR)
	{
	  if (__error_type2 != NO_ERROR)
	    __err_ext += __correc;

	  if (__error_type == NO_ERROR)
	    __error_type = SINGULAR_ERROR;

	  if (__result != _Tp{0} && __area != _Tp{0})
	    {
	      if (__err_ext / std::abs(__res_ext) > __errsum / std::abs(__area))
		{
		  const auto __result = __workspace.total_integral();
		  __check_error(__func__, __error_type, __result, __errsum);
		  return std::make_tuple(__result, __errsum);
		}
	    }
	  else if (__err_ext > __errsum)
	    {
	      const auto __result = __workspace.total_integral();
	      __check_error(__func__, __error_type, __result, __errsum);
	      return std::make_tuple(__result, __errsum);
	    }
	  else if (__area == _Tp{0})
	    {
	      __check_error<_Tp>(__func__, __error_type);
	      std::__throw_runtime_error("qagp_integrate: Unknown error.");
	    }
	}

      // Test on divergence.
      auto __positive_integrand = __test_positivity(__result0, __resabs0);
      auto __max_area = std::max(std::abs(__res_ext), std::abs(__area));
      if (!__positive_integrand && __max_area < 0.01 * __resabs0)
	{
	  __check_error<_Tp>(__func__, __error_type);
	  std::__throw_runtime_error("qagp_integrate: Unknown error.");
	}

      auto __ratio = __res_ext / __area;
      if (__ratio < 0.01 || __ratio > 100 || __errsum > std::abs(__area))
	__error_type = UNKNOWN_ERROR;

      if (__error_type == NO_ERROR)
	return std::make_tuple(__result, __abserr);

      __check_error<_Tp>(__func__, __error_type);
      __throw__IntegrationError("qagp_integrate: Unknown error.",
				UNKNOWN_ERROR, __result, __abserr);
    }

} // namespace __gnu_test

#endif // QAGP_INTEGRATE_H
