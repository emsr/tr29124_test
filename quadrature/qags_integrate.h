// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2016 Free Software Foundation, Inc.
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
//This file implements Gauss-Kronrod integration with singularities
//Based upon gsl-1.9/integration/qags.c

#ifndef QAGS_INTEGRATE_H
#define QAGS_INTEGRATE_H

#include <utility>
#include <limits>
#include <tuple>
#include <stdexcept>
#include <functional>

#include "qk_integrate.h"
#include "integration_workspace.h"
#include "extrapolation_table.h"

namespace __gnu_test
{

  //Throws appropriate error if errcode nonzero
  void __check_error(int errcode);

  // Integrate potentially singular function from a to b using recursive
  // Gauss-Kronrod algorithm
  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qags_integrate(integration_workspace<_VecTp>& __workspace,
		   const _FuncTp& func,
		   const _VecTp a, const _VecTp b,
		   const _VecTp epsabs,
		   const _VecTp epsrel,
		   const std::size_t limit);

  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qags_integrate(const _FuncTp& func,
		   const _VecTp a, const _VecTp b,
		   const _VecTp epsabs,
		   const _VecTp epsrel,
		   const std::size_t limit)
    {
      integration_workspace<_VecTp> __workspace(limit);
      return qags_integrate(__workspace, func, a, b, epsabs, epsrel, limit);
    }

  // Integrate function from -infinity to +infinity
  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagi_integrate(integration_workspace<_VecTp>& __workspace,
		   const _FuncTp& func,
		   const _VecTp epsabs,
		   const _VecTp epsrel,
		   const std::size_t limit);

  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagi_integrate(const _FuncTp& func,
		   const _VecTp epsabs,
		   const _VecTp epsrel,
		   const std::size_t limit)
    {
      integration_workspace<_VecTp> __workspace(limit);
      return qagi_integrate(__workspace, func, epsabs, epsrel, limit);
    }

  // Integrate function from -infinity to b
  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagil_integrate(integration_workspace<_VecTp>& __workspace,
		    const _FuncTp& func, const _VecTp b,
		    const _VecTp epsabs,
		    const _VecTp epsrel,
		    const std::size_t limit);

  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagil_integrate(const _FuncTp& func, const _VecTp b,
		    const _VecTp epsabs,
		    const _VecTp epsrel,
		    const std::size_t limit)
    {
      integration_workspace<_VecTp> __workspace(limit);
      return qagil_integrate(__workspace, func, b, epsabs, epsrel, limit);
    }

  // Integrate function from a to +infinity
  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagiu_integrate(integration_workspace<_VecTp>& __workspace,
		    const _FuncTp& func, const _VecTp a,
		    const _VecTp epsabs,
		    const _VecTp epsrel,
		    const std::size_t limit);

  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagiu_integrate(const _FuncTp& func, const _VecTp a,
		    const _VecTp epsabs,
		    const _VecTp epsrel,
		    const std::size_t limit)
    {
      integration_workspace<_VecTp> __workspace(limit);
      return qagiu_integrate(__workspace, func, a, epsabs, epsrel, limit);
    }

  /**
   *
   */
  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qags_integrate(integration_workspace<_VecTp>& __workspace,
		   const _FuncTp& __func,
		   const _VecTp __a, const _VecTp __b,
		   _VecTp __epsabs,
		   _VecTp __epsrel,
		   const std::size_t __limit)
    {
      const qk_intrule __qkintrule = QK_21;

      _VecTp __area, __errsum;
      _VecTp __res_ext, __err_ext;
      _VecTp __result0, __abserr0, __resabs0, __resasc0;
      _VecTp __tolerance;

      _VecTp __ertest = 0;
      _VecTp __error_over_large_intervals = 0;
      _VecTp __reseps = 0, __abseps = 0, __correc = 0;
      std::size_t __ktmin = 0;
      int __roundoff_type1 = 0, __roundoff_type2 = 0, __roundoff_type3 = 0;
      int __error_type = 0, __error_type2 = 0;

      const auto _S_eps = std::numeric_limits<_VecTp>::epsilon();

      std::size_t __iteration = 0;

      bool __positive_integrand = false;
      int __extrapolate = 0;
      int __disallow_extrapolation = 0;

      /* Test on accuracy */

      if (__epsabs <= 0
	  && (__epsrel < 50 * _S_eps || __epsrel < 0.5e-28))
	std::__throw_logic_error("tolerance cannot be acheived"
			  " with given tolerances in qags_integrate()");

      /* Perform the first integration */

      typedef std::tuple<_VecTp&, _VecTp&, _VecTp&, _VecTp&> __ret_type;

      __ret_type{__result0, __abserr0, __resabs0, __resasc0}
	  = qk_integrate(__func, __a, __b, __qkintrule);

      __workspace.set_initial(__a, __b, __result0, __abserr0);

      __tolerance = std::max(__epsabs, __epsrel * std::abs(__result0));

      if (__abserr0 <= 100 * _S_eps * __resabs0
	  && __abserr0 > __tolerance)
	std::__throw_runtime_error("cannot reach tolerance because of roundoff"
			    " error on first attempt in qags_integrate()");
      else if ((__abserr0 <= __tolerance && __abserr0 != __resasc0)
		|| __abserr0 == 0.0)
	return std::make_pair(__result0, __abserr0);
      else if (__limit == 1)
	std::__throw_runtime_error("a maximum of one iteration was insufficient"
			    " in qags_integrate()");

      extrapolation_table<_VecTp> __table;
      __table.append(__result0);

      __area = __result0;
      __errsum = __abserr0;

      __res_ext = __result0;
      __err_ext = std::numeric_limits<_VecTp>::max();

      __positive_integrand = __test_positivity(__result0, __resabs0);

      __iteration = 1;

      do
	{
	  std::size_t __current_level;
	  _VecTp __a1, __b1, __a2, __b2;
	  _VecTp __a_i, __b_i, __r_i, __e_i;
	  _VecTp __area1 = 0, __area2 = 0, __area12 = 0;
	  _VecTp __error1 = 0, __error2 = 0, __error12 = 0;
	  _VecTp __resasc1, __resasc2;
	  _VecTp __resabs1, __resabs2;
	  _VecTp __last_e_i;

	  /* Bisect the subinterval with the largest error estimate */

	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  __current_level = __workspace.current_level();

	  __a1 = __a_i;
	  __b1 = 0.5 * (__a_i + __b_i);
	  __a2 = __b1;
	  __b2 = __b_i;

	  ++__iteration;

	  __ret_type{__area1, __error1, __resabs1, __resasc1}
	      = qk_integrate(__func, __a1, __b1, __qkintrule);

	  __ret_type{__area2, __error2, __resabs2, __resasc2}
	      = qk_integrate(__func, __a2, __b2, __qkintrule);

	  __area12 = __area1 + __area2;
	  __error12 = __error1 + __error2;
	  __last_e_i = __e_i;

	  /* Improve previous approximations to the integral and test for
	    accuracy.

	    We write these expressions in the same way as the original
	    QUADPACK code so that the rounding errors are the same, which
	    makes testing easier. */

	  __errsum += __error12 - __e_i;
	  __area += __area12 - __r_i;

	  __tolerance = std::max(__epsabs, __epsrel * std::abs(__area));

	  if (__resasc1 != __error1 && __resasc2 != __error2)
	    {
	      _VecTp __delta = __r_i - __area12;

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

	  /* Test for roundoff and eventually set error flag */

	  if (__roundoff_type1 + __roundoff_type2 >= 10
	      || __roundoff_type3 >= 20)
	    __error_type = 2;       /* round off error */

	  if (__roundoff_type2 >= 5)
	    __error_type2 = 1;

	  /* set error flag in the case of bad integrand behaviour at
	    a point of the integration range */

	  if (__workspace.subinterval_too_small (__a1, __a2, __b2))
	    __error_type = 4;

	  /* append the newly-created intervals to the list */
	  __workspace.update(__a1, __b1, __area1, __error1,
			     __a2, __b2, __area2, __error2);

	  if (__errsum <= __tolerance)
	    {
	      __check_error(__error_type);
	      return std::make_pair(__workspace.sum_results(), __errsum);
	    }

	  if (__error_type)
	    break;

	  if (__iteration >= __limit - 1)
	    {
	      __error_type = 1;
	      break;
	    }

	  if (__iteration == 2)       /* set up variables on first iteration */
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
	      /* test whether the interval to be bisected next is the
		smallest interval. */

	      if (__workspace.large_interval())
		continue;

	      __extrapolate = 1;
	      __workspace.set_nrmax(1);
	    }

	  if (!__error_type2 && __error_over_large_intervals > __ertest)
	    if (__workspace.increase_nrmax())
	      continue;

	  /* Perform extrapolation */

	  __table.append(__area);

	  std::pair<_VecTp&, _VecTp&> __qelg_res(__reseps, __abseps);
	  __qelg_res = __table.qelg();

	  ++__ktmin;

	  if (__ktmin > 5 && __err_ext < 0.001 * __errsum)
	    __error_type = 5;

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

	  /* Prepare bisection of the smallest interval. */

	  if (__table.get_nn() == 1)
	    __disallow_extrapolation = 1;

	  if (__error_type == 5)
	    break;

	  /* work on interval with largest error */

	  __workspace.reset_nrmax();
	  __extrapolate = 0;
	  __error_over_large_intervals = __errsum;

	}
      while (__iteration < __limit);

      //__result = __res_ext;
      //__abserr = __err_ext;

      if (__err_ext == std::numeric_limits<_VecTp>::max())
      {
	__check_error(__error_type);
	return std::make_pair(__workspace.sum_results(), __errsum);
      }

      if (__error_type || __error_type2)
	{
	  if (__error_type2)
	    __err_ext += __correc;

	  if (__error_type == 0)
	    __error_type = 3;

	  if (__res_ext != 0.0 && __area != 0.0)
	    {
	      if (__err_ext / std::abs(__res_ext) > __errsum / std::abs(__area))
		{
		  __check_error(__error_type);
		  return std::make_pair(__workspace.sum_results(), __errsum);
		}
	    }
	  else if (__err_ext > __errsum)
	    {
	      __check_error(__error_type);
	      return std::make_pair(__workspace.sum_results(), __errsum);
	    }
	  else if (__area == 0.0)
	    {
	      __check_error(__error_type);
	      std::__throw_runtime_error("qags_integrate: Unknown error");
	    }
	}

      /*  Test on divergence. */

      {
	_VecTp __max_area = std::max(std::abs(__res_ext), std::abs(__area));

	if (!__positive_integrand && __max_area < 0.01 * __resabs0)
	{
	  __check_error(__error_type);
	  std::__throw_runtime_error("qags_integrate: Unknown error");
	}
      }

      {
	_VecTp __ratio = __res_ext / __area;

	if (__ratio < 0.01 || __ratio > 100.0 || __errsum > std::abs(__area))
	  __error_type = 6;
      }

      __check_error(__error_type);
      std::__throw_runtime_error("qags_integrate: Unknown error");
    }

  // Throws appropriate error if errcode nonzero
  void
  __check_error(int __errcode)
  {
    if (__errcode > 2)
      --__errcode;
    switch(__errcode)
      {
      case 0: break;
      case 1:
	std::__throw_runtime_error("qags_integrate: "
				   "Number of iterations was insufficient");
      case 2:
	std::__throw_runtime_error("qags_integrate: "
				   "Cannot reach tolerance "
				   "because of roundoff");
      case 3:
	std::__throw_runtime_error("qags_integrate: "
				   "Bad integrand behavior found "
				   "in the integration interval");
      case 4:
	std::__throw_runtime_error("qags_integrate: "
				   "Roundoff error detected "
				   "in the extrapolation");
      case 5:
	std::__throw_runtime_error("qags_integrate: "
				   "Integral is divergent, "
				   "or slowly convergent");
      default:
	std::__throw_runtime_error("qags_integrate: "
				   "Could not integrate function");
      }
  }

  template<typename _FuncTp, typename _VecTp>
    _VecTp i_transform(const _FuncTp& __func, _VecTp __t)
    {
      _VecTp __x = (1 - __t) / __t;
      _VecTp __y = __func(__x) + __func(-__x);
      return (__y / __t) / __t;
    }

  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagi_integrate(integration_workspace<_VecTp>& __workspace,
		   const _FuncTp& __func,
		   _VecTp __epsabs, _VecTp __epsrel,
		   const std::size_t __limit)
    {
      return qags_integrate(__workspace,
	  std::bind(i_transform<_FuncTp, _VecTp>, __func, std::placeholders::_1),
	      _VecTp{0}, _VecTp{1}, __epsabs, __epsrel, __limit);
    }

  template<typename _FuncTp, typename _VecTp>
    _VecTp
    il_transform(const _FuncTp& __func, _VecTp __b, _VecTp __t)
    {
      _VecTp __x = __b - (1 - __t) / __t;
      _VecTp __y = __func(__x);
      return (__y / __t) / __t;
    }

  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagil_integrate(integration_workspace<_VecTp>& __workspace,
		    const _FuncTp& __func, _VecTp __b,
		    _VecTp __epsabs, _VecTp __epsrel,
		    const std::size_t __limit)
    {
      return qags_integrate(__workspace,
	std::bind(il_transform<_FuncTp, _VecTp>, __func, __b, std::placeholders::_1),
	    _VecTp{0}, _VecTp{1}, __epsabs, __epsrel, __limit);
    }

  template<typename _FuncTp, typename _VecTp>
    _VecTp
    iu_transform(const _FuncTp& __func, _VecTp __a, _VecTp __t)
    {
      _VecTp __x = __a + (1 - __t) / __t;
      _VecTp __y = __func(__x);
      return (__y / __t) / __t;
    }

  template<typename _FuncTp, typename _VecTp>
    std::pair<_VecTp, _VecTp>
    qagiu_integrate(integration_workspace<_VecTp>& __workspace,
		    const _FuncTp& __func, _VecTp __a,
		    _VecTp __epsabs, _VecTp __epsrel,
		    const std::size_t __limit)
    {
      return qags_integrate(__workspace,
	std::bind(iu_transform<_FuncTp, _VecTp>, __func, __a, std::placeholders::_1),
	    _VecTp{0}, _VecTp{1}, __epsabs, __epsrel, __limit);
    }
} // namespace

#endif
