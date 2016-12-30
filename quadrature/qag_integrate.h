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
// Implements integration using a recursive Gauss-Kronrod algorithm
// Based upon gsl-2.3/integration/qag.c

#ifndef QAG_INTEGRATE_H
#define QAG_INTEGRATE_H 1

#include <utility>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>

#include "qk_integrate.h"
#include "integration_workspace.h"

namespace __gnu_test
{

  /**
   * Integrates a function from a to b using the Gauss-Kronrod
   * quadrature rule.  The integration limits are finite.
   *
   * Once either the absolute or relative error limit is reached,
   * qag_integrate() returns
   * @param[in] __func The single-variable function to be integrated
   * @param[in] __a The lower limit of integration
   * @param[in] __b The upper limit of integration
   * @param[in] __max_iter The maximum number of integration steps allowed
   * @param[in] __epsabs The limit on absolute error
   * @param[in] __epsrel The limit on relative error
   * @param[in] __qksz The size of the Gauss-Kronrod integration scheme
   * @return A tuple with the first value being the integration result,
   *	     and the second value being the estimated error.
   */
  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp>
    qag_integrate(integration_workspace<_Tp>& __workspace,
		  const _FuncTp& __func, _Tp __a, _Tp __b,
		  _Tp __epsabs, _Tp __epsrel, const std::size_t __max_iter,
		  const qk_intrule __qkintrule)
    {
      _Tp __area, __errsum;
      _Tp __result0, __abserr0, __resabs0, __resasc0;
      _Tp __tolerance;
      std::size_t __iteration = 0;
      int __roundoff_type1 = 0, __roundoff_type2 = 0, __error_type = 0;

      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      _Tp __round_off;

      _Tp __result = 0;
      _Tp __abserr = 0;

      std::vector<_Tp> __rlist(__max_iter);
      std::vector<_Tp> __elist(__max_iter);

      if (__epsabs <= 0 && (__epsrel < 50 * _S_eps))
	std::__throw_logic_error("tolerance cannot be achieved"
			  " in qag_integrate() with given absolute"
			  " and relative error limits");

      typedef std::tuple<_Tp&,_Tp&,_Tp&,_Tp&> __ret_type;

      __ret_type{__result0, __abserr0,__resabs0,__resasc0}
	  = qk_integrate(__func, __a, __b, __qkintrule);

      __rlist[0] = __result0;
      __elist[0] = __abserr0;

      //Test on accuracy
      __tolerance = std::max(__epsabs, __epsrel * std::abs(__result0));

      //Compute roundoff
      __round_off = 50 * _S_eps * __resabs0;

      if (__abserr0 <= __round_off && __abserr0 > __tolerance)
	{
	  __result = __result0;
	  __abserr = __abserr0;

	  std::__throw_runtime_error("qag_integrate: "
				     "Cannot reach tolerance because"
				     " of roundoff error on first attempt");
	}
      else if ((__abserr0 <= __tolerance && __abserr0 != __resasc0)
		|| __abserr0 == 0.0)
	{
	  __result = __result0;
	  __abserr = __abserr0;

	  return std::make_tuple(__result, __abserr);
	}
      else if (__max_iter == 1)
	{
	  __result = __result0;
	  __abserr = __abserr0;

	  std::__throw_runtime_error("qag_integrate: "
				     "a maximum of one iteration"
				     " was insufficient");
	}

      __area = __result0;
      __errsum = __abserr0;

      __iteration = 1;

      __workspace.set_initial(__a, __b, __result0, __abserr0);
      __result = __workspace.sum_results();

      do
	{
	  _Tp __a1, __b1, __a2, __b2;
	  _Tp __a_i, __b_i, __r_i, __e_i;
	  _Tp __area1 = 0, __area2 = 0, __area12 = 0;
	  _Tp __error1 = 0, __error2 = 0, __error12 = 0;
	  _Tp __resasc1, __resasc2;
	  _Tp __resabs1, __resabs2;

	  // Bisect the subinterval with the largest error estimate

	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  __a1 = __a_i;
	  __b1 = 0.5 * (__a_i + __b_i);
	  __a2 = __b1;
	  __b2 = __b_i;

	  __ret_type{__area1,__error1,__resabs1,__resasc1}
	      = qk_integrate(__func, __a1, __b1, __qkintrule);

	  __ret_type{__area2,__error2,__resabs2,__resasc2}
	      = qk_integrate(__func, __a2, __b2, __qkintrule);

	  __area12 = __area1 + __area2;
	  __error12 = __error1 + __error2;

	  __errsum += (__error12 - __e_i);
	  __area += __area12 - __r_i;

	  if (__resasc1 != __error1 && __resasc2 != __error2)
	    {
	      _Tp __delta = __r_i - __area12;

	      if (std::abs(__delta) <= 1.0e-5 * std::abs(__area12)
		 && __error12 >= 0.99 * __e_i)
		++__roundoff_type1;
	      if (__iteration >= 10 && __error12 > __e_i)
		++__roundoff_type2;
	    }

	  __tolerance = std::max(__epsabs, __epsrel * std::abs(__area));

	  if (__errsum > __tolerance)
	    {
	      if (__roundoff_type1 >= 6 || __roundoff_type2 >= 20)
		__error_type = 2;   /* round off error */

	      /* set error flag in the case of bad integrand behaviour at
		a point of the integration range */

	      if (__workspace.subinterval_too_small(__a1, __a2, __b2))
		__error_type = 3;
	    }

	  __workspace.update(__a1, __b1, __area1, __error1,
			     __a2, __b2, __area2, __error2);

	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  __result = __workspace.sum_results();

	  ++__iteration;
	}
      while (__iteration < __max_iter && !__error_type
	     && __errsum > __tolerance);

      __result = __workspace.sum_results();
      __abserr = __errsum;

      if (__errsum <= __tolerance)
	return std::make_tuple(__result, __abserr);
      else if (__error_type == 2)
	std::__throw_runtime_error("qag_integrate: "
				   "Cannot reach tolerance "
				   "because of roundoff error");
      else if (__error_type == 3)
	std::__throw_runtime_error("qag_integrate: "
				   "Bad integrand behavior found "
				   "in the integrand inteveral");
      else if (__iteration == __max_iter)
	std::__throw_runtime_error("qag_integrate: "
				   "Maximum number of iterations reached");
      else
	std::__throw_runtime_error("qag_integrate: "
				   "Could not integrate function");
    }

  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp>
    qag_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		  _Tp __epsabs, _Tp __epsrel, const std::size_t __max_iter,
		  const qk_intrule __qkintrule)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qag_integrate(__workspace, __func, __a, __b,
			   __epsabs, __epsrel, __max_iter,
			   __qkintrule);
    }

} // namespace __gnu_test

#endif // QAG_INTEGRATE_H
