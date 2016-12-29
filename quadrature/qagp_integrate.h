/* integration/qagp_integrate.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2016 Edward Smith-Rowland
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef QAGP_INTEGRATE_H
#define QAGP_INTEGRATE_H 1

#include "qk_integrate.h"
#include "integration_workspace.h"
#include "extrapolation_table.h"

namespace __gnu_test
{

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
      _Tp __area, __errsum;
      _Tp __res_ext, __err_ext;
      _Tp __tolerance;

      _Tp __ertest = 0;
      _Tp __error_over_large_intervals = 0;
      _Tp __reseps = 0, __abseps = 0, __correc = 0;
      std::size_t __ktmin = 0;
      int __roundoff_type1 = 0, __roundoff_type2 = 0, __roundoff_type3 = 0;
      int __error_type = 0, __error_type2 = 0;

      std::size_t __iteration = 0;

      int __positive_integrand = 0;
      int __extrapolate = 0;
      int __disallow_extrapolation = 0;

      const auto _S_max = std::numeric_limits<_Tp>::max();
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto __limit = __workspace.capacity();

      if (__pts.size() > __workspace.capacity())
	std::__throw_runtime_error("qagp_integrate: "
				   "number of pts exceeds size of workspace");
      if (__epsabs <= 0 && (__epsrel < 50 * _S_eps
			   || __epsrel < 0.5e-28))
	std::__throw_runtime_error("qagp_integrate: "
				   "tolerance cannot be achieved "
				   "with given epsabs and epsrel");

      /* Check that the integration range and break points are an
	 ascending sequence */

      // number of intervals
      const std::size_t __nint = __pts.size() - 1;
      for (std::size_t __i = 0; __i < __nint; ++__i)
	if (__pts[__i + 1] < __pts[__i])
          std::__throw_runtime_error("qagp_integrate: "
				     "points are not in an ascending sequence");

      // Perform the first integration.

      __workspace.set_initial_limits(_Tp{0}, _Tp{0});

      auto __result0 = _Tp{0};
      auto __abserr0 = _Tp{0};
      auto __resabs0 = _Tp{0};
      for (std::size_t __i = 0; __i < __nint; ++__i)
	{
	  _Tp __area1, __error1, __resabs1, __resasc1;
	  const auto __a1 = __pts[__i];
	  const auto __b1 = __pts[__i + 1];

	  std::tie(__area1, __error1, __resabs1, __resasc1)
	    = qk_integrate(__func, __a1, __b1, __qk_rule);

	  __result0 += __area1;
	  __abserr0 += __error1;
	  __resabs0 += __resabs1;
	  __workspace.append(__a1, __b1, __area1, __error1);

	  if (__error1 == __resasc1 && __error1 != 0.0)
            __workspace.set_level(__i, 1);
	  else
            __workspace.set_level(__i, 0);
	}

      // Compute the initial error estimate.
      __errsum = 0.0;
      for (std::size_t __i = 0; __i < __nint; ++__i)
	{
	  if (__workspace.level(__i))
            __workspace.set_abs_error(__i, __abserr0);
	  __errsum += __workspace.abs_error(__i);
	}

      for (std::size_t __i = 0; __i < __nint; ++__i)
	__workspace.set_level(__i, 0);

      // Sort results into order of decreasing error via the indirection
      // array order[]

      __workspace.sort_error();

      // Test on accuracy.
      __tolerance = std::max(__epsabs, __epsrel * std::abs(__result0));

      if (__abserr0 <= 100 * _S_eps * __resabs0
	   && __abserr0 > __tolerance)
	std::__throw_runtime_error("qagp_integrate: "
				   "cannot reach tolerance because "
				   "of roundoff error on first attempt");
      else if (__abserr0 <= __tolerance)
	return std::make_tuple(__result0, __abserr0);
      else if (__limit == 1)
	std::__throw_runtime_error("qagp_integrate: "
				   "a maximum of one iteration was insufficient");

      /* Initialization */

      extrapolation_table<_Tp> __table(__result0);

      __area = __result0;

      __res_ext = __result0;
      __err_ext = _S_max;

      __error_over_large_intervals = __errsum;
      __ertest = __tolerance;

      __positive_integrand = __test_positivity(__result0, __resabs0);

      __iteration = __nint - 1; 

      do
	{
	  // Bisect the subinterval with the largest error estimate.

	  _Tp __a_i, __b_i, __r_i, __e_i;
	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  std::size_t __current_level = __workspace.current_level() + 1;

	  auto __a1 = __a_i;
	  auto __b1 = 0.5 * (__a_i + __b_i);
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

	  // Improve previous approximations to the integral and test for accuracy.

	  __errsum += __error12 - __e_i;
	  __area += __area12 - __r_i;

	  __tolerance = std::max(__epsabs, __epsrel * std::abs(__area));

	  if (__resasc1 != __error1 && __resasc2 != __error2)
            {
              _Tp __delta = __r_i - __area12;

              if (std::abs(__delta) <= 1.0e-5 * std::abs(__area12)
		   && __error12 >= 0.99 * __e_i)
        	{
        	  if (!__extrapolate)
                    ++__roundoff_type1;
        	  else
                    ++__roundoff_type2;
        	}
              // This "i > 10" is just a test on nint.
              if (/*i > 10 &&*/ __error12 > __e_i)
        	++__roundoff_type3;
            }

	  // Test for roundoff and eventually set error flag.

	  if (__roundoff_type1 + __roundoff_type2 >= 10 || __roundoff_type3 >= 20)
            __error_type = 2; // round off error

	  if (__roundoff_type2 >= 5)
            __error_type2 = 1;

	  // Set error flag in the case of bad integrand behaviour at
          // a point of the integration range

	  if (integration_workspace<_Tp>::subinterval_too_small(__a1, __a2, __b2))
            __error_type = 4;

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
              __error_type = 1;
              break;
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

	  /* The smallest interval has the largest error.  Before
             bisecting decrease the sum of the errors over the larger
             intervals (error_over_large_intervals) and perform
             extrapolation. */

	  if (!__error_type2 && __error_over_large_intervals > __ertest)
            if (__workspace.increase_nrmax())
              continue;

	  // Perform extrapolation.

	  __table.append(__area);

	  if (__table.get_nn() < 3) 
            goto skip_extrapolation;

	  std::tie(__reseps, __abseps) = __table.qelg();

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

	  // Prepare bisection of the smallest interval.

	  if (__table.get_nn() == 1)
            __disallow_extrapolation = 1;

	  if (__error_type == 5)
            break;

	skip_extrapolation:

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
            __error_type = 3;
	  if (__result != _Tp{0} && __area != _Tp{0})
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
	      std::__throw_runtime_error("qagp_integrate: Unknown error");
	    }
	}

      // Test on divergence.
      auto __max_area = std::max(std::abs(__res_ext), std::abs(__area));
      if (!__positive_integrand && __max_area < 0.01 * __resabs0)
	{
	  __check_error<_Tp>(__func__, __error_type);
	  std::__throw_runtime_error("qagp_integrate: Unknown error");
	}

      auto __ratio = __res_ext / __area;
      if (__ratio < 0.01 || __ratio > 100 || __errsum > std::abs(__area))
	__error_type = 7;

      if (__error_type == 0)
	return std::make_tuple(__result, __abserr);

      __check_error<_Tp>(__func__, __error_type);
      std::__throw_runtime_error("qagp_integrate: Unknown error");
    }

} // namespace __gnu_test

#endif // QAGP_INTEGRATE_H
