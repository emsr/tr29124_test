/* integration/qagp.c
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

#ifndef QAWC_INTEGRATE_H
#define QAWC_INTEGRATE_H 1

#include <tuple>

#include "integration_workspace.h"
#include "qk_integrate.h"

namespace __gnu_test
{

template<typename _FuncTp, typename _Tp>
  std::tuple<_Tp, _Tp, bool>
  qc25c(const _FuncTp& __func, _Tp __a, _Tp __b, _Tp __c);

template<typename _Tp>
  std::vector<_Tp>
  compute_moments(std::size_t __N, _Tp __cc);

template<typename _FuncTp, typename _Tp>
  std::pair<_Tp, _Tp>
  qawc_integrate(integration_workspace<_Tp>& __workspace,
                 const _FuncTp& __func,
                 const _Tp __a, const _Tp __b, const _Tp __c,
                 const _Tp __epsabs, const _Tp __epsrel,
                 const size_t __limit)
  {
    _Tp __area, __errsum;
    _Tp __result0, __abserr0;
    _Tp __tolerance;
    size_t __iteration = 0;
    int __roundoff_type1 = 0, __roundoff_type2 = 0, __error_type = 0;
    int __err_reliable;
    int __sign = 1;
    _Tp __lower, __higher;

    auto __result = _Tp{};
    auto __abserr = _Tp{};

    if (__limit > __workspace.capacity())
      std::__throw_runtime_error ("iteration limit exceeds available workspace") ;

    if (__b < __a) 
      {
	__lower = __b;
	__higher = __a;
	__sign = -1;
      }
    else
      {
	__lower = __a;
	__higher = __b;
      }

    if (__epsabs <= 0 && (__epsrel < 50 * std::numeric_limits<_Tp>::epsilon()
			 || __epsrel < 0.5e-28))
      std::__throw_runtime_error ("tolerance cannot be achieved with given tolerances");

    if (__c == __a || __c == __b) 
      std::__throw_runtime_error ("cannot integrate with singularity on endpoint");

    /* perform the first integration */

    std::tie(__result0, __abserr0, __err_reliable)
      = qc25c(__func, __lower, __higher, __c);

    __workspace.set_initial(__lower, __higher, __result0, __abserr0);

    /* Test on accuracy, use 0.01 relative error as an extra safety
       margin on the first iteration (ignored for subsequent iterations) */

    __tolerance = std::max(__epsabs, __epsrel * std::abs( __result0));

    if (__abserr0 < __tolerance && __abserr0 < 0.01 * std::abs(__result0)) 
      {
	__result = __sign * __result0;
	__abserr = __abserr0;

	return std::make_pair(__result, __abserr);
      }
    else if (__limit == 1)
      {
	__result = __sign * __result0;
	__abserr = __abserr0;

	std::__throw_runtime_error("a maximum of one iteration was insufficient");
      }

    __area = __result0;
    __errsum = __abserr0;

    __iteration = 1;

    do
      {
	_Tp __a1, __b1, __a2, __b2;
	_Tp __a_i, __b_i, __r_i, __e_i;
	_Tp __area1 = 0, __area2 = 0, __area12 = 0;
	_Tp __error1 = 0, __error2 = 0, __error12 = 0;
	int __err_reliable1, __err_reliable2;

	/* Bisect the subinterval with the largest error estimate */

	__workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	__a1 = __a_i; 
	__b1 = 0.5 * (__a_i + __b_i);
	__a2 = __b1;
	__b2 = __b_i;

	if (__c > __a1 && __c <= __b1) 
          {
            __b1 = 0.5 * (__c + __b2) ;
            __a2 = __b1;
          }
	else if (__c > __b1 && __c < __b2)
          {
            __b1 = 0.5 * (__a1 + __c) ;
            __a2 = __b1;
          }

	std::tie(__area1, __error1, __err_reliable1) = qc25c(__func, __a1, __b1, __c);
	std::tie(__area2, __error2, __err_reliable2) = qc25c(__func, __a2, __b2, __c);

	__area12 = __area1 + __area2;
	__error12 = __error1 + __error2;

	__errsum += (__error12 - __e_i);
	__area += __area12 - __r_i;

	if (__err_reliable1 && __err_reliable2)
          {
            _Tp __delta = __r_i - __area12;

            if (std::abs (__delta) <= 1.0e-5 * std::abs (__area12) && __error12 >= 0.99 * __e_i)
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

            if (integration_workspace<_Tp>::subinterval_too_small(__a1, __a2, __b2))
              __error_type = 3;
          }

	__workspace.update(__a1, __b1, __area1, __error1,
			   __a2, __b2, __area2, __error2);

	__workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	++__iteration;
      }
    while (__iteration < __limit && !__error_type && __errsum > __tolerance);

    __result = __sign * __workspace.sum_results();
    __abserr = __errsum;

    if (__errsum <= __tolerance)
      return std::make_pair(__result, __abserr);
    else if (__error_type == 2)
      std::__throw_runtime_error ("roundoff error prevents tolerance from being achieved");
    else if (__error_type == 3)
      std::__throw_runtime_error ("bad integrand behavior found in the integration interval");
    else if (__iteration == __limit)
      std::__throw_runtime_error ("maximum number of subdivisions reached");
    else
      std::__throw_runtime_error ("could not integrate function");
  }

  /**
   * 
   */
  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp, bool>
    qc25c(const _FuncTp& __func, _Tp __a, _Tp __b, _Tp __c)
    {
      const _Tp __cc = (_Tp{2} * __c - __b - __a) / (__b - __a);

      _Tp __result, __abserr;
      bool __err_reliable;

      if (std::abs(__cc) > 1.1)
	{
          using __qk_ret = std::tuple<_Tp&, _Tp&, _Tp&, _Tp&>;

	  auto __func_cauchy = [__func, __c](_Tp __x)
				-> _Tp
				{ return __func(__x) / (__x - __c); };

	  _Tp __resabs, __resasc;
          __qk_ret{__result, __abserr, __resabs, __resasc}
	    = qk_integrate(__func_cauchy, __a, __b, QK_15);

	  if (__abserr == __resasc)
            __err_reliable = false;
	  else 
            __err_reliable = true;

	  return std::make_tuple(__result, __abserr, __err_reliable);
	}
      else
	{
	  const std::size_t __N = 25;
	  const std::size_t __M = 1 + __N / 2;
	  _Tp __cheb12[__M], __cheb24[__N];
	  _Tp __res12 = 0, __res24 = 0;
	  qcheb_integrate(__func, __a, __b, __cheb12, __cheb24);
	  auto __moment = compute_moments(__N, __cc);

	  for (size_t __i = 0u; __i < __M; ++__i)
            __res12 += __cheb12[__i] * __moment[__i];

	  for (size_t __i = 0u; __i < __N; ++__i)
            __res24 += __cheb24[__i] * __moment[__i];

	  __result = __res24;
	  __abserr = std::abs(__res24 - __res12);
	  __err_reliable = false;

	  return std::make_tuple(__result, __abserr, __err_reliable);
	}
    }

  /**
   * 
   */
  template<typename _Tp>
    std::vector<_Tp>
    compute_moments(std::size_t __N, _Tp __cc)
    {
      std::vector<_Tp> __moment(__N);

      auto __a0 = std::log(std::abs((_Tp{1} - __cc) / (_Tp{1} + __cc)));
      auto __a1 = _Tp{2} + __a0 * __cc;

      __moment[0] = __a0;
      __moment[1] = __a1;

      for (size_t __k = 2; __k < __N; ++__k)
	{
	  _Tp __a2;

	  if ((__k % 2) == 0)
            __a2 = _Tp{2} * __cc * __a1 - __a0;
	  else
            {
              const auto __km1 = _Tp(__k - 1);
              __a2 = _Tp{2} * __cc * __a1
		   - __a0
		   - _Tp{4} / (__km1 * __km1 - _Tp{1});
            }

	  __moment[__k] = __a2;

	  __a0 = __a1;
	  __a1 = __a2;
	}

      return __moment;
    }

} // namespace

#endif // QAWC_INTEGRATE_H
