/* integration/qagp.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2016-2017 Free Software Foundation, Inc.
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
// Ported from GSL by Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements qawc integration algorithm
// Based on gsl/integration/qawc.c

#ifndef QAWC_INTEGRATE_TCC
#define QAWC_INTEGRATE_TCC 1

#include <array>
#include <tuple>

#include "integration_workspace.h"
#include "qk_integrate.tcc"

namespace __gnu_ext
{

  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp, bool>
    qc25c(const _FuncTp& __func, _Tp __a, _Tp __b, _Tp __c);

  template<typename _Tp>
    std::vector<_Tp>
    compute_moments(std::size_t __N, _Tp __cc);

  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qawc_integrate(integration_workspace<_Tp>& __workspace,
		   const _FuncTp& __func,
		   const _Tp __a, const _Tp __b, const _Tp __c,
		   const _Tp __epsabs, const _Tp __epsrel)
    {
      auto __result = _Tp{};
      auto __abserr = _Tp{};

      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto __limit = __workspace.capacity();

      int __sign = 1;
      _Tp __lower, __higher;
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

      if (__epsabs <= 0 && (__epsrel < 50 * _S_eps
			   || __epsrel < 0.5e-28))
	std::__throw_runtime_error ("qawc_integrate: "
				    "Tolerance cannot be achieved "
				    "with given absolute "
				    "and relative error limits.");

      if (__c == __a || __c == __b)
	std::__throw_runtime_error ("qawc_integrate: "
				    "Cannot integrate with singularity "
				    "on endpoint.");

      __workspace.clear();

      // Perform the first integration.
      _Tp __result0, __abserr0;
      bool __err_reliable;
      std::tie(__result0, __abserr0, __err_reliable)
	= qc25c(__func, __lower, __higher, __c);

      __workspace.append(__lower, __higher, __result0, __abserr0);

      // Test on accuracy; Use 0.01 relative error as an extra safety
      // margin on the first iteration (ignored for subsequent iterations).
      auto __tolerance = std::max(__epsabs, __epsrel * std::abs( __result0));
      if (__abserr0 < __tolerance && __abserr0 < 0.01 * std::abs(__result0))
	return std::make_tuple(__sign * __result0, __abserr0);
      else if (__limit == 1)
	__throw__IntegrationError("qawc_integrate: "
				  "a maximum of one iteration was insufficient",
				  MAX_ITER_ERROR, __sign * __result0, __abserr0);

      auto __area = __result0;
      auto __errsum = __abserr0;
      auto __iteration = 1u;
      int __error_type = NO_ERROR;
      int __roundoff_type1 = 0, __roundoff_type2 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  _Tp __a_i, __b_i, __r_i, __e_i;
	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  const auto __a1 = __a_i;
	  auto __b1 = (__a_i + __b_i) / _Tp{2};
	  const auto __b2 = __b_i;
	  if (__c > __a1 && __c <= __b1)
	    __b1 = (__c + __b2) / _Tp{2};
	  else if (__c > __b1 && __c < __b2)
	    __b1 = (__a1 + __c) / _Tp{2};
	  const auto __a2 = __b1;

	  _Tp __area1, __error1;
	  bool __err_reliable1;
	  std::tie(__area1, __error1, __err_reliable1)
	    = qc25c(__func, __a1, __b1, __c);

	  _Tp __area2, __error2;
	  bool __err_reliable2;
	  std::tie(__area2, __error2, __err_reliable2)
	    = qc25c(__func, __a2, __b2, __c);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;

	  __errsum += (__error12 - __e_i);
	  __area += __area12 - __r_i;

	  if (__err_reliable1 && __err_reliable2)
	    {
	      _Tp __delta = __r_i - __area12;

	      if (std::abs (__delta) <= 1.0e-5 * std::abs (__area12)
	    	   && __error12 >= 0.99 * __e_i)
		++__roundoff_type1;
	      if (__iteration >= 10 && __error12 > __e_i)
		++__roundoff_type2;
	    }

	  __tolerance = std::max(__epsabs, __epsrel * std::abs(__area));
	  if (__errsum > __tolerance)
	    {
	      if (__roundoff_type1 >= 6 || __roundoff_type2 >= 20)
		__error_type = ROUNDOFF_ERROR;

	      // set error flag in the case of bad integrand behaviour at
	      // a point of the integration range
	      if (__workspace.subinterval_too_small(__a1, __a2, __b2))
		__error_type = SINGULAR_ERROR;
	    }

	  __workspace.split(__b1, __area1, __error1, __area2, __error2);

	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  ++__iteration;
	}
      while (__iteration < __limit && !__error_type && __errsum > __tolerance);

      __result = __sign * __workspace.total_integral();
      __abserr = __errsum;

      if (__iteration == __limit)
	__error_type = MAX_SUBDIV_ERROR;

      if (__errsum <= __tolerance)
	return std::make_tuple(__result, __abserr);

      if (__error_type == NO_ERROR)
	return std::make_tuple(__result, __abserr);

      __check_error<_Tp>(__func__, __error_type);
      __throw__IntegrationError("qawc_integrate: Unknown error.",
				UNKNOWN_ERROR, __result, __abserr);
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
	  std::array<_Tp, 13> __cheb12;
	  std::array<_Tp, 25> __cheb24;
	  qcheb_integrate(__func, __a, __b, __cheb12, __cheb24);
	  const auto __moment = compute_moments(__cheb24.size(), __cc);

	  auto __res12 = _Tp{0};
	  for (size_t __i = 0u; __i < __cheb12.size(); ++__i)
	    __res12 += __cheb12[__i] * __moment[__i];

	  auto __res24 = _Tp{0};
	  for (size_t __i = 0u; __i < __cheb24.size(); ++__i)
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

#endif // QAWC_INTEGRATE_TCC
