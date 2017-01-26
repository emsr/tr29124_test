/* quadrature/qaws_integrate.h
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
// Implements qawo integration algorithm
// Based on gsl/integration/qawo.c

#ifndef QAWS_INTEGRATE_H
#define QAWS_INTEGRATE_H 1

#include <array>
#include <tuple>

#include "qk_integrate.h"
#include "integration_workspace.h"
#include "qaws_integration_table.h"

namespace __gnu_test
{

  template<typename _Tp>
    struct fn_qaws;

  template<typename _Tp>
    std::tuple<_Tp, _Tp>
    compute_result(const std::array<_Tp, 25>& r,
		   const std::array<_Tp, 13>& cheb12,
		   const std::array<_Tp, 25>& cheb24);

  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp>
    qaws_integrate(integration_workspace<_Tp>& __workspace,
		   qaws_integration_table<_Tp>& __table,
		   const _FuncTp& __func,
		   const _Tp __a, const _Tp __b,
		   const _Tp __epsabs, const _Tp __epsrel)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      if (__b <= __a)
	std::__throw_runtime_error("qaws_integrate: "
				   "Limits must form an ascending sequence");
      if (__epsabs <= 0 && (__epsrel < 50 * _S_eps || __epsrel < 0.5e-28))
	std::__throw_runtime_error("qaws_integrate: "
				   "Tolerance cannot be achieved "
				   "with given absolute "
				   "and relative error limits.");

      const auto __limit = __workspace.capacity();
      __workspace.clear();

      // Perform the first integration.
      _Tp __result0, __abserr0;
      {
	const auto __a1 = __a;
	const auto __b1 = (__a + __b) / _Tp{2};
	const auto __a2 = __b1;
	const auto __b2 = __b;

	_Tp __area1, __error1;
	bool __err_reliable1;
	std::tie(__area1, __error1, __err_reliable1)
	  = qc25s(__table, __func, __a, __b, __a1, __b1);
	__workspace.append(__a1, __b1, __area1, __error1);

	_Tp __area2, __error2;
	bool __err_reliable2;
	std::tie(__area2, __error2, __err_reliable2)
	  = qc25s(__table, __func, __a, __b, __a2, __b2);
	__workspace.append(__a2, __b2, __area2, __error2);

	__result0 = __area1 + __area2;
	__abserr0 = __error1 + __error2;
      }

      // Test on accuracy; Use 0.01 relative error as an extra safety
      // margin on the first iteration (ignored for subsequent iterations).
      auto __tolerance = std::max(__epsabs, __epsrel * std::abs(__result0));
      if (__abserr0 < __tolerance && __abserr0 < 0.01 * std::abs(__result0))
	return std::make_tuple(__result0, __abserr0);
      else if (__limit == 1)
	__throw__IntegrationError("qaws_integrate: "
				  "a maximum of one iteration was insufficient",
				  MAX_ITER_ERROR, __result0, __abserr0);

      auto __area = __result0;
      auto __errsum = __abserr0;
      auto __iteration = 2u;
      int __error_type = NO_ERROR;
      int __roundoff_type1 = 0, __roundoff_type2 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  _Tp __a_i, __b_i, __r_i, __e_i;
	  __workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	  const auto __a1 = __a_i;
	  const auto __b1 = (__a_i + __b_i) / _Tp{2};
	  const auto __a2 = __b1;
	  const auto __b2 = __b_i;

	  _Tp __area1, __error1;
	  bool __err_reliable1;
	  std::tie(__area1, __error1, __err_reliable1)
	    = qc25s(__table, __func, __a, __b, __a1, __b1);

	  _Tp __area2, __error2;
	  bool __err_reliable2;
	  std::tie(__area2, __error2, __err_reliable2)
	    = qc25s(__table, __func, __a, __b, __a2, __b2);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;

	  __errsum += __error12 - __e_i;
	  __area += __area12 - __r_i;

	  if (__err_reliable1 && __err_reliable2)
	    {
	      const auto __delta = __r_i - __area12;

	      if (std::abs (__delta) <= 1.0e-5 * std::abs(__area12)
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

      const auto __result = __workspace.total_integral();
      const auto __abserr = __errsum;

      if (__iteration == __limit)
	__error_type = MAX_SUBDIV_ERROR;

      if (__errsum <= __tolerance)
	return std::make_tuple(__result, __abserr);

      if (__error_type == NO_ERROR)
	return std::make_tuple(__result, __abserr);

      __check_error<_Tp>(__func__, __error_type);
      __throw__IntegrationError("qaws_integrate: Unknown error.",
				UNKNOWN_ERROR, __result, __abserr);
    }

  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp, bool>
    qc25s(qaws_integration_table<_Tp>& __t,
	  const _FuncTp& __func, _Tp __a, _Tp __b, _Tp __a1, _Tp __b1)
    {
      fn_qaws<_Tp> __fqaws(&__t, __func, __a, __b);

      if (__a1 == __a && (__t.alpha != _Tp{0} || __t.mu != 0))
	{
	  const auto __factor
	    = std::pow(0.5 * (__b1 - __a1), __t.alpha + _Tp{1});

	  auto __f = [__fqaws](_Tp __x)->_Tp{ return __fqaws.eval_right(__x); };
	  std::array<_Tp, 13> __cheb12;
	  std::array<_Tp, 25> __cheb24;
	  qcheb_integrate(__f, __a1, __b1, __cheb12, __cheb24);

	  if (__t.mu == 0)
	    {
	      const auto __u = __factor;

	      _Tp __res12 = 0, __res24 = 0;
	      std::tie(__res12, __res24)
		= compute_result(__t.ri, __cheb12, __cheb24);

	      const auto __result = __u * __res24;
	      const auto __abserr = std::abs(__u * (__res24 - __res12));
	      return std::make_tuple(__result, __abserr, false);
	    }
	  else
	    {
	      const auto __u = __factor * std::log(__b1 - __a1);
	      const auto __v = __factor;

	      _Tp __res12a = 0, __res24a = 0;
	      std::tie(__res12a, __res24a)
		= compute_result(__t.ri, __cheb12, __cheb24);
	      _Tp __res12b = 0, __res24b = 0;
	      std::tie(__res12b, __res24b)
		= compute_result(__t.rg, __cheb12, __cheb24);

	      const auto __result = __u * __res24a + __v * __res24b;
	      const auto __abserr = std::abs(__u * (__res24a - __res12a))
				  + std::abs(__v * (__res24b - __res12b));
	      return std::make_tuple(__result, __abserr, false);
	    }
	}
      else if (__b1 == __b && (__t.beta != _Tp{0} || __t.nu != 0))
	{
	  auto __factor = std::pow(0.5 * (__b1 - __a1), __t.beta + _Tp{1});

	  auto __f = [__fqaws](_Tp __x)->_Tp{ return __fqaws.eval_left(__x); };
	  std::array<_Tp, 13> __cheb12;
	  std::array<_Tp, 25> __cheb24;
	  qcheb_integrate(__f, __a1, __b1, __cheb12, __cheb24);

	  if (__t.nu == 0)
	    {
	      const auto __u = __factor;

	      _Tp __res12, __res24;
	      std::tie(__res12, __res24)
		= compute_result(__t.rj, __cheb12, __cheb24);

	      const auto __result = __u * __res24;
	      const auto __abserr = std::abs(__u * (__res24 - __res12));
	      return std::make_tuple(__result, __abserr, false);
	    }
	  else
	    {
	      const auto __u = __factor * std::log(__b1 - __a1);
	      const auto __v = __factor;

	      _Tp __res12a, __res24a;
	      std::tie(__res12a, __res24a)
		= compute_result(__t.rj, __cheb12, __cheb24);
	      _Tp __res12b, __res24b;
	      std::tie(__res12b, __res24b)
		= compute_result(__t.rh, __cheb12, __cheb24);

	      const auto __result = __u * __res24a + __v * __res24b;
	      const auto __abserr = std::abs(__u * (__res24a - __res12a))
				  + std::abs(__v * (__res24b - __res12b));
	      return std::make_tuple(__result, __abserr, false);
	    }
	}
      else
	{
	  auto __f = [__fqaws](_Tp __x)
		     ->_Tp
		     { return __fqaws.eval_middle(__x); };
	  _Tp __result, __abserr, __resabs, __resasc;
	  std::tie(__result, __abserr, __resabs, __resasc)
	    = qk_integrate(__f, __a1, __b1, QK_15);

	  bool __err_reliable;
	  if (__abserr == __resasc)
	    __err_reliable = false;
	  else
	    __err_reliable = true;

	  return std::make_tuple(__result, __abserr, __err_reliable);
	}
    }

  template<typename _Tp>
    struct fn_qaws
    {
      const qaws_integration_table<_Tp> *table;
      std::function<_Tp(_Tp)> func;
      _Tp a;
      _Tp b;

      fn_qaws(const qaws_integration_table<_Tp> *__t,
	      std::function<_Tp(_Tp)> __f, _Tp __a_in, _Tp __b_in)
      : table(__t),
	func(__f), a(__a_in), b(__b_in)
      { }

      _Tp eval_middle(_Tp) const;
      _Tp eval_left(_Tp) const;
      _Tp eval_right(_Tp) const;
    };

  template<typename _Tp>
    _Tp
    fn_qaws<_Tp>::eval_middle(_Tp __x) const
    {
      auto __factor = _Tp{1};

      if (this->table->alpha != _Tp{0})
	__factor *= std::pow(__x - this->a, this->table->alpha);

      if (table->mu == 1)
	__factor *= std::log(__x - this->a);

      if (this->table->beta != _Tp{0})
	__factor *= std::pow(this->b - __x, this->table->beta);

      if (table->nu == 1)
	__factor *= std::log(this->b - __x);

      return __factor * this->func(__x);
    }

  template<typename _Tp>
    _Tp
    fn_qaws<_Tp>::eval_left(_Tp __x) const
    {
      auto __factor = _Tp{1};

      if (this->table->alpha != _Tp{0})
	__factor *= std::pow(__x - this->a, this->table->alpha);

      if (this->table->mu == 1)
	__factor *= std::log(__x - this->a);

      return __factor * this->func(__x);
    }

  template<typename _Tp>
    _Tp
    fn_qaws<_Tp>::eval_right(_Tp __x) const
    {
      auto __factor = _Tp{1};

      if (this->table->beta != _Tp{0})
	__factor *= std::pow(this->b - __x, this->table->beta);

      if (this->table->nu == 1)
	__factor *= std::log(this->b - __x);

      return __factor * this->func(__x);
    }

  template<typename _Tp>
    std::tuple<_Tp, _Tp>
    compute_result(const std::array<_Tp, 25>& __r,
		   const std::array<_Tp, 13>& __cheb12,
		   const std::array<_Tp, 25>& __cheb24)
    {
      auto __res12 = _Tp{0};
      for (size_t __i = 0; __i < __cheb12.size(); ++__i)
	__res12 += __r[__i] * __cheb12[__i];

      auto __res24 = _Tp{0};
      for (size_t __i = 0; __i < __cheb24.size(); ++__i)
	__res24 += __r[__i] * __cheb24[__i];

      return std::make_tuple(__res12, __res24);
    }

} // namespace __gnu_test

#endif // QAWS_INTEGRATE_H
