/* quadrature/qaws_integrate.h
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
    std::pair<_Tp, _Tp>
    compute_result(const std::array<_Tp, 25>& r,
		   const std::array<_Tp, 13>& cheb12,
		   const std::array<_Tp, 25>& cheb24);

template<typename _FuncTp, typename _Tp>
  std::pair<_Tp, _Tp>
  qaws_integrate(integration_workspace<_Tp>& __workspace,
        	 qaws_integration_table<_Tp>& __table,
        	 const _FuncTp& __func,
        	 const _Tp __a, const _Tp __b,
        	 const _Tp __epsabs, const _Tp __epsrel,
        	 const size_t __limit)
  {
    _Tp __result0, __abserr0;

    if (__limit > __workspace.capacity())
      std::__throw_runtime_error("iteration limit exceeds available workspace") ;
    if (__b <= __a) 
      std::__throw_runtime_error("limits must form an ascending sequence, a < b") ;
    if (__epsabs <= 0 && (__epsrel < 50 * std::numeric_limits<_Tp>::epsilon() || __epsrel < 0.5e-28))
      std::__throw_runtime_error("tolerance cannot be achieved with given epsabs and epsrel");

    // Perform the first integration.
    {
      _Tp __area1, __area2;
      _Tp __error1, __error2;
      bool __err_reliable1, __err_reliable2;
      const auto __a1 = __a;
      const auto __b1 = 0.5 * (__a + __b);
      const auto __a2 = __b1;
      const auto __b2 = __b;

      std::tie(__area1, __error1, __err_reliable1)
	= qc25s(__table, __func, __a, __b, __a1, __b1);
      std::tie(__area2, __error2, __err_reliable2)
	= qc25s(__table, __func, __a, __b, __a2, __b2);

      if (__error1 > __error2)
	{
          __workspace.append(__a1, __b1, __area1, __error1);
          __workspace.append(__a2, __b2, __area2, __error2);
	}
      else
	{
          __workspace.append(__a2, __b2, __area2, __error2);
          __workspace.append(__a1, __b1, __area1, __error1);
	}

      __result0 = __area1 + __area2;
      __abserr0 = __error1 + __error2;
      __workspace.set_initial(__a, __b, __result0, __abserr0);
    }

    /* Test on accuracy */

    auto __tolerance = std::max(__epsabs, __epsrel * std::abs(__result0));

    /* Test on accuracy, use 0.01 relative error as an extra safety
       margin on the first iteration (ignored for subsequent iterations) */

    if (__abserr0 < __tolerance && __abserr0 < 0.01 * std::abs(__result0))
      return std::make_pair(__result0, __abserr0);
    else if (__limit == 1)
      std::__throw_runtime_error("a maximum of one iteration was insufficient");

    auto __area = __result0;
    auto __errsum = __abserr0;
    auto __iteration = 2u;
    int __error_type = 0;
    do
      {
	/* Bisect the subinterval with the largest error estimate */

	_Tp __a_i, __b_i, __r_i, __e_i;
	__workspace.retrieve(__a_i, __b_i, __r_i, __e_i);

	const auto __a1 = __a_i; 
	const auto __b1 = 0.5 * (__a_i + __b_i);
	const auto __a2 = __b1;
	const auto __b2 = __b_i;

	_Tp __area1, __error1;
	_Tp __area2, __error2;
	bool __err_reliable1, __err_reliable2;
	std::tie(__area1, __error1, __err_reliable1)
	  = qc25s(__table, __func, __a, __b, __a1, __b1);
	std::tie(__area2, __error2, __err_reliable2)
	  = qc25s(__table, __func, __a, __b, __a2, __b2);

	const auto __area12 = __area1 + __area2;
	const auto __error12 = __error1 + __error2;

	__errsum += __error12 - __e_i;
	__area += __area12 - __r_i;

	int __roundoff_type1 = 0, __roundoff_type2 = 0;
	if (__err_reliable1 && __err_reliable2)
          {
            auto __delta = __r_i - __area12;

            if (std::abs (__delta) <= 1.0e-5 * std::abs(__area12) && __error12 >= 0.99 * __e_i)
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

    const auto __result = __workspace.sum_results();
    const auto __abserr = __errsum;

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

  template<typename _FuncTp, typename _Tp>
    std::tuple<_Tp, _Tp, bool>
    qc25s(qaws_integration_table<_Tp>& t,
	  const _FuncTp& __func, _Tp a, _Tp b, _Tp a1, _Tp b1)
    {
      struct fn_qaws<_Tp> fqaws(&t, __func, a, b);

      if (a1 == a && (t.alpha != _Tp{0} || t.mu != 0))
	{
	  const auto factor = std::pow(0.5 * (b1 - a1), t.alpha + _Tp{1});

	  auto f = [fqaws](_Tp x)->_Tp{ return fqaws.eval_right(x); };
	  std::array<_Tp, 13> cheb12;
	  std::array<_Tp, 25> cheb24;
	  qcheb_integrate(f, a1, b1, cheb12, cheb24);

	  if (t.mu == 0)
            {
              const auto u = factor;

              _Tp res12 = 0, res24 = 0;
              std::tie(res12, res24) = compute_result(t.ri, cheb12, cheb24);

              const auto result = u * res24;
              const auto abserr = std::abs(u * (res24 - res12));
	      return std::make_tuple(result, abserr, false);
            }
	  else 
            {
              const auto u = factor * std::log(b1 - a1);
              const auto v = factor;

              _Tp res12a = 0, res24a = 0;
              std::tie(res12a, res24a) = compute_result(t.ri, cheb12, cheb24);
              _Tp res12b = 0, res24b = 0;
              std::tie(res12b, res24b) = compute_result(t.rg, cheb12, cheb24);

              const auto result = u * res24a + v * res24b;
              const auto abserr = std::abs(u * (res24a - res12a)) + std::abs(v * (res24b - res12b));
	      return std::make_tuple(result, abserr, false);
            }
	}
      else if (b1 == b && (t.beta != _Tp{0} || t.nu != 0))
	{
	  _Tp factor = std::pow(0.5 * (b1 - a1), t.beta + _Tp{1});

	  auto f = [fqaws](_Tp x)->_Tp{ return fqaws.eval_left(x); };
	  std::array<_Tp, 13> cheb12;
	  std::array<_Tp, 25> cheb24;
	  qcheb_integrate(f, a1, b1, cheb12, cheb24);

	  if (t.nu == 0)
            {
              const auto u = factor;

              _Tp res12, res24;
              std::tie(res12, res24) = compute_result(t.rj, cheb12, cheb24);

              const auto result = u * res24;
              const auto abserr = std::abs(u * (res24 - res12));
	      return std::make_tuple(result, abserr, false);
            }
	  else 
            {
              const auto u = factor * std::log(b1 - a1);
              const auto v = factor;

              _Tp res12a, res24a;
              std::tie(res12a, res24a) = compute_result(t.rj, cheb12, cheb24);
              _Tp res12b, res24b;
              std::tie(res12b, res24b) = compute_result(t.rh, cheb12, cheb24);

              const auto result = u * res24a + v * res24b;
              const auto abserr = std::abs(u * (res24a - res12a))
				+ std::abs(v * (res24b - res12b));
	      return std::make_tuple(result, abserr, false);
            }
	}
      else
	{
	  auto f = [fqaws](_Tp x)->_Tp{ return fqaws.eval_middle(x); };
	  _Tp result, abserr, resabs, resasc;
	  std::tie(result, abserr, resabs, resasc)
	    = qk_integrate(f, a1, b1, QK_15);

	  bool err_reliable;
	  if (abserr == resasc)
            err_reliable = false;
	  else
            err_reliable = true;

	  return std::make_tuple(result, abserr, err_reliable);
	}
    }

  template<typename _Tp>
    struct fn_qaws
    {
      const qaws_integration_table<_Tp> *table;
      std::function<_Tp(_Tp)> func;
      _Tp a;
      _Tp b;

      fn_qaws(const qaws_integration_table<_Tp> *t,
	      std::function<_Tp(_Tp)> f, _Tp a_in, _Tp b_in)
      : table(t),
        func(f), a(a_in), b(b_in)
      { }

      _Tp eval_middle(_Tp) const;
      _Tp eval_left(_Tp) const;
      _Tp eval_right(_Tp) const;
    };

  template<typename _Tp>
    _Tp
    fn_qaws<_Tp>::eval_middle(_Tp x) const
    {
      auto factor = _Tp{1};

      if (this->table->alpha != _Tp{0})
	factor *= std::pow(x - this->a, this->table->alpha);

      if (table->mu == 1)
	factor *= std::log(x - this->a);

      if (this->table->beta != _Tp{0})
	factor *= std::pow(this->b - x, this->table->beta);

      if (table->nu == 1)
	factor *= std::log(this->b - x);

      return factor * this->func(x);
    }

  template<typename _Tp>
    _Tp
    fn_qaws<_Tp>::eval_left(_Tp x) const
    {
      auto factor = _Tp{1};

      if (this->table->alpha != _Tp{0})
	factor *= std::pow(x - this->a, this->table->alpha);

      if (this->table->mu == 1)
	factor *= std::log(x - this->a);

      return factor * this->func(x);
    }

  template<typename _Tp>
    _Tp
    fn_qaws<_Tp>::eval_right(_Tp x) const
    {
      auto factor = _Tp{1};

      if (this->table->beta != _Tp{0})
	factor *= std::pow(this->b - x, this->table->beta);

      if (this->table->nu == 1)
	factor *= std::log(this->b - x);

      return factor * this->func(x);
    }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    compute_result(const std::array<_Tp, 25>& r,
		   const std::array<_Tp, 13>& cheb12,
		   const std::array<_Tp, 25>& cheb24)
    {  
      auto res12 = _Tp{0};
      for (size_t i = 0; i < cheb12.size(); ++i)
	res12 += r[i] * cheb12[i];

      auto res24 = _Tp{0};
      for (size_t i = 0; i < cheb24.size(); ++i)
	res24 += r[i] * cheb24[i];

      return std::make_pair(res12, res24);
    }

} // namespace __gnu_test

#endif // QAWS_INTEGRATE_H
