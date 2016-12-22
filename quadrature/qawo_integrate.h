/* integration/qawo.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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

#ifndef QAWO_INTEGRATE_H
#define QAWO_INTEGRATE_H 1

#include <cmath>

//#include "initialise.c"
//#include "set_initial.c"
//#include "reset.c"
//#include "qpsrt.c"
//#include "util.c"
//#include "qpsrt2.c"
//#include "qelg.c"
//#include "positivity.c"

//#include "qc25f.c"

#include "oscillatory_integration_table.h"

namespace __gnu_test
{

  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp, _Tp, _Tp>
    qc25f(oscillatory_integration_table<_Tp>& wf,
	  const _FuncTp& __func, _Tp a, _Tp b, 
	  size_t level);

  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp>
    qawo_integrate(integration_workspace<_Tp>& workspace,
        	   oscillatory_integration_table<_Tp>& wf,
		   const _FuncTp& __func,
        	   const _Tp a,
        	   const _Tp epsabs, const _Tp epsrel,
        	   const size_t limit)
    {
      _Tp area, errsum;
      _Tp res_ext, err_ext;
      _Tp result0, abserr0, resabs0, resasc0;
      _Tp tolerance;

      _Tp ertest = 0;
      _Tp error_over_large_intervals = 0;
      _Tp reseps = 0, abseps = 0, correc = 0;
      size_t ktmin = 0;
      int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0;
      int error_type = 0, error_type2 = 0;

      size_t iteration = 0;

      int positive_integrand = 0;
      int extrapolate = 0;
      int extall = 0;
      int disallow_extrapolation = 0;

      extrapolation_table<_Tp> table;

      auto b = a + wf.get_length();
      auto abs_omega = std::abs(wf.omega);

      workspace.set_initial_limits(a, b);

      auto result = _Tp{0};
      auto abserr = _Tp{0};

      if (limit > workspace.capacity())
	std::__throw_runtime_error("iteration limit exceeds available workspace");

      /* Test on accuracy */

      if (epsabs <= 0 && (epsrel < 50 * std::numeric_limits<_Tp>::epsilon() || epsrel < 0.5e-28))
	std::__throw_runtime_error("tolerance cannot be achieved with given epsabs and epsrel");

      // Perform the first integration.

      std::tie(result0, abserr0, resabs0, resasc0)
	= qc25f(wf, __func, a, b, 0);

      workspace.set_initial_results(result0, abserr0);

      tolerance = std::max(epsabs, epsrel * std::abs(result0));

      if (abserr0 <= 100 * std::numeric_limits<_Tp>::epsilon() * resabs0 && abserr0 > tolerance)
	{
	  result = result0;
	  abserr = abserr0;

	  std::__throw_runtime_error("cannot reach tolerance because of roundoff error"
                     "on first attempt");
	}
      else if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0)
	{
	  result = result0;
	  abserr = abserr0;

	  return std::make_tuple(result, abserr);
	}
      else if (limit == 1)
	{
	  result = result0;
	  abserr = abserr0;

	  std::__throw_runtime_error("a maximum of one iteration was insufficient");
	}

      //initialise_table(&table);

      if (0.5 * abs_omega * std::abs(b - a) <= 2)
	{
	  table.append(result0);
	  extall = 1;
	}

      area = result0;
      errsum = abserr0;

      res_ext = result0;
      err_ext = std::numeric_limits<_Tp>::max();

      positive_integrand = __test_positivity(result0, resabs0);

      iteration = 1;

      do
	{
	  size_t current_level;
	  _Tp a1, b1, a2, b2;
	  _Tp a_i, b_i, r_i, e_i;
	  _Tp area1 = 0, area2 = 0, area12 = 0;
	  _Tp error1 = 0, error2 = 0, error12 = 0;
	  _Tp resasc1, resasc2;
	  _Tp resabs1, resabs2;
	  _Tp last_e_i;

	  /* Bisect the subinterval with the largest error estimate */

	  workspace.retrieve(a_i, b_i, r_i, e_i);

	  current_level = workspace.current_level() + 1;

	  if (current_level >= wf.n) 
            {
              error_type = -1; /* exceeded limit of table */
              break;
            }

	  a1 = a_i;
	  b1 = 0.5 * (a_i + b_i);
	  a2 = b1;
	  b2 = b_i;

	  iteration++;

	  std::tie(area1, error1, resabs1, resasc1)
	    = qc25f(wf, __func, a1, b1, current_level);

	  std::tie(area2, error2, resabs2, resasc2)
	    = qc25f(wf, __func, a2, b2, current_level);

	  area12 = area1 + area2;
	  error12 = error1 + error2;
	  last_e_i = e_i;

	  /* Improve previous approximations to the integral and test for
             accuracy.

             We write these expressions in the same way as the original
             QUADPACK code so that the rounding errors are the same, which
             makes testing easier. */

	  errsum = errsum + error12 - e_i;
	  area = area + area12 - r_i;

	  tolerance = std::max(epsabs, epsrel * std::abs(area));

	  if (resasc1 != error1 && resasc2 != error2)
            {
              _Tp delta = r_i - area12;

              if (std::abs(delta) <= 1.0e-5 * std::abs(area12) && error12 >= 0.99 * e_i)
        	{
        	  if (!extrapolate)
                    ++roundoff_type1;
        	  else
                    ++roundoff_type2;
        	}
              if (iteration > 10 && error12 > e_i)
        	++roundoff_type3;
            }

	  /* Test for roundoff and eventually set error flag */

	  if (roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20)
            error_type = 2; // round off error

	  if (roundoff_type2 >= 5)
            error_type2 = 1;

	  /* set error flag in the case of bad integrand behaviour at
             a point of the integration range */

	  if (integration_workspace<_Tp>::subinterval_too_small (a1, a2, b2))
            error_type = 4;

	  /* append the newly-created intervals to the list */

	  workspace.update(a1, b1, area1, error1, a2, b2, area2, error2);

	  if (errsum <= tolerance)
            goto compute_result;

	  if (error_type)
            break;

	  if (iteration >= limit - 1)
            {
              error_type = 1;
              break;
            }

	  /* set up variables on first iteration */

	  if (iteration == 2 && extall)     
            {
              error_over_large_intervals = errsum;
              ertest = tolerance;
              table.append(area);
              continue;
            }

	  if (disallow_extrapolation)
            continue;

	  if (extall)
            {
              error_over_large_intervals += -last_e_i;

              if (current_level < workspace.max_level())
        	error_over_large_intervals += error12;

              if (extrapolate)
        	goto label70;
            }

	  if (workspace.large_interval())
            continue;

	  if (extall)
            {
              extrapolate = 1;
              workspace.set_nrmax(1);
            }
	  else
            {
              /* test whether the interval to be bisected next is the
        	 smallest interval. */
              size_t i = workspace.current_index(); // ???
              auto width = workspace.upper_lim(i) - workspace.lower_lim(i);

              if (0.25 * std::abs(width) * abs_omega > _Tp{2})
        	continue;

              extall = 1;
              error_over_large_intervals = errsum;
              ertest = tolerance;
              continue;
            }

	label70:
	  if (!error_type2 && error_over_large_intervals > ertest)
            {
              if (workspace.increase_nrmax())
        	continue;
            }

	  /* Perform extrapolation */

	  table.append(area);

	  if (table.get_nn() < 3)
            {
              workspace.reset_nrmax();
              extrapolate = 0;
              error_over_large_intervals = errsum;
              continue;
            }

	  std::tie(reseps, abseps) = table.qelg();

	  ++ktmin;

	  if (ktmin > 5 && err_ext < 0.001 * errsum)
            {
              error_type = 5;
            }

	  if (abseps < err_ext)
            {
              ktmin = 0;
              err_ext = abseps;
              res_ext = reseps;
              correc = error_over_large_intervals;
              ertest = std::max(epsabs, epsrel * std::abs(reseps));
              if (err_ext <= ertest)
        	break;
            }

	  /* Prepare bisection of the smallest interval. */

	  if (table.get_nn() == 1)
            disallow_extrapolation = 1;

	  if (error_type == 5)
            break;

	  /* work on interval with largest error */

	  workspace.reset_nrmax();
	  extrapolate = 0;
	  error_over_large_intervals = errsum;

	}
      while (iteration < limit);

      result = res_ext;
      abserr = err_ext;

      if (err_ext == std::numeric_limits<_Tp>::max())
	goto compute_result;

      if (error_type || error_type2)
	{
	  if (error_type2)
            err_ext += correc;

	  if (error_type == 0)
            error_type = 3;

	  if (result != 0 && area != 0)
            {
              if (err_ext / std::abs(res_ext) > errsum / std::abs(area))
        	goto compute_result;
            }
	  else if (err_ext > errsum)
            {
              goto compute_result;
            }
	  else if (area == 0.0)
            {
              goto return_error;
            }
	}

      /*  Test on divergence. */

      {
	_Tp max_area = std::max(std::abs(res_ext), std::abs(area));

	if (!positive_integrand && max_area < 0.01 * resabs0)
	  goto return_error;
      }

      {
	_Tp ratio = res_ext / area;

	if (ratio < 0.01 || ratio > 100 || errsum > std::abs(area))
	  error_type = 6;
      }

      goto return_error;

    compute_result:

      result = workspace.sum_results();
      abserr = errsum;

    return_error:

      if (error_type > 2)
	--error_type;

      if (error_type == 0)
       return std::make_pair(result, abserr);
      else if (error_type == 1)
	std::__throw_runtime_error("number of iterations was insufficient");
      else if (error_type == 2)
	std::__throw_runtime_error("cannot reach tolerance because of roundoff error");
      else if (error_type == 3)
	std::__throw_runtime_error("bad integrand behavior found in the integration interval");
      else if (error_type == 4)
	std::__throw_runtime_error("roundoff error detected in extrapolation table");
      else if (error_type == 5)
	std::__throw_runtime_error("integral is divergent, or slowly convergent");
      else if (error_type == -1) 
	std::__throw_runtime_error("exceeded limit of trigonometric table");
      else
	std::__throw_runtime_error("could not integrate function");

    }


  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp, _Tp, _Tp>
    qc25f(oscillatory_integration_table<_Tp>& wf,
	  const _FuncTp& __func, _Tp a, _Tp b, 
	  size_t level)
    {
      const auto center = 0.5 * (a + b);
      const auto half_length = 0.5 * (b - a);
      const auto omega = wf.omega ;

      const auto par = omega * half_length;

      if (std::abs(par) < _Tp{2})
	{
	  if (wf.circfun == oscillatory_integration_table<_Tp>::INTEG_SINE)
	    {
	      auto wfun = [__func, omega](_Tp x)
			  ->_Tp
			  { return std::sin(omega * x) * __func(x); };
	      return qk_integrate(wfun, a, b, QK_15);
	    }
	  else
	    {
	      auto wfun = [__func, omega](_Tp x)
			  ->_Tp
			  { return std::cos(omega * x) * __func(x); };
	      return qk_integrate(wfun, a, b, QK_15);
	    }
	}
      else
	{
	  std::array<_Tp, 13> cheb12;
	  std::array<_Tp, 25> cheb24;

	  qcheb_integrate(__func, a, b, cheb12, cheb24);

	  if (level >= wf.n)
            {
              // table overflow should not happen, check before calling.
              std::__throw_runtime_error("table overflow in internal function");
            }

	  // obtain moments from the table..

	  const auto& moment = wf.get_moments(level);

	  auto res12_cos = cheb12[12] * moment[12];
	  auto res12_sin = _Tp{0};
	  for (int i = 0; i < 6; ++i)
            {
              size_t k = 10 - 2 * i;
              res12_cos += cheb12[k] * moment[k];
              res12_sin += cheb12[k + 1] * moment[k + 1];
            }

	  auto res24_cos = cheb24[24] * moment[24];
	  auto res24_sin = _Tp{0};
	  auto result_abs = std::abs(cheb24[24]);
	  for (int i = 0; i < 12; ++i)
            {
              size_t k = 22 - 2 * i;
              res24_cos += cheb24[k] * moment[k];
              res24_sin += cheb24[k + 1] * moment[k + 1];
              result_abs += std::abs(cheb24[k]) + std::abs(cheb24[k+1]);
            }

	  const auto est_cos = std::abs(res24_cos - res12_cos);
	  const auto est_sin = std::abs(res24_sin - res12_sin);

	  const auto c = half_length * std::cos(center * omega);
	  const auto s = half_length * std::sin(center * omega);

	  _Tp result, abserr, resabs, resasc;
	  if (wf.circfun == oscillatory_integration_table<_Tp>::INTEG_SINE)
            {
              result = c * res24_sin + s * res24_cos;
              abserr = std::abs(c * est_sin) + std::abs(s * est_cos);
            }
	  else
            {
              result = c * res24_cos - s * res24_sin;
              abserr = std::abs(c * est_cos) + std::abs(s * est_sin);
            }

	  resabs = result_abs * half_length;
	  resasc = std::numeric_limits<_Tp>::max();

	  return std::make_tuple(result, abserr, resabs, resasc);
	}
    }

} // namespace __gnu_test

#endif // QAWO_INTEGRATE_H

