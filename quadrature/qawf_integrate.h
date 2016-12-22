/* integration/qawf.c
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

#ifndef QAWF_INTEGRATE_H
#define QAWF_INTEGRATE_H 1

#include <cmath>

#include "integration_workspace.h"
#include "oscillatory_integration_table.h"

namespace __gnu_test
{

  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp>
    qawf_integrate(integration_workspace<_Tp>& workspace,
                   integration_workspace<_Tp>& cycle_workspace,
                   oscillatory_integration_table<_Tp>& wf,
		   const _FuncTp& f,
                   const _Tp a,
                   const _Tp epsabs,
                   const size_t limit)
    {
      _Tp area, errsum;
      _Tp res_ext, err_ext;
      _Tp correc, total_error = 0.0, truncation_error;

      size_t ktmin = 0;
      size_t iteration = 0;

      extrapolation_table<_Tp> table;

      _Tp cycle;
      auto omega = wf.omega;

      const _Tp p = 0.9;
      _Tp factor = 1;
      _Tp initial_eps, eps;
      int error_type = 0;

      workspace.set_initial_limits(a, a);

      int status = 0;
      auto result = _Tp{0};
      auto abserr = _Tp{0};

      if (limit > workspace.capacity())
	std::__throw_runtime_error("iteration limit exceeds available workspace") ;

      /* Test on accuracy */

      if (epsabs <= 0)
	std::__throw_runtime_error("absolute tolerance epsabs must be positive") ;

      if (omega == 0.0)
	{
	  if (wf.circfun == oscillatory_integration_table<_Tp>::INTEG_SINE)
            // The function sin(w x) f(x) is always zero for w = 0.
            return std::make_tuple(_Tp{0}, _Tp{0});
	  else
            // The function cos(w x) f(x) is always f(x) for w = 0.
	    return qagiu_integrate(cycle_workspace, f, a, epsabs, 0.0,
                                   cycle_workspace.capacity());
	}

      if (epsabs > std::numeric_limits<_Tp>::min() / (1 - p)) // seems small
	eps = epsabs * (1 - p);
      else
	eps = epsabs;

      initial_eps = eps;

      area = 0;
      errsum = 0;

      res_ext = 0;
      err_ext = std::numeric_limits<_Tp>::max();
      correc = 0;

      cycle = (2 * floor (std::abs(omega)) + 1) * M_PI / std::abs(omega);

      wf.set_length(cycle);

      for (iteration = 0; iteration < limit; ++iteration)
	{
	  _Tp area1, error1, reseps, erreps;

	  auto a1 = a + iteration * cycle;
	  auto b1 = a1 + cycle;

	  auto epsabs1 = eps * factor;

	  std::tie(area1, error1) = qawo_integrate(cycle_workspace, wf, f, a1, epsabs1, 0.0, limit);

	  workspace.append(a1, b1, area1, error1);

	  factor *= p;

	  area += area1;
	  errsum += error1;

	  // estimate the truncation error as 50 times the final term.

	  truncation_error = 50 * std::abs(area1);

	  total_error = errsum + truncation_error;

	  if (total_error < epsabs && iteration > 4)
            goto compute_result;

	  if (error1 > correc)
            correc = error1;

	  if (status)
            eps = std::max(initial_eps, correc * (1.0 - p));

	  if (status && total_error < 10 * correc && iteration > 3)
            goto compute_result;

	  table.append(area);

	  if (table.get_nn() < 2)
            continue;

	  std::tie(reseps, erreps) = table.qelg();

	  ++ktmin;

	  if (ktmin >= 15 && err_ext < 0.001 * total_error)
            error_type = 4;

	  if (erreps < err_ext)
            {
              ktmin = 0;
              err_ext = erreps;
              res_ext = reseps;

              if (err_ext + 10 * correc <= epsabs)
        	break;
              if (err_ext <= epsabs && 10 * correc >= epsabs)
        	break;
            }
	}

      if (iteration == limit)
	error_type = 1;

      if (err_ext == std::numeric_limits<_Tp>::max())
	goto compute_result;

      err_ext = err_ext + 10 * correc;

      result = res_ext;
      abserr = err_ext;

      if (error_type == 0)
	return std::make_tuple(result, abserr);

      if (res_ext != _Tp{0} && area != _Tp{0})
	{
	  if (err_ext / std::abs(res_ext) > errsum / std::abs(area))
	    goto compute_result;
	}
      else if (err_ext > errsum)
	goto compute_result;
      else if (area == _Tp{0})
	goto return_error;

      if (error_type == 4)
	err_ext += truncation_error;

      goto return_error;

    compute_result:

      result = area;
      abserr = total_error;

    return_error:

      if (error_type > 2)
	--error_type;

      if (error_type == 0)
	return std::make_tuple(result, abserr);
      else if (error_type == 1)
	std::__throw_runtime_error ("number of iterations was insufficient");
      else if (error_type == 2)
	std::__throw_runtime_error ("cannot reach tolerance because of roundoff error");
      else if (error_type == 3)
	std::__throw_runtime_error ("bad integrand behavior found in the integration interval");
      else if (error_type == 4)
	std::__throw_runtime_error ("roundoff error detected in the extrapolation table");
      else if (error_type == 5)
	std::__throw_runtime_error ("integral is divergent, or slowly convergent");
      else
	std::__throw_runtime_error ("could not integrate function");

    }

} // namespace __gn_test

#endif // QAWF_INTEGRATE_H
