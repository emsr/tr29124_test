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

#include <tuple>

#include "qk_integrate.h"
#include "integration_workspace.h"
#include "qaws_integration_table.h"

namespace __gnu_test
{

template<typename _FuncTp, typename _Tp>
  std::pair<_Tp, _Tp>
  qaws_integrate(integration_workspace<_Tp>& workspace,
        	 const _FuncTp& __func,
        	 const _Tp a, const _Tp b,
        	 qaws_integration_table<_Tp>& t,
        	 const _Tp epsabs, const _Tp epsrel,
        	 const size_t limit)
  {
    _Tp area, errsum;
    _Tp result0, abserr0;
    _Tp tolerance;
    size_t iteration = 0;
    int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

    /* Initialize results */

    workspace.initialise(a, b);

    _Tp result = _Tp{};
    _Tp abserr = _Tp{};

    if (limit > workspace->limit)
      std::__throw_runtime_error ("iteration limit exceeds available workspace") ;
    if (b <= a) 
      std::__throw_runtime_error ("limits must form an ascending sequence, a < b") ;
    if (epsabs <= 0 && (epsrel < 50 * std::numeric_limits<_Tp>::epsilon() || epsrel < 0.5e-28))
      std::__throw_runtime_error ("tolerance cannot be achieved with given epsabs and epsrel");

    /* perform the first integration */

    {
      _Tp area1, area2;
      _Tp error1, error2;
      bool err_reliable1, err_reliable2;
      auto a1 = a;
      auto b1 = 0.5 * (a + b);
      auto a2 = b1;
      auto b2 = b;

      std::tie(area1, error1, err_reliable1) = qc25s(__func, a, b, a1, b1, t);
      std::tie(area2, error2, err_reliable2) = qc25s(__func, a, b, a2, b2, t);

      if (error1 > error2)
	{
          workspace.append_interval(a1, b1, area1, error1);
          workspace.append_interval(a2, b2, area2, error2);
	}
      else
	{
          workspace.append_interval(a2, b2, area2, error2);
          workspace.append_interval(a1, b1, area1, error1);
	}

      result0 = area1 + area2;
      abserr0 = error1 + error2;
    }

    /* Test on accuracy */

    tolerance = std::max(epsabs, epsrel * std::abs (result0));

    /* Test on accuracy, use 0.01 relative error as an extra safety
       margin on the first iteration (ignored for subsequent iterations) */

    if (abserr0 < tolerance && abserr0 < 0.01 * std::abs(result0))
      std::make_pair(result0, abserr0);
    else if (limit == 1)
      std::__throw_runtime_error("a maximum of one iteration was insufficient");

    area = result0;
    errsum = abserr0;

    iteration = 2;

    do
      {
	_Tp area1 = 0, area2 = 0;
	_Tp error1 = 0, error2 = 0;
	int err_reliable1, err_reliable2;

	/* Bisect the subinterval with the largest error estimate */

	_Tp a_i, b_i, r_i, e_i;
	workspace.retrieve(a_i, b_i, r_i, e_i);

	auto a1 = a_i; 
	auto b1 = 0.5 * (a_i + b_i);
	auto a2 = b1;
	auto b2 = b_i;

	std::tie(area1, error1, err_reliable1) = qc25s(__func, a, b, a1, b1, t);
	std::tie(area2, error2, err_reliable2) = qc25s(__func, a, b, a2, b2, t);

	auto area12 = area1 + area2;
	auto error12 = error1 + error2;

	errsum += (error12 - e_i);
	area += area12 - r_i;

	if (err_reliable1 && err_reliable2)
          {
            auto delta = r_i - area12;

            if (std::abs (delta) <= 1.0e-5 * std::abs (area12) && error12 >= 0.99 * e_i)
              ++roundoff_type1;
            if (iteration >= 10 && error12 > e_i)
              ++roundoff_type2;
          }

	tolerance = std::max (epsabs, epsrel * std::abs (area));

	if (errsum > tolerance)
          {
            if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
              error_type = 2;   /* round off error */

            /* set error flag in the case of bad integrand behaviour at
               a point of the integration range */

            if (integration_workspace<_Tp>::subinterval_too_small(a1, a2, b2))
              error_type = 3;
          }

	workspace.update(a1, b1, area1, error1,
			 a2, b2, area2, error2);

	workspace.retrieve(a_i, b_i, r_i, e_i);

	++iteration;

      }
    while (iteration < limit && !error_type && errsum > tolerance);

    result = sum_results (workspace);
    abserr = errsum;

    if (errsum <= tolerance)
      std::make_pair(result, abserr);
    else if (error_type == 2)
      std::__throw_runtime_error ("roundoff error prevents tolerance from being achieved");
    else if (error_type == 3)
      std::__throw_runtime_error ("bad integrand behavior found in the integration interval");
    else if (iteration == limit)
      std::__throw_runtime_error ("maximum number of subdivisions reached");
    else
      std::__throw_runtime_error ("could not integrate function");

  }

} // namespace __gnu_test

#endif // QAWS_INTEGRATE_H
