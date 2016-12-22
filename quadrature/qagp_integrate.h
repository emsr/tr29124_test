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

#ifndef QAGP_INTEGRATE_H
#define QAGP_INTEGRATE_H 1

#include "qk_integrate.h"
#include "integration_workspace.h"
#include "extrapolation_table.h"

namespace __gnu_test
{

template<typename _FuncTp, typename _Tp>
  std::pair<_Tp, _Tp>
  qagp_integrate(integration_workspace<_Tp>& workspace,
                 const _FuncTp& __func,
                 std::vector<_Tp> pts,
                 _Tp epsabs, _Tp epsrel, size_t limit)
  {
    return qagp_integrate(workspace, __func, pts,
			  epsabs, epsrel, limit, QK_21);
  }

template<typename _FuncTp, typename _Tp>
  std::pair<_Tp, _Tp>
  qagp_integrate(integration_workspace<_Tp>& workspace,
                 const _FuncTp& __func,
                 std::vector<_Tp> pts,
                 _Tp epsabs, _Tp epsrel, size_t limit,
		 const qk_intrule __qk_rule)
  {
    using qk_return = std::tuple<_Tp&, _Tp&, _Tp&, _Tp&>;

    _Tp area, errsum;
    _Tp res_ext, err_ext;
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
    int disallow_extrapolation = 0;

    const size_t nint = pts.size() - 1; /* number of intervals */

    size_t i;

    /* Initialize results */
    auto result = _Tp{0};
    auto abserr = _Tp{0};

    /* Test on validity of parameters */
    if (limit > workspace.capacity())
      std::__throw_runtime_error("qagp_integrate: "
				 "iteration limit exceeds available workspace");
    if (pts.size() > workspace.capacity())
      std::__throw_runtime_error("qagp_integrate: "
				 "number of pts exceeds size of workspace");
    if (epsabs <= 0 && (epsrel < 50 * std::numeric_limits<_Tp>::epsilon()
			 || epsrel < 0.5e-28))
      std::__throw_runtime_error("qagp_integrate: "
				 "tolerance cannot be achieved "
				 "with given epsabs and epsrel");

    /* Check that the integration range and break points are an
       ascending sequence */

    for (i = 0; i < nint; ++i)
      if (pts[i + 1] < pts[i])
        std::__throw_runtime_error("qagp_integrate: "
				   "points are not in an ascending sequence");

    /* Perform the first integration */

    auto result0 = _Tp{0};
    auto abserr0 = _Tp{0};
    auto resabs0 = _Tp{0};

    //workspace.set_initial(_Tp{0}, _Tp{0}, _Tp{0}, _Tp{0});

    for (i = 0; i < nint; ++i)
      {
	_Tp area1, error1, resabs1, resasc1;
	const _Tp a1 = pts[i];
	const _Tp b1 = pts[i + 1];

	qk_return{area1, error1, resabs1, resasc1}
	  = qk_integrate(__func, a1, b1, __qk_rule);

	result0 += area1;
	abserr0 += error1;
	resabs0 += resabs1;

	workspace.append(a1, b1, area1, error1);

	if (error1 == resasc1 && error1 != 0.0)
          workspace.set_level(i, 1);
	else
          workspace.set_level(i, 0);
      }

    // Compute the initial error estimate.
    errsum = 0.0;
    for (i = 0; i < nint; ++i)
      {
	if (workspace.level(i))
          workspace.set_abs_error(i, abserr0);
	errsum += workspace.abs_error(i);
      }

    for (i = 0; i < nint; ++i)
      workspace.set_level(i, 0);

    // Sort results into order of decreasing error via the indirection
    // array order[]

    workspace.sort_error();

    // Test on accuracy.
    tolerance = std::max(epsabs, epsrel * std::abs(result0));

    if (abserr0 <= 100 * std::numeric_limits<_Tp>::epsilon() * resabs0
	 && abserr0 > tolerance)
      std::__throw_runtime_error("qagp_integrate: "
				 "cannot reach tolerance because "
				 "of roundoff error on first attempt");
    else if (abserr0 <= tolerance)
      return std::make_pair(result0, abserr0);
    else if (limit == 1)
      std::__throw_runtime_error("qagp_integrate: "
				 "a maximum of one iteration was insufficient");

    /* Initialization */

    extrapolation_table<_Tp> table(result0);

    area = result0;

    res_ext = result0;
    err_ext = std::numeric_limits<_Tp>::max();

    error_over_large_intervals = errsum;
    ertest = tolerance;

    positive_integrand = __test_positivity(result0, resabs0);

    iteration = nint - 1; 

    do
      {
	// Bisect the subinterval with the largest error estimate.

	_Tp a_i, b_i, r_i, e_i;
	workspace.retrieve(a_i, b_i, r_i, e_i);

	size_t current_level = workspace.current_level() + 1;

	auto a1 = a_i;
	auto b1 = 0.5 * (a_i + b_i);
	auto a2 = b1;
	auto b2 = b_i;

	++iteration;

	_Tp area1, error1, resabs1, resasc1;
	qk_return{area1, error1, resabs1, resasc1}
	  = qk_integrate(__func, a1, b1, __qk_rule);

	_Tp area2, error2, resabs2, resasc2;
	qk_return{area2, error2, resabs2, resasc2}
	  = qk_integrate(__func, a2, b2, __qk_rule);

	auto area12 = area1 + area2;
	auto error12 = error1 + error2;
	auto last_e_i = e_i;

	/* Improve previous approximations to the integral and test for
           accuracy.

           We write these expressions in the same way as the original
           QUADPACK code so that the rounding errors are the same, which
           makes testing easier. */

	errsum += error12 - e_i;
	area += area12 - r_i;

	tolerance = std::max(epsabs, epsrel * std::abs (area));

	if (resasc1 != error1 && resasc2 != error2)
          {
            _Tp delta = r_i - area12;

            if (std::abs(delta) <= 1.0e-5 * std::abs(area12)
		 && error12 >= 0.99 * e_i)
              {
        	if (!extrapolate)
                  ++roundoff_type1;
        	else
                  ++roundoff_type2;
              }

            if (i > 10 && error12 > e_i)
              ++roundoff_type3;
          }

	// Test for roundoff and eventually set error flag.

	if (roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20)
          error_type = 2; // round off error

	if (roundoff_type2 >= 5)
          error_type2 = 1;

	// Set error flag in the case of bad integrand behaviour at
        // a point of the integration range

	if (integration_workspace<_Tp>::subinterval_too_small(a1, a2, b2))
          error_type = 4;

	// Append the newly-created intervals to the list.
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

	if (disallow_extrapolation)
          continue;

	error_over_large_intervals += -last_e_i;

	if (current_level < workspace.max_level())
          error_over_large_intervals += error12;

	if (!extrapolate)
          {
            // Test whether the interval to be bisected next is the
            // smallest interval.
            if (workspace.large_interval())
              continue;

            extrapolate = 1;
            workspace.set_nrmax(1);
          }

	/* The smallest interval has the largest error.  Before
           bisecting decrease the sum of the errors over the larger
           intervals (error_over_large_intervals) and perform
           extrapolation. */

	if (!error_type2 && error_over_large_intervals > ertest)
          if (workspace.increase_nrmax())
            continue;

	// Perform extrapolation.

	table.append(area);

	if (table.get_nn() < 3) 
          goto skip_extrapolation;

	std::tie(reseps, abseps) = table.qelg();

	++ktmin;
	if (ktmin > 5 && err_ext < 0.001 * errsum)
          error_type = 5;

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

	// Prepare bisection of the smallest interval.

	if (table.get_nn() == 1)
          disallow_extrapolation = 1;

	if (error_type == 5)
          break;

      skip_extrapolation:

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
	if (result != _Tp{0} && area != _Tp{0})
	  {
	    if (err_ext / std::abs(res_ext) > errsum / std::abs(area))
	      goto compute_result;
	  }
	else if (err_ext > errsum)
          goto compute_result;
	else if (area == _Tp{0})
          goto return_error;
      }

    // Test on divergence.
    {
      auto max_area = std::max(std::abs(res_ext), std::abs(area));
      if (!positive_integrand && max_area < 0.01 * resabs0)
	goto return_error;
    }

    {
      auto ratio = res_ext / area;

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
      std::__throw_runtime_error("qagp_integrate: "
				 "Number of iterations was insufficient");
    else if (error_type == 2)
      std::__throw_runtime_error("qagp_integrate: "
				 "Cannot reach tolerance "
				 "because of roundoff error");
    else if (error_type == 3)
      std::__throw_runtime_error("qagp_integrate: "
				 "Bad integrand behavior found "
				 "in the integration interval");
    else if (error_type == 4)
      std::__throw_runtime_error("qagp_integrate: "
				 "Roundoff error detected "
				 "in the extrapolation table");
    else if (error_type == 5)
      std::__throw_runtime_error("qagp_integrate: "
				 "Integral is divergent, or slowly convergent");
    else
      std::__throw_runtime_error("qagp_integrate: "
				 "Could not integrate function");
  }

} // namespace __gnu_test

#endif // QAGP_INTEGRATE_H
