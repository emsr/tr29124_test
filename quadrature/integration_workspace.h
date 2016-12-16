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
// Implements the integration_workspace class which stores temporary data
// for performing integrals
// Based upon gsl-1.9/integration/workspace.c

#ifndef INTEGRATION_WORKSPACE_H
#define INTEGRATION_WORKSPACE_H

#include <vector>
#include <limits>
#include <cmath>

namespace __gnu_test
{

  template<typename _VecTp>
    class integration_workspace
    {
    private:

      std::size_t _M_limit;
      std::size_t _M_size;
      std::size_t _M_nrmax;
      std::size_t _M_ii;
      std::size_t _M_maximum_level;
      std::vector<_VecTp> _M_lower_lim, _M_upper_lim, _M_result, _M_abs_error;
      std::vector<std::size_t> _M_order, _M_level;

      void qpsrt();

    public:

      integration_workspace(std::size_t __lim)
      : _M_limit(__lim),
	_M_size(0),
	_M_nrmax(0),
	_M_ii(0),
	_M_maximum_level(0),
	_M_lower_lim(__lim),
	_M_upper_lim(__lim),
	_M_result(__lim),
	_M_abs_error(__lim),
	_M_order(__lim),
	_M_level(__lim)
      { }

      void
      set_initial(_VecTp __a0, _VecTp __b0, _VecTp __result0, _VecTp __error0)
      {
	_M_lower_lim[0] = __a0;
	_M_upper_lim[0] = __b0;
	_M_result[0] = __result0;
	_M_abs_error[0] = __error0;
	_M_order[0] = 0;
	_M_level[0] = 0;
	_M_size = 1;
      }

      void update(_VecTp __a1, _VecTp __b1, _VecTp __area1, _VecTp __error1,
		  _VecTp __a2, _VecTp __b2, _VecTp __area2, _VecTp __error2);

      void
      retrieve(_VecTp& __lolim, _VecTp& __uplim, _VecTp& __res, _VecTp& __err) const
      {
	__lolim = _M_lower_lim[_M_ii];
	__uplim = _M_upper_lim[_M_ii];
	__res = _M_result[_M_ii];
	__err = _M_abs_error[_M_ii];
      }

      size_t
      size() const
      { return _M_size; }

      _VecTp
      lower_lim(size_t ii) const
      { return _M_lower_lim[ii]; }

      _VecTp
      upper_lim(size_t ii) const
      { return _M_upper_lim[ii]; }

      _VecTp
      result(size_t ii) const
      { return _M_result[ii]; }

      _VecTp
      abs_error(size_t ii) const
      { return _M_abs_error[ii]; }

      size_t
      order(size_t ii) const
      { return _M_order[ii]; }

      size_t
      level(size_t ii) const
      { return _M_level[ii]; }

      bool
      increase_nrmax()
      {
	int __k;
	int __id = _M_nrmax;
	int __jupbnd;

	std::size_t __last = _M_size - 1;

	if (__last > (1 + _M_limit / 2))
	  __jupbnd = _M_limit + 1 - __last;
	else
	  __jupbnd = __last;

	for (__k = __id; __k <= __jupbnd; ++__k)
	  {
	    std::size_t __i_max = _M_order[_M_nrmax];

	    _M_ii = __i_max;

	    if (_M_level[__i_max] < _M_maximum_level)
	      return true;

	    ++_M_nrmax;

	  }
	return false;
      }

      void
      reset_nrmax()
      {
	_M_nrmax = 0;
	_M_ii = _M_order[0];
      }

      std::size_t
      get_nrmax() const
      {
	return _M_nrmax;
      }

      void
      set_nrmax(std::size_t NRMAX)
      {
	_M_nrmax = NRMAX;
      }

      _VecTp
      sum_results() const
      {
	_VecTp __result_sum = 0;

	for (std::size_t __kk = 0; __kk < _M_size; ++__kk)
	  __result_sum += _M_result[__kk];

	return __result_sum;
      }

      static bool
      subinterval_too_small(_VecTp __a1, _VecTp __a2, _VecTp __b2)
      {
	const _VecTp __e = std::numeric_limits<_VecTp>::epsilon();
	const _VecTp __u = std::numeric_limits<_VecTp>::min();

	_VecTp __tmp = (1 + 100 * __e) * (std::abs(__a2) + 1000 * __u);

	return std::abs(__a1) <= __tmp && std::abs(__b2) <= __tmp;
      }

      std::size_t
      current_level() const
      {
	return _M_level[_M_ii];
      }

      std::size_t
      max_level() const
      {
	return _M_maximum_level;
      }

      bool
      large_interval() const
      {
	if (_M_level[_M_ii] < _M_maximum_level)
	  return true;
	else
	  return false;
      }
    };

  template<typename _VecTp>
    void
    integration_workspace<_VecTp>::qpsrt()
    {
      const std::size_t __last = _M_size - 1;

      double __errmax;
      double __errmin;
      int __jj, __kk, __top;

      std::size_t __i_nrmax = _M_nrmax;
      std::size_t __i_maxerr = _M_order[__i_nrmax];

      // Check whether the list contains more than two error estimates
      if (__last < 2)
	{
	  _M_order[0] = 0;
	  _M_order[1] = 1;
	  _M_ii = __i_maxerr;
	  return;
	}

      __errmax = _M_abs_error[__i_maxerr];

      // This part of the routine is only executed if, due to a difficult
      // integrand, subdivision increased the error estimate. In the normal
      // case the insert procedure should start after the nrmax-th largest
      // error estimate.

      while (__i_nrmax > 0 && __errmax > _M_abs_error[_M_order[__i_nrmax - 1]])
	{
	  _M_order[__i_nrmax] = _M_order[__i_nrmax - 1];
	  __i_nrmax--;
	}

      // Compute the number of elements in the list to be maintained in
      // descending order. This number depends on the number of
      // subdivisions still allowed.

      if(__last < (_M_limit/2 + 2))
	__top = __last;
      else
	__top = _M_limit - __last + 1;

      // Insert errmax by traversing the list top-down, starting
      // comparison from the element elist(order(i_nrmax+1)).


      // The order of the tests in the following line is important to
      // prevent a segmentation fault

      __jj = __i_nrmax + 1;
      while (__jj < __top && __errmax < _M_abs_error[_M_order[__jj]])
	{
	  _M_order[__jj - 1] = _M_order[__jj];
	  ++__jj;
	}
      _M_order[__jj - 1] = __i_maxerr;

      // Insert errmin by traversing the list bottom-up

      __errmin = _M_abs_error[__last];

      __kk = __top - 1;
      while (__kk > __jj - 2 && __errmin >= _M_abs_error[_M_order[__kk]])
	{
	  _M_order[__kk + 1] = _M_order[__kk];
	  --__kk;
	}
      _M_order[__kk + 1] = __last;

      // Set i_max and e_max

      __i_maxerr = _M_order[__i_nrmax];

      _M_ii = __i_maxerr;
      _M_nrmax = __i_nrmax;
    }

  template<typename _VecTp>
    void
    integration_workspace<_VecTp>::update(_VecTp __a1, _VecTp __b1,
					  _VecTp __area1, _VecTp __error1,
					  _VecTp __a2, _VecTp __b2,
					  _VecTp __area2, _VecTp __error2)
    {
      const std::size_t __i_max = _M_ii;
      const std::size_t __i_new = _M_size;

      const std::size_t __new_level = _M_level[__i_max] + 1;

      // append the newly-created intervals to the list

      if (__error2 > __error1)
	{
	  _M_lower_lim[__i_max] = __a2;	// blist[maxerr] is already == b2
	  _M_result[__i_max] = __area2;
	  _M_abs_error[__i_max] = __error2;
	  _M_level[__i_max] = __new_level;

	  _M_lower_lim[__i_new] = __a1;
	  _M_upper_lim[__i_new] = __b1;
	  _M_result[__i_new] = __area1;
	  _M_abs_error[__i_new] = __error1;
	  _M_level[__i_new] = __new_level;
	}
      else
	{
	  _M_upper_lim[__i_max] = __b1;	// alist[maxerr] is already == a1
	  _M_result[__i_max] = __area1;
	  _M_abs_error[__i_max] = __error1;
	  _M_level[__i_max] = __new_level;

	  _M_lower_lim[__i_new] = __a2;
	  _M_upper_lim[__i_new] = __b2;
	  _M_result[__i_new] = __area2;
	  _M_abs_error[__i_new] = __error2;
	  _M_level[__i_new] = __new_level;
	}

      ++_M_size;

      if (__new_level > _M_maximum_level)
	_M_maximum_level = __new_level;

      qpsrt();
    }

} // namespace __gnu_test

#endif // INTEGRATION_WORKSPACE_H
