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

      std::size_t __limit;
      std::size_t __size;
      std::size_t __nrmax;
      std::size_t __ii;
      std::size_t __maximum_level;
      std::vector<_VecTp> __alist, __blist, __rlist, __elist;
      std::vector<std::size_t> __order, __level;

      void qpsrt();

    public:

      integration_workspace(std::size_t __lim)
      : __limit(__lim),
	__size(0),
	__nrmax(0),
	__ii(0),
	__maximum_level(0),
	__alist(__lim),
	__blist(__lim),
	__rlist(__lim),
	__elist(__lim),
	__order(__lim),
	__level(__lim)
      { }

      void
      set_initial(_VecTp __a0, _VecTp __b0, _VecTp __result0, _VecTp __error0)
      {
	__alist[0] = __a0;
	__blist[0] = __b0;
	__rlist[0] = __result0;
	__elist[0] = __error0;
	__order[0] = 0;
	__level[0] = 0;
	__size = 1;
      }

      void update (_VecTp __a1, _VecTp __b1, _VecTp __area1, _VecTp __error1,
		  _VecTp __a2, _VecTp __b2, _VecTp __area2, _VecTp __error2);

      void
      retrieve(_VecTp& __a, _VecTp& __b, _VecTp& __r, _VecTp& __e) const
      {
	__a = __alist[__ii];
	__b = __blist[__ii];
	__r = __rlist[__ii];
	__e = __elist[__ii];
      }

      bool
      increase_nrmax()
      {
	int __k;
	int __id = __nrmax;
	int __jupbnd;

	std::size_t __last = __size - 1;

	if (__last > (1 + __limit / 2))
	  __jupbnd = __limit + 1 - __last;
	else
	  __jupbnd = __last;

	for (__k = __id; __k <= __jupbnd; ++__k)
	  {
	    std::size_t __i_max = __order[__nrmax];

	    __ii = __i_max;

	    if (__level[__i_max] < __maximum_level)
	      return true;

	    ++__nrmax;

	  }
	return false;
      }

      void
      reset_nrmax()
      {
	__nrmax = 0;
	__ii = __order[0];
      }

      std::size_t
      get_nrmax() const
      {
	return __nrmax;
      }

      void
      set_nrmax(std::size_t NRMAX)
      {
	__nrmax = NRMAX;
      }

      _VecTp
      sum_results () const
      {
	std::size_t __kk;
	_VecTp __result_sum = 0;

	for (__kk = 0; __kk < __size; ++__kk)
	  {
	    __result_sum += __rlist[__kk];
	  }

	return __result_sum;
      }

      static bool
      subinterval_too_small (_VecTp __a1, _VecTp __a2, _VecTp __b2)
      {
	const _VecTp __e = std::numeric_limits<_VecTp>::epsilon();
	const _VecTp __u = std::numeric_limits<_VecTp>::min();

	_VecTp __tmp = (1 + 100 * __e) * (std::abs(__a2) + 1000 * __u);

	return std::abs(__a1) <= __tmp && std::abs(__b2) <= __tmp;
      }

      std::size_t
      current_level() const
      {
	return __level[__ii];
      }

      std::size_t
      max_level() const
      {
	return __maximum_level;
      }

      bool
      large_interval() const
      {
	if (__level[__ii] < __maximum_level)
	  return true;
	else
	  return false;
      }
    };

  template<typename _VecTp>
    void
    integration_workspace<_VecTp>::qpsrt()
    {
      const std::size_t __last = __size - 1;

      double __errmax;
      double __errmin;
      int __jj, __kk, __top;

      std::size_t __i_nrmax = __nrmax;
      std::size_t __i_maxerr = __order[__i_nrmax];

      // Check whether the list contains more than two error estimates
      if (__last < 2)
	{
	  __order[0] = 0;
	  __order[1] = 1;
	  __ii = __i_maxerr;
	  return;
	}

      __errmax = __elist[__i_maxerr];

      /* This part of the routine is only executed if, due to a difficult
	integrand, subdivision increased the error estimate. In the normal
	case the insert procedure should start after the nrmax-th largest
	error estimate. */

      while (__i_nrmax > 0 && __errmax > __elist[__order[__i_nrmax - 1]])
	{
	  __order[__i_nrmax] = __order[__i_nrmax - 1];
	  __i_nrmax--;
	}

      /* Compute the number of elements in the list to be maintained in
	descending order. This number depends on the number of
	subdivisions still allowed. */

      if(__last < (__limit/2 + 2))
	{
	  __top = __last;
	}
      else
	{
	  __top = __limit - __last + 1;
	}

      /* Insert errmax by traversing the list top-down, starting
	comparison from the element elist(order(i_nrmax+1)). */

      __jj = __i_nrmax + 1;

      /* The order of the tests in the following line is important to
	prevent a segmentation fault */

      while (__jj < __top && __errmax < __elist[__order[__jj]])
	{
	  __order[__jj-1] = __order[__jj];
	  ++__jj;
	}

      __order[__jj-1] = __i_maxerr;

      // Insert errmin by traversing the list bottom-up

      __errmin = __elist[__last];

      __kk = __top - 1;

      while (__kk > __jj - 2 && __errmin >= __elist[__order[__kk]])
	{
	  __order[__kk+1] = __order[__kk];
	  __kk--;
	}

      __order[__kk+1] = __last;

      // Set i_max and e_max

      __i_maxerr = __order[__i_nrmax];

      __ii = __i_maxerr;
      __nrmax = __i_nrmax;
    }

  template<typename _VecTp>
    void
    integration_workspace<_VecTp>::update(_VecTp __a1, _VecTp __b1,
					 _VecTp __area1, _VecTp __error1,
					 _VecTp __a2, _VecTp __b2,
					 _VecTp __area2, _VecTp __error2)
    {
      const std::size_t __i_max = __ii;
      const std::size_t __i_new = __size;

      const std::size_t __new_level = __level[__i_max] + 1;

      // append the newly-created intervals to the list

      if (__error2 > __error1)
	{
	  __alist[__i_max] = __a2;	// blist[maxerr] is already == b2
	  __rlist[__i_max] = __area2;
	  __elist[__i_max] = __error2;
	  __level[__i_max] = __new_level;

	  __alist[__i_new] = __a1;
	  __blist[__i_new] = __b1;
	  __rlist[__i_new] = __area1;
	  __elist[__i_new] = __error1;
	  __level[__i_new] = __new_level;
	}
      else
	{
	  __blist[__i_max] = __b1;	// alist[maxerr] is already == a1
	  __rlist[__i_max] = __area1;
	  __elist[__i_max] = __error1;
	  __level[__i_max] = __new_level;

	  __alist[__i_new] = __a2;
	  __blist[__i_new] = __b2;
	  __rlist[__i_new] = __area2;
	  __elist[__i_new] = __error2;
	  __level[__i_new] = __new_level;
	}

      ++__size;

      if (__new_level > __maximum_level)
	{
	  __maximum_level = __new_level;
	}

      qpsrt();
    }

} // namespace

#endif
