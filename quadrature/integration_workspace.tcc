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

#ifndef INTEGRATION_WORKSPACE_TCC
#define INTEGRATION_WORKSPACE_TCC 1

namespace __gnu_test
{

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::sort_error()
    {
      const std::size_t __last = _M_size - 1;

      double __errmax;
      double __errmin;
      int __jj, __kk, __top;

      std::size_t __i_nrmax = _M_nrmax;
      std::size_t __i_maxerr = _M_order[__i_nrmax];

      // Check whether the list contains more than two error estimates.
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

      if(__last < (_M_capacity/2 + 2))
	__top = __last;
      else
	__top = _M_capacity - __last + 1;

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

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::append(_Tp __a, _Tp __b,
					  _Tp __area, _Tp __error)
    {
      const std::size_t __i_max = this->_M_ii;
      const std::size_t __i_new = this->_M_size;

      const std::size_t __new_level = this->_M_level[__i_max] + 1;

      // append the newly-created interval to the list

      this->_M_lower_lim[__i_new] = __a;
      this->_M_upper_lim[__i_new] = __b;
      this->_M_result[__i_new] = __area;
      this->_M_abs_error[__i_new] = __error;
      this->_M_level[__i_new] = __new_level;

      ++this->_M_size;

      if (__new_level > this->_M_maximum_level)
	this->_M_maximum_level = __new_level;

      this->sort_error();
    }

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::update(_Tp __a1, _Tp __b1,
					  _Tp __area1, _Tp __error1,
					  _Tp __a2, _Tp __b2,
					  _Tp __area2, _Tp __error2)
    {
      const std::size_t __i_max = this->_M_ii;
      const std::size_t __i_new = this->_M_size;

      const std::size_t __new_level = this->_M_level[__i_max] + 1;

      // append the newly-created intervals to the list

      if (__error2 > __error1)
	{
	  this->_M_lower_lim[__i_max] = __a2;	// blist[maxerr] is already == b2
	  this->_M_result[__i_max] = __area2;
	  this->_M_abs_error[__i_max] = __error2;
	  this->_M_level[__i_max] = __new_level;

	  this->_M_lower_lim[__i_new] = __a1;
	  this->_M_upper_lim[__i_new] = __b1;
	  this->_M_result[__i_new] = __area1;
	  this->_M_abs_error[__i_new] = __error1;
	  this->_M_level[__i_new] = __new_level;
	}
      else
	{
	  this->_M_upper_lim[__i_max] = __b1;	// alist[maxerr] is already == a1
	  this->_M_result[__i_max] = __area1;
	  this->_M_abs_error[__i_max] = __error1;
	  this->_M_level[__i_max] = __new_level;

	  this->_M_lower_lim[__i_new] = __a2;
	  this->_M_upper_lim[__i_new] = __b2;
	  this->_M_result[__i_new] = __area2;
	  this->_M_abs_error[__i_new] = __error2;
	  this->_M_level[__i_new] = __new_level;
	}

      ++this->_M_size;

      if (__new_level > this->_M_maximum_level)
	this->_M_maximum_level = __new_level;

      this->sort_error();
    }

} // namespace __gnu_test

#endif // INTEGRATION_WORKSPACE_TCC
