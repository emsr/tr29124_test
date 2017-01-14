// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2017 Free Software Foundation, Inc.
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
// Implements the integration_workspace class which stores temporary data
// for performing integrals
// Based on gsl/integration/workspace.c

#ifndef INTEGRATION_WORKSPACE_TCC
#define INTEGRATION_WORKSPACE_TCC 1

#include "integration_error.h"

namespace __gnu_test
{

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::sort_error()
    {
/*
      interval_comp __cmp;
      //std::make_heap(std::begin(this->_M_ival),
	//	     std::begin(this->_M_ival) + this->_M_size, __cmp);
      std::sort(std::begin(this->_M_ival),
		std::begin(this->_M_ival) + this->_M_size, __cmp);
      std::reverse(std::begin(this->_M_ival),
		   std::begin(this->_M_ival) + this->_M_size);
      for (auto __k = 0u; __k < this->_M_size; ++__k)
	this->_M_ival[__k]._M_order = __k;
      this->_M_current_index = 0;
      this->_M_nrmax = 0;
      return;
*/
      if (this->_M_size < 2)
	return;

      const auto __last = this->_M_size - 1;

      auto __i_nrmax = this->_M_nrmax;
      auto __i_maxerr = this->_M_ival[__i_nrmax]._M_order;

      // Check whether the list contains more than two error estimates.
      if (__last < 2)
	{
	  this->_M_ival[0]._M_order = 0;
	  this->_M_ival[1]._M_order = 1;
	  this->_M_current_index = __i_maxerr;
	  return;
	}

      auto __errmax = this->_M_ival[__i_maxerr]._M_abs_error;

      // This part of the routine is only executed if, due to a difficult
      // integrand, subdivision increased the error estimate. In the normal
      // case the insert procedure should start after the nrmax-th largest
      // error estimate.
      while (__i_nrmax > 0
	  && __errmax > this->_M_ival[this->_M_ival[__i_nrmax - 1]._M_order]._M_abs_error)
	{
	  this->_M_ival[__i_nrmax]._M_order = this->_M_ival[__i_nrmax - 1]._M_order;
	  --__i_nrmax;
	}

      // Compute the number of elements in the list to be maintained in
      // descending order. This number depends on the number of
      // subdivisions still allowed.
      std::size_t __top;
      if(__last < (2 + this->_M_capacity / 2))
	__top = __last;
      else
	__top = this->_M_capacity - __last + 1;

      // Insert errmax by traversing the list top-down, starting
      // comparison from the element abs_error(order(i_nrmax + 1)).
      auto __jj = __i_nrmax + 1;

      // The order of the tests in the following line is important to
      // prevent a segmentation fault
      while (__jj < __top
	 && __errmax < this->_M_ival[this->_M_ival[__jj]._M_order]._M_abs_error)
	{
	  this->_M_ival[__jj - 1]._M_order = this->_M_ival[__jj]._M_order;
	  ++__jj;
	}
      this->_M_ival[__jj - 1]._M_order = __i_maxerr;

      // Insert errmin by traversing the list bottom-up
      const auto __errmin = this->_M_ival[__last]._M_abs_error;
      auto __kk = __top - 1;
      while (__kk > __jj - 2
	  && __errmin >= this->_M_ival[this->_M_ival[__kk]._M_order]._M_abs_error)
	{
	  this->_M_ival[__kk + 1]._M_order = this->_M_ival[__kk]._M_order;
	  --__kk;
	}
      this->_M_ival[__kk + 1]._M_order = __last;

      // Set i_max and e_max
      __i_maxerr = this->_M_ival[__i_nrmax]._M_order;
      this->_M_current_index = __i_maxerr;
      this->_M_nrmax = __i_nrmax;
    }

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::sort_results()
    {
      for (std::size_t __i = 0; __i < this->_M_size; ++__i)
	{
	  auto __i1 = this->_M_ival[__i]._M_order;
	  auto __e1 = this->_M_ival[__i1]._M_abs_error;
	  auto __i_max = __i1;

	  for (auto __j = __i + 1; __j < this->_M_size; ++__j)
	    {
	      auto __i2 = this->_M_ival[__j]._M_order;
	      auto __e2 = this->_M_ival[__i2]._M_abs_error;

	      if (__e2 >= __e1)
		{
		  __i_max = __i2;
		  __e1 = __e2;
		}
	    }

	  if (__i_max != __i1)
	    {
	      this->_M_ival[__i]._M_order = this->_M_ival[__i_max]._M_order;
	      this->_M_ival[__i_max]._M_order = __i1;
	    }
	}

      this->_M_current_index = this->_M_ival[0]._M_order;
    }

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::append(_Tp __a, _Tp __b,
				       _Tp __area, _Tp __error)
    {
      const std::size_t __i_new = this->_M_size;
      if (__i_new >= this->capacity())
	{
	  if (this->_M_try_resize)
	    this->resize();
	  else
	    __throw__IntegrationError<_Tp>("integration_workspace: "
					   "Exhausted work space.",
					   MAX_ITER_ERROR);
	}

      // Append the newly-created interval to the list.
      this->_M_ival[__i_new]._M_lower_lim = __a;
      this->_M_ival[__i_new]._M_upper_lim = __b;
      this->_M_ival[__i_new]._M_result = __area;
      this->_M_ival[__i_new]._M_abs_error = __error;
      this->_M_ival[__i_new]._M_order = __i_new;
      this->_M_ival[__i_new]._M_level = 0;

      ++this->_M_size;
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
      const auto __i_max = this->_M_current_index;
      const auto __i_new = this->_M_size;
      if (__i_new >= this->capacity())
	{
	  if (this->_M_try_resize)
	    this->resize();
	  else
	    __throw__IntegrationError<_Tp>("integration_workspace: "
					   "Exhausted work space.",
					   MAX_ITER_ERROR);
	}

      const auto __new_level = this->_M_ival[__i_max]._M_level + 1;

      // Append the newly-created intervals to the list.
      if (__error2 > __error1)
	{
	  this->_M_ival[__i_max]._M_lower_lim = __a2;
	  // upper_lim[i_max] is already == b2
	  this->_M_ival[__i_max]._M_result = __area2;
	  this->_M_ival[__i_max]._M_abs_error = __error2;
	  this->_M_ival[__i_max]._M_level = __new_level;

	  this->_M_ival[__i_new]._M_lower_lim = __a1;
	  this->_M_ival[__i_new]._M_upper_lim = __b1;
	  this->_M_ival[__i_new]._M_result = __area1;
	  this->_M_ival[__i_new]._M_abs_error = __error1;
	  this->_M_ival[__i_new]._M_level = __new_level;
	}
      else
	{
	  // lower_lim[i_max] is already == a1
	  this->_M_ival[__i_max]._M_upper_lim = __b1;
	  this->_M_ival[__i_max]._M_result = __area1;
	  this->_M_ival[__i_max]._M_abs_error = __error1;
	  this->_M_ival[__i_max]._M_level = __new_level;

	  this->_M_ival[__i_new]._M_lower_lim = __a2;
	  this->_M_ival[__i_new]._M_upper_lim = __b2;
	  this->_M_ival[__i_new]._M_result = __area2;
	  this->_M_ival[__i_new]._M_abs_error = __error2;
	  this->_M_ival[__i_new]._M_level = __new_level;
	}

      ++this->_M_size;

      if (__new_level > this->_M_maximum_level)
	this->_M_maximum_level = __new_level;

      this->sort_error();
    }

} // namespace __gnu_test

#endif // INTEGRATION_WORKSPACE_TCC
