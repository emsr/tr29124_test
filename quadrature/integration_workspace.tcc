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
// Based upon gsl-2.3/integration/workspace.c

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
      const std::size_t __last = this->_M_size - 1;

      std::size_t __i_nrmax = this->_M_nrmax;
      std::size_t __i_maxerr = this->_M_order[__i_nrmax];

      // Check whether the list contains more than two error estimates.
      if (__last < 2)
	{
	  this->_M_order[0] = 0;
	  this->_M_order[1] = 1;
	  this->_M_current_index = __i_maxerr;
	  return;
	}

      auto __errmax = this->_M_abs_error[__i_maxerr];

      // This part of the routine is only executed if, due to a difficult
      // integrand, subdivision increased the error estimate. In the normal
      // case the insert procedure should start after the nrmax-th largest
      // error estimate.

      while (__i_nrmax > 0
	  && __errmax > this->_M_abs_error[this->_M_order[__i_nrmax - 1]])
	{
	  this->_M_order[__i_nrmax] = this->_M_order[__i_nrmax - 1];
	  --__i_nrmax;
	}

      // Compute the number of elements in the list to be maintained in
      // descending order. This number depends on the number of
      // subdivisions still allowed.

      int __top;
      if(__last < (this->_M_capacity/2 + 2))
	__top = __last;
      else
	__top = this->_M_capacity - __last + 1;

      // Insert errmax by traversing the list top-down, starting
      // comparison from the element abs_error(order(i_nrmax+1)).

      // The order of the tests in the following line is important to
      // prevent a segmentation fault

      int __jj = __i_nrmax + 1;
      while (__jj < __top
	  && __errmax < this->_M_abs_error[this->_M_order[__jj]])
	{
	  this->_M_order[__jj - 1] = this->_M_order[__jj];
	  ++__jj;
	}
      this->_M_order[__jj - 1] = __i_maxerr;

      // Insert errmin by traversing the list bottom-up

      auto __errmin = this->_M_abs_error[__last];

      int __kk = __top - 1;
      while (__kk > __jj - 2
	  && __errmin >= this->_M_abs_error[this->_M_order[__kk]])
	{
	  this->_M_order[__kk + 1] = this->_M_order[__kk];
	  --__kk;
	}
      this->_M_order[__kk + 1] = __last;

      // Set i_max and e_max

      __i_maxerr = this->_M_order[__i_nrmax];

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
	  auto __i1 = this->_M_order[__i];
	  auto __e1 = this->_M_abs_error[__i1];
	  auto __i_max = __i1;

	  for (auto __j = __i + 1; __j < this->_M_size; ++__j)
            {
              auto __i2 = this->_M_order[__j];
              auto __e2 = this->_M_abs_error[__i2];

              if (__e2 >= __e1)
        	{
        	  __i_max = __i2;
        	  __e1 = __e2;
        	}
            }

	  if (__i_max != __i1)
            {
              this->_M_order[__i] = this->_M_order[__i_max];
              this->_M_order[__i_max] = __i1;
            }
	}

      this->_M_current_index = this->_M_order[0];
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

      // append the newly-created interval to the list

      this->_M_lower_lim[__i_new] = __a;
      this->_M_upper_lim[__i_new] = __b;
      this->_M_result[__i_new] = __area;
      this->_M_abs_error[__i_new] = __error;
      this->_M_order[__i_new] = __i_new;
      this->_M_level[__i_new] = 0;

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
      const std::size_t __i_max = this->_M_current_index;
      const std::size_t __i_new = this->_M_size;

      const std::size_t __new_level = this->_M_level[__i_max] + 1;

      // append the newly-created intervals to the list

      if (__error2 > __error1)
	{
	  this->_M_lower_lim[__i_max] = __a2;
	  // upper_lim[i_max] is already == b2
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
	  // lower_lim[i_max] is already == a1
	  this->_M_upper_lim[__i_max] = __b1;
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
