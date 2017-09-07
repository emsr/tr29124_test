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

namespace __gnu_ext
{

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::sort_error()
    {
      std::make_heap(this->begin(), this->end(), interval_comp{});
      return;
    }

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::append(_Tp __a, _Tp __b,
				       _Tp __area, _Tp __error,
				       std::size_t __depth)
    {
      interval __iv;
      __iv._M_lower_lim = __a;
      __iv._M_upper_lim = __b;
      __iv._M_result = __area;
      __iv._M_abs_error = __error;
      __iv._M_depth = __depth;
      this->push(__iv);
    }

  /**
   *
   */
  template<typename _Tp>
    void
    integration_workspace<_Tp>::split(_Tp __ab,
				      _Tp __area1, _Tp __error1,
				      _Tp __area2, _Tp __error2)
    {
      auto __iv = this->top();
      const auto __a1 = __iv._M_lower_lim;
      const auto __b1 = __ab;
      const auto __a2 = __ab;
      const auto __b2 = __iv._M_upper_lim;
      const auto __depth = __iv._M_depth + 1;
      this->pop();

      interval __iv1;
      __iv1._M_lower_lim = __a1;
      __iv1._M_upper_lim = __b1;
      __iv1._M_result = __area1;
      __iv1._M_abs_error = __error1;
      __iv1._M_depth = __depth;
      this->push(__iv1);

      interval __iv2;
      __iv2._M_lower_lim = __a2;
      __iv2._M_upper_lim = __b2;
      __iv2._M_result = __area2;
      __iv2._M_abs_error = __error2;
      __iv2._M_depth = __depth;
      this->push(__iv2);

      if (__depth > this->_M_max_depth)
	this->_M_max_depth = __depth;
    }

  /**
   * 
   */
  template<typename _Tp>
    bool
    integration_workspace<_Tp>::increment_start()
    {
      const auto __i_max = this->_M_start;
      for (auto __k = __i_max; __k < this->size(); ++__k)
	{
	  if (this->_M_ival[__k]._M_depth < this->_M_max_depth)
	    return true;
	  else if (this->_M_start + 1 < this->size())
	    {
	      ++this->_M_start;
	      this->sort_error();
	    }
	}
      return false;
    }

} // namespace __gnu_ext

#endif // INTEGRATION_WORKSPACE_TCC
