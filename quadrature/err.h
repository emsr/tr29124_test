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
//Used to interpret the error estimate of the itnegration routines
//Based upon gsl-1.9/integration/err.c

#ifndef ERR_H
#define ERR_H

#include <cmath>
#include <limits>

namespace __gnu_test
{
  template<typename _VecTp>
    _VecTp __rescale_error(_VecTp __err,
			 const _VecTp __result_abs, const _VecTp __result_asc);

  template<typename _VecTp>
    _VecTp
    __rescale_error(_VecTp __err,
		    const _VecTp __result_abs, const _VecTp __result_asc)
    {
      __err = std::abs(__err);
      if (__result_asc != 0 && __err != 0)
	{
	  auto __scale = std::pow((200 * __err / __result_asc), 1.5);

	  if (__scale < 1)
	    __err = __result_asc * __scale;
	  else
	    __err = __result_asc;
	}
      if (__result_abs >
	    std::numeric_limits<_VecTp>::min()
		/ (50 * std::numeric_limits<_VecTp>::epsilon()))
	{
	  auto __min_err = 50 * std::numeric_limits<_VecTp>::epsilon() * __result_abs;

	  if (__min_err > __err)
	    __err = __min_err;
	}

      return __err;
    }

} //namespace

#endif
