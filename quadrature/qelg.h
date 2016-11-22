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
//This file implements an extrapolation table for use in integration schemes
//Based upon gsl-1.9/integration/qelg.c

#include <array>
#include <utility>
#include <limits>
#include <cmath>
#include <algorithm>

#ifndef QELG_H
#define QELG_H
namespace __gnu_test
{
  template<typename _VecTp>
    class extrapolation_table
    {

    private:

      size_t __nn;
      std::array<_VecTp, 52> __rlist2;
      size_t __nres;
      std::array<_VecTp, 3> __res3la;

    public:

      extrapolation_table()
      : __nn(0),
	__nres(0)
      {}

      explicit extrapolation_table(_VecTp __y)
      : __nn(0),
	__nres(0)
      { __rlist2[0] = __y; }

      void
      append(_VecTp __y)
      {
	__rlist2[__nn] = __y;
	++__nn;
      }

      std::pair<_VecTp, _VecTp>
      qelg();

      size_t
      get_nn() const
      { return __nn; }
    };

  template<typename _VecTp>
    std::pair<_VecTp, _VecTp>
    extrapolation_table<_VecTp>::qelg()
    {
      const size_t __cur_n = __nn - 1;
      const _VecTp __current = __rlist2[__cur_n];

      _VecTp __absolute = std::numeric_limits<_VecTp>::max();
      _VecTp __relative = 5 * std::numeric_limits<_VecTp>::epsilon() * std::abs(__current);

      const size_t __newelm = __cur_n / 2;
      const size_t __n_orig = __cur_n;
      size_t __n_final = __cur_n;
      size_t __ii;

      const size_t __nres_orig = __nres;

      _VecTp __result = __current;
      _VecTp __abserr = std::numeric_limits<_VecTp>::max();

      if (__cur_n < 2)
	{
	  __result = __current;
	  __abserr = std::max(__absolute, __relative);
	  return std::make_pair(__result, __abserr);
	}

      __rlist2[__cur_n + 2] = __rlist2[__cur_n];
      __rlist2[__cur_n] = std::numeric_limits<_VecTp>::max();

      for (size_t __ii = 0; __ii < __newelm; ++__ii)
	{
	  _VecTp __res = __rlist2[__cur_n - 2 * __ii + 2];
	  _VecTp __e0 = __rlist2[__cur_n - 2 * __ii - 2];
	  _VecTp __e1 = __rlist2[__cur_n - 2 * __ii - 1];
	  _VecTp __e2 = __res;

	  _VecTp __e1abs = std::abs(__e1);
	  _VecTp __delta2 = __e2 - __e1;
	  _VecTp __err2 = std::abs(__delta2);
	  _VecTp __tol2 = std::max(std::abs(__e2), __e1abs) * std::numeric_limits<_VecTp>::epsilon();
	  _VecTp __delta3 = __e1 - __e0;
	  _VecTp __err3 = std::abs(__delta3);
	  _VecTp __tol3 = std::max(__e1abs, std::abs(__e0)) * std::numeric_limits<_VecTp>::epsilon();

	  _VecTp __e3, __delta1, __err1, __tol1, __ss;

	  if (__err2 <= __tol2 && __err3 <= __tol3)
	    {
	      /* If e0, e1 and e2 are equal to within machine accuracy,
		convergence is assumed.  */

	      __result = __res;
	      __absolute = __err2 + __err3;
	      __relative = 5 * std::numeric_limits<_VecTp>::epsilon() * std::abs(__res);
	      __abserr = std::max(__absolute, __relative);
	      return std::make_pair(__result, __abserr);
	    }

	  __e3 = __rlist2[__cur_n - 2 * __ii];
	  __rlist2[__cur_n - 2 * __ii] = __e1;
	  __delta1 = __e1 - __e3;
	  __err1 = std::abs(__delta1);
	  __tol1 = std::max(__e1abs, std::abs(__e3)) * std::numeric_limits<_VecTp>::epsilon();

	  /* If two elements are very close to each other, omit a part of
	    the table by adjusting the value of n */

	  if (__err1 <= __tol1 || __err2 <= __tol2 || __err3 <= __tol3)
	    {
	      __n_final = 2 * __ii;
	      break;
	    }

	  __ss = (1 / __delta1 + 1 / __delta2) - 1 / __delta3;

	  /* Test to detect irregular behaviour in the table, and
	    eventually omit a part of the table by adjusting the value of
	    n. */

	  if (std::abs(__ss * __e1) <= 0.0001)
	    {
	      __n_final = 2 * __ii;
	      break;
	    }

	  /* Compute a new element and eventually adjust the value of
	    result. */

	  __res = __e1 + 1 / __ss;
	  __rlist2[__cur_n - 2 * __ii] = __res;

	  {
	    const _VecTp __error = __err2 + std::abs(__res - __e2) + __err3;

	    if (__error <= __abserr)
	      {
		__abserr = __error;
		__result = __res;
	      }
	  }
	}

      /* Shift the table */

      {
	const size_t __limexp = 50 - 1;

	if (__n_final == __limexp)
	  __n_final = 2 * (__limexp / 2);
      }

      if (__n_orig % 2 == 1)
	{
	  for (__ii = 0; __ii <= __newelm; ++__ii)
	    __rlist2[1 + __ii * 2] = __rlist2[__ii * 2 + 3];
	}
      else
	{
	  for (__ii = 0; __ii <= __newelm; ++__ii)
	    __rlist2[__ii * 2] = __rlist2[__ii * 2 + 2];
	}

      if (__n_orig != __n_final)
	{
	  for (__ii = 0; __ii <= __n_final; ++__ii)
	    __rlist2[__ii] = __rlist2[__n_orig - __n_final + __ii];
	}

      __nn = __n_final + 1;

      if (__nres_orig < 3)
	{
	  __res3la[__nres_orig] = __result;
	  __abserr = std::numeric_limits<_VecTp>::max();
	}
      else
	{			   /* Compute error estimate */
	  __abserr = (std::abs(__result - __res3la[2]) + std::abs(__result - __res3la[1])
		    + std::abs(__result - __res3la[0]));

	  __res3la[0] = __res3la[1];
	  __res3la[1] = __res3la[2];
	  __res3la[2] = __result;
	}

      /* In QUADPACK the variable table->nres is incremented at the top of
	qelg, so it increases on every call. This leads to the array
	res3la being accessed when its elements are still undefined, so I
	have moved the update to this point so that its value more
	useful. */

      __nres = __nres_orig + 1;

      __abserr = std::max(__abserr, 5 * std::numeric_limits<_VecTp>::epsilon() * std::abs(__result));

      return std::make_pair(__result, __abserr);
    }

}//namespace

#endif
