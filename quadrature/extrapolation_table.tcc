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
// Ported from GSL by Ed Smith-Rowland
// Originally written by Brian Gaugh
//
//This file implements an extrapolation table for use in integration schemes
//Based on gsl/integration/qelg.c

#ifndef EXTRAPOLATION_TABLE_TCC
#define EXTRAPOLATION_TABLE_TCC 1

namespace __gnu_test
{

  template<typename _Tp>
    std::tuple<_Tp, _Tp>
    extrapolation_table<_Tp>::qelg()
    {
      const auto __cur_n = this->_M_nn - 1;
      const auto __current = this->_M_rlist2[__cur_n];

      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_max = std::numeric_limits<_Tp>::max();

      auto __absolute = _S_max;
      auto __relative = 5 * _S_eps * std::abs(__current);

      const auto __newelm = __cur_n / 2;
      const auto __n_orig = __cur_n;
      auto __n_final = __cur_n;

      const auto __nres_orig = this->_M_nres;

      auto __result = __current;
      auto __abserr = _S_max;

      if (__cur_n < 2)
	{
	  __result = __current;
	  __abserr = std::max(__absolute, __relative);
	  return std::make_tuple(__result, __abserr);
	}

      this->_M_rlist2[__cur_n + 2] = this->_M_rlist2[__cur_n];
      this->_M_rlist2[__cur_n] = _S_max;

      for (size_t __ii = 0; __ii < __newelm; ++__ii)
	{
	  auto __res = this->_M_rlist2[__cur_n - 2 * __ii + 2];
	  const auto __e0 = this->_M_rlist2[__cur_n - 2 * __ii - 2];
	  const auto __e1 = this->_M_rlist2[__cur_n - 2 * __ii - 1];
	  const auto __e2 = __res;

	  const auto __e1abs = std::abs(__e1);
	  const auto __delta2 = __e2 - __e1;
	  auto __err2 = std::abs(__delta2);
	  auto __tol2 = std::max(std::abs(__e2), __e1abs) * _S_eps;
	  auto __delta3 = __e1 - __e0;
	  auto __err3 = std::abs(__delta3);
	  auto __tol3 = std::max(__e1abs, std::abs(__e0)) * _S_eps;

	  if (__err2 <= __tol2 && __err3 <= __tol3)
	    {
	      // If e0, e1 and e2 are equal to within machine accuracy,
	      // convergence is assumed.
	      __result = __res;
	      __absolute = __err2 + __err3;
	      __relative = 5 * _S_eps * std::abs(__res);
	      __abserr = std::max(__absolute, __relative);
	      return std::make_tuple(__result, __abserr);
	    }

	  auto __e3 = this->_M_rlist2[__cur_n - 2 * __ii];
	  this->_M_rlist2[__cur_n - 2 * __ii] = __e1;
	  auto __delta1 = __e1 - __e3;
	  auto __err1 = std::abs(__delta1);
	  auto __tol1 = std::max(__e1abs, std::abs(__e3)) * _S_eps;

	  // If two elements are very close to each other, omit a part of
	  // the table by adjusting the value of n.
	  if (__err1 <= __tol1 || __err2 <= __tol2 || __err3 <= __tol3)
	    {
	      __n_final = 2 * __ii;
	      break;
	    }

	  auto __ss = (_Tp{1} / __delta1 + _Tp{1} / __delta2)
		    - _Tp{1} / __delta3;

	  // Test to detect irregular behaviour in the table,
	  // and eventually omit a part of the table by adjusting
	  // the value of n.
	  if (std::abs(__ss * __e1) <= 0.0001)
	    {
	      __n_final = 2 * __ii;
	      break;
	    }

	  // Compute a new element and eventually adjust the value of result.
	  __res = __e1 + _Tp{1} / __ss;
	  this->_M_rlist2[__cur_n - 2 * __ii] = __res;
	  {
	    const auto __error = __err2 + std::abs(__res - __e2) + __err3;

	    if (__error <= __abserr)
	      {
		__abserr = __error;
		__result = __res;
	      }
	  }
	}

      // Shift the table.

      {
	const size_t __limexp = 50 - 1;

	if (__n_final == __limexp)
	  __n_final = 2 * (__limexp / 2);
      }

      if (__n_orig % 2 == 1)
	{
	  for (size_t __ii = 0; __ii <= __newelm; ++__ii)
	    this->_M_rlist2[1 + __ii * 2] = this->_M_rlist2[__ii * 2 + 3];
	}
      else
	{
	  for (size_t __ii = 0; __ii <= __newelm; ++__ii)
	    this->_M_rlist2[__ii * 2] = this->_M_rlist2[__ii * 2 + 2];
	}

      if (__n_orig != __n_final)
	{
	  for (size_t __ii = 0; __ii <= __n_final; ++__ii)
	    this->_M_rlist2[__ii] = this->_M_rlist2[__n_orig - __n_final + __ii];
	}

      this->_M_nn = __n_final + 1;

      if (__nres_orig < 3)
	{
	  this->_M_res3la[__nres_orig] = __result;
	  __abserr = _S_max;
	}
      else
	{ // Compute error estimate.
	  __abserr = (std::abs(__result - this->_M_res3la[2])
		    + std::abs(__result - this->_M_res3la[1])
		    + std::abs(__result - this->_M_res3la[0]));

	  this->_M_res3la[0] = this->_M_res3la[1];
	  this->_M_res3la[1] = this->_M_res3la[2];
	  this->_M_res3la[2] = __result;
	}

      /* In QUADPACK the variable table->nres is incremented at the top of
	qelg, so it increases on every call. This leads to the array
	res3la being accessed when its elements are still undefined, so I
	have moved the update to this point so that its value more
	useful. */

      this->_M_nres = __nres_orig + 1;

      __abserr = std::max(__abserr, 5 * _S_eps * std::abs(__result));

      return std::make_tuple(__result, __abserr);
    }

} // namespace __gnu_test

#endif // EXTRAPOLATION_TABLE_TCC
