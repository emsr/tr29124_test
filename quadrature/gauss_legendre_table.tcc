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
// Ported from GSL by Jason Dick and Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// This file implements an Gauss-Legendre table for use in integration schemes
// Based upon gsl-2.3/integration/qelg.c

#ifndef GAUSS_LEGENDRE_TABLE_TCC
#define GAUSS_LEGENDRE_TABLE_TCC 1

namespace __gnu_test
{

  template<typename _Tp>
      gauss_legendre_table<_Tp>::gauss_legendre_table(std::size_t n)
      : order(n),
	point(nullptr),
	weight(nullptr),
	precomputed(false),
	i_precomp(-1),
	rule()
      {
	using __prec_t = decltype(gauss_legendre_precomp[0]);
	auto prec_beg = gauss_legendre_precomp;
	auto prec_end = gauss_legendre_precomp + num_gauss_legendre_precomp;
	auto prec = std::find_if(prec_beg, prec_end,
				 [this, n](__prec_t __tab)
				 { return __tab.order == this->order; });
	if (prec == prec_end)
	  this->rule = std::__detail::__legendre_zeros<_Tp>(this->order);
	else
	  {
	    this->precomputed = true;
	    this->i_precomp = prec - prec_beg;
	  }
      }

  /**
   * Routine to retrieve the i-th Gauss-Legendre point and weight.
   * Useful when the caller wishes to access the information stored in
   * the high-precision gsl_integration_glfixed_table struct.  Points
   * are indexed and presented in increasing order to the caller.
   */
  template<typename _Tp>
    std::tuple<_Tp, _Tp>
    gauss_legendre_table<_Tp>::get_point(_Tp a, _Tp b, size_t i) const
    {
      const auto A = (b - a) / _Tp{2};  /* Length of [a,b] */
      const auto B = (a + b) / _Tp{2};  /* Midpoint of [a,b] */

      if (i >= this->order)
	std::__throw_domain_error("i must be less than n");

      // See comments above gsl_integration_glfixed for struct's x, w layout.
      // Simply unpack that layout into a sorted set of points, weights.
      _Tp xi, wi;
      if (this->order & 1) // n is odd
	{
	  const auto k = int(i) - int(this->order) / 2;
	  const auto sign = (k < 0 ? -1 : +1);

	  xi = B + sign * A * this->pt(sign * k);
	  wi =		  A * this->wt(sign * k);
	}
      else if (/* n is even && */ i < this->order / 2)
	{
	  i = int(this->order) / 2 - 1 - int(i);
	  xi = B - A * this->pt(i);
	  wi =     A * this->wt(i);
	}
      else // n is even and i >= n / 2
	{
	  i  -= this->order / 2;
	  xi = B + A * this->pt(i);
	  wi =     A * this->wt(i);
	}

      return std::make_tuple(xi, wi);
    }

  template<typename _Tp>
    _Tp
    gauss_legendre_table<_Tp>::pt(size_t i) const
    {
      if (this->precomputed)
	if (this->i_precomp == std::size_t(-1))
	  return this->point[i];
	else
	  return _Tp(gauss_legendre_precomp[this->i_precomp].point[i]);
      else
	return this->rule[i + this->order / 2].__zero;
    }

  template<typename _Tp>
    _Tp
    gauss_legendre_table<_Tp>::wt(size_t i) const
    {
      if (this->precomputed)
	if (this->i_precomp == std::size_t(-1))
	  return this->weight[i];
	else
	  return _Tp(gauss_legendre_precomp[this->i_precomp].weight[i]);
      else
	return this->rule[i + this->order / 2].__weight;
    }

} // namespace __gnu_test

#endif // GAUSS_LEGENDRE_TABLE_TCC
