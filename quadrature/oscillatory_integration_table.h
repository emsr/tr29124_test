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
// Ported from GSL by Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// This file implements an oscillatory integrand table for use
// in integration schemes
// Based upon gsl-2.3/integration/qmomof.c

#ifndef OSCILLATORY_INTEGRATION_TABLE_H
#define OSCILLATORY_INTEGRATION_TABLE_H 1

namespace __gnu_test
{

  template<typename _Tp>
    struct oscillatory_integration_table
    {
      enum circular_function
      {
	INTEG_COSINE,
	INTEG_SINE
      };

      std::size_t n;
      _Tp omega;
      _Tp length;
      _Tp par;
      enum circular_function circfun;
      std::vector<_Tp> chebmo;

      oscillatory_integration_table(_Tp omega_in, _Tp length_in,
				    circular_function circfun_in,
				    std::size_t n_in)
      : n(n_in),
	omega(omega_in),
	length(length_in),
	par(0.5 * omega_in * length_in),
	circfun(circfun_in),
	chebmo(25 * n_in)
      {
	auto __scale = _Tp{1};
	for (auto __i = 0u; __i < this->n; ++__i)
	  {
	    this->compute_moments(this->par * __scale, __i);
	    __scale *= 0.5;
	  }
      }

      _Tp
      get_length() const
      { return this->length; }

      oscillatory_integration_table<_Tp>&
      set_length(_Tp length_in)
      {
	this->length = length_in;
	this->par = 0.5 * this->omega * this->length;
	this->reset(this->omega, this->length, this->circfun);
	return *this;
      }

      void
      reset(_Tp omega_in, _Tp length_in,
	    circular_function circfun_in)
      {
	this->omega = omega_in;
	this->length = length_in;
	this->par = 0.5 * omega_in * length_in;
	this->circfun = circfun_in;
	auto __scale = _Tp{1};
	for (auto __i = 0u; __i < this->n; ++__i)
	  {
	    this->compute_moments(this->par * __scale, __i);
	    __scale *= 0.5;
	  }
      }

      inline const _Tp*
      get_moments(std::size_t __level) const
      { return this->chebmo.data() + 25 * __level; }

    private:

      void compute_moments(_Tp par, std::size_t __level);

    };

} // namespace __gnu_test

#include "oscillatory_integration_table.tcc"

#endif // OSCILLATORY_INTEGRATION_TABLE_H
