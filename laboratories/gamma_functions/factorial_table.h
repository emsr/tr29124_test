// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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
// Written by Jason Dick.
//
// Includes classes for handling factorial tables.
// Both plain factorials and the logarithms of factorials are supported.

#ifndef FACTORIAL_TABLE_H
#define FACTORIAL_TABLE_H 1

#include <vector>
#include <cmath>

namespace emsr
{

  // Look-up table for storing and obtaining factorials
  template<typename _Tp>
    class factorial_table
    {

    private:

      static std::vector<_Tp> s_ft;
      factorial_table() = delete;

    public:

      static _Tp
      get(std::size_t n)
      {
	if (n >= s_ft.size())
	  add_factorials(n);
	return s_ft[n];
      }

      static void
      add_factorials(std::size_t max_fact)
      {
	if (s_ft.size() == 0)
	  s_ft.push_back(_Tp(1));

	for (_Tp ii = s_ft.size(); ii <= _Tp(max_fact); ++ii)
	  s_ft.push_back(ii * s_ft[s_ft.size() - 1]);
      }
    };

  template<typename _Tp>
    std::vector<_Tp>
    factorial_table<_Tp>::s_ft = std::vector<_Tp>();

  // Look-up table for storing and obtaining logarithms of factorials
  template<typename _Tp>
    class lnfactorial_table
    {

    private:

      static std::vector<_Tp> s_lnft;
      lnfactorial_table() = delete;

    public:

      static _Tp
      get(std::size_t n)
      {
	if (n >= s_lnft.size())
	  add_factorials(n);
	return s_lnft[n];
      }

      static void
      add_factorials(std::size_t max_fact)
      {
	if (s_lnft.size() == 0)
	  s_lnft.push_back(_Tp(0));

	for (_Tp ii = s_lnft.size(); ii <= _Tp(max_fact); ++ii)
	  s_lnft.push_back(std::log(ii) + s_lnft[s_lnft.size() - 1]);
      }
    };

  template<typename _Tp>
    std::vector<_Tp>
    lnfactorial_table<_Tp>::s_lnft = std::vector<_Tp>();

  template<typename _Tp>
    inline _Tp
    factorial(std::size_t n)
    { return factorial_table<_Tp>::get(n); }

  inline float
  factorialf(std::size_t n)
  { return factorial<float>(n); }

  inline double
  factoriald(std::size_t n)
  { return factorial<double>(n); }

  inline long double
  factorialld(std::size_t n)
  { return factorial<long double>(n); }

  inline int
  factoriali(std::size_t n)
  { return factorial<int>(n); }

  inline long
  factoriall(std::size_t n)
  { return factorial<long>(n); }

  inline long long
  factorialll(std::size_t n)
  { return factorial<long long>(n); }

  template<typename _Tp>
    inline _Tp
    lnfactorial(std::size_t n)
    { return lnfactorial_table<_Tp>::get(n); }

  inline float
  lnfactorialf(std::size_t n)
  { return lnfactorial<float>(n); }

  inline double
  lnfactoriald(std::size_t n)
  { return lnfactorial<double>(n); }

  inline long double
  lnfactorialld(std::size_t n)
  { return lnfactorial<long double>(n); }

} // namespace

#endif
