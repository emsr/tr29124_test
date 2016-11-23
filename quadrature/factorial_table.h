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
// Written by Jason Dick.
//
// Includes classes for handling factorial tables.
// Both plain factorials and the logarithms of factorials are supported.

#ifndef FACTORIAL_TABLE_H
#define FACTORIAL_TABLE_H

#include <vector>
#include <cmath>

namespace __gnu_test
{

  // Look-up table for storing and obtaining factorials
  template<typename _Tp>
    class factorial_table
    {

    private:

      static std::vector<_Tp> _S_ft;
      factorial_table(); //Disallow construction

    public:

      static _Tp
      get(std::size_t __n)
      {
	if (__n >= _S_ft.size())
	  add_factorials(__n);
	return _S_ft[__n];
      }

      static void
      add_factorials(std::size_t __max_fact)
      {
	if (_S_ft.size() == 0)
	  _S_ft.push_back(_Tp(1));

	for (_Tp __ii = _S_ft.size(); __ii <= _Tp(__max_fact); ++__ii)
	  _S_ft.push_back(__ii * _S_ft[_S_ft.size() - 1]);
      }
    };

  template<typename _Tp>
    std::vector<_Tp>
    factorial_table<_Tp>::_S_ft = std::vector<_Tp>();

  // Look-up table for storing and obtaining logarithms of factorials
  template<typename _Tp>
    class lnfactorial_table
    {

    private:

      static std::vector<_Tp> _S_lnft;
      lnfactorial_table(); //Disallow construction

    public:

      static _Tp
      get(std::size_t __n)
      {
	if (__n >= _S_lnft.size())
	  add_factorials(__n);
	return _S_lnft[__n];
      }

      static void
      add_factorials(std::size_t __max_fact)
      {
	if (_S_lnft.size() == 0)
	  _S_lnft.push_back(_Tp(0));

	for (_Tp __ii = _S_lnft.size(); __ii <= _Tp(__max_fact); ++__ii)
	  _S_lnft.push_back(std::log(__ii) + _S_lnft[_S_lnft.size() - 1]);
      }
    };

  template<typename _Tp>
    std::vector<_Tp>
    lnfactorial_table<_Tp>::_S_lnft = std::vector<_Tp>();

  template<typename _Tp>
    inline _Tp
    factorial(std::size_t __n)
    { return factorial_table<_Tp>::get(__n); }

  inline float
  factorialf(std::size_t __n)
  { return factorial<float>(__n); }

  inline double
  factoriald(std::size_t __n)
  { return factorial<double>(__n); }

  inline long double
  factorialld(std::size_t __n)
  { return factorial<long double>(__n); }

  inline int
  factoriali(std::size_t __n)
  { return factorial<int>(__n); }

  inline long
  factoriall(std::size_t __n)
  { return factorial<long>(__n); }

  inline long long
  factorialll(std::size_t __n)
  { return factorial<long long>(__n); }

  template<typename _Tp>
    inline _Tp
    lnfactorial(std::size_t __n)
    { return lnfactorial_table<_Tp>::get(__n); }

  inline float
  lnfactorialf(std::size_t __n)
  { return lnfactorial<float>(__n); }

  inline double
  lnfactoriald(std::size_t __n)
  { return lnfactorial<double>(__n); }

  inline long double
  lnfactorialld(std::size_t __n)
  { return lnfactorial<long double>(__n); }

} // namespace

#endif
