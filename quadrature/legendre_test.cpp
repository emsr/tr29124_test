// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016 Free Software Foundation, Inc.
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

#include <iostream>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <sstream>
#include <string>

//#include "simple_integrate.h"
#include "factorial_table.h"
#include "integration.h"

using namespace __gnu_test;

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_legendre(int l1, int l2, _Tp x)
  {
    return (_Tp(l1 + l2 + 1) / _Tp{2})
	 * std::legendre(l1,x)
	 * std::legendre(l2,x);
  }

template<typename _Tp>
  _Tp
  delta(int l1, int l2)
  { return l1 == l2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_legendre()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();
    for (int l1 = 0; l1 <= 720; ++l1)
      {
	for (int l2 = 0; l2 <= l1; ++l2)
	  {
	    std::function<_Tp(_Tp)> func(std::bind(&normalized_legendre<_Tp>, l1, l2,
                                         std::placeholders::_1));
	    _Tp integ_precision = _Tp{1000} * eps;
	    _Tp comp_precision = _Tp{10} * integ_precision;
	    _Tp integration_result, integration_error;

	    typedef std::pair<_Tp&,_Tp&> ret_type;
	    ret_type{integration_result, integration_error}
        	= integrate(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});

            if (std::abs(delta<_Tp>(l1, l2) - integration_result) > comp_precision)
              {
        	std::stringstream ss;
        	ss.precision(-int(log10(eps)));
        	ss << "Integration failed at l1=" << l1 << ", l2=" << l2
        	   << ", returning result " << integration_result
        	   << " instead of the expected " << delta<_Tp>(l1, l2) << '\n';
        	throw std::logic_error(ss.str());
              }
	  }
	std::cout << "Integration successful for legendre polynomials up to l = " << l1
             << '\n';
      }
  }

int
main()
{
  try
    {
      test_legendre<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  try
    {
      test_legendre<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
