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

/**
 * Neumann's number
 */
template<typename _Tp>
  _Tp
  epsilon(int m)
  { return m == 0 ? _Tp{2} : _Tp{1}; }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_radpoly(int n1, int m1, int n2, int m2, _Tp rho)
  {
    auto norm = _Tp{1} / std::sqrt(_Tp(2 * n1 + 2) * _Tp(2 * n2 + 2));
    return rho
	 * __gnu_cxx::radpoly(n1, m1, rho)
	 * __gnu_cxx::radpoly(n2, m2, rho) / norm;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_radpoly()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();

    // Neverending loop: runs until integration fails
    for (int n1 = 0; ; ++n1)
      {
	for (int m1 = 0; m1 <= n1; ++m1)
	  {
	    if ((n1 - m1) & 1)
	      continue;
	    for (int n2 = 0; n2 <= n1; ++n2)
	      {
		for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    if ((n2 - m2) & 1)
		      continue;
		    std::function<_Tp(_Tp)> func(std::bind(&normalized_radpoly<_Tp>,
						 n1, m1, n2, m2,
                                        	 std::placeholders::_1));
		    _Tp integ_precision = _Tp{1000} * eps;
		    _Tp comp_precision = _Tp{10} * integ_precision;

		    auto [result, error]
			= integrate(func, _Tp{0}, _Tp{1}, integ_precision, _Tp{0});

		    if (std::abs(delta<_Tp>(n1, m1, n2, m2) - result) > comp_precision)
		      {
        		std::stringstream ss;
        		ss.precision(std::numeric_limits<_Tp>::digits10);
			auto w = 8 + ss.precision();
			ss << std::showpoint << std::scientific;
        		ss << "Integration failed at n1=" << n1 << ", m1=" << m1
        		   << ", n2=" << n2 << ", m2=" << m2
        		   << ", returning result = " << std::setw(w) << result
        		   << ", with error = " << std::setw(w) << error
        		   << " instead of the expected " << delta<_Tp>(n1, m1, n2, m2) << '\n';
			std::cerr << ss.str();
        		//throw std::logic_error(ss.str());
		      }
		  }
	      }
	  }
	std::cout << "Integration successful for radpoly polynomials up to n = " << n1
	     << '\n';
      }
  }

int
main()
{
  try
    {
      test_radpoly<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  try
    {
      test_radpoly<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}