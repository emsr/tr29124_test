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

#include <iostream>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <sstream>
#include <string>

#include "integration.h"

using namespace __gnu_test;

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_hermite(int n1, int n2, _Tp x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(x);
    _Tp lnorm = _Tp{0.5} * (log(_S_pi) + (n1 + n2) * log(_Tp{2})
			  + __gnu_cxx::lfactorial<_Tp>(n1)
			  + __gnu_cxx::lfactorial<_Tp>(n2));
    return std::hermite(n2, x)
	 * std::exp(-x * x - lnorm)
	 * std::hermite(n1, x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_hermite()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp infty = std::numeric_limits<_Tp>::infinity();

    // Neverending loop: runs until integration fails
    for (int n1 = 0; ; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    std::function<_Tp(_Tp)> func(std::bind(&normalized_hermite<_Tp>, n1, n2,
                                         std::placeholders::_1));
	    _Tp integ_precision = _Tp{1000} * eps;
	    _Tp comp_precision = _Tp{10} * integ_precision;

	    auto [result, error]
        	= integrate(func, -infty, infty, integ_precision, _Tp{0});

            if (std::abs(delta<_Tp>(n1, n2) - result) > comp_precision)
              {
        	std::stringstream ss;
        	ss.precision(std::numeric_limits<_Tp>::digits10);
		ss << std::showpoint << std::scientific;
        	ss << "Integration failed at n1=" << n1 << ", n2=" << n2
        	   << ", returning result " << result
        	   << ", with error " << error
        	   << " instead of the expected " << delta<_Tp>(n1, n2) << '\n';
        	throw std::logic_error(ss.str());
              }
	  }
	std::cout << "Integration successful for hermite polynomials up to n = " << n1
             << '\n';
      }
  }

int
main()
{
  try
    {
      test_hermite<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  try
    {
      test_hermite<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
