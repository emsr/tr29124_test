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
//Provides test example of using integration suite with Hermite polynomials
//Written by Jason Dick

#include <iostream>
#include <cmath>
#include <tr1/cmath>
#include <functional>
#include <stdexcept>
#include <sstream>
#include <string>

//#include "simple_integrate.h"
#include "factorial_table.h"
#include "integration.h"

using namespace __gnu_test;

//Function which should integrate to 1 for n1=n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_hermite(int n1, int n2, _Tp x)
  {
    _Tp lnnorm_fact = _Tp{0.5} * (log(M_PI) + (n1 + n2) * log(_Tp{2})
                       + lnfactorialld(n1) + lnfactorialld(n2));
    return std::tr1::hermite(n2,x) * std::exp(-x*x - lnnorm_fact) * std::tr1::hermite(n1,x);
  }

template<typename _Tp>
  void
  test_hermite()
  {
    const _Tp infty = std::numeric_limits<_Tp>::infinity();

    // Neverending loop: runs until integration fails
    for (int n1 = 0; ; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    std::function<_Tp(_Tp)> func(std::bind(&normalized_hermite<_Tp>, n1, n2,
                                         std::placeholders::_1));
	    _Tp integ_precision = _Tp{1000} * std::numeric_limits<_Tp>::epsilon();
	    _Tp comp_precision = _Tp{10} * integ_precision;
	    _Tp integration_result, integration_error;

	    typedef std::pair<_Tp&,_Tp&> ret_type;
	    ret_type{integration_result, integration_error}
        	= integrate(func, -infty, infty, integ_precision, _Tp{0});

	    if (n1 == n2)
	      {
        	//integration_result should be 1, throw error if differs from one by
        	//more than integration precision (with an additional fudge factor in
        	//case integration isn't quite that accurate)
        	if (fabs(_Tp{1} - integration_result) > comp_precision)
        	  {
        	    std::stringstream ss;
        	    ss.precision(-int(log10(std::numeric_limits<_Tp>::epsilon())));
        	    ss << "Integration failed at n1=" << n1 << ", n2=" << n2
        	       << ", returning result " << integration_result
        	       << " instead of the expected 1" << std::endl;
        	    throw std::logic_error(ss.str());
        	  }
	      }
	    else
	      {
        	//integration_result should be 0, throw error if differs from zero by
        	//more than integration precision (with an additional fudge factor in
        	//case integration isn't quite that accurate)
        	if (fabs(integration_result) > comp_precision)
        	  {
        	    std::stringstream ss;
        	    ss.precision(14);
        	    ss << "Integration failed at n1=" << n1 << ", n2=" << n2
        	       << ", returning result " << integration_result
        	       << " instead of the expected 0" << std::endl;
        	    throw std::logic_error(ss.str());
        	  }
	      }
	  }
	std::cout << "Integration successful for hermite polynomials up to n=" << n1
             << std::endl;
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
