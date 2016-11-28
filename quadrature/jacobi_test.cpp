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
#include "integration.h"

using namespace __gnu_test;

// Try to manage the four-gamma ratio.
// alpha > -1, beta > -1.
template<typename _Tp>
  _Tp
  gamma_ratio(int n, _Tp alpha, _Tp beta)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(alpha);
    if (std::abs(_Tp(1) + alpha + beta) < _S_eps)
      return _Tp(0);
    else
      {
	auto gaman1 = std::tgamma(_Tp(1) + alpha);
	auto gambn1 = std::tgamma(_Tp(1) + beta);
	auto gamabn1 = std::tgamma(_Tp(1) + alpha + beta);
	auto fact = gaman1 * gambn1 / gamabn1;
	for (int k = 1; k <= n; ++k)
	  fact *= (_Tp(k) + alpha) * (_Tp(k) + beta)
		/ (_Tp(k) + alpha + beta) / _Tp(k);
	return fact;
      }
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_jacobi(int n1, int n2, _Tp alpha, _Tp beta, _Tp x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(x);
    if (std::abs(x - _Tp{1}) < _S_eps)
      return _Tp{0};
    else if (std::abs(x + _Tp{1}) < _S_eps)
      return _Tp{0};
    else
      {
	auto gam = gamma_ratio(n1, alpha, beta);
	auto norm = std::pow(_Tp{2}, _Tp{1} + alpha + beta)
		  * gam / (_Tp(2 * n1 + 1) + alpha + beta);
	return std::pow(_Tp{1} - x, alpha) * std::pow(_Tp{1} + x, beta)
	     * __gnu_cxx::jacobi(n1, alpha, beta, x)
	     * __gnu_cxx::jacobi(n2, alpha, beta, x) / norm;
      }
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_jacobi(_Tp alpha, _Tp beta)
  {
    const auto eps = __gnu_cxx::__epsilon(alpha + beta);

    for (int n1 = 0; n1 <= 720; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    std::function<_Tp(_Tp)>
	      func([n1, n2, alpha, beta](_Tp x)
		   -> _Tp
		   { return normalized_jacobi<_Tp>(n1, n2, alpha, beta, x); });
	    _Tp integ_precision = _Tp{1000} * eps;
	    _Tp comp_precision = _Tp{10} * integ_precision;
	    _Tp integration_result, integration_error;

	    typedef std::pair<_Tp&,_Tp&> ret_type;
	    ret_type{integration_result, integration_error}
        	= integrate_smooth(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});
//        	= integrate(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});

            if (std::abs(delta<_Tp>(n1, n2) - integration_result) > comp_precision)
              {
        	std::stringstream ss;
        	ss.precision(-int(log10(eps)));
        	ss << "Integration failed at n1=" << n1 << ", n2=" << n2
        	   << ", returning result " << integration_result
        	   << " instead of the expected " << delta<_Tp>(n1, n2) << '\n';
        	throw std::logic_error(ss.str());
              }
	  }
	std::cout << "Integration successful for jacobi polynomials up to n = " << n1
             << '\n';
      }
  }

int
main()
{
  try
    {
      test_jacobi<double>(0.5, 1.5);
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  try
    {
      test_jacobi<long double>(0.5L, 1.5L);
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
