// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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
#include <stdexcept>
#include <sstream>
#include <string>

#include <ext/integration.h>
#include <ext/math_constants.h>

using namespace __gnu_cxx;

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_gegenbauer(int n1, int n2, _Tp lambda, _Tp x)
  {
    const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Tp>;
    auto gama = std::tgamma(lambda);
    auto gamn2a = std::tgamma(n1 + _Tp{2} * lambda);
    auto norm = _S_pi * std::pow(_Tp{2}, _Tp{1} - _Tp{2} * lambda) * gamn2a
	      / __gnu_cxx::factorial<_Tp>(n1) / (_Tp(n1) + lambda) / gama / gama;
    return std::pow(_Tp{1} - x * x, lambda - _Tp{0.5})
	 * __gnu_cxx::gegenbauer(n1, lambda, x)
	 * __gnu_cxx::gegenbauer(n2, lambda, x) / norm;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_gegenbauer(_Tp lambda)
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = _Tp{10} * rel_precision;

    const bool singular = (lambda < _Tp{0.5});

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    auto func = [n1, n2, lambda](_Tp x)
			-> _Tp
			{ return normalized_gegenbauer<_Tp>(n1, n2, lambda, x); };

	    // Using integrate_singular works pretty well.
	    auto [result, error]
		= singular
		? integrate_singular_endpoints(func,
					       _Tp{-1}, _Tp{1},
					       lambda - _Tp{0.5}, lambda - _Tp{0.5}, 0, 0,
					       abs_precision, rel_precision)
		: integrate(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision);
		//= integrate_singular(func, _Tp{-1}, _Tp{1}, rel_precision, _Tp{0});

	    if (std::abs(delta<_Tp>(n1, n2) - result) > cmp_precision)
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
	std::cout << "Integration successful for gegenbauer polynomials up to n = " << n1
		  << '\n' << std::flush;
      }

    int ibot = n1 - 1;
    int itop = 2 * ibot;
    int del = 2;
    bool breakout = false;
    while (itop != ibot)
      {
	RESTART:
	for (int n2 = 0; n2 <= itop; n2 += del)
	  {
	    auto func = [n1 = itop, n2, lambda](_Tp x)
			-> _Tp
			{ return normalized_gegenbauer<_Tp>(n1, n2, lambda, x); };

	    // Using integrate_singular works pretty well.
	    auto [result, error]
		= singular
		? integrate_singular_endpoints(func,
					       _Tp{-1}, _Tp{1},
					       lambda - _Tp{0.5}, lambda - _Tp{0.5}, 0, 0,
					       abs_precision, rel_precision)
		: integrate(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision);
//		= integrate_singular(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision);

	    if (std::abs(delta<_Tp>(itop, n2) - result) > cmp_precision)
	      {
		if ((ibot + itop) / 2 < itop)
		  {
		    itop = (ibot + itop) / 2;
		    goto RESTART;
		  }
		else
		  {
		    breakout = true;
		    break;
		  }
	      }
	  }

	std::cout << "Integration successful for gegenbauer polynomials up to n = " << itop
		  << '\n' << std::flush;

	if (breakout)
	  break;

	ibot = itop;
	if (itop > 1000)
	  {
	    std::cout << "\nGood enough!\n" << std::flush;
	    break;
	  }
	else if (itop <= std::numeric_limits<int>::max() / 2)
	  itop *= 2;
	else
	  break;
        del *= 2;
      }
  }

int
main()
{
  std::cout << "\n\nOrthonormality tests for float\n";
  try
    {
      test_gegenbauer<float>(0.25F);
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for float\n";
  try
    {
      test_gegenbauer<float>(0.5F);
    }
  catch (__gnu_cxx::__integration_error<float>& ierr)
    {
      std::cerr << ierr.what() << '\n';
      std::cerr << " result = " << ierr.result() << " abserr = " << ierr.abserr() << '\n';
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for double\n";
  try
    {
      test_gegenbauer<double>(0.5);
    }
  catch (__gnu_cxx::__integration_error<double>& ierr)
    {
      std::cerr << ierr.what() << '\n';
      std::cerr << " result = " << ierr.result() << " abserr = " << ierr.abserr() << '\n';
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for long double\n";
  try
    {
      test_gegenbauer<long double>(0.5L);
    }
  catch (__gnu_cxx::__integration_error<long double>& ierr)
    {
      std::cerr << ierr.what() << '\n';
      std::cerr << " result = " << ierr.result() << " abserr = " << ierr.abserr() << '\n';
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
