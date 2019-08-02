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

using namespace __gnu_cxx;

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_assoc_legendre(int l1, int m1, int l2, int m2, _Tp x)
  {
    return (_Tp(l1 + l2 + 1) / _Tp{2})
	 * std::assoc_legendre(l1, m1, x)
	 * std::assoc_legendre(l2, m2, x);
  }

template<typename _Tp>
  _Tp
  delta(int l1, int l2)
  { return l1 == l2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_assoc_legendre(int m1, int m2)
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = _Tp{10} * rel_precision;

    int l1 = 0;
    for (; l1 <= 128; ++l1)
      {
	if (m1 > l1)
	  continue;

	for (int l2 = 0; l2 <= l1; ++l2)
	  {
	    if (m2 > l2)
	      continue;

	    auto func = [l1, m1, l2, m2](_Tp x)
			-> _Tp
			{ return normalized_assoc_legendre(l1, m1, l2, m2, x); };

	    auto [result, error]
//		= integrate(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision);
		= integrate_tanh_sinh(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision, 6);

	    if (std::abs(delta<_Tp>(l1, l2) * delta<_Tp>(m1, m2) - result) > cmp_precision)
	      {
		std::stringstream ss;
		ss.precision(std::numeric_limits<_Tp>::digits10);
		ss << std::showpoint << std::scientific;
		ss << "Integration failed at l1=" << l1 << ", l2=" << l2
		   << ", returning result " << result
		   << ", with error " << error
		   << " instead of the expected " << delta<_Tp>(l1, l2) << '\n';
		throw std::logic_error(ss.str());
	      }
	  }
	std::cout << "Integration successful for assoc_legendre polynomials up to l = " << l1
		  << '\n' << std::flush;
      }

    int ibot = l1 - 1;
    int itop = 2 * ibot;
    int del = 2;
    bool breakout = false;
    while (itop != ibot)
      {
	RESTART:
	for (int l2 = 0; l2 <= itop; l2 += del)
	  {
	    auto func = [l1 = itop, m1, l2, m2](_Tp x)
			-> _Tp
			{ return normalized_assoc_legendre(l1, m1, l2, m2, x); };

	    auto [result, error]
//		= integrate(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision);
		= integrate_tanh_sinh(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision, 6);

	    if (std::abs(delta<_Tp>(itop, l2) * delta<_Tp>(m1, m2) - result) > cmp_precision)
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

	std::cout << "Integration successful for assoc_legendre polynomials up to l = " << itop
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
      test_assoc_legendre<float>(0, 0);
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
      test_assoc_legendre<double>(0, 0);
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
      test_assoc_legendre<long double>(0, 0);
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
