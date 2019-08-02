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

// Try to manage the four-gamma ratio.
// alpha > -1, beta > -1.
template<typename _Tp>
  _Tp
  gamma_ratio(int n, _Tp alpha, _Tp beta)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
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
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
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
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = _Tp{10} * rel_precision;

    const bool singular = (alpha < _Tp{0} || beta < _Tp{0});

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    auto func = [n1, n2, alpha, beta](_Tp x)
			-> _Tp
			{ return normalized_jacobi<_Tp>(n1, n2, alpha, beta, x); };

	    auto [result, error]
		= singular
		? integrate_singular_endpoints(func,
					       _Tp{-1}, _Tp{1},
					       alpha, beta, 0, 0,
					       abs_precision, rel_precision)
		: integrate(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision);

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
	std::cout << "Integration successful for jacobi polynomials up to n = " << n1
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
	    auto func = [itop, n2, alpha, beta](_Tp x)
			-> _Tp
			{ return normalized_jacobi<_Tp>(itop, n2, alpha, beta, x); };

	    auto [result, error]
		= singular
		? integrate_singular_endpoints(func,
					       _Tp{-1}, _Tp{1},
					       alpha, beta, 0, 0,
					       abs_precision, rel_precision)
		: integrate(func, _Tp{-1}, _Tp{1}, abs_precision, rel_precision);

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

	std::cout << "Integration successful for jacobi polynomials up to n = " << itop
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
      test_jacobi<float>(0.5F, 1.5F);
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
      test_jacobi<double>(0.5, 1.5);
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
      test_jacobi<long double>(0.5L, 1.5L);
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
