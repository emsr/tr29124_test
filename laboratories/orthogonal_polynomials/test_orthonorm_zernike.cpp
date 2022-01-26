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
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include <emsr/integration.h>
#include <emsr/specfun.h>

// Neumann's number
template<typename _Tp>
  _Tp
  epsilon(int m)
  { return m == 0 ? _Tp{2} : _Tp{1}; }

// Azimuthal integral of zernike product.


// Function which should integrate to 1 for n1 == n2, 0 otherwise.
// This does the angular integral.
// FIXME: Try an FFTish integral here. Fail.
template<typename _Tp>
  _Tp
  normalized_zernike(int n1, int m1, int n2, int m2, _Tp rho)
  {
    const auto _S_eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto _S_eps = _S_eps_factor * std::numeric_limits<_Tp>::epsilon();
    const auto _S_2pi = _Tp{2} * emsr::pi_v<_Tp>;

    auto z1 = [n1, m1, rho](_Tp phi)
	      -> _Tp
	      { return emsr::zernike(n1, m1, rho, phi); };
    auto z2 = [n2, m2, rho](_Tp phi)
	      -> _Tp
	      { return emsr::zernike(n2, m2, rho, phi); };

    auto norm = _Tp{1} / std::sqrt(_Tp(2 * n1 + 2) * _Tp(2 * n2 + 2));
    auto fun = [n1, m1, rho, z1, z2, norm](_Tp phi)
		-> _Tp
		{ return rho * z1(phi) * z2(phi) / norm; };

    auto val
	= emsr::integrate_tanh_sinh(fun, _Tp{0}, _Tp{_S_2pi},
			      _S_eps, _S_eps, 8);
	//= emsr::integrate_oscillatory(fun, _Tp{0}, _Tp{_S_2pi},
	//			_S_eps, _S_eps, 1024);

    return _Tp{2} * val.result / _S_2pi / epsilon<_Tp>(m1);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_zernike()
  {
    bool full_range = false;

    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_precision = eps_factor * eps;
    const auto rel_precision = eps_factor * eps;
    const auto cmp_precision = _Tp{10} * rel_precision;

    std::vector<int> degree{0, 1, 2, 3, 4, 5, 8, 9, 16, 17};

    int n1_save = 0;
    for (int n1 : degree)
      {
	for (int m1 : degree)
	  {
	    if (m1 > n1 || (n1 - m1) & 1)
	      continue;
	    for (int n2 : degree)
	      {
		if (n2 > n1)
		  continue;
		for (int m2 : degree)
		  {
		    if (m2 > n2 || (n2 - m2) & 1)
		      continue;
		    auto func = [n1, m1, n2, m2](_Tp x)
				-> _Tp
				{ return normalized_zernike(n1, m1, n2, m2, x); };

		    auto [result, error]
			= emsr::integrate_tanh_sinh(func, _Tp{0}, _Tp{1},
					      abs_precision, rel_precision, 8);

		    if (std::abs(delta<_Tp>(n1, m1, n2, m2) - result) > cmp_precision)
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
		      }
		  }
	      }
	  }
	n1_save = n1;
	std::cout << "Integration successful for zernike polynomials up to n = " << n1_save
		  << '\n' << std::flush;
      }

    if (!full_range)
      return;

    int n1_lower = n1_save;
    int n1_upper = 2 * n1_lower;
    int del = 2;
    bool breakout = false;
    while (n1_upper != n1_lower)
      {
	RESTART:
	for (int m1 = 0; m1 <= n1_upper; ++m1)
	  {
	    if (m1 > n1_upper || (n1_upper - m1) & 1)
	      continue;
	    for (int n2 = 0; n2 <= n1_upper; n2 += del)
	      {
		for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    if (m2 > n2 || (n2 - m2) & 1)
		      continue;
		    auto func = [n1 = n1_upper, m1, n2, m2](_Tp x)
				-> _Tp
				{ return normalized_zernike(n1, m1, n2, m2, x); };

		    auto [result, error]
			= emsr::integrate_tanh_sinh(func, _Tp{0}, _Tp{1},
					      abs_precision, rel_precision, 8);

		    if (std::abs(delta<_Tp>(n1_upper, m1, n2, m2) - result) > cmp_precision)
		      {
			if ((n1_lower + n1_upper) / 2 < n1_upper)
			  {
			    n1_upper = (n1_lower + n1_upper) / 2;
			    goto RESTART;
			  }
			else
			  {
			    breakout = true;
			    break;
			  }
		      }
		  }
	      }
	    if (breakout)
	      break;
	  }

	std::cout << "Integration successful for zernike polynomials up to n = " << n1_upper
		  << '\n' << std::flush;

	if (breakout)
	  break;

	n1_lower = n1_upper;
	if (n1_upper > 1000)
	  {
	    std::cout << "\nGood enough!\n" << std::flush;
	    break;
	  }
	else if (n1_upper <= std::numeric_limits<int>::max() / 2)
	  n1_upper *= 2;
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
      test_zernike<float>();
    }
  catch (emsr::integration_error<float, float>& ierr)
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
      test_zernike<double>();
    }
  catch (emsr::integration_error<double, double>& ierr)
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
      test_zernike<long double>();
    }
  catch (emsr::integration_error<long double, long double>& ierr)
    {
      std::cerr << ierr.what() << '\n';
      std::cerr << " result = " << ierr.result() << " abserr = " << ierr.abserr() << '\n';
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
