// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
// { dg-skip-if "no extensions in strict dialects" { *-*-* } { "-std=c++*" } }

// Copyright (C) 2019 Free Software Foundation, Inc.
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

#include <cmath>

#include <ext/integration.h>

#if defined(__TEST_DEBUG)
#  include <iostream>
#  define VERIFY(A) \
  if (!(A)) \
    { \
      std::cout << "line " << __LINE__ << '\n'; \
    }
#else
#  include <testsuite_hooks.h>
#endif
#include <specfun_testcase.h>

// Neumann's number
template<typename _Tp>
  _Tp
  epsilon(int m)
  { return m == 0 ? _Tp{2} : _Tp{1}; }

// Azimuthal integral of zernike product.


// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_zernike(int n1, int m1, int n2, int m2, _Tp rho)
  {
    const auto _S_eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto _S_eps = _S_eps_factor * std::numeric_limits<_Tp>::epsilon();
    const auto _S_2pi = _Tp{2} * __gnu_cxx::numbers::__pi_v<_Tp>;

    auto z1 = [n1, m1, rho](_Tp phi)
	      -> _Tp
	      { return __gnu_cxx::zernike(n1, m1, rho, phi); };
    auto z2 = [n2, m2, rho](_Tp phi)
	      -> _Tp
	      { return __gnu_cxx::zernike(n2, m2, rho, phi); };

    auto norm = _Tp{1} / std::sqrt(_Tp(2 * n1 + 2) * _Tp(2 * n2 + 2));
    auto fun = [n1, m1, rho, z1, z2, norm](_Tp phi)
		-> _Tp
		{ return rho * z1(phi) * z2(phi) / norm; };

    auto val
	= __gnu_cxx::integrate_tanh_sinh(fun, _Tp{0}, _Tp{_S_2pi},
			      _S_eps, _S_eps, 8);

    return _Tp{2} * val.__result / _S_2pi / epsilon<_Tp>(m1);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  int
  test_zernike(int m1, int m2)
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = _Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 3, 4, 5, 8, 9, 16, 17}};
    int fail = 0;

    for (const auto n1 : degree)
      {
	if ((n1 - m1) & 1)
	  continue;
	for (const auto n2 : degree)
	  {
	    //int m2 = m1;
	    if (m2 > n2)
	      continue;
	    if ((n2 - m2) & 1)
	      continue;

	    auto func = [n1, m1, n2, m2](_Tp x)
			-> _Tp
			{ return norm_zernike(n1, m1, n2, m2, x); };

	    auto [result, error]
		= __gnu_cxx::integrate_tanh_sinh(func, _Tp{0}, _Tp{1},
					abs_prec, rel_prec, 8);

	    if (std::abs(delta<_Tp>(n1, m1, n2, m2) - result) > cmp_prec)
	      ++fail;
	  }
      }
    return fail;
  }

int
main()
{
  VERIFY(0 == test_zernike<float>(1, 2));
  VERIFY(0 == test_zernike<double>(1, 2));
  VERIFY(0 == test_zernike<long double>(1, 2));
}