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

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_hermite(int n1, int n2, _Tp x)
  {
    const auto
    _S_pi = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    auto lnorm = _Tp{0.5} * (log(_S_pi) + _Tp(n1 + n2) * log(_Tp{2})
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
  int
  test_hermite()
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = _Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int fail = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
	    auto func = [n1, n2](_Tp x)
			-> _Tp
			{ return norm_hermite(n1, n2, x); };

	    auto [result, error]
		  = __gnu_cxx::integrate_sinh_sinh(func, abs_prec, rel_prec);

	    if (std::abs(delta<_Tp>(n1, n2) - result) > cmp_prec)
	      ++fail;
	  }
      }
    return fail;
  }

int
main()
{
  VERIFY(0 == test_hermite<float>());
  VERIFY(0 == test_hermite<double>());
  VERIFY(0 == test_hermite<long double>());
}
