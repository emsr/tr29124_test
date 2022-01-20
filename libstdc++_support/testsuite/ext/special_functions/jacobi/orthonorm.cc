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

#include <emsr/integration.h>

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

// Try to manage the four-gamma ratio.
// alpha > -1, beta > -1.
template<typename _Tp>
  _Tp
  gamma_ratio(int n, _Tp alpha, _Tp beta)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    if (std::abs(_Tp{1} + alpha + beta) < _S_eps)
      return _Tp{0};
    else
      {
	auto gaman1 = std::tgamma(_Tp{1} + alpha);
	auto gambn1 = std::tgamma(_Tp{1} + beta);
	auto gamabn1 = std::tgamma(_Tp{1} + alpha + beta);
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
  norm_jacobi(int n1, int n2, _Tp alpha, _Tp beta, _Tp x)
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
  int
  test_jacobi(_Tp alpha, _Tp beta)
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = _Tp{10} * rel_prec;

    const bool singular = (alpha < _Tp{0} || beta < _Tp{0});

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int fail = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
	    auto func = [n1, n2, alpha, beta](_Tp x)
			-> _Tp
			{ return norm_jacobi<_Tp>(n1, n2, alpha, beta, x); };

	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func,
					       _Tp{-1}, _Tp{1},
					       alpha, beta, 0, 0,
					       abs_prec, rel_prec)
		: emsr::integrate_tanh_sinh(func, _Tp{-1}, _Tp{1},
						 abs_prec, rel_prec);

	    if (std::abs(delta<_Tp>(n1, n2) - result) > cmp_prec)
	      ++fail;
	  }
      }
    return fail;
  }

int
main()
{
  VERIFY(0 == test_jacobi<float>(0.5F, 1.5F));
  VERIFY(0 == test_jacobi<double>(0.5, 1.5));
  VERIFY(0 == test_jacobi<long double>(0.5L, 1.5L));
}
