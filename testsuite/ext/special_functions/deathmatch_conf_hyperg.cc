// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
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

//  conf_hyperg
//  Test relations between special functions.

#include <limits>
#include <ext/cmath>

#if defined(__TEST_DEBUG)
#  include <iostream>
#  define VERIFY(A) \
  if (!(A)) \
    { \
      std::cout << "line " << __LINE__ \
	<< "  max_abs_frac = " << stats.max_abs_frac \
	<< std::endl; \
    } \
  else \
    { \
      std::cout << "Success: " \
	<< "  samples = " << stats.num \
	<< "  mean = " << stats.mean_diff() \
	<< "  variance = " << stats.variance_diff() \
	<< std::endl; \
    }
#else
#  include <testsuite_hooks.h>
#endif

#include <specfun_stats.h>

template<typename Tp>
  void
  test_erf(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    const auto _S_sqrt_pi = __gnu_cxx::__math_constants<Tp>::__root_pi;
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = Tp{2} * z * __gnu_cxx::conf_hyperg(Tp{0.5L}, Tp{1.5L}, -z * z) / _S_sqrt_pi;
	auto erf = std::erf(z);
	stats << std::make_pair(hyp, erf);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test_kummer_xform_m(Tp a, Tp c, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto M = __gnu_cxx::conf_hyperg(a, c, z);
	auto Mt = std::exp(z) * __gnu_cxx::conf_hyperg(c - a, c, -z);
	stats << std::make_pair(M, Mt);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test_kummer_xform_u(Tp a, Tp c, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto U = __gnu_cxx::tricomi_u(a, c, z);
	auto Ut = std::pow(z, Tp{1} - c)
		* __gnu_cxx::tricomi_u(Tp{1} + a - c, Tp{1} - c, z);
	stats << std::make_pair(U, Ut);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    test_erf<Tp>(toler);

    test_kummer_xform_m(Tp{1.5L}, Tp{0.5L}, toler);
    test_kummer_xform_m(Tp{2.0L}, Tp{2.5L}, toler);
    test_kummer_xform_m(Tp{3.5L}, Tp{0.5L}, toler);
    test_kummer_xform_m(Tp{5.0L}, Tp{2.5L}, toler);

    test_kummer_xform_u(Tp{1.5L}, Tp{0.5L}, toler);
    test_kummer_xform_u(Tp{2.0L}, Tp{2.5L}, toler);
    test_kummer_xform_u(Tp{3.5L}, Tp{0.5L}, toler);
    test_kummer_xform_u(Tp{5.0L}, Tp{2.5L}, toler);
  }

int
main()
{
  test<double>();
}
