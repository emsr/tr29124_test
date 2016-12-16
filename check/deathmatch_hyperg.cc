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

//  hyperg
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

#include <specfun_testcase.h>
#include <specfun_stats.h>

template<typename Tp>
  void
  test_log(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = z * __gnu_cxx::hyperg(Tp{1}, Tp{1}, Tp{2}, -z);
	auto ln1p = std::log1p(z);
	stats << std::make_pair(hyp, ln1p);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test_pow(Tp a, Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = __gnu_cxx::hyperg(a, Tp{1}, Tp{1}, z);
	auto binom = std::pow(Tp{1} - z, -a);
	stats << std::make_pair(hyp, binom);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test_asin(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = z * __gnu_cxx::hyperg(Tp{0.5L}, Tp{0.5L}, Tp{1.5L}, z * z);
	auto arcsin = std::asin(z);
	stats << std::make_pair(hyp, arcsin);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test_comp_ellint_1(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    const auto _S_pi2 = __gnu_cxx::__math_constants<Tp>::__pi_half;
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto k = Tp{0.1L} * i;
	auto hyp = _S_pi2 * __gnu_cxx::hyperg(Tp{0.5L}, Tp{0.5L}, Tp{1}, k * k);
	auto K = std::comp_ellint_1(k);
	stats << std::make_pair(hyp, K);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test_comp_ellint_2(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    const auto _S_pi2 = __gnu_cxx::__math_constants<Tp>::__pi_half;
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = 1; i < 10; ++i)
      {
	auto k = Tp{0.1L} * i;
	auto hyp = _S_pi2 * __gnu_cxx::hyperg(Tp{-0.5L}, Tp{0.5L}, Tp{1}, k * k);
	auto E = std::comp_ellint_2(k);
	stats << std::make_pair(hyp, E);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test_jacobi(unsigned int n, Tp alpha, Tp beta,
	      Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    bool test __attribute__((unused)) = true;
    Stats<Tp> stats(toler);
    for (auto i = -9; i < 10; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto hyp = __gnu_cxx::hyperg(-Tp(n), Tp(1 + n) + alpha + beta, Tp{1} + alpha, z);
	auto P = __gnu_cxx::factorial<Tp>(n)
	       * __gnu_cxx::jacobi(n, alpha, beta, Tp{1} - Tp{2} * z)
	       / __gnu_cxx::pochhammer(alpha + Tp{1}, int(n));
	stats << std::make_pair(hyp, P);
      }
    VERIFY(stats.max_abs_frac < toler);
  }

template<typename Tp>
  void
  test(Tp toler = 100 * std::numeric_limits<Tp>::epsilon())
  {
    test_log<Tp>(toler);

    test_pow(Tp{0.5L}, toler);
    test_pow(Tp{1.5L}, toler);
    test_pow(Tp{5.0L}, toler);

    test_asin<Tp>(toler);

    test_comp_ellint_1<Tp>(toler);
    test_comp_ellint_2<Tp>(toler);

    test_jacobi(0, Tp{0.5L}, Tp{0.5L}, toler);
    test_jacobi(3, Tp{0.5L}, Tp{0.5L}, toler);
    test_jacobi(3, Tp{1.5L}, Tp{4.0L}, toler);
  }

int
main()
{
  test<double>();
}
