/* integration/test.c
 *
 * Copyright (C) 1996-2000, 2007 Brian Gough
 * Copyright (C) 2016-2017 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
$HOME/bin/bin/g++ -std=gnu++17 -fconcepts -g -Wall -Wextra -Wno-psabi -I.. -c -o obj/test_quadrature..o test_quadrature.cpp
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>

#include "integration.h"
#include "testcase.h"

/**
 *
 */
template<typename _Tp, typename... _Parms>
  std::function<_Tp(_Tp)>
  make_function(_Tp(f)(_Tp, _Parms...), _Parms... p)
  {
    return [f, p...](_Tp x)->_Tp{ return f(x, p...); };
  }

/**
 *
 */
template<typename _Tp>
  struct counted_function
  {
    counted_function(std::function<_Tp(_Tp)> f)
    : func(f), neval(new int(0))
    { }

    _Tp
    operator()(_Tp x) const
    {
      ++(*this->neval);
      return this->func(x);
    }

    int
    num_evals() const
    { return *this->neval; }

    void
    num_evals(int num)
    { *this->neval = num; }

    std::function<_Tp(_Tp)> func;
    mutable std::shared_ptr<int> neval;
  };

template<typename _Tp>
  struct quadrature_test
  {
    unsigned int num_tests = 0;
    unsigned int num_passed = 0;
    unsigned int num_failed = 0;
    bool verbose = false;

    enum STATUS
    {
      SUCCESS,
      FAILURE
    };

    quadrature_test()
    : num_tests{},
      num_passed{},
      num_failed{},
      verbose{true}
    { }

    ~quadrature_test()
    {
      tot_num_tests += this->num_tests;
      tot_num_passed += this->num_passed;
      tot_num_failed += this->num_failed;
    }

    void test_update(int status);
    void test_update(int status, const char* test_desc);
    void test_relative(_Tp result, _Tp expected, _Tp rel_error, const char* test_desc);
    void test_absolute(_Tp result, _Tp expected, _Tp abs_error, const char* test_desc);
    void test_factor(_Tp result, _Tp expected, _Tp factor, const char* test_desc);
    void test_integer(int result, int expected, const char* test_desc);

    static unsigned int tot_num_tests;
    static unsigned int tot_num_passed;
    static unsigned int tot_num_failed;
    static bool tot_verbose;

    static int test_summary();
  };

template<typename _Tp>
  unsigned int
  quadrature_test<_Tp>::tot_num_tests = 0;

template<typename _Tp>
  unsigned int
  quadrature_test<_Tp>::tot_num_passed = 0;

template<typename _Tp>
  unsigned int
  quadrature_test<_Tp>::tot_num_failed = 0;

template<typename _Tp>
  bool
  quadrature_test<_Tp>::tot_verbose = true;


template<typename _Tp>
  void
  quadrature_test<_Tp>::test_update(int status)
  {
    ++this->num_tests;

    if (status == 0)
      ++this->num_passed;
    else
      ++this->num_failed;
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_update(int status, const char* test_desc)
  {
    this->test_update(status);

    if (status || this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status && !this->verbose)
          printf(" [%u]", this->num_tests);

	printf(" %s", test_desc);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_relative(_Tp result, _Tp expected, _Tp rel_error,
				      const char* test_desc)
  {
    // Check for NaN vs inf vs number.
    int status;
    if (std::isnan(result) || std::isnan(expected))
      status = std::isnan(result) != std::isnan(expected);
    else if (std::isinf(result) || std::isinf(expected))
      status = std::isinf(result) != std::isinf(expected);
    else if ((expected > 0 && expected < std::numeric_limits<_Tp>::min())
	  || (expected < 0 && expected > -(std::numeric_limits<_Tp>::min())))
      status = -1;
    else if (expected != 0)
      status = (std::abs(result - expected) > rel_error * std::abs(expected));
    else
      status = (std::abs(result) > rel_error);

    this->test_update(status);

    if (status || this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status == 0)
          printf(" (%g observed vs %g expected)", result, expected);
	else
          printf(" (%.18g observed vs %.18g expected)  \"%s\"",
		 result, expected, test_desc);

	if (status == -1)
          printf(" [test uses subnormal value]");

	if (status && !this->verbose)
          printf(" [%u]", this->num_tests);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_absolute(_Tp result, _Tp expected, _Tp abs_error,
				      const char* test_desc)
  {
    // Check for NaN vs inf vs number.
    int status;
    if (std::isnan(result) || std::isnan(expected))
      status = std::isnan(result) != std::isnan(expected);
    else if (std::isinf(result) || std::isinf(expected))
      status = std::isinf(result) != std::isinf(expected);
    else if ((expected > 0 && expected < std::numeric_limits<_Tp>::min())
	  || (expected < 0 && expected > -(std::numeric_limits<_Tp>::min())))
      status = -1;
    else
      status = std::abs(result - expected) > abs_error;

    this->test_update(status);

    if (status || this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status == 0)
          printf(" (%g observed vs %g expected)", result, expected);
	else
          printf(" (%.18g observed vs %.18g expected)  \"%s\"",
		 result, expected, test_desc);

	if (status == -1)
          printf(" [test uses subnormal value]");

	if (status && !this->verbose)
          printf(" [%u]", this->num_tests);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_factor(_Tp result, _Tp expected, _Tp factor,
				    const char* test_desc)
  {
    int status;
    if ((expected > 0 && expected < std::numeric_limits<_Tp>::min())
	|| (expected < 0 && expected > -(std::numeric_limits<_Tp>::min())))
      status = -1;
    else if (result == expected)
      status = 0;
    else if (expected == _Tp{0})
      status = (result > expected || result < expected);
    else
      {
	_Tp u = result / expected;
	status = (u > factor || u < _Tp{1} / factor);
      }

    this->test_update (status);

    if (status || this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status == 0)
          printf(" (%g observed vs %g expected)", result, expected);
	else
          printf(" (%.18g observed vs %.18g expected)  \"%s\"",
		 result, expected, test_desc);

	if (status == -1)
          printf(" [test uses subnormal value]");

	if (status && !this->verbose)
          printf(" [%u]", this->num_tests);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_integer(int result, int expected,
				     const char* test_desc)
  {
    int status = (result != expected);

    this->test_update(status);

    if (status || this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status == 0)
          printf(" (%d observed vs %d expected)", result, expected);
	else
          printf(" (%d observed vs %d expected)  \"%s\"",
		 result, expected, test_desc);

	if (status && !this->verbose)
          printf(" [%u]", this->num_tests);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  int
  quadrature_test<_Tp>::test_summary()
  {
    if (tot_verbose)
      printf("%d tests, passed %d, failed %d.\n",
	     tot_num_tests, tot_num_passed, tot_num_failed);

    if (tot_num_failed != 0)
      return FAILURE;

    if (tot_num_tests !=  tot_num_passed +  tot_num_failed)
      {
	if (tot_verbose)
          printf("TEST RESULTS DO NOT ADD UP %d != %d + %d\n",
                   tot_num_tests, tot_num_passed, tot_num_failed);
	return FAILURE;
      }

    if (tot_num_passed ==  tot_num_tests)
      {
	if (! tot_verbose)
          printf("Completed [%d/%d]\n", tot_num_passed, tot_num_tests);

	return SUCCESS;
      }

    return FAILURE;
  }


template<typename _Tp>
  void
  belch(const __gnu_test::_IntegrationError<double>& iex)
  {
    std::cout << "ERROR: " << iex.what()
	      << "       status = " << iex.error_code()
	      << "       result = " << iex.result()
	      << "       abserr = " << iex.abserr()
	      << std::endl;
    std::cerr << "ERROR: " << iex.what()
	      << "       status = " << iex.error_code()
	      << "       result = " << iex.result()
	      << "       abserr = " << iex.abserr()
	      << '\n';
  }

template<typename _Tp>
  struct
  test_ival
  {
    _Tp a, b, r, e;
  };


void
test_quadrature()
{
  using _Tp = double;

  const auto _S_pi = __gnu_cxx::__const_pi<_Tp>();

  // Test the basic Gauss-Kronrod rules with a smooth positive function.
  const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 15 with a smooth positive function..." << std::endl;

      double exp_result = 7.716049357767090777e-02;
      double exp_abserr = 2.990224871000550874e-06;
      double exp_resabs = 7.716049357767090777e-02;
      double exp_resasc = 4.434273814139995384e-02;
      quadrature_test<double> qtest;

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_15);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk15(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk15(f1) smooth abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk15(f1) smooth resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk15(f1) smooth resasc");

      std::tie(result, abserr, resabs, resasc)
	= qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_15);

      qtest.test_relative(result, -exp_result, 1.0e-15, "qk15(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk15(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk15(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk15(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 21 with a smooth positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 21 with a smooth positive function..." << std::endl;

      double exp_result = 7.716049379303084599e-02;
      double exp_abserr = 9.424302194248481445e-08;
      double exp_resabs = 7.716049379303084599e-02;
      double exp_resasc = 4.434311425038358484e-02;
      quadrature_test<double> qtest;

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_21);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk21(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk21(f1) smooth abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk21(f1) smooth resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk21(f1) smooth resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_21);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk21(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk21(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk21(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk21(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 31 with a smooth positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 31 with a smooth positive function..." << std::endl;

      double exp_result = 7.716049382494900855e-02;
      double exp_abserr = 1.713503193600029893e-09;
      double exp_resabs = 7.716049382494900855e-02;
      double exp_resasc = 4.427995051868838933e-02;
      quadrature_test<double> qtest;

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_31);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk31(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk31(f1) smooth abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk31(f1) smooth resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk31(f1) smooth resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_31);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk31(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk31(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk31(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk31(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 41 with a smooth positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 41 with a smooth positive function..." << std::endl;

      double exp_result = 7.716049382681375302e-02;
      double exp_abserr = 9.576386660975511224e-11;
      double exp_resabs = 7.716049382681375302e-02;
      double exp_resasc = 4.421521169637691873e-02;
      quadrature_test<double> qtest;

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_41);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk41(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk41(f1) smooth abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk41(f1) smooth resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk41(f1) smooth resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_41);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk41(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk41(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk41(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk41(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 51 with a smooth positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 51 with a smooth positive function..." << std::endl;

      double exp_result = 7.716049382708510540e-02;
      double exp_abserr = 1.002079980317363772e-11;
      double exp_resabs = 7.716049382708510540e-02;
      double exp_resasc = 4.416474291216854892e-02;
      quadrature_test<double> qtest;

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_51);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk51(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qk51(f1) smooth abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk51(f1) smooth resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk51(f1) smooth resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_51);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk51(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qk51(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk51(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk51(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 61 with a smooth positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 61 with a smooth positive function..." << std::endl;

      double exp_result = 7.716049382713800753e-02;
      double exp_abserr = 1.566060362296155616e-12;
      double exp_resabs = 7.716049382713800753e-02;
      double exp_resasc = 4.419287685934316506e-02;
      quadrature_test<double> qtest;

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_61);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk61(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qk61(f1) smooth abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk61(f1) smooth resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk61(f1) smooth resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_61);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk61(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qk61(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk61(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk61(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  /* Now test the basic rules with a positive function that has a
     singularity. This should give large values of abserr which would
     find discrepancies in the abserr calculation. */

  // Test Gauss-Kronrod 15 with a singular positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 15 with a singular positive function..." << std::endl;

      double exp_result = 1.555688196612745777e+01;
      double exp_abserr = 2.350164577239293706e+01;
      double exp_resabs = 1.555688196612745777e+01;
      double exp_resasc = 2.350164577239293706e+01;
      quadrature_test<double> qtest;

      double alpha = -0.9;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_15);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk15(f1) singular result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk15(f1) singular abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk15(f1) singular resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk15(f1) singular resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_15);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk15(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk15(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk15(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk15(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 21 with a singular positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 21 with a singular positive function..." << std::endl;

      double exp_result = 1.799045317938126232e+01;
      double exp_abserr = 2.782360287710622515e+01;
      double exp_resabs = 1.799045317938126232e+01;
      double exp_resasc = 2.782360287710622515e+01;
      quadrature_test<double> qtest;

      double alpha = -0.9;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_21);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk21(f1) singular result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk21(f1) singular abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk21(f1) singular resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk21(f1) singular resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_21);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk21(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk21(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk21(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk21(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 31 with a singular positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 31 with a singular positive function..." << std::endl;

      double exp_result = 2.081873305159121657e+01;
      double exp_abserr = 3.296500137482590276e+01;
      double exp_resabs = 2.081873305159121301e+01;
      double exp_resasc = 3.296500137482590276e+01;
      quadrature_test<double> qtest;

      double alpha = -0.9;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_31);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk31(f1) singular result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk31(f1) singular abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk31(f1) singular resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk31(f1) singular resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_31);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk31(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk31(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk31(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk31(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 41 with a singular positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 41 with a singular positive function..." << std::endl;

      double exp_result = 2.288677623903126701e+01;
      double exp_abserr = 3.671538820274916048e+01;
      double exp_resabs = 2.288677623903126701e+01;
      double exp_resasc = 3.671538820274916048e+01;
      quadrature_test<double> qtest;

      double alpha = -0.9;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_41);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk41(f1) singular result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk41(f1) singular abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk41(f1) singular resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk41(f1) singular resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_41);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk41(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk41(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk41(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk41(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 51 with a singular positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 51 with a singular positive function..." << std::endl;

      double exp_result = 2.449953612016972215e+01;
      double exp_abserr = 3.967771249391228849e+01;
      double exp_resabs = 2.449953612016972215e+01;
      double exp_resasc = 3.967771249391228849e+01;
      quadrature_test<double> qtest;

      double alpha = -0.9;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_51);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk51(f1) singular result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk51(f1) singular abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk51(f1) singular resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk51(f1) singular resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_51);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk51(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk51(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk51(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk51(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 61 with a singular positive function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 61 with a singular positive function..." << std::endl;

      double exp_result = 2.583030240976628988e+01;
      double exp_abserr = 4.213750493076978643e+01;
      double exp_resabs = 2.583030240976628988e+01;
      double exp_resasc = 4.213750493076978643e+01;
      quadrature_test<double> qtest;

      double alpha = -0.9;
      auto f = make_function<double>(f1, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_test::QK_61);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk61(f1) singular result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk61(f1) singular abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk61(f1) singular resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk61(f1) singular resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_test::QK_61);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk61(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk61(f1) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk61(f1) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk61(f1) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  /* Test the basic Gauss-Kronrod rules with a smooth oscillating
     function, over an unsymmetric range. This should find any
     discrepancies in the abscissae. */

  // Test Gauss-Kronrod 15 with a smooth oscillating function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 15 with a smooth oscillating function..." << std::endl;

      double exp_result = -7.238969575483799046e-01;
      double exp_abserr =  8.760080200939757174e-06;
      double exp_resabs =  1.165564172429140788e+00;
      double exp_resasc =  9.334560307787327371e-01;
      quadrature_test<double> qtest;

      double alpha = 1.3;
      auto f = make_function<double>(f3, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_15);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk15(f3) oscill result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk15(f3) oscill abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk15(f3) oscill resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk15(f3) oscill resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_15);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk15(f3) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk15(f3) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk15(f3) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk15(f3) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 21 with a smooth oscillating function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 21 with a smooth oscillating function..." << std::endl;

      double exp_result =-7.238969575482959717e-01;
      double exp_abserr = 7.999213141433641888e-11;
      double exp_resabs = 1.150829032708484023e+00;
      double exp_resasc = 9.297591249133687619e-01;
      quadrature_test<double> qtest;

      double alpha = 1.3;
      auto f = make_function<double>(f3, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_21);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk21(f3) oscill result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qk21(f3) oscill abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk21(f3) oscill resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk21(f3) oscill resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_21);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk21(f3) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qk21(f3) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk21(f3) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk21(f3) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 31 with a smooth oscillating function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 31 with a smooth oscillating function..." << std::endl;

      double exp_result =-7.238969575482959717e-01;
      double exp_abserr = 1.285805464427459261e-14;
      double exp_resabs = 1.158150602093290571e+00;
      double exp_resasc = 9.277828092501518853e-01;
      quadrature_test<double> qtest;

      double alpha = 1.3;
      auto f = make_function<double>(f3, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_31);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk31(f3) oscill result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk31(f3) oscill abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk31(f3) oscill resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk31(f3) oscill resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_31);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk31(f3) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk31(f3) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk31(f3) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk31(f3) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 41 with a smooth oscillating function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 41 with a smooth oscillating function..." << std::endl;

      double exp_result =-7.238969575482959717e-01;
      double exp_abserr = 1.286535726271015626e-14;
      double exp_resabs = 1.158808363486595328e+00;
      double exp_resasc = 9.264382258645686985e-01;
      quadrature_test<double> qtest;

      double alpha = 1.3;
      auto f = make_function<double>(f3, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_41);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk41(f3) oscill result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk41(f3) oscill abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk41(f3) oscill resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk41(f3) oscill resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_41);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk41(f3) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk41(f3) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk41(f3) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk41(f3) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 51 with a smooth oscillating function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 51 with a smooth oscillating function..." << std::endl;

      double exp_result =-7.238969575482961938e-01;
      double exp_abserr = 1.285290995039385778e-14;
      double exp_resabs = 1.157687209264406381e+00;
      double exp_resasc = 9.264666884071264263e-01;
      quadrature_test<double> qtest;

      double alpha = 1.3;
      auto f = make_function<double>(f3, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_51);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk51(f3) oscill result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk51(f3) oscill abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk51(f3) oscill resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk51(f3) oscill resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_51);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk51(f3) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk51(f3) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk51(f3) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk51(f3) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Gauss-Kronrod 61 with a smooth oscillating function.
  try
    {
      std::cout << ">>>> Test Gauss-Kronrod 61 with a smooth oscillating function..." << std::endl;

      double exp_result =-7.238969575482959717e-01;
      double exp_abserr = 1.286438572027470736e-14;
      double exp_resabs = 1.158720854723590099e+00;
      double exp_resasc = 9.270469641771273972e-01;
      quadrature_test<double> qtest;

      double alpha = 1.3;
      auto f = make_function<double>(f3, alpha);

      auto [result, abserr, resabs, resasc]
	= __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_61);
      qtest.test_relative(result, exp_result, 1.0e-15, "qk61(f3) oscill result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk61(f3) oscill abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk61(f3) oscill resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk61(f3) oscill resasc");

      std::tie(result, abserr, resabs, resasc)
	= __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_61);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qk61(f3) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qk61(f3) reverse abserr");
      qtest.test_relative(resabs, exp_resabs, 1.0e-15, "qk61(f3) reverse resabs");
      qtest.test_relative(resasc, exp_resasc, 1.0e-15, "qk61(f3) reverse resasc");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      int status = __gnu_test::NO_ERROR;
      _Tp exp_result = 7.716049379303083211e-02;
      _Tp exp_abserr = 9.424302199601294244e-08;
      int exp_neval  =  21;
      int exp_ier    =   __gnu_test::NO_ERROR;
      quadrature_test<_Tp> qtest;

      _Tp alpha = _Tp{2.6};
      auto f = make_function<_Tp>(f1, alpha);
      auto fc = counted_function<double>(f);

      auto [result, abserr]
	= __gnu_test::qng_integrate(fc, _Tp{0}, _Tp{1}, _Tp{1.0e-1}, _Tp{0});
      qtest.test_relative(result, exp_result, 1.0e-15, "qng(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qng(f1) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) smooth neval");
      qtest.test_integer(status, exp_ier, "qng(f1) smooth status");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qng_integrate(fc, _Tp{1}, _Tp{0}, _Tp{1.0e-1}, _Tp{0});
      qtest.test_relative(result, -exp_result, 1.0e-15, "qng(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qng(f1) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) reverse neval");
      qtest.test_integer(status, exp_ier, "qng(f1) reverse status");
    }
  catch (__gnu_test::_IntegrationError<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      quadrature_test<double> qtest;

      int status = 0;
      double exp_result = 7.716049382706505200e-02;
      double exp_abserr = 2.666893044866214501e-12;
      int exp_neval  =  43;
      int exp_ier    =   __gnu_test::NO_ERROR;

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);
      auto fc = counted_function<double>(f);

      auto [result, abserr]
	= __gnu_test::qng_integrate(fc, _Tp{0}, _Tp{1}, _Tp{0}, 1.0e-9);
      qtest.test_relative(result, exp_result, 1.0e-15, "qng(f1) smooth 43pt result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qng(f1) smooth 43pt abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) smooth 43pt neval");
      qtest.test_integer(status, exp_ier, "qng(f1) smooth 43pt status");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qng_integrate(fc, _Tp{1}, _Tp{0}, _Tp{0}, 1.0e-9);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qng(f1) reverse 43pt result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qng(f1) reverse 43pt abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) reverse 43pt neval");
      qtest.test_integer(status, exp_ier, "qng(f1) reverse 43pt status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      int status = 0;
      double exp_result =-7.238969575482961938e-01;
      double exp_abserr = 1.277676889520056369e-14;
      int exp_neval  =  43;
      int exp_ier    =   __gnu_test::NO_ERROR;
      quadrature_test<double> qtest;

      double alpha = 1.3;
      auto f = make_function<double>(f3, alpha);
      auto fc = counted_function<double>(f);

      auto [result, abserr]
	= __gnu_test::qng_integrate(fc, 0.3, 2.71, _Tp{0}, 1.0e-12);
      qtest.test_relative(result, exp_result, 1.0e-15, "qnq(f3) oscill result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qng(f3) oscill abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f3) oscill neval");
      qtest.test_integer(status, exp_ier, "qng(f3) oscill status");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qng_integrate(fc, 2.71, 0.3, _Tp{0}, 1.0e-12);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qnq(f3) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qng(f3) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f3) reverse neval");
      qtest.test_integer(status, exp_ier, "qng(f3) reverse status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      quadrature_test<double> qtest;

      int status = 0;
      double exp_result = 7.716049382716029525e-02;
      double exp_abserr = 8.566535680046930668e-16;
      int exp_neval  =  87;
      int exp_ier    =   __gnu_test::NO_ERROR;

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);
      auto fc = counted_function<double>(f);

      auto [result, abserr]
	= __gnu_test::qng_integrate(fc, _Tp{0}, _Tp{1}, _Tp{0}, 1.0e-13);
      qtest.test_relative(result, exp_result, 1.0e-15, "qng(f1) 87pt smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qng(f1) 87pt smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) 87pt smooth neval");
      qtest.test_integer(status, exp_ier, "qng(f1) 87pt smooth status");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qng_integrate(fc, _Tp{1}, _Tp{0}, _Tp{0}, 1.0e-13);
      qtest.test_relative(result, -exp_result, 1.0e-15, "qng(f1) 87pt reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qng(f1) 87pt reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) 87pt reverse neval");
      qtest.test_integer(status, exp_ier, "qng(f1) 87pt reverse status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      quadrature_test<double> qtest;

      int status = 0;
      double exp_result = 3.222948711817264211e+01;
      double exp_abserr = 2.782360287710622870e+01;
      int exp_neval  =  87;
      int exp_ier    =  __gnu_test::TOLERANCE_ERROR;

      double alpha = -0.9;
      auto f = make_function<double>(f1, alpha);
      auto fc = counted_function<double>(f);

      double result, abserr;
      try
	{
	  std::tie(result, abserr)
	    = __gnu_test::qng_integrate(fc, _Tp{0}, _Tp{1}, _Tp{0}, 1.0e-3);
	}
      catch (__gnu_test::_IntegrationError<double>& iex)
	{
	  status = iex.error_code();
	  result = iex.result();
	  abserr = iex.abserr();
	}
      qtest.test_relative(result, exp_result, 1.0e-15, "qng(f1) sing beyond 87pt result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qng(f1) sing beyond 87pt abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) sing beyond 87pt neval");
      qtest.test_integer(status, exp_ier, "qng(f1) sing beyond 87pt status");

      fc.num_evals(0);
      try
	{
	  std::tie(result, abserr)
	    = __gnu_test::qng_integrate(fc, _Tp{1}, _Tp{0}, _Tp{0}, 1.0e-3);
	}
      catch (__gnu_test::_IntegrationError<double>& iex)
	{
	  status = iex.error_code();
	  result = iex.result();
	  abserr = iex.abserr();
	}
      qtest.test_relative(result, -exp_result, 1.0e-15, "qng(f1) reverse beyond 87pt result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-7, "qng(f1) rev beyond 87pt abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) rev beyond 87pt neval");
      qtest.test_integer(status, exp_ier, "qng(f1) rev beyond 87pt status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the adaptive Gaussian integrator.
  try
    {
      std::cout << ">>>> Test adaptive Gaussian integrator..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = 7.716049382715854665e-02;
      double exp_abserr = 6.679384885865053037e-12;
      int exp_neval  =     165;
      int exp_ier    =       __gnu_test::NO_ERROR;
      int exp_last   =       6;

      constexpr std::size_t num_test = 6;
      test_ival<double> test[num_test]
      {
	{0.0,     0.03125, 3.966769831709074375e-06, 6.678528276336181873e-12},
	{0.5,     1.0,     5.491842501998222409e-02, 6.097169993333454062e-16},
	{0.125,   0.25,    2.776531175604360531e-03, 3.082568839745514608e-17},
	{0.0625,  0.125,   3.280661030752063693e-04, 3.642265412331439511e-18},
	{0.25,    0.5,     1.909827770934243926e-02, 2.120334764359736934e-16},
	{0.03125, 0.0625,  3.522704932261797744e-05, 3.910988124757650942e-19},
      };

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);
      auto fc = counted_function<double>(f);

      auto [result, abserr]
	= __gnu_test::qag_integrate(w, fc, _Tp{0}, _Tp{1}, _Tp{0}, 1.0e-10, 1000,
				    __gnu_test::QK_15);

      qtest.test_relative(result, exp_result, 1.0e-15, "qag(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qag(f1) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f1) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qag(f1) smooth last");
      qtest.test_integer(status, exp_ier, "qag(f1) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qag(f1) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qag(f1) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, 1.0e-15, "qag(f1) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-6, "qag(f1) smooth abs error");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qag_integrate(w, fc, _Tp{1}, _Tp{0}, _Tp{0}, 1.0e-10, 1000,
				    __gnu_test::QK_15);

      qtest.test_relative(result, -exp_result, 1.0e-15, "qag(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qag(f1) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f1) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qag(f1) reverse last");
      qtest.test_integer(status, exp_ier, "qag(f1) reverse status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  try
    {
      std::cout << ">>>> Test adaptive Gaussian integrator with absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = 7.716049382716050342e-02;
      double exp_abserr = 2.227969521869139532e-15;
      int exp_neval  =     315;
      int exp_ier    =       __gnu_test::NO_ERROR;
      int exp_last   =       8;

      constexpr std::size_t num_test = 8;
      test_ival<double> test[num_test]
      {
	{0.0,       0.0078125, 3.696942726831556522e-08, 1.371316364034059572e-15},
	{0.25,      0.5,	     1.909827770934243579e-02, 2.120334764359736441e-16},
	{0.5,       1.0,	     5.491842501998223103e-02, 6.097169993333454062e-16},
	{0.0625,    0.125,     3.280661030752062609e-04, 3.642265412331439511e-18},
	{0.03125,   0.0625,    3.522704932261797744e-05, 3.910988124757650460e-19},
	{0.015625,  0.03125,   3.579060884684503576e-06, 3.973555800712018091e-20},
	{0.125,     0.25,	     2.776531175604360097e-03, 3.082568839745514608e-17},
	{0.0078125, 0.015625,  3.507395216921808047e-07, 3.893990926286736620e-21},
      };

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);
      auto fc = counted_function<double>(f);

      auto [result, abserr]
	= __gnu_test::qag_integrate(w, fc, _Tp{0}, _Tp{1}, 1.0e-14, _Tp{0}, 1000, __gnu_test::QK_21);

      qtest.test_relative(result, exp_result, 1.0e-15, "qag(f1, 21pt) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qag(f1, 21pt) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f1, 21pt) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qag(f1, 21pt) smooth last");
      qtest.test_integer(status, exp_ier, "qag(f1, 21pt) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qag(f1, 21pt) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qag(f1, 21pt) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, 1.0e-15, "qag(f1, 21pt) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-6, "qag(f1, 21pt) smooth abs error");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qag_integrate(w, fc, _Tp{1}, _Tp{0}, 1.0e-14, _Tp{0}, 1000, __gnu_test::QK_21);

      qtest.test_relative(result, -exp_result, 1.0e-15, "qag(f1, 21pt) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qag(f1, 21pt) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f1, 21pt) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qag(f1, 21pt) reverse last");
      qtest.test_integer(status, exp_ier, "qag(f1, 21pt) reverse status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Adaptive integration of an oscillatory function which terminates because
  // of roundoff error, uses the 31-pt rule.
  try
    {
      std::cout << ">>>> Test adaptive integration of an oscillatory function\n"
		   ">>>> which terminates because of roundoff error, uses the 31-pt rule..." << std::endl;

      int status = __gnu_test::NO_ERROR;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = -7.238969575482959717e-01;
      double exp_abserr =  1.285805464427459261e-14;
      int exp_neval   =     31;
      int exp_ier     =     __gnu_test::ROUNDOFF_ERROR;
      int exp_last    =     1;

      double alpha = 1.3;
      auto f = make_function<double>(f3, alpha);

      auto fc = counted_function<double>(f);

      double result, abserr;
      try
	{
	  std::tie(result, abserr)
	    = __gnu_test::qag_integrate(w, fc, 0.3, 2.71, 1.0e-14, _Tp{0}, 1000,
				       __gnu_test::QK_31);
	}
      catch (__gnu_test::_IntegrationError<double>& iex)
	{
	  result = iex.result();
	  abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_relative(result, exp_result, 1.0e-15, "qag(f3, 31pt) oscill result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qag(f3, 31pt) oscill abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f3, 31pt) oscill neval");
      qtest.test_integer(w.size(), exp_last, "qag(f3, 31pt) oscill last");
      qtest.test_integer(status, exp_ier, "qag(f3, 31pt) oscill status");

      fc.num_evals(0);
      try
	{
	  std::tie(result, abserr)
	    = __gnu_test::qag_integrate(w, fc, 2.71, 0.3, 1.0e-14, _Tp{0}, 1000,
				       __gnu_test::QK_31);
	}
      catch (__gnu_test::_IntegrationError<double>& iex)
	{
	  result = iex.result();
	  abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_relative(result, -exp_result, 1.0e-15, "qag(f3, 31pt) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qag(f3, 31pt) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f3, 31pt) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qag(f3, 31pt) reverse last");
      qtest.test_integer(status, exp_ier, "qag(f3, 31pt) reverse status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Check the singularity detection (singularity at x=-0.1 in this example).
  try
    {
      std::cout << ">>>> Test singularity detection (singularity at x=-0.1 in this example)..." << std::endl;

      int status = __gnu_test::NO_ERROR;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      int exp_neval  =     5151;
      int exp_ier    =     __gnu_test::SINGULAR_ERROR;
      int exp_last   =     51;

      double alpha = _Tp{2};
      auto f = make_function<double>(f16, alpha);
      auto fc = counted_function<double>(f);

      double result, abserr;
      try
	{
	  std::tie(result, abserr)
	    = __gnu_test::qag_integrate(w, fc, _Tp{-1}, _Tp{1}, 1.0e-14, _Tp{0}, 1000,
					__gnu_test::QK_51);
	}
      catch (__gnu_test::_IntegrationError<double>& iex)
	{
	  result = iex.result();
	  abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f16, 51pt) sing neval");
      qtest.test_integer(w.size(), exp_last, "qag(f16, 51pt) sing last");
      qtest.test_integer(status, exp_ier, "qag(f16, 51pt) sing status");

      fc.num_evals(0);
      try
	{
	  std::tie(result, abserr)
	    = __gnu_test::qag_integrate(w, fc, _Tp{1}, _Tp{-1}, 1.0e-14, _Tp{0}, 1000,
					__gnu_test::QK_51);
	}
      catch (__gnu_test::_IntegrationError<double>& iex)
	{
	  result = iex.result();
	  abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f16, 51pt) rev neval");
      qtest.test_integer(w.size(), exp_last, "qag(f16, 51pt) rev last");
      qtest.test_integer(status, exp_ier, "qag(f16, 51pt) rev status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Check for hitting the iteration limit.
  try
    {
      std::cout << ">>>> Test hitting the iteration limit..." << std::endl;

      int status = __gnu_test::NO_ERROR;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(3);

      double exp_result =  9.565151449233894709;
      double exp_abserr =  1.570369823891028460e+01;
      int exp_neval  =     305;
      int exp_ier    =     __gnu_test::MAX_ITER_ERROR;
      int exp_last   =     3;

      constexpr std::size_t num_test = 3;
      test_ival<double> test[num_test]
      {
	{-5.000000000000000000e-01,  0.000000000000000000,     9.460353469435913709,     1.570369823891028460e+01},
	{-1.000000000000000000,     -5.000000000000000000e-01, 1.388888888888888812e-02, 1.541976423090495140e-16},
	{ 0.000000000000000000,      1.000000000000000000,     9.090909090909091161e-02, 1.009293658750142399e-15},
      };

      double alpha = _Tp{1};
      auto f = make_function<double>(f16, alpha);
      auto fc = counted_function<double>(f);

      double result, abserr;
      try
	{
	  std::tie(result, abserr)
	    = __gnu_test::qag_integrate(w, fc, _Tp{-1}, _Tp{1}, 1.0e-14, _Tp{0}, 3,
					__gnu_test::QK_61);
	}
      catch (__gnu_test::_IntegrationError<double>& iex)
	{
	  result = iex.result();
	  abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_relative(result, exp_result, 1.0e-15, "qag(f16, 61pt) limit result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qag(f16, 61pt) limit abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f16, 61pt) limit neval");
      qtest.test_integer(w.size(), exp_last, "qag(f16, 61pt) limit last");
      qtest.test_integer(status, exp_ier, "qag(f16, 61pt) limit status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qag(f16, 61pt) limit lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qag(f16, 61pt) limit upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, 1.0e-15, "qag(f16, 61pt) limit integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-6, "qag(f16, 61pt) limit abs error");

      fc.num_evals(0);
      try
	{
	  std::tie(result, abserr)
	    = __gnu_test::qag_integrate(w, fc, _Tp{1}, _Tp{-1}, 1.0e-14, _Tp{0}, 3,
					__gnu_test::QK_61);
	}
      catch (__gnu_test::_IntegrationError<double>& iex)
	{
	  result = iex.result();
	  abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_relative(result, -exp_result, 1.0e-15, "qag(f16, 61pt) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qag(f16, 61pt) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f16, 61pt) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qag(f16, 61pt) reverse last");
      qtest.test_integer(status, exp_ier, "qag(f16, 61pt) reverse status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the adaptive integrator with extrapolation.
  try
    {
      std::cout << ">>>> Test the adaptive integrator with extrapolation..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = 7.716049382715789440e-02;
      double exp_abserr = 2.216394961010438404e-12;
      int exp_neval  =     189;
      int exp_ier    =       __gnu_test::NO_ERROR;
      int exp_last   =       5;

      constexpr std::size_t num_test = 5;
      test_ival<double> test[num_test]
      {
	0.0,    0.0625, 3.919381915366914693e-05, 2.215538742580964735e-12,
	0.5,    1.0,    5.491842501998223103e-02, 6.097169993333454062e-16,
	0.25,   0.5,    1.909827770934243579e-02, 2.120334764359736441e-16,
	0.125,  0.25,   2.776531175604360097e-03, 3.082568839745514608e-17,
	0.0625, 0.125,  3.280661030752062609e-04, 3.642265412331439511e-18,
      };

      double alpha = 2.6;
      auto f = make_function<double>(f1, alpha);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto [result, abserr]
	= __gnu_test::qags_integrate(w, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1e-15, "qags(f1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1e-6, "qags(f1) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qags(f1) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qags(f1) smooth last");
      qtest.test_integer(status, exp_ier, "qags(f1) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qags(f1) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qags(f1) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, 1.0e-15, "qags(f1) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-6, "qags(f1) smooth abs error");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qags_integrate(w, fc, _Tp{1}, _Tp{0}, _Tp{0}, 1e-10);

      qtest.test_relative(result, -exp_result, 1.0e-15, "qags(f1) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qags(f1) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qags(f1) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qags(f1) reverse last");
      qtest.test_integer(status, exp_ier, "qags(f1) reverse status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test f11 using an absolute error bound.
  try
    {
      std::cout << ">>>> Test f11 using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = -5.908755278982136588e+03;
      double exp_abserr = 1.299646281053874554e-10;
      int exp_neval  =     357;
      int exp_ier    =       __gnu_test::NO_ERROR;
      int exp_last   =       9;

      constexpr std::size_t num_test = 9;
      test_ival<double> test[num_test]
      {
	{1.000000000000000000e+00, 4.902343750000000000e+00, -3.890977835520834649e+00, 6.448276035006137169e-11},
	{5.005000000000000000e+02, 1.000000000000000000e+03, -3.297343675805121620e+03, 3.660786868980994028e-11},
	{1.258750000000000000e+02, 2.507500000000000000e+02, -6.517404019686431411e+02, 7.235772003440423011e-12},
	{2.507500000000000000e+02, 5.005000000000000000e+02, -1.475904154146372775e+03, 1.638582774073219226e-11},
	{3.221875000000000000e+01, 6.343750000000000000e+01, -1.201692001973227519e+02, 1.334146129098576244e-12},
	{6.343750000000000000e+01, 1.258750000000000000e+02, -2.829354222635842007e+02, 3.141214202790722909e-12},
	{1.660937500000000000e+01, 3.221875000000000000e+01, -4.959999906099650246e+01, 5.506706097890446534e-13},
	{8.804687500000000000e+00, 1.660937500000000000e+01, -1.971441499411640308e+01, 2.188739744348345039e-13},
	{4.902343750000000000e+00, 8.804687500000000000e+00, -7.457032710459004399e+00, 8.278969410534525339e-14},
      };

      double alpha = _Tp{2};
      auto f = make_function<double>(f11, alpha);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto [result, abserr]
	= __gnu_test::qags_integrate(w, fc, _Tp{1}, _Tp{1000}, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1e-15, "qags(f11) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1e-3, "qags(f11) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qags(f11) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qags(f11) smooth last");
      qtest.test_integer(status, exp_ier, "qags(f11) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qags(f11) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qags(f11) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, 1.0e-15, "qags(f11) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-5, "qags(f11) smooth abs error");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qags_integrate(w, fc, _Tp{1000}, _Tp{1}, 1e-7, _Tp{0});

      qtest.test_relative(result, -exp_result, 1.0e-15, "qags(f11) reverse result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-3, "qags(f11) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qags(f11) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qags(f11) reverse last");
      qtest.test_integer(status, exp_ier, "qags(f11) reverse status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test infinite range integral f455 using a relative error bound.
  try
    {
      std::cout << ">>>> Test infinite range integral f455 using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = -3.616892186127022568e-01;
      double exp_abserr = 3.016716913328831851e-06;
      int exp_neval  =      285;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =       10;

      constexpr std::size_t num_test = 10;
      test_ival<double> test[num_test]
      {
	{0.000000000000000000e+00, 3.125000000000000000e-02,  1.429785306003466313e-03, 2.161214992172538524e-04},
	{8.750000000000000000e-01, 9.375000000000000000e-01, -8.653752279614615461e-02, 9.607595030230581153e-16},
	{2.500000000000000000e-01, 5.000000000000000000e-01,  2.995321156568048898e-03, 3.325474514168701167e-17},
	{9.375000000000000000e-01, 9.687500000000000000e-01, -8.398745675010892142e-02, 9.324480826368044019e-16},
	{5.000000000000000000e-01, 7.500000000000000000e-01, -1.229943369113085765e-02, 5.720644840858777846e-14},
	{1.250000000000000000e-01, 2.500000000000000000e-01,  2.785385934678596704e-03, 3.092399597147240624e-17},
	{6.250000000000000000e-02, 1.250000000000000000e-01,  1.736218164975512294e-03, 1.927589382528252344e-17},
	{7.500000000000000000e-01, 8.750000000000000000e-01, -4.980050133751051655e-02, 3.147380432198176412e-14},
	{3.125000000000000000e-02, 6.250000000000000000e-02,  1.041689192004495576e-03, 1.156507325466566521e-17},
	{9.687500000000000000e-01, 1.000000000000000000e+00, -1.390003415539725340e-01, 2.395037249893453013e-02},
      };

      auto f = make_function<double>(f455);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-3};
      auto [result, abserr]
	= __gnu_test::qagiu_integrate(w, fc, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(result, exp_result, /*1.0e-14*/epsrel, "qagiu(f455) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qagiu(f455) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagiu(f455) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagiu(f455) smooth last");
      qtest.test_integer(status, exp_ier, "qagiu(f455) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qagiu(f455) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qagiu(f455) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-15*/epsrel, "qagiu(f455) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-4, "qagiu(f455) smooth abs error");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test infinite range integral f15 using a relative error bound.
  try
    {
      std::cout << ">>>> Test infinite range integral f15 using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = 6.553600000000024738e+04;
      double exp_abserr = 7.121667111456009280e-04;
      int exp_neval  =      285;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =       10;

      constexpr std::size_t num_test = 10;
      test_ival<double> test[num_test]
      {
	{0.000000000000000000e+00, 1.953125000000000000e-03, 1.099297665754340292e+00, 7.101865971621337814e-04},
	{7.812500000000000000e-03, 1.562500000000000000e-02, 2.899418134793237914e+04, 1.934652413547325474e-07},
	{1.562500000000000000e-02, 3.125000000000000000e-02, 1.574317583220441520e+04, 1.380003928453846583e-07},
	{3.906250000000000000e-03, 7.812500000000000000e-03, 1.498314766425578091e+04, 3.408933028357320364e-07},
	{2.500000000000000000e-01, 5.000000000000000000e-01, 8.064694554185326325e+00, 9.167763417119923333e-08},
	{1.250000000000000000e-01, 2.500000000000000000e-01, 8.873128656118993263e+01, 3.769501719163865578e-07},
	{6.250000000000000000e-02, 1.250000000000000000e-01, 6.977679035845269482e+02, 6.973493131275552509e-07},
	{5.000000000000000000e-01, 1.000000000000000000e+00, 3.256176475185617591e-01, 1.912660677170175771e-08},
	{3.125000000000000000e-02, 6.250000000000000000e-02, 4.096981198511257389e+03, 1.205653952340679711e-07},
	{1.953125000000000000e-03, 3.906250000000000000e-03, 9.225251570832365360e+02, 2.132473175465897029e-09},
      };

      double alpha = _Tp{5};

      auto f = make_function<double>(f15, alpha);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto [result, abserr]
	= __gnu_test::qagiu_integrate(w, fc, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(result, exp_result, /*1.0e-14*/epsrel, "qagiu(f15) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qagiu(f15) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagiu(f15) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagiu(f15) smooth last");
      qtest.test_integer(status, exp_ier, "qagiu(f15) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qagiu(f15) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qagiu(f15) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-15*/epsrel, "qagiu(f15) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-4, "qagiu(f15) smooth abs error");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test infinite range integral f16 using an absolute error bound.
  try
    {
      std::cout << ">>>> Test infinite range integral f16 using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = 1.000000000006713292e-04;
      double exp_abserr = 3.084062020905636316e-09;
      int exp_neval  =      165;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =        6;

      constexpr std::size_t num_test = 6;
      test_ival<double> test[num_test]
      {
	{0.000000000000000000e+00, 3.125000000000000000e-02, 7.633587786326674618e-05, 3.084061858351569051e-09},
	{6.250000000000000000e-02, 1.250000000000000000e-01, 6.501422186103209199e-06, 3.014338672269481784e-17},
	{2.500000000000000000e-01, 5.000000000000000000e-01, 1.922522349322310737e-06, 4.543453652226561245e-17},
	{5.000000000000000000e-01, 1.000000000000000000e+00, 9.900990099009899620e-07, 3.112064814755089674e-17},
	{1.250000000000000000e-01, 2.500000000000000000e-01, 3.629434715543053753e-06, 4.908618166361344548e-17},
	{3.125000000000000000e-02, 6.250000000000000000e-02, 1.062064387653501389e-05, 6.795996738013555461e-18},
      };

      double alpha = _Tp{1};

      auto f = make_function<double>(f16, alpha);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto [result, abserr]
	= __gnu_test::qagiu_integrate(w, fc, 99.9, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1.0e-14, "qagiu(f16) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qagiu(f16) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagiu(f16) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagiu(f16) smooth last");
      qtest.test_integer(status, exp_ier, "qagiu(f16) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qagiu(f16) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qagiu(f16) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-15*/epsabs, "qagiu(f16) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-4, "qagiu(f16) smooth abs error");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test infinite range integral myfn1 using an absolute error bound.
  try
    {
      std::cout << ">>>> Test infinite range integral myfn1 using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = 2.275875794468747770e+00;
      double exp_abserr = 7.436490118267390744e-09;
      int exp_neval  =      270;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =        5;

      constexpr std::size_t num_test = 5;
      test_ival<double> test[num_test]
      {
	{0.000000000000000000e+00, 1.250000000000000000e-01, 4.379392477350953574e-20, 8.360902986775307673e-20},
	{5.000000000000000000e-01, 1.000000000000000000e+00, 1.691664195356748834e+00, 4.265988974874425043e-09},
	{2.500000000000000000e-01, 3.750000000000000000e-01, 1.146307471900291086e-01, 1.231954072964969637e-12},
	{1.250000000000000000e-01, 2.500000000000000000e-01, 4.639317228058405717e-04, 3.169263960393051137e-09},
	{3.750000000000000000e-01, 5.000000000000000000e-01, 4.691169201991640669e-01, 5.208244060463541433e-15},
      };

      auto f = make_function<double>(myfn1);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto [result, abserr] = __gnu_test::qagi_integrate(w, fc, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1e-14, "qagi(myfn1) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1e-5, "qagi(myfn1) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagi(myfn1) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagi(myfn1) smooth last");
      qtest.test_integer(status, exp_ier, "qagi(myfn1) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qagi(myfn1) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qagi(myfn1) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-14*/epsabs, "qagi(myfn1) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-4, "qagi(myfn1) smooth abs error");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test infinite range integral myfn2 using an absolute error bound.
  try
    {
      std::cout << ">>>> Test infinite range integral myfn2 using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = 2.718281828459044647e+00;
      double exp_abserr = 1.588185109253204805e-10;
      int exp_neval  =      135;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =        5;

      constexpr std::size_t num_test = 5;
      test_ival<double> test[num_test]
      {
	{0.000000000000000000e+00, 6.250000000000000000e-02, 8.315287189746029816e-07, 1.533437090413525935e-10},
	{2.500000000000000000e-01, 5.000000000000000000e-01, 8.646647167633871867e-01, 7.802455785301941044e-13},
	{5.000000000000000000e-01, 1.000000000000000000e+00, 1.718281828459045091e+00, 4.117868247943567505e-12},
	{1.250000000000000000e-01, 2.500000000000000000e-01, 1.328565310599463256e-01, 5.395586026138397182e-13},
	{6.250000000000000000e-02, 1.250000000000000000e-01, 2.477920647947255521e-03, 3.713312434866150125e-14},
      };

      double alpha = _Tp{1};
      auto f = make_function<double>(myfn2, alpha);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto [result, abserr]
	= __gnu_test::qagil_integrate(w, fc, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1.0e-14, "qagil(myfn2) smooth result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qagil(myfn2) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagil(myfn2) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagil(myfn2) smooth last");
      qtest.test_integer(status, exp_ier, "qagil(myfn2) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qagil(myfn2) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qagil(myfn2) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-14*/epsabs, "qagil(myfn2) smooth integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-4, "qagil(myfn2) smooth abs error");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test integral f454 with integrable singular points.
  try
    {
      std::cout << ">>>> Test integral f454 with integrable singular points..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = 5.274080611672716401e+01;
      double exp_abserr = 1.755703848687062418e-04;
      int exp_neval  = 777;
      int exp_ier    = __gnu_test::NO_ERROR;
      int exp_last   = 20;

      constexpr std::size_t num_test = 20;
      test_ival<double> test[num_test]
      {
	{1.000000000000000000e+00, 1.051776695296636976e+00, -1.830392049835374568e-01, 3.252808038935910834e-02},
	{1.401269388548935790e+00, 1.414213562373095145e+00, -1.565132123531515207e-01, 2.730454695485963826e-02},
	{2.207106781186547462e+00, 3.000000000000000000e+00,  4.873920540843067783e+01, 5.411138804637469780e-13},
	{9.687500000000000000e-01, 1.000000000000000000e+00, -1.125078814079027711e-01, 2.506431410088378817e-02},
	{1.612436867076458391e+00, 1.810660171779821415e+00,  5.911661670635662835e-01, 6.573952690524728748e-15},
	{1.207106781186547462e+00, 1.310660171779821415e+00, -2.991531901645863023e-01, 3.321267596107916554e-15},
	{1.362436867076458391e+00, 1.388325214724776657e+00, -1.589345454585119055e-01, 1.764527917763735212e-15},
	{1.463769388548935790e+00, 1.513325214724776657e+00, -2.208730668000830344e-01, 2.452183642810224359e-15},
	{1.810660171779821415e+00, 2.207106781186547462e+00,  6.032891565603589079e+00, 6.697855121200013106e-14},
	{0.000000000000000000e+00, 5.000000000000000000e-01,  6.575875041899758092e-03, 7.300687878575027348e-17},
	{7.500000000000000000e-01, 8.750000000000000000e-01, -5.647871991778510847e-02, 6.270397525408045936e-16},
	{1.103553390593273731e+00, 1.207106781186547462e+00, -2.431894410706912923e-01, 2.699945168224041491e-15},
	{1.051776695296636976e+00, 1.103553390593273731e+00, -1.305470403178642658e-01, 1.449363299575615261e-15},
	{5.000000000000000000e-01, 7.500000000000000000e-01, -7.326282608704996063e-03, 1.417509685426979386e-16},
	{8.750000000000000000e-01, 9.375000000000000000e-01, -7.406626263352669715e-02, 8.223007012367522077e-16},
	{1.513325214724776657e+00, 1.612436867076458391e+00, -1.721363984401322045e-01, 1.911097929242846383e-15},
	{1.388325214724776657e+00, 1.401269388548935790e+00, -1.048692749517999567e-01, 1.164282836272345215e-15},
	{9.375000000000000000e-01, 9.687500000000000000e-01, -6.302287584527696551e-02, 6.996944784151910810e-16},
	{1.310660171779821415e+00, 1.362436867076458391e+00, -2.236786562536174916e-01, 2.483331942899818875e-15},
	{1.414213562373095145e+00, 1.463769388548935790e+00, -4.225328513207429193e-01, 1.017446081816190118e-01},
      };

      auto f = make_function<double>(f454);
      auto fc = counted_function<double>(f);

      std::vector<double> pts{_Tp{0}, _Tp{1}, std::sqrt(_Tp{2}), _Tp{3}};

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-3};
      auto [result, abserr] = __gnu_test::qagp_integrate(w, fc, pts, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1.0e-14, "qagp(f454) singular result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-5, "qagp(f454) singular abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagp(f454) singular neval");
      qtest.test_integer(w.size(), exp_last, "qagp(f454) singular last");
      qtest.test_integer(status, exp_ier, "qagp(f454) singular status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qagp(f454) singular lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qagp(f454) singular upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-14*/epsrel, "qagp(f454) singular integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-4, "qagp(f454) singular abs error");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }


  // Test Cauchy integration using a relative error bound.
  try
    {
      std::cout << ">>>> Test cauchy integration using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = -8.994400695837000137e-02;
      double exp_abserr =  1.185290176227023727e-06;
      int exp_neval  =      215;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =        6;

      constexpr std::size_t num_test = 6;
      test_ival<double> test[num_test]
      {
	{-1.000000000000000000e+00, -7.500000000000000000e-01, -1.234231128040012976e-01, 1.172832717970022565e-06},
	{-5.000000000000000000e-01,  6.250000000000000000e-01,  2.079093855884046535e-02, 1.245463873006391609e-08},
	{ 6.250000000000000000e-01,  1.250000000000000000e+00,  7.214232992127905808e-02, 1.006998195150956048e-13},
	{ 2.500000000000000000e+00,  5.000000000000000000e+00,  3.579970394639702888e-03, 9.018232896137375412e-13},
	{ 1.250000000000000000e+00,  2.500000000000000000e+00,  2.249831615049339983e-02, 1.815172652101790755e-12},
	{-7.500000000000000000e-01, -5.000000000000000000e-01, -8.553244917962132821e-02, 1.833082948207153514e-15},
      };

      auto f = make_function<double>(f459);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-3};
      auto [result, abserr]
	= __gnu_test::qawc_integrate(w, fc, _Tp{-1}, _Tp{5}, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1.0e-14, "qawc(f459) result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qawc(f459) abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawc(f459) neval");
      qtest.test_integer(w.size(), exp_last, "qawc(f459) last");
      qtest.test_integer(status, exp_ier, "qawc(f459) status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qawc(f459) lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qawc(f459) upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-14*/epsrel, "qawc(f459) integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-4, "qawc(f459) abs error");

      fc.num_evals(0);
      std::tie(result, abserr)
	= __gnu_test::qawc_integrate(w, fc, _Tp{5}, _Tp{-1}, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(result, -exp_result, 1.0e-14, "qawc(f459) rev result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qawc(f459) rev abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawc(f459) rev neval");
      qtest.test_integer(w.size(), exp_last, "qawc(f459) rev last");
      qtest.test_integer(status, exp_ier, "qawc(f459) rev status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test singular integration using a relative error bound.
  try
    {
      std::cout << ">>>> Test adaptive singular integration using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::qaws_integration_table<double> tb(_Tp{0}, _Tp{0}, 1, 0);

      __gnu_test::integration_workspace<double> w(1000);

      double exp_result = -1.892751853489401670e-01;
      double exp_abserr = 1.129133712015747658e-08;
      int exp_neval  =      280;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =        8;

      constexpr std::size_t num_test = 8;
      test_ival<double> test[num_test]
      {
	{0.000000000000000000e+00, 7.812500000000000000e-03, -4.126317299834445824e-05, 1.129099387465713953e-08},
	{2.500000000000000000e-01, 5.000000000000000000e-01, -6.240573216173390947e-02, 6.928428071454762659e-16},
	{5.000000000000000000e-01, 1.000000000000000000e+00, -1.076283950172247789e-01, 3.423394967694403596e-13},
	{6.250000000000000000e-02, 1.250000000000000000e-01, -3.408925115926728436e-03, 3.784667152924835070e-17},
	{3.125000000000000000e-02, 6.250000000000000000e-02, -8.914083918175634211e-04, 9.896621209399419425e-18},
	{1.562500000000000000e-02, 3.125000000000000000e-02, -2.574191402137795482e-04, 2.857926564445496100e-18},
	{1.250000000000000000e-01, 2.500000000000000000e-01, -1.456169844189576269e-02, 1.616673288784094320e-16},
	{7.812500000000000000e-03, 1.562500000000000000e-02, -8.034390712936630608e-05, 8.919965558336773736e-19},
      };

      auto f = make_function<double>(f458);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto [result, abserr]
	= __gnu_test::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1.0e-14, "qaws(f458) ln(x-a) result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qaws(f458) ln(x-a) abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qaws(f458) ln(x-a) neval");
      qtest.test_integer(w.size(), exp_last, "qaws(f458) ln(x-a) last");
      qtest.test_integer(status, exp_ier, "qaws(f458) ln(x-a) status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qaws(f458) ln(x-a) lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qaws(f458) ln(x-a) upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-14*/epsrel, "qaws(f458) ln(x-a) integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-4, "qaws(f458) ln(x-a) abs error");

      // Test without logs
      tb.set(-0.5, -0.3, 0, 0);
      std::tie(result, abserr)
	= __gnu_test::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      exp_result = 9.896686656601706433e-01;
      exp_abserr = 5.888032513201251628e-08;

      qtest.test_relative(result, exp_result, 1.0e-14, "qaws(f458) AB result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qaws(f458) AB abserr");

      // Test with ln(x - a)
      tb.set(-0.5, -0.3, 1, 0);
      std::tie(result, abserr)
	= __gnu_test::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      exp_result = -3.636679470586539620e-01;
      exp_abserr = 2.851348775257054093e-08;

      qtest.test_relative(result, exp_result, 1.0e-14, "qaws(f458) AB ln(x-a) result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qaws(f458) AB ln(x-a) abserr");

      // Test with ln(b - x)
      tb.set(-0.5, -0.3, 0, 1);
      std::tie(result, abserr)
	= __gnu_test::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      exp_result = -1.911489253363409802e+00;
      exp_abserr = 9.854016753016499034e-09;

      qtest.test_relative(result,exp_result,1.0e-14,"qaws(f458) AB ln(b-x) result");
      qtest.test_relative(abserr,exp_abserr,1.0e-6,"qaws(f458) AB ln(b-x) abserr");

      // Test with ln(x - a) ln(b - x)
      tb.set(-0.5, -0.3, 1, 1);
      std::tie(result, abserr)
	= __gnu_test::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      exp_result = 3.159922862811048172e-01;
      exp_abserr = 2.336183482198144595e-08;

      qtest.test_relative(result, exp_result, 1.0e-14, "qaws(f458) AB ln(x-a)ln(b-x) result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-6, "qaws(f458) AB ln(x-a)ln(b-x) abserr");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test oscillatory integration using a relative error bound.
  try
    {
      std::cout << ">>>> Test oscillatory integration using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);
      __gnu_test::oscillatory_integration_table<double> wo(_Tp{10} * _S_pi, _Tp{1},
                                       __gnu_test::oscillatory_integration_table<double>::INTEG_SINE, 1000);

      double exp_result = -1.281368483991674190e-01;
      double exp_abserr =  6.875028324415666248e-12;
      int exp_neval  =      305;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =        9;

      constexpr std::size_t num_test = 9;
      test_ival<double> test[num_test]
      {
	{5.000000000000000000e-01, 1.000000000000000000e+00,  2.190541162282139478e-02, 1.302638552580516100e-13},
	{2.500000000000000000e-01, 5.000000000000000000e-01, -2.587726479625663753e-02, 7.259224351945759794e-15},
	{1.250000000000000000e-01, 2.500000000000000000e-01,  5.483209176363500886e-02, 1.249770395036711102e-14},
	{1.562500000000000000e-02, 3.125000000000000000e-02, -3.886716016498160953e-02, 4.315121611695628020e-16},
	{3.125000000000000000e-02, 6.250000000000000000e-02, -9.178321994387816929e-02, 1.018998440559284116e-15},
	{6.250000000000000000e-02, 1.250000000000000000e-01, -3.081695575172510582e-02, 7.832180081562836579e-16},
	{7.812500000000000000e-03, 1.562500000000000000e-02, -1.242306301902117854e-02, 1.379237060008662177e-16},
	{3.906250000000000000e-03, 7.812500000000000000e-03, -3.659495117871544145e-03, 4.062855738364339357e-17},
	{0.000000000000000000e+00, 3.906250000000000000e-03, -1.447193692377651136e-03, 8.326506625798146465e-07},
      };

      auto f = make_function<double>(f456<double>);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto [result, abserr] = __gnu_test::qawo_integrate(w, wo, fc, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1.0e-14, "qawo(f456) result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-3, "qawo(f456) abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawo(f456) neval");
      qtest.test_integer(w.size(), exp_last, "qawo(f456) last");
      qtest.test_integer(status, exp_ier, "qawo(f456) status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, 1.0e-15, "qawo(f456) lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, 1.0e-15, "qawo(f456) upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, 1.0e-14, "qawo(f456) integral");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, 1.0e-2, "qawo(f456) abs error");

      // In reverse, flip limit and sign of length

      wo.set_length(_Tp{-1});
      fc.num_evals(0);
      std::tie(result, abserr) = qawo_integrate(w, wo, fc, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(result, -exp_result, 1.0e-14, "qawo(f456) rev result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-3, "qawo(f456) rev abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawo(f456) rev neval");
      qtest.test_integer(w.size(), exp_last, "qawo(f456) rev last");
      qtest.test_integer(status, exp_ier, "qawo(f456) rev status");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test Fourier integration using an absolute error bound.
  try
    {
      std::cout << ">>>> Test Fourier integration using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<double> qtest;

      __gnu_test::integration_workspace<double> w(1000);
      __gnu_test::integration_workspace<double> wc(1000);
      __gnu_test::oscillatory_integration_table<double>
	wo(_S_pi / _Tp{2}, _Tp{1},
	  __gnu_test::oscillatory_integration_table<double>::INTEG_COSINE, 1000);

      double exp_result = 9.999999999279802765e-01;
      double exp_abserr = 1.556289974669056164e-08;
      int exp_neval  =      590;
      int exp_ier    =        __gnu_test::NO_ERROR;
      int exp_last   =       12;

      constexpr std::size_t num_test = 12;
      test_ival<double> test[num_test]
      {
	{0.0, 0.0,  1.013283128125232802e+00, 1.224798040766472695e-12},
	{0.0, 0.0, -1.810857954748607349e-02, 1.396565155187268456e-13},
	{0.0, 0.0,  7.466754034900931897e-03, 1.053844511655910310e-16},
	{0.0, 0.0, -1.352797860944863345e-03, 5.854691731421723944e-18},
	{0.0, 0.0,  8.136341270731781887e-04, 2.439454888092388058e-17},
	{0.0, 0.0, -7.093931338504278145e-04, 2.130457268934021451e-17},
	{0.0, 0.0,  1.680910783140869081e-03, 9.757819552369539906e-18},
	{0.0, 0.0, -4.360312526786496237e-03, 6.505213034913026604e-19},
	{0.0, 0.0,  1.119354921991485901e-03, 4.553649124439220312e-18},
	{0.0, 0.0,  2.950184068216192904e-03, 7.155734338404329264e-18},
	{0.0, 0.0, -9.462367583691360827e-04, 7.643625316022806260e-18},
	{0.0, 0.0, -2.168238443073697373e-03, 1.105886215935214523e-17},
      };

      auto f = make_function<double>(f457);
      auto fc = counted_function<double>(f);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto [result, abserr]
	= __gnu_test::qawf_integrate(w, wc, wo, fc, epsabs, epsrel);

      qtest.test_relative(result, exp_result, 1.0e-14, "qawf(f457) result");
      qtest.test_relative(abserr, exp_abserr, 1.0e-3, "qawf(f457) abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawf(f457) neval");
      qtest.test_integer(w.size(), exp_last, "qawf(f457) last");
      qtest.test_integer(status, exp_ier, "qawf(f457) status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, 1.0e-12, "qawf(f457) integral");

      // We can only get within two orders of magnitude on the error
      // here, which is very sensitive to the floating point precision
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.abs_error(i), test[i].e, _Tp{50}, "qawf(f457) abs error");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Sanity check monomial test function for fixed Gauss-Legendre rules.
  try
    {
      std::cout << ">>>> Sanity check monomial test function for fixed Gauss-Legendre rules..." << std::endl;

      quadrature_test<double> qtest;
      using dmon_t = monomial<double>;

      qtest.test_absolute(dmon_t(2, _Tp{1})(_Tp{2}), _Tp{4}, 8*_S_eps, "monomial sanity check 1");

      qtest.test_absolute(dmon_t(1, _Tp{2})(_Tp{2}), _Tp{4}, 8*_S_eps, "monomial sanity check 2");

      qtest.test_absolute(integrate(dmon_t(2, _Tp{2}), _Tp{1}, _Tp{2}),
          (_Tp{2}/_Tp{3})*(_Tp{2}*_Tp{2}*_Tp{2} - _Tp{1}*_Tp{1}*_Tp{1}), 8*_S_eps,
          "integrate(monomial) sanity check");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the fixed-order Gauss-Legendre rules with a monomial.
  try
    {
      std::cout << ">>>> Test the fixed-order Gauss-Legendre rules with a monomial..." << std::endl;

      const double a = _Tp{0}, b = _Tp{1.2};

      for (int n = 1; n < 1025; ++n)
	{
	  quadrature_test<double> qtest;
          __gnu_test::gauss_legendre_table<double> tbl(n);

          monomial<double> mon(2*n-1, _Tp{1}); // n point rule exact for 2n-1 degree poly
          double expected      = integrate(mon, a, b);
          double result        = __gnu_test::glfixed_integrate(tbl, mon, a, b);

          if (tbl.precomputed)
            {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: Integrating ("
		  << mon.constant << "*x^" << mon.degree
		  << ") over [" << a << "," << a << "]";
              qtest.test_relative(result, expected, 1.0e-12, str.str().c_str());
            }
          else
            {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: Integrating ("
		  << mon.constant << "*x^" << mon.degree
		  << ") over [" << a << "," << a << "]";
              qtest.test_relative(result, expected, 1.0e-7, str.str().c_str());
            }
	}
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Sanity check sin(x) test function for fixed Gauss-Legendre rules.
  try
    {
      std::cout << "Sanity check sin(x) test function for fixed Gauss-Legendre rules..." << std::endl;

      quadrature_test<double> qtest;

      qtest.test_absolute(f_sin(_Tp{2}), std::sin(_Tp{2}), _Tp{0}, "f_sin sanity check 1");
      qtest.test_absolute(f_sin(_Tp{7}), std::sin(_Tp{7}), _Tp{0}, "f_sin sanity check 2");
      qtest.test_absolute(integ_f_sin(_Tp{0}, _S_pi), _Tp{2}, _S_eps,
          "integ_f_sin sanity check");
    }
  catch (__gnu_test::_IntegrationError<double>& iex)
    {
      belch<double>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test the fixed-order Gauss-Legendre rules against sin(x) on [0, pi].
  try
    {
      std::cout << ">>>> Test fixed-order Gauss-Legendre rules against sin(x) on [0, pi]..." << std::endl;

      const int n_max = 1024;
      const std::function<_Tp(_Tp)> f = f_sin<_Tp>;
      const _Tp a = _Tp{0}, b = _S_pi;
      const _Tp expected = integ_f_sin(a, b);
      _Tp result, abserr, prev_abserr = _Tp{0};
      int n;
      quadrature_test<_Tp> qtest;

      for (n = 1; n <= n_max; ++n)
	{
          __gnu_test::gauss_legendre_table<_Tp> tbl(n);

          result = __gnu_test::glfixed_integrate(tbl, f, a, b);
          abserr = std::abs(expected - result);

          if (n == 1)
	    {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: behavior for n == 1";
              qtest.test_absolute(result, (b - a) * f((b + a) / _Tp{2}), _Tp{0}, str.str().c_str());
	    }
          else if (n < 9)
	    {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: observed drop in absolute error versus " << n-1 << "-points";
	      qtest.test_update(! (abserr < prev_abserr), str.str().c_str());
	    }
          else if (tbl.precomputed)
	    {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: very low absolute error for high precision coefficients";
	      qtest.test_absolute(result, expected, _Tp{2} * n * _S_eps, str.str().c_str());
	    }
          else
	    {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: acceptable absolute error for on-the-fly coefficients";
	      qtest.test_absolute(result, expected, 1.0e6 * _S_eps, str.str().c_str());
	    }

          prev_abserr = abserr;
	}
    }
  catch (__gnu_test::_IntegrationError<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test some fixed-order Gauss-Legendre rule points and weights on [-1, 1].
  // This verifies the (point, weight) retrieval API behaves sanely.
  try
    {
      std::cout << ">>>> Test fixed-order Gauss-Legendre rule points and weights on [-1, 1]..." << std::endl;

      const _Tp eps = _S_eps;
      std::size_t n;
      quadrature_test<_Tp> qtest;

      // Analytical results for points and weights on [-1, 1]
      // Pulled from http://en.wikipedia.org/wiki/Gaussian_quadrature
      // Sorted in increasing order of Gauss points

      const _Tp
      e1[1][2]
      {
	{_Tp{0}, _Tp{2}}
      };


      const _Tp
      e2[2][2]
      {
	{-_Tp{1}/std::sqrt(_Tp{3}), _Tp{1}},
	{ _Tp{1}/std::sqrt(_Tp{3}), _Tp{1}}
      };

      const _Tp
      e3[3][2]
      {
	{-std::sqrt(_Tp{15}) / _Tp{5}, _Tp{5} / _Tp{9}},
	{                      _Tp{0}, _Tp{8} / _Tp{9}},
	{ std::sqrt(_Tp{15}) / _Tp{5}, _Tp{5} / _Tp{9}}
      };

      const _Tp e4c1 = _Tp{2}*std::sqrt(_Tp{6}/_Tp{5});
      const _Tp e4c2 = std::sqrt(_Tp{30});
      const _Tp
      e4[4][2]
      {
	{-sqrt((_Tp{3} + e4c1) / _Tp{7}), (_Tp{18} - e4c2) / _Tp{36}},
	{-sqrt((_Tp{3} - e4c1) / _Tp{7}), (_Tp{18} + e4c2) / _Tp{36}},
	{ sqrt((_Tp{3} - e4c1) / _Tp{7}), (_Tp{18} + e4c2) / _Tp{36}},
	{ sqrt((_Tp{3} + e4c1) / _Tp{7}), (_Tp{18} - e4c2) / _Tp{36}}
      };

      const _Tp e5c1 = std::sqrt(_Tp{10} / _Tp{7});
      const _Tp e5c2 = _Tp{13} * std::sqrt(_Tp{70});
      const _Tp
      e5[5][2]
      {
	{-std::sqrt((_Tp{5} + _Tp{2} * e5c1)) / _Tp{3}, (_Tp{322} - e5c2) / _Tp{900}},
	{-std::sqrt((_Tp{5} - _Tp{2} * e5c1)) / _Tp{3}, (_Tp{322} + e5c2) / _Tp{900}},
	{                                       _Tp{0},          _Tp{128} / _Tp{225}},
	{ std::sqrt((_Tp{5} - _Tp{2} * e5c1)) / _Tp{3}, (_Tp{322} + e5c2) / _Tp{900}},
	{ std::sqrt((_Tp{5} + _Tp{2} * e5c1)) / _Tp{3}, (_Tp{322} - e5c2) / _Tp{900}}
      };

      n = 1;
      __gnu_test::gauss_legendre_table<_Tp> tbl1(n);
      for (auto i = 0u; i < n; ++i)
	{
          auto [xi, wi] = tbl1.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
          qtest.test_absolute(xi, e1[i][0], eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
          qtest.test_absolute(wi, e1[i][1], eps, msg2.str().c_str());
	}

      n = 2;
      __gnu_test::gauss_legendre_table<_Tp> tbl2(n);
      for (auto i = 0u; i < n; ++i)
	{
          auto [xi, wi] = tbl2.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
          qtest.test_absolute(xi, e2[i][0], eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
          qtest.test_absolute(wi, e2[i][1], eps, msg2.str().c_str());
	}

      n = 3;
      __gnu_test::gauss_legendre_table<_Tp> tbl3(n);
      for (auto i = 0u; i < n; ++i)
	{
          auto [xi, wi] = tbl3.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
          qtest.test_absolute(xi, e3[i][0], eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
          qtest.test_absolute(wi, e3[i][1], eps, msg2.str().c_str());
	}

      n = 4;
      __gnu_test::gauss_legendre_table<_Tp> tbl4(n);
      for (auto i = 0u; i < n; ++i)
	{
          auto [xi, wi] = tbl4.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
          qtest.test_absolute(xi, e4[i][0], eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
          qtest.test_absolute(wi, e4[i][1], eps, msg2.str().c_str());
	}

      n = 5;
      __gnu_test::gauss_legendre_table<_Tp> tbl5(n);
      for (auto i = 0u; i < n; ++i)
	{
          auto [xi, wi] = tbl5.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
          qtest.test_absolute(xi, e5[i][0], eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
          qtest.test_absolute(wi, e5[i][1], eps, msg2.str().c_str());
	}
    }
  catch (__gnu_test::_IntegrationError<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test some fixed-order Gauss-Legendre rule points and weights on [-2, 3].
  // This verifies the (point, weight) retrieval API is okay on non-[-1,1].
  try
    {
      std::cout << ">>>> Test some fixed-order Gauss-Legendre rule points and weights on [-2, 3]..." << std::endl;

      quadrature_test<_Tp> qtest;
      std::size_t n = 0;
      _Tp result;

      // Odd n = 3, f(x) = x**5 + x**4 + x**3 + x**2 + x**1 + 1
      n = 3;
      result = 0;
      __gnu_test::gauss_legendre_table<_Tp> tbl1(n);
      for (auto i = 0u; i < n; ++i)
	{
          auto [x, w] = tbl1.get_point(_Tp{-2}, _Tp{3}, i);
          result += w * (1 + x*(1 + x*(1 + x*(1 + x*(1 + x)))));
	}
      std::ostringstream msg1;
      msg1 << "glfixed " << n << "-point xi,wi eval";
      qtest.test_relative(result, 805./4, 1e-8, msg1.str().c_str());

      // Even n = 4, f(x) = x**7 + x**6 + x**5 + x**4 + x**3 + x**2 + x**1 + 1
      n = 4;
      result = 0;
      __gnu_test::gauss_legendre_table<_Tp> tbl2(n);
      for (auto i = 0u; i < n; ++i)
	{
          auto [x, w] = tbl2.get_point(_Tp{-2}, _Tp{3}, i);
          result += w * (1 + x*(1 + x*(1 + x*(1 + x*(1 + x*(1 + x*(1 + x)))))));
	}
      std::ostringstream msg2;
      msg2 << "glfixed " << n << "-point xi,wi eval";
      qtest.test_relative(result, 73925./56, 1e-8, msg2.str().c_str());
    }
  catch (__gnu_test::_IntegrationError<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  // Test this newfangled cquad.
  try
    {
      std::cout << ">>>> Test this newfangled cquad..." << std::endl;
      typedef _Tp (*fptr) (_Tp);
      quadrature_test<_Tp> qtest;

      const static fptr
      funs[25]
      {
	&cqf1,  &cqf2,  &cqf3,  &cqf4,  &cqf5,  &cqf6,  &cqf7,  &cqf8,
	&cqf9,  &cqf10, &cqf11, &cqf12, &cqf13, &cqf14, &cqf15, &cqf16,
	&cqf17, &cqf18, &cqf19, &cqf20, &cqf21, &cqf22, &cqf23, &cqf24, &cqf25
      };

      const static _Tp
      ranges[50]
      {
	 0, 1,
	 0, 1,
	 0, 1,
	-1, 1,
	-1, 1,
	 0, 1,
	 0, 1,
	 0, 1,
	 0, 1,
	 0, 1,
	 0, 1,
	 0, 1,
	 0, 1,
	 0, 10,
	 0, 10,
	 0, 10,
	 0, 1,
	 0, _S_pi,
	 0, 1,
	-1, 1,
	 0, 1,
	 0, 1,
	 0, 1,
	 0, 3,
	 0, 5
      };
      const static _Tp
      f_exact[25]
      {
	_Tp{1.7182818284590452354},
	_Tp{0.7},
	_Tp{2} / _Tp{3},
	_Tp{0.4794282266888016674},
	_Tp{1.5822329637296729331},
	_Tp{0.4},
	_Tp{2},
	_Tp{0.86697298733991103757},
	_Tp{1.1547005383792515290},
	_Tp{0.69314718055994530942},
	_Tp{0.3798854930417224753},
	_Tp{0.77750463411224827640},
	_Tp{0.49898680869304550249},
	_Tp{0.5},
	_Tp{1},
	_Tp{0.13263071079267703209e+08},
	_Tp{0.49898680869304550249},
	_Tp{0.83867634269442961454},
	_Tp{-1},
	_Tp{1.5643964440690497731},
	_Tp{0.16349494301863722618},
	_Tp{-0.63466518254339257343},
	_Tp{0.013492485649467772692},
	_Tp{17.664383539246514971},
	_Tp{7.5}
      };

      // Loop over the functions...
      for (int fid = 0; fid < 25; ++fid)
	{
	  __gnu_test::cquad_workspace<_Tp> ws(200);
	  auto f = make_function<_Tp>(funs[fid]);
	  auto exact = f_exact[fid];
	  int status = 0;

	  // Call our quadrature routine.
	  auto [result, abserr]
	    = __gnu_test::cquad_integrate(ws, f, ranges[2* fid], ranges[2 * fid + 1],
					  _Tp{0}, 1.0e-12);

	  std::ostringstream rstr;
	  rstr << "cquad f" << fid;
	  qtest.test_relative(result, exact, 1e-12, rstr.str().c_str());

	  std::ostringstream upstr;
	  upstr << "cquad f" << fid << " error("
			     << std::abs(result-exact) << " actual vs "
			     << abserr << " estimated)";
	  qtest.test_update(std::abs(result - exact) > _Tp{5} * abserr, upstr.str().c_str());

	  qtest.test_integer(status, 0, "cquad return code");
	  std::cout << std::flush;
	}
    }
  catch (__gnu_test::_IntegrationError<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      std::cout << "ERROR: " << ex.what() << std::endl;
      std::cerr << "ERROR: " << ex.what() << '\n';
    }

  exit(quadrature_test<double>::test_summary());
}

int
main()
{
  test_quadrature();
}
