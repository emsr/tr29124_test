/* integration/test.c
 *
 * Copyright (C) 1996-2000, 2007 Brian Gough
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
    : func(f), neval(0)
    { }

    _Tp
    operator()(_Tp x) const
    {
      ++this->neval;
      return this->func(x);
    }

    std::function<_Tp(_Tp)> func;
    mutable int neval;
  };

enum
{
  ROUND_ERROR,
  SINGULAR_ERROR,
  MAX_ITER_ERROR,
  TOLERANCE_ERROR
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
    void test_rel(_Tp result, _Tp expected, _Tp rel_error, const char* test_desc);
    void test_abs(_Tp result, _Tp expected, _Tp abs_error, const char* test_desc);
    void test_factor(_Tp result, _Tp expected, _Tp factor, const char* test_desc);
    void test_int(int result, int expected, const char* test_desc);

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

    if (status ||  this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status && ! this->verbose)
          printf(" [%u]", this->num_tests);

	 printf(" %s", test_desc);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_rel(_Tp result, _Tp expected, _Tp rel_error,
		      const char* test_desc)
  {
    int status;

    /* Check for NaN vs inf vs number */

    if (std::isnan(result) || std::isnan(expected))
      status = std::isnan(result) != std::isnan(expected);
    else if (std::isinf(result) || std::isinf(expected))
      status = std::isinf(result) != std::isinf(expected);
    else if ((expected > 0 && expected < std::numeric_limits<_Tp>::min())
             || (expected < 0 && expected > -(std::numeric_limits<_Tp>::min())))
      status = -1;
    else if (expected != 0)
      status = (std::abs(result-expected)/std::abs(expected) > rel_error);
    else
      status = (std::abs(result) > rel_error);

    this->test_update(status);

    if (status ||  this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status == 0)
          {
            if (strlen(test_desc) < 45)
              printf(" (%g observed vs %g expected)", result, expected);
            else
              printf(" (%g obs vs %g exp)", result, expected);
          }
	else
          printf(" (%.18g observed vs %.18g expected)", result, expected);

	if (status == -1)
          printf(" [test uses subnormal value]");

	if (status && ! this->verbose)
          printf(" [%u]",  this->num_tests);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_abs(_Tp result, _Tp expected, _Tp abs_error,
		      const char* test_desc)
  {
    int status;

    /* Check for NaN vs inf vs number */

    if (std::isnan(result) || std::isnan(expected))
      status = std::isnan(result) != std::isnan(expected);
    else if (std::isinf(result) || std::isinf(expected))
      status = std::isinf(result) != std::isinf(expected);
    else if ((expected > 0 && expected < std::numeric_limits<_Tp>::min())
             || (expected < 0 && expected > -(std::numeric_limits<_Tp>::min())))
      status = -1;
    else
      status = std::abs(result-expected) > abs_error;

    this->test_update(status);

    if (status ||  this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status == 0)
          {
            if (strlen(test_desc) < 45)
              printf(" (%g observed vs %g expected)", result, expected);
            else
              printf(" (%g obs vs %g exp)", result, expected);
          }
	else
          printf(" (%.18g observed vs %.18g expected)", result, expected);

	if (status == -1)
          printf(" [test uses subnormal value]");

	if (status && ! this->verbose)
          printf(" [%u]",  this->num_tests);

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
    else if (expected == 0.0)
      status = (result > expected || result < expected);
    else
      {
	_Tp u = result / expected;
	status = (u > factor || u < 1.0 / factor);
      }

    this->test_update (status);

    if (status ||  this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status == 0)
          {
            if (strlen(test_desc) < 45)
              printf(" (%g observed vs %g expected)", result, expected);
            else
              printf(" (%g obs vs %g exp)", result, expected);
          }
	else
          printf(" (%.18g observed vs %.18g expected)", result, expected);

	if (status == -1)
          printf(" [test uses subnormal value]");

	if (status && ! this->verbose)
          printf(" [%u]",  this->num_tests);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_int(int result, int expected, const char*/*test_desc*/)
  {
    int status = (result != expected);

    this->test_update(status);

    if (status ||  this->verbose)
      {
	printf(status ? "FAIL: " : "PASS: ");

	if (status == 0)
          printf(" (%d observed vs %d expected)", result, expected);
	else
          printf(" (%d observed vs %d expected)", result, expected);

	if (status && ! this->verbose)
          printf(" [%u]",  this->num_tests);

	printf("\n");
	fflush(stdout);
      }
  }

template<typename _Tp>
  int
  quadrature_test<_Tp>::test_summary()
  {
    if (tot_verbose)
      printf("%d tests, passed %d, failed %d.\n",  tot_num_tests,  tot_num_passed,  tot_num_failed);

    if (tot_num_failed != 0)
      {
	return FAILURE;
      }

    if (tot_num_tests !=  tot_num_passed +  tot_num_failed)
      {
	if (tot_verbose)
          printf("TEST RESULTS DO NOT ADD UP %d != %d + %d\n",
                   tot_num_tests,  tot_num_passed,  tot_num_failed);
	return FAILURE;
      }

    if (tot_num_passed ==  tot_num_tests)
      {
	if (! tot_verbose)
          printf("Completed [%d/%d]\n",  tot_num_passed,  tot_num_tests);

	return SUCCESS;
      }

    return FAILURE;
  }


template<typename _Tp>
  void
  belch(const __gnu_test::_IntegrationError<double>& iex)
  {
    std::cout << "ERROR: " << iex.what()
	      << "       status = " << iex.status()
	      << "       result = " << iex.result()
	      << "       abserr = " << iex.abserr()
	      << std::endl;
    std::cerr << "ERROR: " << iex.what()
	      << "       status = " << iex.status()
	      << "       result = " << iex.result()
	      << "       abserr = " << iex.abserr()
	      << '\n';
  }


int
main()
{
  // Test the basic Gauss-Kronrod rules with a smooth positive function.
  const auto _S_eps = std::numeric_limits<double>::epsilon();
  try
  {
    std::cout << ">>>> Test Gauss-Kronrod 15 with a smooth positive function..." << std::endl;

    double exp_result = 7.716049357767090777E-02;
    double exp_abserr = 2.990224871000550874E-06;
    double exp_resabs = 7.716049357767090777E-02;
    double exp_resasc = 4.434273814139995384E-02;
    quadrature_test<double> qtest;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = qk_integrate(f, 0.0, 1.0, __gnu_test::QK_15);
    qtest.test_rel(result,exp_result,1e-15,"qk15(f1) smooth result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk15(f1) smooth abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk15(f1) smooth resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk15(f1) smooth resasc");

    std::tie(result, abserr, resabs, resasc)
      = qk_integrate(f, 1.0, 0.0, __gnu_test::QK_15);

    qtest.test_rel(result,-exp_result,1e-15,"qk15(f1) reverse result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk15(f1) reverse abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk15(f1) reverse resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk15(f1) reverse resasc");
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

    double exp_result = 7.716049379303084599E-02;
    double exp_abserr = 9.424302194248481445E-08;
    double exp_resabs = 7.716049379303084599E-02;
    double exp_resasc = 4.434311425038358484E-02;
    quadrature_test<double> qtest;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_21);
    qtest.test_rel(result,exp_result,1e-15,"qk21(f1) smooth result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk21(f1) smooth abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk21(f1) smooth resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk21(f1) smooth resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_21);
    qtest.test_rel(result,-exp_result,1e-15,"qk21(f1) reverse result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk21(f1) reverse abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk21(f1) reverse resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk21(f1) reverse resasc");
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

    double exp_result = 7.716049382494900855E-02;
    double exp_abserr = 1.713503193600029893E-09;
    double exp_resabs = 7.716049382494900855E-02;
    double exp_resasc = 4.427995051868838933E-02;
    quadrature_test<double> qtest;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_31);
    qtest.test_rel(result,exp_result,1e-15,"qk31(f1) smooth result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk31(f1) smooth abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk31(f1) smooth resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk31(f1) smooth resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_31);
    qtest.test_rel(result,-exp_result,1e-15,"qk31(f1) reverse result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk31(f1) reverse abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk31(f1) reverse resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk31(f1) reverse resasc");
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

    double exp_result = 7.716049382681375302E-02;
    double exp_abserr = 9.576386660975511224E-11;
    double exp_resabs = 7.716049382681375302E-02;
    double exp_resasc = 4.421521169637691873E-02;
    quadrature_test<double> qtest;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_41);
    qtest.test_rel(result,exp_result,1e-15,"qk41(f1) smooth result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk41(f1) smooth abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk41(f1) smooth resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk41(f1) smooth resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_41);
    qtest.test_rel(result,-exp_result,1e-15,"qk41(f1) reverse result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk41(f1) reverse abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk41(f1) reverse resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk41(f1) reverse resasc");
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

    double exp_result = 7.716049382708510540E-02;
    double exp_abserr = 1.002079980317363772E-11;
    double exp_resabs = 7.716049382708510540E-02;
    double exp_resasc = 4.416474291216854892E-02;
    quadrature_test<double> qtest;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_51);
    qtest.test_rel(result,exp_result,1e-15,"qk51(f1) smooth result");
    qtest.test_rel(abserr,exp_abserr,1e-5,"qk51(f1) smooth abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk51(f1) smooth resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk51(f1) smooth resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_51);
    qtest.test_rel(result,-exp_result,1e-15,"qk51(f1) reverse result");
    qtest.test_rel(abserr,exp_abserr,1e-5,"qk51(f1) reverse abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk51(f1) reverse resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk51(f1) reverse resasc");
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

    double exp_result = 7.716049382713800753E-02;
    double exp_abserr = 1.566060362296155616E-12;
    double exp_resabs = 7.716049382713800753E-02;
    double exp_resasc = 4.419287685934316506E-02;
    quadrature_test<double> qtest;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_61);
    qtest.test_rel(result,exp_result,1e-15,"qk61(f1) smooth result");
    qtest.test_rel(abserr,exp_abserr,1e-5,"qk61(f1) smooth abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk61(f1) smooth resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk61(f1) smooth resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_61);
    qtest.test_rel(result,-exp_result,1e-15,"qk61(f1) reverse result");
    qtest.test_rel(abserr,exp_abserr,1e-5,"qk61(f1) reverse abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk61(f1) reverse resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk61(f1) reverse resasc");
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

    double exp_result = 1.555688196612745777E+01;
    double exp_abserr = 2.350164577239293706E+01;
    double exp_resabs = 1.555688196612745777E+01;
    double exp_resasc = 2.350164577239293706E+01;
    quadrature_test<double> qtest;

    double alpha = -0.9;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_15);
    qtest.test_rel(result,exp_result,1e-15,"qk15(f1) singular result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk15(f1) singular abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk15(f1) singular resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk15(f1) singular resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_15);
    qtest.test_rel(result,-exp_result,1e-15,"qk15(f1) reverse result");
    qtest.test_rel(abserr,exp_abserr,1e-7,"qk15(f1) reverse abserr");
    qtest.test_rel(resabs,exp_resabs,1e-15,"qk15(f1) reverse resabs");
    qtest.test_rel(resasc,exp_resasc,1e-15,"qk15(f1) reverse resasc");
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

    double exp_result = 1.799045317938126232E+01;
    double exp_abserr = 2.782360287710622515E+01;
    double exp_resabs = 1.799045317938126232E+01;
    double exp_resasc = 2.782360287710622515E+01;
    quadrature_test<double> qtest;

    double alpha = -0.9;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_21);
    qtest.test_rel(result, exp_result, 1e-15, "qk21(f1) singular result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk21(f1) singular abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk21(f1) singular resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk21(f1) singular resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_21);
    qtest.test_rel(result, -exp_result, 1e-15, "qk21(f1) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk21(f1) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk21(f1) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk21(f1) reverse resasc");
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

    double exp_result = 2.081873305159121657E+01;
    double exp_abserr = 3.296500137482590276E+01;
    double exp_resabs = 2.081873305159121301E+01;
    double exp_resasc = 3.296500137482590276E+01;
    quadrature_test<double> qtest;

    double alpha = -0.9;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_31);
    qtest.test_rel(result, exp_result, 1e-15, "qk31(f1) singular result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk31(f1) singular abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk31(f1) singular resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk31(f1) singular resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_31);
    qtest.test_rel(result, -exp_result, 1e-15, "qk31(f1) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk31(f1) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk31(f1) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk31(f1) reverse resasc");
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

    double exp_result = 2.288677623903126701E+01;
    double exp_abserr = 3.671538820274916048E+01;
    double exp_resabs = 2.288677623903126701E+01;
    double exp_resasc = 3.671538820274916048E+01;
    quadrature_test<double> qtest;

    double alpha = -0.9;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_41);
    qtest.test_rel(result, exp_result, 1e-15, "qk41(f1) singular result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk41(f1) singular abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk41(f1) singular resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk41(f1) singular resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_41);
    qtest.test_rel(result, -exp_result, 1e-15, "qk41(f1) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk41(f1) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk41(f1) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk41(f1) reverse resasc");
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

    double exp_result = 2.449953612016972215E+01;
    double exp_abserr = 3.967771249391228849E+01;
    double exp_resabs = 2.449953612016972215E+01;
    double exp_resasc = 3.967771249391228849E+01;
    quadrature_test<double> qtest;

    double alpha = -0.9;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_51);
    qtest.test_rel(result, exp_result, 1e-15, "qk51(f1) singular result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk51(f1) singular abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk51(f1) singular resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk51(f1) singular resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_51);
    qtest.test_rel(result, -exp_result, 1e-15, "qk51(f1) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk51(f1) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk51(f1) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk51(f1) reverse resasc");
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

    double exp_result = 2.583030240976628988E+01;
    double exp_abserr = 4.213750493076978643E+01;
    double exp_resabs = 2.583030240976628988E+01;
    double exp_resasc = 4.213750493076978643E+01;
    quadrature_test<double> qtest;

    double alpha = -0.9;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.0, 1.0, __gnu_test::QK_61);
    qtest.test_rel(result, exp_result, 1e-15, "qk61(f1) singular result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk61(f1) singular abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk61(f1) singular resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk61(f1) singular resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 1.0, 0.0, __gnu_test::QK_61);
    qtest.test_rel(result, -exp_result, 1e-15, "qk61(f1) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk61(f1) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk61(f1) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk61(f1) reverse resasc");
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

    double exp_result =-7.238969575483799046E-01;
    double exp_abserr = 8.760080200939757174E-06;
    double exp_resabs = 1.165564172429140788E+00;
    double exp_resasc = 9.334560307787327371E-01;
    quadrature_test<double> qtest;

    double alpha = 1.3;
    auto f = make_function<double>(f3, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_15);
    qtest.test_rel(result, exp_result, 1e-15, "qk15(f3) oscill result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk15(f3) oscill abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk15(f3) oscill resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk15(f3) oscill resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_15);
    qtest.test_rel(result, -exp_result, 1e-15, "qk15(f3) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk15(f3) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk15(f3) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk15(f3) reverse resasc");
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

    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 7.999213141433641888E-11;
    double exp_resabs = 1.150829032708484023E+00;
    double exp_resasc = 9.297591249133687619E-01;
    quadrature_test<double> qtest;

    double alpha = 1.3;
    auto f = make_function<double>(f3, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_21);
    qtest.test_rel(result, exp_result, 1e-15, "qk21(f3) oscill result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qk21(f3) oscill abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk21(f3) oscill resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk21(f3) oscill resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_21);
    qtest.test_rel(result, -exp_result, 1e-15, "qk21(f3) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qk21(f3) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk21(f3) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk21(f3) reverse resasc");
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

    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.285805464427459261E-14;
    double exp_resabs = 1.158150602093290571E+00;
    double exp_resasc = 9.277828092501518853E-01;
    quadrature_test<double> qtest;

    double alpha = 1.3;
    auto f = make_function<double>(f3, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_31);
    qtest.test_rel(result, exp_result, 1e-15, "qk31(f3) oscill result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk31(f3) oscill abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk31(f3) oscill resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk31(f3) oscill resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_31);
    qtest.test_rel(result, -exp_result, 1e-15, "qk31(f3) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk31(f3) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk31(f3) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk31(f3) reverse resasc");
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

    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.286535726271015626E-14;
    double exp_resabs = 1.158808363486595328E+00;
    double exp_resasc = 9.264382258645686985E-01;
    quadrature_test<double> qtest;

    double alpha = 1.3;
    auto f = make_function<double>(f3, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_41);
    qtest.test_rel(result, exp_result, 1e-15, "qk41(f3) oscill result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk41(f3) oscill abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk41(f3) oscill resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk41(f3) oscill resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_41);
    qtest.test_rel(result, -exp_result, 1e-15, "qk41(f3) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk41(f3) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk41(f3) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk41(f3) reverse resasc");
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

    double exp_result =-7.238969575482961938E-01;
    double exp_abserr = 1.285290995039385778E-14;
    double exp_resabs = 1.157687209264406381E+00;
    double exp_resasc = 9.264666884071264263E-01;
    quadrature_test<double> qtest;

    double alpha = 1.3;
    auto f = make_function<double>(f3, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_51);
    qtest.test_rel(result, exp_result, 1e-15, "qk51(f3) oscill result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk51(f3) oscill abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk51(f3) oscill resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk51(f3) oscill resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_51);
    qtest.test_rel(result, -exp_result, 1e-15, "qk51(f3) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk51(f3) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk51(f3) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk51(f3) reverse resasc");
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

    double exp_result =-7.238969575482959717E-01;
    double exp_abserr = 1.286438572027470736E-14;
    double exp_resabs = 1.158720854723590099E+00;
    double exp_resasc = 9.270469641771273972E-01;
    quadrature_test<double> qtest;

    double alpha = 1.3;
    auto f = make_function<double>(f3, alpha);

    auto [result, abserr, resabs, resasc]
      = __gnu_test::qk_integrate(f, 0.3, 2.71, __gnu_test::QK_61);
    qtest.test_rel(result, exp_result, 1e-15, "qk61(f3) oscill result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk61(f3) oscill abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk61(f3) oscill resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk61(f3) oscill resasc");

    std::tie(result, abserr, resabs, resasc)
      = __gnu_test::qk_integrate(f, 2.71, 0.3, __gnu_test::QK_61);
    qtest.test_rel(result, -exp_result, 1e-15, "qk61(f3) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qk61(f3) reverse abserr");
    qtest.test_rel(resabs, exp_resabs, 1e-15, "qk61(f3) reverse resabs");
    qtest.test_rel(resasc, exp_resasc, 1e-15, "qk61(f3) reverse resasc");
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
    double exp_result = 7.716049379303083211E-02;
    double exp_abserr = 9.424302199601294244E-08;
    int exp_neval  =  21;
    int exp_ier    =   0;
    quadrature_test<double> qtest;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, neval]
      = __gnu_test::qng_integrate(f, 0.0, 1.0, 1e-1, 0.0);
    qtest.test_rel(result, exp_result, 1e-15, "qng(f1) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qng(f1) smooth abserr");
    qtest.test_int(neval, exp_neval, "qng(f1) smooth neval");
    qtest.test_int(status, exp_ier, "qng(f1) smooth status");

    std::tie(result, abserr, neval)
      = __gnu_test::qng_integrate(f, 1.0, 0.0, 1e-1, 0.0);
    qtest.test_rel(result, -exp_result, 1e-15, "qng(f1) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qng(f1) reverse abserr");
    qtest.test_int(neval, exp_neval, "qng(f1) reverse neval");
    qtest.test_int(status, exp_ier, "qng(f1) reverse status");
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
    double exp_result = 7.716049382706505200E-02;
    double exp_abserr = 2.666893044866214501E-12;
    int exp_neval  =  43;
    int exp_ier    =   0;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, neval]
      = __gnu_test::qng_integrate(f, 0.0, 1.0, 0.0, 1e-9);
    qtest.test_rel(result, exp_result, 1e-15, "qng(f1) smooth 43pt result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qng(f1) smooth 43pt abserr");
    qtest.test_int(neval, exp_neval, "qng(f1) smooth 43pt neval");
    qtest.test_int(status, exp_ier, "qng(f1) smooth 43pt status");

    std::tie(result, abserr, neval)
      = __gnu_test::qng_integrate(f, 1.0, 0.0, 0.0, 1e-9);
    qtest.test_rel(result, -exp_result, 1e-15, "qng(f1) reverse 43pt result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qng(f1) reverse 43pt abserr");
    qtest.test_int(neval, exp_neval, "qng(f1) reverse 43pt neval");
    qtest.test_int(status, exp_ier, "qng(f1) reverse 43pt status");
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
    double exp_result =-7.238969575482961938E-01;
    double exp_abserr = 1.277676889520056369E-14;
    int exp_neval  =  43;
    int exp_ier    =   0;
    quadrature_test<double> qtest;

    double alpha = 1.3;
    auto f = make_function<double>(f3, alpha);

    auto [result, abserr, neval]
      = __gnu_test::qng_integrate(f, 0.3, 2.71, 0.0, 1e-12);
    qtest.test_rel(result, exp_result, 1e-15, "qnq(f3) oscill result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qng(f3) oscill abserr");
    qtest.test_int(neval, exp_neval, "qng(f3) oscill neval");
    qtest.test_int(status, exp_ier, "qng(f3) oscill status");

    std::tie(result, abserr, neval)
      = __gnu_test::qng_integrate(f, 2.71, 0.3, 0.0, 1e-12);
    qtest.test_rel(result, -exp_result, 1e-15, "qnq(f3) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qng(f3) reverse abserr");
    qtest.test_int(neval, exp_neval, "qng(f3) reverse neval");
    qtest.test_int(status, exp_ier, "qng(f3) reverse status");
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
    double exp_result = 7.716049382716029525E-02;
    double exp_abserr = 8.566535680046930668E-16;
    int exp_neval  =  87;
    int exp_ier    =   0;

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, neval]
      = __gnu_test::qng_integrate(f, 0.0, 1.0, 0.0, 1e-13);
    qtest.test_rel(result, exp_result, 1e-15, "qng(f1) 87pt smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qng(f1) 87pt smooth abserr");
    qtest.test_int(neval, exp_neval, "qng(f1) 87pt smooth neval");
    qtest.test_int(status, exp_ier, "qng(f1) 87pt smooth status");

    std::tie(result, abserr, neval)
      = __gnu_test::qng_integrate(f, 1.0, 0.0, 0.0, 1e-13);
    qtest.test_rel(result, -exp_result, 1e-15, "qng(f1) 87pt reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qng(f1) 87pt reverse abserr");
    qtest.test_int(neval, exp_neval, "qng(f1) 87pt reverse neval");
    qtest.test_int(status, exp_ier, "qng(f1) 87pt reverse status");
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
    double exp_result = 3.222948711817264211E+01;
    double exp_abserr = 2.782360287710622870E+01;
    int exp_neval  =  87;
    int exp_ier    =  TOLERANCE_ERROR;

    double alpha = -0.9;
    auto f = make_function<double>(f1, alpha);

    auto [result, abserr, neval]
      = __gnu_test::qng_integrate(f, 0.0, 1.0, 0.0, 1e-3);
    qtest.test_rel(result, exp_result, 1e-15, "qng(f1) sing beyond 87pt result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qng(f1) sing beyond 87pt abserr");
    qtest.test_int(neval, exp_neval, "qng(f1) sing beyond 87pt neval");
    qtest.test_int(status, exp_ier, "qng(f1) sing beyond 87pt status");

    std::tie(result, abserr, neval)
      = __gnu_test::qng_integrate(f, 1.0, 0.0, 0.0, 1e-3);
    qtest.test_rel(result, -exp_result, 1e-15, "qng(f1) reverse beyond 87pt result");
    qtest.test_rel(abserr, exp_abserr, 1e-7, "qng(f1) rev beyond 87pt abserr");
    qtest.test_int(neval, exp_neval, "qng(f1) rev beyond 87pt neval");
    qtest.test_int(status, exp_ier, "qng(f1) rev beyond 87pt status");
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

    double exp_result = 7.716049382715854665E-02;
    double exp_abserr = 6.679384885865053037E-12;
    int exp_neval  =     165;
    int exp_ier    =       0;
    int exp_last   =       6;

    double a[6] = { 0, 0.5, 0.25, 0.125, 0.0625, 0.03125 };
    double b[6] = { 0.03125, 1, 0.5, 0.25, 0.125, 0.0625 };
    double r[6] = { 3.966769831709074375E-06, 5.491842501998222409E-02,
                    1.909827770934243926E-02, 2.776531175604360531E-03,
                    3.280661030752063693E-04, 3.522704932261797744E-05 };
    double e[6] = { 6.678528276336181873E-12, 6.097169993333454062E-16,
                    2.120334764359736934E-16, 3.082568839745514608E-17,
                    3.642265412331439511E-18, 3.910988124757650942E-19 };
    int order[6] = { 1, 2, 3, 4, 5, 6 };

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);

    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qag_integrate(w, fc, 0.0, 1.0, 0.0, 1e-10, 1000,
				  __gnu_test::QK_15);

    qtest.test_rel(result, exp_result, 1e-15, "qag(f1) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qag(f1) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qag(f1) smooth neval");
    qtest.test_int(w.size(), exp_last, "qag(f1) smooth last");
    qtest.test_int(status, exp_ier, "qag(f1) smooth status");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qag(f1) smooth alist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qag(f1) smooth blist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-15, "qag(f1) smooth rlist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-6, "qag(f1) smooth elist");

    for (int i = 0; i < 6; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qag(f1) smooth order");

    fc.neval = 0;

    std::tie(result, abserr)
      = __gnu_test::qag_integrate(w, fc, 1.0, 0.0, 0.0, 1e-10, 1000,
				  __gnu_test::QK_15);

    qtest.test_rel(result, -exp_result, 1e-15, "qag(f1) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qag(f1) reverse abserr");
    qtest.test_int(fc.neval, exp_neval, "qag(f1) reverse neval");
    qtest.test_int(w.size(), exp_last, "qag(f1) reverse last");
    qtest.test_int(status, exp_ier, "qag(f1) reverse status");
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

    double exp_result = 7.716049382716050342E-02;
    double exp_abserr = 2.227969521869139532E-15;
    int exp_neval  =     315;
    int exp_ier    =       0;
    int exp_last   =       8;

    double a[8] = { 0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625,
                    0.0078125 };
    double b[8] = { 0.0078125, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125,
                    0.015625 };
    double r[8] = { 3.696942726831556522E-08, 5.491842501998223103E-02,
                    1.909827770934243579E-02, 2.776531175604360097E-03,
                    3.280661030752062609E-04, 3.522704932261797744E-05,
                    3.579060884684503576E-06, 3.507395216921808047E-07 };
    double e[8] = { 1.371316364034059572E-15, 6.097169993333454062E-16,
                    2.120334764359736441E-16, 3.082568839745514608E-17,
                    3.642265412331439511E-18, 3.910988124757650460E-19,
                    3.973555800712018091E-20, 3.893990926286736620E-21 };
    int order[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qag_integrate(w, fc, 0.0, 1.0, 1e-14, 0.0, 1000, __gnu_test::QK_21);

    qtest.test_rel(result, exp_result, 1e-15, "qag(f1, 21pt) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qag(f1, 21pt) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qag(f1, 21pt) smooth neval");
    qtest.test_int(w.size(), exp_last, "qag(f1, 21pt) smooth last");
    qtest.test_int(status, exp_ier, "qag(f1, 21pt) smooth status");

    for (int i = 0; i < 8; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qag(f1, 21pt) smooth alist");

    for (int i = 0; i < 8; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qag(f1, 21pt) smooth blist");

    for (int i = 0; i < 8; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-15, "qag(f1, 21pt) smooth rlist");

    for (int i = 0; i < 8; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-6, "qag(f1, 21pt) smooth elist");

    for (int i = 0; i < 8; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qag(f1, 21pt) smooth order");


    fc.neval = 0;
    std::tie(result, abserr) = __gnu_test::qag_integrate(w, fc, 1.0, 0.0, 1e-14, 0.0, 1000, __gnu_test::QK_21);

    qtest.test_rel(result, -exp_result, 1e-15, "qag(f1, 21pt) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qag(f1, 21pt) reverse abserr");
    qtest.test_int(fc.neval, exp_neval, "qag(f1, 21pt) reverse neval");
    qtest.test_int(w.size(), exp_last, "qag(f1, 21pt) reverse last");
    qtest.test_int(status, exp_ier, "qag(f1, 21pt) reverse status");
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

    int status = 0;
    quadrature_test<double> qtest;

    __gnu_test::integration_workspace<double> w(1000);

    double exp_result = -7.238969575482959717E-01;
    double exp_abserr =  1.285805464427459261E-14;
    int exp_neval   =     31;
    int exp_ier     =     ROUND_ERROR;
    int exp_last    =     1;

    double alpha = 1.3;
    auto f = make_function<double>(f3, alpha);

    counted_function<double> fc(f);

    auto [result, abserr]
     = __gnu_test::qag_integrate(w, fc, 0.3, 2.71, 1e-14, 0.0, 1000,
				 __gnu_test::QK_31);

    qtest.test_rel(result, exp_result, 1e-15, "qag(f3, 31pt) oscill result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qag(f3, 31pt) oscill abserr");
    qtest.test_int(fc.neval, exp_neval, "qag(f3, 31pt) oscill neval");
    qtest.test_int(w.size(), exp_last, "qag(f3, 31pt) oscill last");
    qtest.test_int(status, exp_ier, "qag(f3, 31pt) oscill status");

    fc.neval = 0;
    std::tie(result, abserr)
     = __gnu_test::qag_integrate(w, fc, 2.71, 0.3, 1e-14, 0.0, 1000,
				 __gnu_test::QK_31);

    qtest.test_rel(result, -exp_result, 1e-15, "qag(f3, 31pt) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qag(f3, 31pt) reverse abserr");
    qtest.test_int(fc.neval, exp_neval, "qag(f3, 31pt) reverse neval");
    qtest.test_int(w.size(), exp_last, "qag(f3, 31pt) reverse last");
    qtest.test_int(status, exp_ier, "qag(f3, 31pt) reverse status");
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

    int status = 0;
    quadrature_test<double> qtest;

    __gnu_test::integration_workspace<double> w(1000);

    int exp_neval  =     5151;
    int exp_ier    =     SINGULAR_ERROR;
    int exp_last   =     51;

    double alpha = 2.0;
    auto f = make_function<double>(f16, alpha);

    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qag_integrate(w, fc, -1.0, 1.0, 1e-14, 0.0, 1000,
				  __gnu_test::QK_51);

    qtest.test_int(fc.neval, exp_neval, "qag(f16, 51pt) sing neval");
    qtest.test_int(w.size(), exp_last, "qag(f16, 51pt) sing last");
    qtest.test_int(status, exp_ier, "qag(f16, 51pt) sing status");

    fc.neval = 0;
    std::tie(result, abserr)
      = __gnu_test::qag_integrate(w, fc, 1.0, -1.0, 1e-14, 0.0, 1000,
				  __gnu_test::QK_51);

    qtest.test_int(fc.neval, exp_neval, "qag(f16, 51pt) rev neval");
    qtest.test_int(w.size(), exp_last, "qag(f16, 51pt) rev last");
    qtest.test_int(status, exp_ier, "qag(f16, 51pt) rev status");
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

    int status = 0;
    quadrature_test<double> qtest;

    __gnu_test::integration_workspace<double> w(3);

    double exp_result =  9.565151449233894709;
    double exp_abserr =  1.570369823891028460E+01;
    int exp_neval  =     305;
    int exp_ier    =     MAX_ITER_ERROR;
    int exp_last   =     3;

    double a[3] = { -5.000000000000000000E-01,
                    0.000000000000000000,
                    -1.000000000000000000 };

    double b[3] = { 0.000000000000000000,
                    1.000000000000000000,
                    -5.000000000000000000E-01 };

    double r[3] = { 9.460353469435913709,
                    9.090909090909091161E-02,
                    1.388888888888888812E-02 };

    double e[3] = { 1.570369823891028460E+01,
                    1.009293658750142399E-15,
                    1.541976423090495140E-16 };

    int order[3] = { 1, 2, 3 };

    double alpha = 1.0;
    auto f = make_function<double>(f16, alpha);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qag_integrate(w, fc, -1.0, 1.0, 1e-14, 0.0, 3,
				  __gnu_test::QK_61);

    qtest.test_rel(result, exp_result, 1e-15, "qag(f16, 61pt) limit result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qag(f16, 61pt) limit abserr");
    qtest.test_int(fc.neval, exp_neval, "qag(f16, 61pt) limit neval");
    qtest.test_int(w.size(), exp_last, "qag(f16, 61pt) limit last");
    qtest.test_int(status, exp_ier, "qag(f16, 61pt) limit status");

    for (int i = 0; i < 3; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qag(f16, 61pt) limit alist");

    for (int i = 0; i < 3; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qag(f16, 61pt) limit blist");

    for (int i = 0; i < 3; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-15, "qag(f16, 61pt) limit rlist");

    for (int i = 0; i < 3; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-6, "qag(f16, 61pt) limit elist");

    for (int i = 0; i < 3; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qag(f16, 61pt) limit order");

    fc.neval = 0;
    std::tie(result, abserr)
      = __gnu_test::qag_integrate(w, fc, 1.0, -1.0, 1e-14, 0.0, 1000,
				  __gnu_test::QK_61);

    qtest.test_rel(result, -exp_result, 1e-15, "qag(f16, 61pt) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qag(f16, 61pt) reverse abserr");
    qtest.test_int(fc.neval, exp_neval, "qag(f16, 61pt) reverse neval");
    qtest.test_int(w.size(), exp_last, "qag(f16, 61pt) reverse last");
    qtest.test_int(status, exp_ier, "qag(f16, 61pt) reverse status");
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

    double exp_result = 7.716049382715789440E-02;
    double exp_abserr = 2.216394961010438404E-12;
    int exp_neval  =     189;
    int exp_ier    =       0;
    int exp_last   =       5;

    double a[5] = { 0, 0.5, 0.25, 0.125, 0.0625 };
    double b[5] = { 0.0625, 1, 0.5, 0.25, 0.125 };
    double r[5] = { 3.919381915366914693E-05,
                    5.491842501998223103E-02,
                    1.909827770934243579E-02,
                    2.776531175604360097E-03,
                    3.280661030752062609E-04 };
    double e[5] = { 2.215538742580964735E-12,
                    6.097169993333454062E-16,
                    2.120334764359736441E-16,
                    3.082568839745514608E-17,
                    3.642265412331439511E-18 };
    int order[5] = { 1, 2, 3, 4, 5 };

    double alpha = 2.6;
    auto f = make_function<double>(f1, alpha);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qags_integrate(w, fc, 0.0, 1.0, 0.0, 1e-10);

    qtest.test_rel(result, exp_result, 1e-15, "qags(f1) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qags(f1) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qags(f1) smooth neval");
    qtest.test_int(w.size(), exp_last, "qags(f1) smooth last");
    qtest.test_int(status, exp_ier, "qags(f1) smooth status");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qags(f1) smooth alist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qags(f1) smooth blist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-15, "qags(f1) smooth rlist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-6, "qags(f1) smooth elist");

    for (int i = 0; i < 5; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qags(f1) smooth order");

    fc.neval = 0;
    std::tie(result, abserr)
      = __gnu_test::qags_integrate(w, fc, 1.0, 0.0, 0.0, 1e-10);

    qtest.test_rel(result, -exp_result, 1e-15, "qags(f1) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qags(f1) reverse abserr");
    qtest.test_int(fc.neval, exp_neval, "qags(f1) reverse neval");
    qtest.test_int(w.size(), exp_last, "qags(f1) reverse last");
    qtest.test_int(status, exp_ier, "qags(f1) reverse status");
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

    double exp_result = -5.908755278982136588E+03;
    double exp_abserr = 1.299646281053874554E-10;
    int exp_neval  =     357;
    int exp_ier    =       0;
    int exp_last   =       9;

    double a[9] = { 1.000000000000000000E+00,
                    5.005000000000000000E+02,
                    2.507500000000000000E+02,
                    1.258750000000000000E+02,
                    6.343750000000000000E+01,
                    3.221875000000000000E+01,
                    1.660937500000000000E+01,
                    8.804687500000000000E+00,
                    4.902343750000000000E+00 };
    double b[9] = { 4.902343750000000000E+00,
                    1.000000000000000000E+03,
                    5.005000000000000000E+02,
                    2.507500000000000000E+02,
                    1.258750000000000000E+02,
                    6.343750000000000000E+01,
                    3.221875000000000000E+01,
                    1.660937500000000000E+01,
                    8.804687500000000000E+00 };
    double r[9] = { -3.890977835520834649E+00,
                    -3.297343675805121620E+03,
                    -1.475904154146372775E+03,
                    -6.517404019686431411E+02,
                    -2.829354222635842007E+02,
                    -1.201692001973227519E+02,
                    -4.959999906099650246E+01,
                    -1.971441499411640308E+01,
                    -7.457032710459004399E+00 };
    double e[9] = { 6.448276035006137169E-11,
                    3.660786868980994028E-11,
                    1.638582774073219226E-11,
                    7.235772003440423011E-12,
                    3.141214202790722909E-12,
                    1.334146129098576244E-12,
                    5.506706097890446534E-13,
                    2.188739744348345039E-13,
                    8.278969410534525339E-14 };
    int order[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    double alpha = 2.0;
    auto f = make_function<double>(f11, alpha);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qags_integrate(w, fc, 1.0, 1000.0, 1e-7, 0.0);

    qtest.test_rel(result, exp_result, 1e-15, "qags(f11) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-3, "qags(f11) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qags(f11) smooth neval");
    qtest.test_int(w.size(), exp_last, "qags(f11) smooth last");
    qtest.test_int(status, exp_ier, "qags(f11) smooth status");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qags(f11) smooth alist");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qags(f11) smooth blist");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-15, "qags(f11) smooth rlist");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-5, "qags(f11) smooth elist");

    for (int i = 0; i < 9; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qags(f11) smooth order");

    fc.neval = 0;
    std::tie(result, abserr)
      = __gnu_test::qags_integrate(w, fc, 1000.0, 1.0, 1e-7, 0.0);

    qtest.test_rel(result, -exp_result, 1e-15, "qags(f11) reverse result");
    qtest.test_rel(abserr, exp_abserr, 1e-3, "qags(f11) reverse abserr");
    qtest.test_int(fc.neval, exp_neval, "qags(f11) reverse neval");
    qtest.test_int(w.size(), exp_last, "qags(f11) reverse last");
    qtest.test_int(status, exp_ier, "qags(f11) reverse status");
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

    double exp_result = -3.616892186127022568E-01;
    double exp_abserr = 3.016716913328831851E-06;
    int exp_neval  =      285;
    int exp_ier    =        0;
    int exp_last   =       10;

    double a[10] = { 9.687500000000000000E-01,
                     0.000000000000000000E+00,
                     5.000000000000000000E-01,
                     2.500000000000000000E-01,
                     7.500000000000000000E-01,
                     1.250000000000000000E-01,
                     8.750000000000000000E-01,
                     6.250000000000000000E-02,
                     9.375000000000000000E-01,
                     3.125000000000000000E-02 };
    double b[10] = { 1.000000000000000000E+00,
                     3.125000000000000000E-02,
                     7.500000000000000000E-01,
                     5.000000000000000000E-01,
                     8.750000000000000000E-01,
                     2.500000000000000000E-01,
                     9.375000000000000000E-01,
                     1.250000000000000000E-01,
                     9.687500000000000000E-01,
                     6.250000000000000000E-02 };
    double r[10] = { -1.390003415539725340E-01,
                     1.429785306003466313E-03,
                     -1.229943369113085765E-02,
                     2.995321156568048898E-03,
                     -4.980050133751051655E-02,
                     2.785385934678596704E-03,
                     -8.653752279614615461E-02,
                     1.736218164975512294E-03,
                     -8.398745675010892142E-02,
                     1.041689192004495576E-03 };
    double e[10] = { 2.395037249893453013E-02,
                     2.161214992172538524E-04,
                     5.720644840858777846E-14,
                     3.325474514168701167E-17,
                     3.147380432198176412E-14,
                     3.092399597147240624E-17,
                     9.607595030230581153E-16,
                     1.927589382528252344E-17,
                     9.324480826368044019E-16,
                     1.156507325466566521E-17 };
    int order[10] = { 1, 2, 3, 5, 7, 9, 4, 6, 8, 10 };

    auto f = make_function<double>(f455);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qagiu_integrate(w, fc, 0.0, 0.0, 1.0e-3);

    qtest.test_rel(result, exp_result, 1e-14, "qagiu(f455) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qagiu(f455) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qagiu(f455) smooth neval");
    qtest.test_int(w.size(), exp_last, "qagiu(f455) smooth last");
    qtest.test_int(status, exp_ier, "qagiu(f455) smooth status");

    for (int i = 0; i < 10; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qagiu(f455) smooth alist");

    for (int i = 0; i < 10; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qagiu(f455) smooth blist");

    for (int i = 0; i < 10; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-15, "qagiu(f455) smooth rlist");

    for (int i = 0; i < 10; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-4, "qagiu(f455) smooth elist");

    for (int i = 0; i < 10; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qagiu(f455) smooth order");
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

    double exp_result = 6.553600000000024738E+04;
    double exp_abserr = 7.121667111456009280E-04;
    int exp_neval  =      285;
    int exp_ier    =        0;
    int exp_last   =       10;

    double a[10] = { 0.000000000000000000E+00,
                     5.000000000000000000E-01,
                     2.500000000000000000E-01,
                     1.250000000000000000E-01,
                     6.250000000000000000E-02,
                     3.125000000000000000E-02,
                     1.562500000000000000E-02,
                     7.812500000000000000E-03,
                     3.906250000000000000E-03,
                     1.953125000000000000E-03 };
    double b[10] = { 1.953125000000000000E-03,
                     1.000000000000000000E+00,
                     5.000000000000000000E-01,
                     2.500000000000000000E-01,
                     1.250000000000000000E-01,
                     6.250000000000000000E-02,
                     3.125000000000000000E-02,
                     1.562500000000000000E-02,
                     7.812500000000000000E-03,
                     3.906250000000000000E-03 };
    double r[10] = { 1.099297665754340292E+00,
                     3.256176475185617591E-01,
                     8.064694554185326325E+00,
                     8.873128656118993263E+01,
                     6.977679035845269482E+02,
                     4.096981198511257389E+03,
                     1.574317583220441520E+04,
                     2.899418134793237914E+04,
                     1.498314766425578091E+04,
                     9.225251570832365360E+02 };
    double e[10] = { 7.101865971621337814E-04,
                     1.912660677170175771E-08,
                     9.167763417119923333E-08,
                     3.769501719163865578E-07,
                     6.973493131275552509E-07,
                     1.205653952340679711E-07,
                     1.380003928453846583E-07,
                     1.934652413547325474E-07,
                     3.408933028357320364E-07,
                     2.132473175465897029E-09 };
    int order[10] = { 1, 5, 4, 9, 8, 7, 6, 3, 2, 10 };

    double alpha = 5.0;

    auto f = make_function<double>(f15, alpha);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qagiu_integrate(w, fc, 0.0, 0.0, 1.0e-7);

    qtest.test_rel(result, exp_result, 1e-14, "qagiu(f15) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qagiu(f15) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qagiu(f15) smooth neval");
    qtest.test_int(w.size(), exp_last, "qagiu(f15) smooth last");
    qtest.test_int(status, exp_ier, "qagiu(f15) smooth status");

    for (int i = 0; i < 10; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qagiu(f15) smooth alist");

    for (int i = 0; i < 10; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qagiu(f15) smooth blist");

    for (int i = 0; i < 10; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-15, "qagiu(f15) smooth rlist");

    for (int i = 0; i < 10; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-4, "qagiu(f15) smooth elist");

    for (int i = 0; i < 10; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qagiu(f15) smooth order");
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

    double exp_result = 1.000000000006713292E-04;
    double exp_abserr = 3.084062020905636316E-09;
    int exp_neval  =      165;
    int exp_ier    =        0;
    int exp_last   =        6;

    double a[6] = { 0.000000000000000000E+00,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    1.250000000000000000E-01,
                    6.250000000000000000E-02,
                    3.125000000000000000E-02 };
    double b[6] = { 3.125000000000000000E-02,
                    1.000000000000000000E+00,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    1.250000000000000000E-01,
                    6.250000000000000000E-02 };
    double r[6] = { 7.633587786326674618E-05,
                    9.900990099009899620E-07,
                    1.922522349322310737E-06,
                    3.629434715543053753E-06,
                    6.501422186103209199E-06,
                    1.062064387653501389E-05 };
    double e[6] = { 3.084061858351569051E-09,
                    3.112064814755089674E-17,
                    4.543453652226561245E-17,
                    4.908618166361344548E-17,
                    3.014338672269481784E-17,
                    6.795996738013555461E-18 };
    int order[6] = { 1, 4, 3, 2, 5, 6 };

    double alpha = 1.0;

    auto f = make_function<double>(f16, alpha);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qagiu_integrate(w, fc, 99.9, 1.0e-7, 0.0);

    qtest.test_rel(result, exp_result, 1e-14, "qagiu(f16) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qagiu(f16) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qagiu(f16) smooth neval");
    qtest.test_int(w.size(), exp_last, "qagiu(f16) smooth last");
    qtest.test_int(status, exp_ier, "qagiu(f16) smooth status");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qagiu(f16) smooth alist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qagiu(f16) smooth blist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-15, "qagiu(f16) smooth rlist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-4, "qagiu(f16) smooth elist");

    for (int i = 0; i < 6; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qagiu(f16) smooth order");
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

    double exp_result = 2.275875794468747770E+00;
    double exp_abserr = 7.436490118267390744E-09;
    int exp_neval  =      270;
    int exp_ier    =        0;
    int exp_last   =        5;

    double a[5] = { 1.250000000000000000E-01,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    0.000000000000000000E+00,
                    3.750000000000000000E-01 };
    double b[5] = { 2.500000000000000000E-01,
                    1.000000000000000000E+00,
                    3.750000000000000000E-01,
                    1.250000000000000000E-01,
                    5.000000000000000000E-01 };
    double r[5] = { 4.639317228058405717E-04,
                    1.691664195356748834E+00,
                    1.146307471900291086E-01,
                    4.379392477350953574E-20,
                    4.691169201991640669E-01 };
    double e[5] = { 3.169263960393051137E-09,
                    4.265988974874425043E-09,
                    1.231954072964969637E-12,
                    8.360902986775307673E-20,
                    5.208244060463541433E-15 };
    int order[5] = { 2, 1, 3, 5, 4 };

    auto f = make_function<double>(myfn1);
    counted_function<double> fc(f);

    auto [result, abserr] = __gnu_test::qagi_integrate(w, fc, 1.0e-7, 0.0);

    qtest.test_rel(result, exp_result, 1e-14, "qagi(myfn1) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qagi(myfn1) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qagi(myfn1) smooth neval");
    qtest.test_int(w.size(), exp_last, "qagi(myfn1) smooth last");
    qtest.test_int(status, exp_ier, "qagi(myfn1) smooth status");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qagi(myfn1) smooth alist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qagi(myfn1) smooth blist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-14, "qagi(myfn1) smooth rlist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-4, "qagi(myfn1) smooth elist");

    for (int i = 0; i < 5; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qagi(myfn1) smooth order");
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

    double exp_result = 2.718281828459044647E+00;
    double exp_abserr = 1.588185109253204805E-10;
    int exp_neval  =      135;
    int exp_ier    =        0;
    int exp_last   =        5;

    double a[5] = { 0.000000000000000000E+00,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    1.250000000000000000E-01,
                    6.250000000000000000E-02 };
    double b[5] = { 6.250000000000000000E-02,
                    1.000000000000000000E+00,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    1.250000000000000000E-01 };
    double r[5] = { 8.315287189746029816E-07,
                    1.718281828459045091E+00,
                    8.646647167633871867E-01,
                    1.328565310599463256E-01,
                    2.477920647947255521E-03 };
    double e[5] = { 1.533437090413525935E-10,
                    4.117868247943567505E-12,
                    7.802455785301941044E-13,
                    5.395586026138397182E-13,
                    3.713312434866150125E-14 };
    int order[5] = { 1, 2, 3, 4, 5 };

    double alpha = 1.0;
    auto f = make_function<double>(myfn2, alpha);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qagil_integrate(w, fc, 1.0, 1.0e-7, 0.0);

    qtest.test_rel(result, exp_result, 1e-14, "qagil(myfn2) smooth result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qagil(myfn2) smooth abserr");
    qtest.test_int(fc.neval, exp_neval, "qagil(myfn2) smooth neval");
    qtest.test_int(w.size(), exp_last, "qagil(myfn2) smooth last");
    qtest.test_int(status, exp_ier, "qagil(myfn2) smooth status");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qagil(myfn2) smooth alist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qagil(myfn2) smooth blist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-14, "qagil(myfn2) smooth rlist");

    for (int i = 0; i < 5; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-4, "qagil(myfn2) smooth elist");

    for (int i = 0; i < 5; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qagil(myfn2) smooth order");
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

    double exp_result = 5.274080611672716401E+01;
    double exp_abserr = 1.755703848687062418E-04;
    int exp_neval  =        777;
    int exp_ier    =          0;
    int exp_last   =         20;

    double a[20] = { 9.687500000000000000E-01,
                    1.401269388548935790E+00,
                    1.414213562373095145E+00,
                    1.000000000000000000E+00,
                    0.000000000000000000E+00,
                    2.207106781186547462E+00,
                    1.810660171779821415E+00,
                    1.207106781186547462E+00,
                    5.000000000000000000E-01,
                    1.103553390593273731E+00,
                    1.612436867076458391E+00,
                    1.310660171779821415E+00,
                    7.500000000000000000E-01,
                    1.051776695296636976E+00,
                    1.513325214724776657E+00,
                    1.362436867076458391E+00,
                    8.750000000000000000E-01,
                    1.463769388548935790E+00,
                    1.388325214724776657E+00,
                    9.375000000000000000E-01};
    double b[20] = { 1.000000000000000000E+00,
                     1.414213562373095145E+00,
                     1.463769388548935790E+00,
                     1.051776695296636976E+00,
                     5.000000000000000000E-01,
                     3.000000000000000000E+00,
                     2.207106781186547462E+00,
                     1.310660171779821415E+00,
                     7.500000000000000000E-01,
                     1.207106781186547462E+00,
                     1.810660171779821415E+00,
                     1.362436867076458391E+00,
                     8.750000000000000000E-01,
                     1.103553390593273731E+00,
                     1.612436867076458391E+00,
                     1.388325214724776657E+00,
                     9.375000000000000000E-01,
                     1.513325214724776657E+00,
                     1.401269388548935790E+00,
                     9.687500000000000000E-01};
    double r[20] = { -1.125078814079027711E-01,
                     -1.565132123531515207E-01,
                     -4.225328513207429193E-01,
                     -1.830392049835374568E-01,
                     6.575875041899758092E-03,
                     4.873920540843067783E+01,
                     6.032891565603589079E+00,
                     -2.991531901645863023E-01,
                     -7.326282608704996063E-03,
                     -2.431894410706912923E-01,
                     5.911661670635662835E-01,
                     -2.236786562536174916E-01,
                     -5.647871991778510847E-02,
                     -1.305470403178642658E-01,
                     -1.721363984401322045E-01,
                     -1.589345454585119055E-01,
                     -7.406626263352669715E-02,
                     -2.208730668000830344E-01,
                     -1.048692749517999567E-01,
                     -6.302287584527696551E-02};
    double e[20] = { 2.506431410088378817E-02,
                     2.730454695485963826E-02,
                     1.017446081816190118E-01,
                     3.252808038935910834E-02,
                     7.300687878575027348E-17,
                     5.411138804637469780E-13,
                     6.697855121200013106E-14,
                     3.321267596107916554E-15,
                     1.417509685426979386E-16,
                     2.699945168224041491E-15,
                     6.573952690524728748E-15,
                     2.483331942899818875E-15,
                     6.270397525408045936E-16,
                     1.449363299575615261E-15,
                     1.911097929242846383E-15,
                     1.764527917763735212E-15,
                     8.223007012367522077E-16,
                     2.452183642810224359E-15,
                     1.164282836272345215E-15,
                     6.996944784151910810E-16};
    int order[20] = { 3, 4, 2, 1, 6, 7, 11, 8, 10, 12, 18,
                     15, 16, 14, 19, 17, 20, 13, 9, 5 };

    auto f = make_function<double>(f454);
    counted_function<double> fc(f);

    std::vector<double> pts{0.0, 1.0, std::sqrt(2.0), 3.0};

    auto [result, abserr] = __gnu_test::qagp_integrate(w, fc, pts, 0.0, 1.0e-3);

    qtest.test_rel(result, exp_result, 1e-14, "qagp(f454) singular result");
    qtest.test_rel(abserr, exp_abserr, 1e-5, "qagp(f454) singular abserr");
    qtest.test_int(fc.neval, exp_neval, "qagp(f454) singular neval");
    qtest.test_int(w.size(), exp_last, "qagp(f454) singular last");
    qtest.test_int(status, exp_ier, "qagp(f454) singular status");

    for (int i = 0; i < 20; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qagp(f454) singular alist");

    for (int i = 0; i < 20; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qagp(f454) singular blist");

    for (int i = 0; i < 20; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-14, "qagp(f454) singular rlist");

    for (int i = 0; i < 20; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-4, "qagp(f454) singular elist");

    for (int i = 0; i < 20; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qagp(f454) singular order");
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

    double exp_result = -8.994400695837000137E-02;
    double exp_abserr =  1.185290176227023727E-06;
    int exp_neval  =      215;
    int exp_ier    =        0;
    int exp_last   =        6;

    double a[6] = { -1.000000000000000000E+00,
                    2.500000000000000000E+00,
                    1.250000000000000000E+00,
                    6.250000000000000000E-01,
                    -5.000000000000000000E-01,
                    -7.500000000000000000E-01};
    double b[6] = { -7.500000000000000000E-01,
                    5.000000000000000000E+00,
                    2.500000000000000000E+00,
                    1.250000000000000000E+00,
                    6.250000000000000000E-01,
                    -5.000000000000000000E-01};
    double r[6] = { -1.234231128040012976E-01,
                    3.579970394639702888E-03,
                    2.249831615049339983E-02,
                    7.214232992127905808E-02,
                    2.079093855884046535E-02,
                    -8.553244917962132821E-02};
    double e[6] = { 1.172832717970022565E-06,
                    9.018232896137375412E-13,
                    1.815172652101790755E-12,
                    1.006998195150956048E-13,
                    1.245463873006391609E-08,
                    1.833082948207153514E-15 };
    int order[6] = { 1, 5, 3, 2, 4, 6 };

    auto f = make_function<double>(f459);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qawc_integrate(w, fc, -1.0, 5.0, 0.0, 0.0, 1.0e-3);

    qtest.test_rel(result, exp_result, 1e-14, "qawc(f459) result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qawc(f459) abserr");
    qtest.test_int(fc.neval, exp_neval, "qawc(f459) neval");
    qtest.test_int(w.size(), exp_last, "qawc(f459) last");
    qtest.test_int(status, exp_ier, "qawc(f459) status");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qawc(f459) alist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qawc(f459) blist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-14, "qawc(f459) rlist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-4, "qawc(f459) elist");

    for (int i = 0; i < 6; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qawc(f459) order");

    fc.neval = 0;
    std::tie(result, abserr)
      = __gnu_test::qawc_integrate(w, fc, 5.0, -1.0, 0.0, 0.0, 1.0e-3);

    qtest.test_rel(result, -exp_result, 1e-14, "qawc(f459) rev result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qawc(f459) rev abserr");
    qtest.test_int(fc.neval, exp_neval, "qawc(f459) rev neval");
    qtest.test_int(w.size(), exp_last, "qawc(f459) rev last");
    qtest.test_int(status, exp_ier, "qawc(f459) rev status");
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

    __gnu_test::qaws_integration_table<double> t(0.0, 0.0, 1, 0);

    __gnu_test::integration_workspace<double> w(1000);

    double exp_result = -1.892751853489401670E-01;
    double exp_abserr = 1.129133712015747658E-08;
    int exp_neval  =      280;
    int exp_ier    =        0;
    int exp_last   =        8;

    double a[8] = { 0.000000000000000000E+00,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    1.250000000000000000E-01,
                    6.250000000000000000E-02,
                    3.125000000000000000E-02,
                    1.562500000000000000E-02,
                    7.812500000000000000E-03};
    double b[8] = { 7.812500000000000000E-03,
                    1.000000000000000000E+00,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    1.250000000000000000E-01,
                    6.250000000000000000E-02,
                    3.125000000000000000E-02,
                    1.562500000000000000E-02};
    double r[8] = { -4.126317299834445824E-05,
                    -1.076283950172247789E-01,
                    -6.240573216173390947E-02,
                    -1.456169844189576269E-02,
                    -3.408925115926728436E-03,
                    -8.914083918175634211E-04,
                    -2.574191402137795482E-04,
                    -8.034390712936630608E-05};
    double e[8] = { 1.129099387465713953E-08,
                    3.423394967694403596E-13,
                    6.928428071454762659E-16,
                    1.616673288784094320E-16,
                    3.784667152924835070E-17,
                    9.896621209399419425E-18,
                    2.857926564445496100E-18,
                    8.919965558336773736E-19};
    int order[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };

    auto f = make_function<double>(f458);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qaws_integrate(w, t, fc, 0.0, 1.0, 0.0, 1.0e-7);

    qtest.test_rel(result, exp_result, 1e-14, "qaws(f458) ln(x-a) result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qaws(f458) ln(x-a) abserr");
    qtest.test_int(fc.neval, exp_neval, "qaws(f458) ln(x-a) neval");
    qtest.test_int(w.size(), exp_last, "qaws(f458) ln(x-a) last");
    qtest.test_int(status, exp_ier, "qaws(f458) ln(x-a) status");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qaws(f458) ln(x-a) alist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qaws(f458) ln(x-a) blist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-14, "qaws(f458) ln(x-a) rlist");

    for (int i = 0; i < 6; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-4, "qaws(f458) ln(x-a) elist");

    for (int i = 0; i < 6; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qaws(f458) ln(x-a) order");

    // Test without logs

    t.set(-0.5, -0.3, 0, 0);

    std::tie(result, abserr)
      = __gnu_test::qaws_integrate(w, t, fc, 0.0, 1.0, 0.0, 1.0e-7);

    exp_result = 9.896686656601706433E-01;
    exp_abserr = 5.888032513201251628E-08;

    qtest.test_rel(result, exp_result, 1e-14, "qaws(f458) AB result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qaws(f458) AB abserr");

    // Test with ln(x - a)

    t.set(-0.5, -0.3, 1, 0);

    std::tie(result, abserr)
      = __gnu_test::qaws_integrate(w, t, fc, 0.0, 1.0, 0.0, 1.0e-7);

    exp_result = -3.636679470586539620E-01;
    exp_abserr = 2.851348775257054093E-08;

    qtest.test_rel(result, exp_result, 1e-14, "qaws(f458) AB ln(x-a) result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qaws(f458) AB ln(x-a) abserr");

    // Test with ln(b - x)

    t.set(-0.5, -0.3, 0, 1);

    std::tie(result, abserr)
      = __gnu_test::qaws_integrate(w, t, fc, 0.0, 1.0, 0.0, 1.0e-7);

    exp_result = -1.911489253363409802E+00;
    exp_abserr = 9.854016753016499034E-09;

    qtest.test_rel(result,exp_result,1e-14,"qaws(f458) AB ln(b-x) result");
    qtest.test_rel(abserr,exp_abserr,1e-6,"qaws(f458) AB ln(b-x) abserr");

    // Test with ln(x - a) ln(b - x)

    t.set (-0.5, -0.3, 1, 1);

    std::tie(result, abserr)
      = __gnu_test::qaws_integrate(w, t, fc, 0.0, 1.0, 0.0, 1.0e-7);

    exp_result = 3.159922862811048172E-01;
    exp_abserr = 2.336183482198144595E-08;

    qtest.test_rel(result, exp_result, 1e-14, "qaws(f458) AB ln(x-a)ln(b-x) result");
    qtest.test_rel(abserr, exp_abserr, 1e-6, "qaws(f458) AB ln(x-a)ln(b-x) abserr");
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
    __gnu_test::oscillatory_integration_table<double> wo(10.0 * M_PI, 1.0,
                                     __gnu_test::oscillatory_integration_table<double>::INTEG_SINE, 1000);

    double exp_result = -1.281368483991674190E-01;
    double exp_abserr =  6.875028324415666248E-12;
    int exp_neval  =      305;
    int exp_ier    =        0;
    int exp_last   =        9;

    double a[9] = { 0.000000000000000000E+00,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    1.250000000000000000E-01,
                    6.250000000000000000E-02,
                    3.125000000000000000E-02,
                    1.562500000000000000E-02,
                    7.812500000000000000E-03,
                    3.906250000000000000E-03 };
    double b[9] = { 3.906250000000000000E-03,
                    1.000000000000000000E+00,
                    5.000000000000000000E-01,
                    2.500000000000000000E-01,
                    1.250000000000000000E-01,
                    6.250000000000000000E-02,
                    3.125000000000000000E-02,
                    1.562500000000000000E-02,
                    7.812500000000000000E-03 };
    double r[9] = { -1.447193692377651136E-03,
                    2.190541162282139478E-02,
                    -2.587726479625663753E-02,
                    5.483209176363500886E-02,
                    -3.081695575172510582E-02,
                    -9.178321994387816929E-02,
                    -3.886716016498160953E-02,
                    -1.242306301902117854E-02,
                    -3.659495117871544145E-03};
    double e[9] = { 8.326506625798146465E-07,
                    1.302638552580516100E-13,
                    7.259224351945759794E-15,
                    1.249770395036711102E-14,
                    7.832180081562836579E-16,
                    1.018998440559284116E-15,
                    4.315121611695628020E-16,
                    1.379237060008662177E-16,
                    4.062855738364339357E-17 };
    int order[9] = { 1, 2, 4, 3, 6, 5, 7, 8, 9 };

    auto f = make_function<double>(f456<double>);
    counted_function<double> fc(f);

    auto [result, abserr] = __gnu_test::qawo_integrate(w, wo, fc, 0.0, 0.0, 1e-7);

    qtest.test_rel(result, exp_result, 1e-14, "qawo(f456) result");
    qtest.test_rel(abserr, exp_abserr, 1e-3, "qawo(f456) abserr");
    qtest.test_int(fc.neval, exp_neval, "qawo(f456) neval");
    qtest.test_int(w.size(), exp_last, "qawo(f456) last");
    qtest.test_int(status, exp_ier, "qawo(f456) status");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.lower_lim(i), a[i], 1e-15, "qawo(f456) alist");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.upper_lim(i), b[i], 1e-15, "qawo(f456) blist");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-14, "qawo(f456) rlist");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 1e-2, "qawo(f456) elist");

    for (int i = 0; i < 9; ++i)
      qtest.test_int(w.order(i), order[i]-1, "qawo(f456) order");


    // In reverse, flip limit and sign of length

    wo.set_length(-1.0);
    fc.neval = 0;
    std::tie(result, abserr) = qawo_integrate(w, wo, fc, 1.0, 0.0, 1e-7);

    qtest.test_rel(result, -exp_result, 1e-14, "qawo(f456) rev result");
    qtest.test_rel(abserr, exp_abserr, 1e-3, "qawo(f456) rev abserr");
    qtest.test_int(fc.neval, exp_neval, "qawo(f456) rev neval");
    qtest.test_int(w.size(), exp_last, "qawo(f456) rev last");
    qtest.test_int(status, exp_ier, "qawo(f456) rev status");
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
      wo(M_PI / 2.0, 1.0,
	__gnu_test::oscillatory_integration_table<double>::INTEG_COSINE, 1000);

    double exp_result = 9.999999999279802765E-01;
    double exp_abserr = 1.556289974669056164E-08;
    int exp_neval  =      590;
    int exp_ier    =        0;
    int exp_last   =       12;

    double r[12] = { 1.013283128125232802E+00,
                    -1.810857954748607349E-02,
                    7.466754034900931897E-03,
                    -4.360312526786496237E-03,
                    2.950184068216192904E-03,
                    -2.168238443073697373E-03,
                    1.680910783140869081E-03,
                    -1.352797860944863345E-03,
                    1.119354921991485901E-03,
                    -9.462367583691360827E-04,
                    8.136341270731781887E-04,
                    -7.093931338504278145E-04 };
    double e[12] = { 1.224798040766472695E-12,
                    1.396565155187268456E-13,
                    1.053844511655910310E-16,
                    6.505213034913026604E-19,
                    7.155734338404329264E-18,
                    1.105886215935214523E-17,
                    9.757819552369539906E-18,
                    5.854691731421723944E-18,
                    4.553649124439220312E-18,
                    7.643625316022806260E-18,
                    2.439454888092388058E-17,
                    2.130457268934021451E-17 };

    auto f = make_function<double>(f457);
    counted_function<double> fc(f);

    auto [result, abserr]
      = __gnu_test::qawf_integrate(w, wc, wo, fc, 0.0, 1e-7);

    qtest.test_rel(result, exp_result, 1e-14, "qawf(f457) result");
    qtest.test_rel(abserr, exp_abserr, 1e-3, "qawf(f457) abserr");
    qtest.test_int(fc.neval, exp_neval, "qawf(f457) neval");
    qtest.test_int(w.size(), exp_last, "qawf(f457) last");
    qtest.test_int(status, exp_ier, "qawf(f457) status");

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.result(i), r[i], 1e-12, "qawf(f457) rlist");

    // We can only get within two orders of magnitude on the error
    // here,  which is very sensitive to the floating point precision

    for (int i = 0; i < 9; ++i)
      qtest.test_rel(w.abs_error(i), e[i], 50.0, "qawf(f457) elist");
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

    qtest.test_abs(dmon_t(2, 1.0)(2.0), 4.0, 8*_S_eps, "monomial sanity check 1");

    qtest.test_abs(dmon_t(1, 2.0)(2.0), 4.0, 8*_S_eps, "monomial sanity check 2");

    qtest.test_abs(integrate(dmon_t(2, 2.0), 1.0, 2.0),
        (2.0/3.0)*(2.0*2.0*2.0 - 1.0*1.0*1.0), 8*_S_eps,
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

    const double a = 0.0, b = 1.2;

    for (int n = 1; n < 1025; ++n)
      {
	quadrature_test<double> qtest;
        __gnu_test::gauss_legendre_table<double> tbl(n);

        monomial<double> mon(2*n-1, 1.0); // n point rule exact for 2n-1 degree poly
        double expected      = integrate(mon, a, b);
        double result        = __gnu_test::glfixed_integrate(tbl, mon, a, b);

        if (tbl.precomputed)
          {
	    std::ostringstream str;
	    str << "glfixed " << n << "-point: Integrating ("
		<< mon.constant << "*x^" << mon.degree
		<< ") over [" << a << "," << a << "]";
            qtest.test_rel(result, expected, 1.0e-12, str.str().c_str());
          }
        else
          {
	    std::ostringstream str;
	    str << "glfixed " << n << "-point: Integrating ("
		<< mon.constant << "*x^" << mon.degree
		<< ") over [" << a << "," << a << "]";
            qtest.test_rel(result, expected, 1.0e-7, str.str().c_str());
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

    qtest.test_abs(f_sin(2.0), std::sin(2.0), 0.0, "f_sin sanity check 1");
    qtest.test_abs(f_sin(7.0), std::sin(7.0), 0.0, "f_sin sanity check 2");
    qtest.test_abs(integ_f_sin(0.0, M_PI), 2.0, _S_eps,
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
    const std::function<double(double)> f = f_sin<double>;
    const double a = 0.0, b = M_PI;
    const double expected = integ_f_sin(a, b);
    double result, abserr, prev_abserr = 0.0;
    int n;
    quadrature_test<double> qtest;

    for (n = 1; n <= n_max; ++n)
      {
        __gnu_test::gauss_legendre_table<double> tbl(n);

        result = __gnu_test::glfixed_integrate(tbl, f, a, b);
        abserr = std::abs(expected - result);

        if (n == 1)
	  {
	    std::ostringstream str;
	    str << "glfixed " << n << "-point: behavior for n == 1";
            qtest.test_abs(result, (b - a) * f((b + a) / 2.0), 0.0, str.str().c_str());
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
	    qtest.test_abs(result, expected, 2.0 * n * _S_eps, str.str().c_str());
	  }
        else
	  {
	    std::ostringstream str;
	    str << "glfixed " << n << "-point: acceptable absolute error for on-the-fly coefficients";
	    qtest.test_abs(result, expected, 1.0e6 * _S_eps, str.str().c_str());
	  }

        prev_abserr = abserr;
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

  // Test some fixed-order Gauss-Legendre rule points and weights on [-1, 1].
  // This verifies the (point, weight) retrieval API behaves sanely.
  try
  {
    std::cout << ">>>> Test fixed-order Gauss-Legendre rule points and weights on [-1, 1]..." << std::endl;

    const double eps = _S_eps;
    std::size_t n;
    quadrature_test<double> qtest;

    // Analytical results for points and weights on [-1, 1]
    // Pulled from http://en.wikipedia.org/wiki/Gaussian_quadrature
    // Sorted in increasing order of Gauss points

    const double
    e1[1][2]
    {
      {0.0, 2.0}
    };

    const double
    e2[2][2]
    {
      {-1.0/std::sqrt(3.0), 1.0},
      { 1.0/std::sqrt(3.0), 1.0}
    };

    const double
    e3[3][2]
    {
      {-std::sqrt(15.0)/5.0, 5.0/9.0},
      {                 0.0, 8.0/9.0},
      { std::sqrt(15.0)/5.0, 5.0/9.0}
    };

    const double
    e4[4][2]
    {
      {-sqrt((3.0+2.0*std::sqrt(6./5.0))/7.0), (18.0-std::sqrt(30.0))/36.0},
      {-sqrt((3.0-2.0*std::sqrt(6./5.0))/7.0), (18.0+std::sqrt(30.0))/36.0},
      { sqrt((3.0-2.0*std::sqrt(6./5.0))/7.0), (18.0+std::sqrt(30.0))/36.0},
      { sqrt((3.0+2.0*std::sqrt(6./5.0))/7.0), (18.0-std::sqrt(30.0))/36.0}
    };

    const double
    e5[5][2]
    {
      {-std::sqrt((5.0+2.0*std::sqrt(10.0/7.0)))/3, (322.0-13.0*std::sqrt(70.0))/900.0},
      {-std::sqrt((5.0-2.0*std::sqrt(10.0/7.0)))/3, (322.0+13.0*std::sqrt(70.0))/900.0},
      {                                        0.0,                        128.0/225.0},
      { std::sqrt((5.0-2.0*std::sqrt(10.0/7.0)))/3, (322.0+13.0*std::sqrt(70.0))/900.0},
      { std::sqrt((5.0+2.0*std::sqrt(10.0/7.0)))/3, (322.0-13.0*std::sqrt(70.0))/900.0}
    };

    n = 1;
    __gnu_test::gauss_legendre_table<double> tbl1(n);
    for (auto i = 0u; i < n; ++i)
      {
        auto [xi, wi] = tbl1.get_point(-1.0, 1.0, i);
	std::ostringstream msg1, msg2;
	msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
        qtest.test_abs(xi, e1[i][0], eps, msg1.str().c_str());
	msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
        qtest.test_abs(wi, e1[i][1], eps, msg2.str().c_str());
      }

    n = 2;
    __gnu_test::gauss_legendre_table<double> tbl2(n);
    for (auto i = 0u; i < n; ++i)
      {
        auto [xi, wi] = tbl2.get_point(-1.0, 1.0, i);
	std::ostringstream msg1, msg2;
	msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
        qtest.test_abs(xi, e2[i][0], eps, msg1.str().c_str());
	msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
        qtest.test_abs(wi, e2[i][1], eps, msg2.str().c_str());
      }

    n = 3;
    __gnu_test::gauss_legendre_table<double> tbl3(n);
    for (auto i = 0u; i < n; ++i)
      {
        auto [xi, wi] = tbl3.get_point(-1.0, 1.0, i);
	std::ostringstream msg1, msg2;
	msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
        qtest.test_abs(xi, e3[i][0], eps, msg1.str().c_str());
	msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
        qtest.test_abs(wi, e3[i][1], eps, msg2.str().c_str());
      }

    n = 4;
    __gnu_test::gauss_legendre_table<double> tbl4(n);
    for (auto i = 0u; i < n; ++i)
      {
        auto [xi, wi] = tbl4.get_point(-1.0, 1.0, i);
	std::ostringstream msg1, msg2;
	msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
        qtest.test_abs(xi, e4[i][0], eps, msg1.str().c_str());
	msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
        qtest.test_abs(wi, e4[i][1], eps, msg2.str().c_str());
      }

    n = 5;
    __gnu_test::gauss_legendre_table<double> tbl5(n);
    for (auto i = 0u; i < n; ++i)
      {
        auto [xi, wi] = tbl5.get_point(-1.0, 1.0, i);
	std::ostringstream msg1, msg2;
	msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
        qtest.test_abs(xi, e5[i][0], eps, msg1.str().c_str());
	msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
        qtest.test_abs(wi, e5[i][1], eps, msg2.str().c_str());
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

  // Test some fixed-order Gauss-Legendre rule points and weights on [-2, 3].
  // This verifies the (point, weight) retrieval API is okay on non-[-1,1].
  try
  {
    std::cout << ">>>> Test some fixed-order Gauss-Legendre rule points and weights on [-2, 3]..." << std::endl;

    quadrature_test<double> qtest;
    std::size_t n = 0;
    double result;

    // Odd n = 3, f(x) = x**5 + x**4 + x**3 + x**2 + x**1 + 1
    n = 3;
    result = 0;
    __gnu_test::gauss_legendre_table<double> tbl1(n);
    for (auto i = 0u; i < n; ++i)
      {
        auto [x, w] = tbl1.get_point(-2.0, 3.0, i);
        result += w * (1 + x*(1 + x*(1 + x*(1 + x*(1 + x)))));
      }
    std::ostringstream msg1;
    msg1 << "glfixed " << n << "-point xi,wi eval";
    qtest.test_rel(result, 805./4, 1e-8, msg1.str().c_str());

    // Even n = 4, f(x) = x**7 + x**6 + x**5 + x**4 + x**3 + x**2 + x**1 + 1
    n = 4;
    result = 0;
    __gnu_test::gauss_legendre_table<double> tbl2(n);
    for (auto i = 0u; i < n; ++i)
      {
        auto [x, w] = tbl2.get_point(-2.0, 3.0, i);
        result += w * (1 + x*(1 + x*(1 + x*(1 + x*(1 + x*(1 + x*(1 + x)))))));
      }
    std::ostringstream msg2;
    msg2 << "glfixed " << n << "-point xi,wi eval";
    qtest.test_rel(result, 73925./56, 1e-8, msg2.str().c_str());
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

  // Test this newfangled cquad.
  try
  {
    std::cout << ">>>> Test this newfangled cquad..." << std::endl;
    typedef double (*fptr) (double);
    quadrature_test<double> qtest;

    const static fptr
    funs[25]
    {
      &cqf1, &cqf2, &cqf3, &cqf4, &cqf5, &cqf6, &cqf7,
      &cqf8, &cqf9, &cqf10, &cqf11, &cqf12, &cqf13, &cqf14, &cqf15, &cqf16, &cqf17,
      &cqf18, &cqf19, &cqf20, &cqf21, &cqf22, &cqf23, &cqf24, &cqf25
    };

    const static double
    ranges[50]
    {
      0, 1, 0, 1, 0, 1, -1, 1, -1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
      0, 1, 0, 1, 0, 1, 0, 1, 0, 10, 0, 10, 0, 10, 0, 1, 0, M_PI,
      0, 1, -1, 1, 0, 1, 0, 1, 0, 1, 0, 3, 0, 5
    };
    const static double
    f_exact[25]
    {
      1.7182818284590452354, 0.7, 2.0/3, 0.4794282266888016674,
      1.5822329637296729331, 0.4, 2, 0.86697298733991103757,
      1.1547005383792515290, 0.69314718055994530942, 0.3798854930417224753,
      0.77750463411224827640, 0.49898680869304550249,
      0.5, 1, 0.13263071079267703209e+08, 0.49898680869304550249,
      0.83867634269442961454, -1, 1.5643964440690497731,
      0.16349494301863722618, -0.63466518254339257343,
      0.013492485649467772692, 17.664383539246514971, 7.5
    };

    // Loop over the functions...
    for (int fid = 0; fid < 25; ++fid)
    {
      __gnu_test::cquad_workspace<double> ws(200);
      auto f = make_function<double>(funs[fid]);
      auto exact = f_exact[fid];
      int status = 0;

      // Call our quadrature routine.
      auto [result, abserr]
	= __gnu_test::cquad_integrate(ws, f, ranges[2* fid], ranges[2 * fid + 1], 0.0, 1.0e-12);

      std::ostringstream rstr;
      rstr << "cquad f" << fid;
      qtest.test_rel(result, exact, 1e-12, rstr.str().c_str());

      std::ostringstream upstr;
      upstr << "cquad f" << fid << " error("
			 << std::abs(result-exact) << " actual vs "
			 << abserr << " estimated)";
      qtest.test_update(std::abs(result - exact) > 5.0 * abserr, upstr.str().c_str());

      qtest.test_int(status, 0, "cquad return code");
      std::cout << std::flush;
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

  exit(quadrature_test<double>::test_summary());
}
