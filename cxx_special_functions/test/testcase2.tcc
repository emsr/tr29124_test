#ifndef TESTCASE2_TCC
#define TESTCASE2_TCC 1

#include <complex>
#include <sstream>
#include <type_traits>
#include <string_view>
#include <algorithm> // For clamp.
#include "complex_compare.h" // For the Statistics min/max (maybe rethink that there)
#include "statistics.h"
#include "spaceship.h"

const std::string_view boilerplate = 
R"(// { dg-do run { target c++11 } }
// { dg-options "-D__STDCPP_WANT_MATH_SPEC_FUNCS__" }
//
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.
)";

const std::string_view header = 
R"(//  Compare against values generated by the GNU Scientific Library.
//  The GSL can be found on the web: http://www.gnu.org/software/gsl/
#include <limits>
#include <cmath>
#if defined(__TEST_DEBUG)
#  include <iostream>
#  define VERIFY(A) \
  if (!(A)) \
    { \
      std::cout << "line " << __LINE__ \
	<< "  max_abs_frac = " << max_abs_frac \
	<< '\n'; \
    }
#else
#  include <testsuite_hooks.h>
#endif
#include <specfun_testcase.h>
)";

const std::string_view riemann_limits = 
R"(// This can take long on simulators, timing out the test.
// { dg-options "-DMAX_ITERATIONS=5" { target simulator } }

#ifndef MAX_ITERATIONS
#define MAX_ITERATIONS (sizeof(data001) / sizeof(testcase_riemann_zeta<double>))
#endif
)";


/// A class to abstract the scalar data type in a generic way.
template<typename Tp>
  struct num_traits
  {
    using __value_type = Tp;
  };

template<>
  template<typename Tp>
    struct num_traits<std::complex<Tp>>
    {
      using __value_type = typename std::complex<Tp>::value_type;
    };

template<typename Tp>
  using num_traits_t = typename num_traits<Tp>::__value_type;

template<typename Tp>
  struct type_strings
  {
    static const std::string_view
    type()
    { return std::string_view(""); }

    static const std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<float>
  {
    static const std::string_view
    type()
    { return std::string("float"); }

    static const std::string_view
    suffix()
    { return std::string_view("F"); }
  };

template<>
  struct type_strings<double>
  {
    static const std::string_view
    type()
    { return std::string_view("double"); }

    static const std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<long double>
  {
    static const std::string_view
    type()
    { return std::string_view("long double"); }

    static const std::string_view
    suffix()
    { return std::string_view("L"); }
  };

template<>
  struct type_strings<__float128>
  {
    static const std::string_view
    type()
    { return std::string_view("__float128"); }

    static const std::string_view
    suffix()
    { return std::string_view("Q"); }
  };

template<>
  struct type_strings<std::complex<float>>
  {
    static const std::string_view
    type()
    { return std::string_view("std::complex<float>"); }

    static const std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<std::complex<double>>
  {
    static const std::string_view
    type()
    { return std::string_view("std::complex<double>"); }

    static const std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<std::complex<long double>>
  {
    static const std::string_view
    type()
    { return std::string_view("std::complex<long double>"); }

    static const std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<std::complex<__float128>>
  {
    static const std::string_view
    type()
    { return std::string_view("std::complex<__float128>"); }

    static const std::string_view
    suffix()
    { return std::string_view(""); }
  };


///
///  @brief  Append two testcase vectors.
///
template<typename Tp>
  std::vector<Tp>
  fill2(std::vector<Tp>&& arg, std::vector<Tp>&& more)
  {
    std::vector<Tp> ret;
    ret.reserve(arg.size() + more.size());
    ret.insert(ret.end(), arg.begin(), arg.end());
    ret.insert(ret.end(), more.begin(), more.end());
    return ret;
  }

///
///  @brief  Fill an array with evenly spaced values between two limits.
///
template<typename Tp>
  std::vector<Tp>
  fill2(Tp start, Tp stop, unsigned int num_points)
  {
    std::vector<Tp> argument;

    if (num_points == 1)
      {
	argument.push_back(start);
	return argument;
      }

    auto delta = (stop - start) / (num_points - 1);
    auto min = std::min(start, stop);
    auto max = std::max(start, stop);
    auto clamp = [min, max](Tp x)
		 -> Tp
		 { return std::clamp(x, min, max); };
    for (unsigned int i = 0; i < num_points; ++i)
      argument.push_back(clamp(start + i * delta));

    return argument;
  }


///
///  @brief  Find a nice round number limit.
///
template<typename Tp>
  Tp
  get_tolerance(Tp delta, Tp min_tol, bool & ok)
  {
    using Val = num_traits_t<Tp>;

    const auto abs_delta = std::abs(delta);
    //  Make this some number larger because you lose some accuracy writing and reading.
    const auto eps = Tp{10} * std::numeric_limits<Val>::epsilon();
    Tp tol = min_tol;
    while (tol > std::abs(delta)
	&& tol > eps)
      {
	if (Tp(0.5L) * tol > abs_delta
	 && Tp(0.5L) * tol > eps)
	  {
	    if (Tp(0.2L) * tol > abs_delta
	     && Tp(0.2L) * tol > eps)
	      {
		if (Tp(0.1L) * tol > abs_delta
		 && Tp(0.1L) * tol > eps)
		  {
		    tol *= Tp(0.1L);
		  }
		else
		  {
		    tol *= Tp(0.2L);
		    break;
		  }
	      }
	    else
	      {
		tol *= Tp(0.5L);
		break;
	      }
	  }
	else
	  break;
      }
    ok = true;
    if (tol < min_tol && tol <= abs_delta)
      {
	ok = false;
	std::cerr << "**** Error in get_tolerance:"
		  << " abs(delta)=" << abs_delta
		  << " tol=" << tol
		  << '\n';
      }
    if (tol == min_tol && tol < abs_delta)
      {
	ok = false;
	std::cerr << "Note in get_tolerance:"
		  << " delta=" << delta
		  << " tol=" << tol
		  << '\n';
      }

    //  Somehow, we seem to need extra space to get the tests to pass.
    //  TODO: Figure this out.
    return Tp(50.0L) * tol;
  }

/**
 * A concept for mask functions.
 */
template<typename Fun, typename... Arg>
  concept bool
  MaskFunction()
  {
    return requires(Fun mask, Arg... args)
    {
      { mask(args...) } -> bool;
    };
  }

/**
 * A concept for special functions.
 */
template<typename Fun, typename Ret, typename... Arg>
  concept bool
  SpecialFunction()
  {
    return requires(Fun specfun, Arg... args)
    {
      { specfun(args...) } -> Ret;
    };
  }

/**
 * A class for function argument sets.
 */
template<typename Arg>
  struct argument
  {
    using value_type = Arg;

    argument(std::string_view arg,
	     const std::vector<Arg> & val)
    : name(arg),
      value(val)
    {}

    std::string_view name;
    std::vector<Arg> value;
  };

template<typename Arg>
  argument<Arg>
  make_argument(std::string_view name,
		const std::vector<Arg> & val)
  { return argument<Arg>(name, val); }

/**
 * A class for test function - the function to be tested.
 */
template<typename TestFun,
	 typename Ret, typename... Arg>
requires SpecialFunction<TestFun, Ret, Arg...>()
  struct test_function
  {
    test_function(std::string_view name,
		  TestFun func)
    : funcname(name),
      function(func)
    {}

    std::string_view funcname;
    TestFun function;
  };

template<typename TestFun,
	 typename Ret, typename... Arg>
requires SpecialFunction<TestFun, Ret, Arg...>()
  test_function<TestFun, Ret, Arg...>
  make_test_function(std::string_view name,
		     TestFun func)
  { return test_function<TestFun, Ret, Arg...>(name, func); }

/**
 * A class for baseline function - the function that is a baseline.
 */
template<typename BaselineFun,
	 typename Ret, typename... Arg>
requires SpecialFunction<BaselineFun, Ret, Arg...>()
  struct baseline_function
  {
    baseline_function(std::string_view src,
		      std::string_view name,
		      BaselineFun func)
    : source(src),
      funcname(name),
      function(func)
    {}

    std::string_view source;
    std::string_view funcname;
    BaselineFun function;
  };

template<typename BaselineFun,
	 typename Ret, typename... Arg>
requires SpecialFunction<BaselineFun, Ret, Arg...>()
  baseline_function<BaselineFun, Ret, Arg...>
  make_baseline_function(std::string_view src,
			 std::string_view name,
			 BaselineFun func)
  { return baseline_function<BaselineFun, Ret, Arg...>(src, name, func); }


/**
 * A class for testcases.
 */
template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
requires MaskFunction<MaskFun, Arg...>()
      && SpecialFunction<TestFun, Ret, Arg...>()
      && SpecialFunction<BaselineFun, Ret, Arg...>()
  class testcase2
  {

  public:

    testcase2(test_function<TestFun, Ret, Arg...> test,
	      baseline_function<BaselineFun, Ret, Arg...> base,
	      MaskFun mask,
	      std::string_view structname,
	      Arg... arg)
    : _M_testfun(test),
      _M_basefun(base),
      _M_maskfun(mask),
      _M_structname(structname),
      _M_range(arg...)
    {
      funcall = get_funcall();

      // Some functions taking only integral arguments must have explicit template args.
      auto _M_all_integral = !(std::is_floating_point_v<Arg> || ...);
      if (_M_all_integral)
	_M_templparm += "<Ret>";
    }

    void
    operator()(std::ostream & output) const
    {
      // Some targets take too long for riemann zeta.
      bool riemann_zeta_limits
	= (_S_arity == 1 && _M_testfun.funcname == "riemann_zeta");

      output << boilerplate << '\n';
      output << "//  " << _M_testfun.funcname << '\n' << '\n';
      if (riemann_zeta_limits)
	output << riemann_limits << '\n';
      output << header;

      write_test_data(output, std::index_sequence_for<Arg...>{});

      write_main(output, riemann_zeta_limits, get_funcall().str());
    }

  private:

    static const auto _S_arity = sizeof...(Arg);

    std::ostringstream get_funcall() const;

    template<std::size_t... Index>
      void get_funcall_help(std::ostringstream &, std::index_sequence<Index...>) const;

    template<std::size_t... Index>
      void write_test_data(std::ostream &, std::index_sequence<Index...>) const;

    template<std::size_t Index, typename Saucer>
      struct Test_Data_Help;

    void write_main(std::ostream & output,
		    bool riemann_zeta_limits,
		    std::string_view funcall) const;

    std::string_view _M_structname;
    mutable unsigned int _M_start_test = 1;
    test_function<TestFun, Ret, Arg...> _M_testfun;
    baseline_function<BaselineFun, Ret, Arg...> _M_basefun;
    MaskFun _M_maskfun;
    std::tuple<argument<Arg>...> _M_range;
    std::ostringstream funcall;
    std::string _M_templparm;
    bool _M_all_integral = true;
  };

/**
 * A testcase maker - ADL to the rescue.
 */
template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
  testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>
  make_testcase2(test_function<TestFun, Ret, Arg...> test,
		 baseline_function<BaselineFun, Ret, Arg...> base,
		 MaskFun mask,
		 std::string_view structname,
		 Ret ret, argument<Arg>... arg)
  {
    return testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>(test, base, mask, structname, arg...);
  }

/**
 * Generate and write the test data.  Accumulate statistics.
 */
template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
  template<std::size_t... Index>
    void
    testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>::
    write_test_data(std::ostream & output, std::index_sequence<Index...>) const
    {
      //using Ret = decltype(TestFun(Arg{}.value[0]...));
      using Val = num_traits_t<Ret>;

      const int old_prec = output.precision(std::numeric_limits<Val>::max_digits10);
      output.flags(std::ios::showpoint);

      constexpr auto Sign = _Spaceship<std::size_t>{}(1, _S_arity);
      using Saucer = _SpaceshipType<std::size_t, Sign>;
      Test_Data_Help<0, Saucer>{}(output, *this);
    }

/**
 * Generate and write the test data.  Accumulate statistics.
 */
template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
  template<std::size_t Index>
    struct testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>::
    Test_Data_Help<Index, _SpaceLess<std::size_t>>
    {
      void
      operator()(std::ostream & output, const testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...> & outer) const
      {
	using Outer = testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>;
	constexpr auto Sign = _Spaceship<std::size_t>{}(Index + 1, Outer::_S_arity);
	using Saucer = _SpaceshipType<std::size_t, Sign>;
        for (const auto x : std::get<Index>(outer._M_range).value)
	  Test_Data_Help<Index + 1, Saucer>{}(output, outer);
      }
    };

/**
 * Generate and write the test data.  Accumulate statistics.
 */
template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
  template<std::size_t Index>
    struct testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>::
    Test_Data_Help<Index, _SpaceEqual<std::size_t>>
    {
      void
      operator()(std::ostream &, const testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...> &) const;
    };

template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
  template<std::size_t Index>
    void
    testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>::
    Test_Data_Help<Index, _SpaceEqual<std::size_t>>::
    operator()(std::ostream & output, const testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...> & outer) const
    {
      //using Ret = decltype(TestFun(Arg{}.value[0]...));
      using Val = num_traits_t<Ret>;
      constexpr auto eps = std::numeric_limits<Val>::epsilon();
      constexpr auto inf = std::numeric_limits<Val>::infinity();
      constexpr auto NaN = std::numeric_limits<Val>::quiet_NaN();
      constexpr auto ret_complex = emsr::is_complex_v<Ret>;

      auto numname = type_strings<Val>::type();
      auto structname = outer._M_structname.to_string() + '<' + numname + '>';
      auto baseline = outer._M_basefun.source;
      auto arg = std::get<Index>(outer._M_range).name;

      std::vector<std::tuple<Arg...>> crud;
      _Statistics<Ret> raw_stats;
      _Statistics<decltype(std::abs(Ret{}))> abs_stats;
      auto num_divergences = 0;
      std::tuple<Ret, Arg...> last_divergence;
      auto max_abs_frac = Val{-1};
      for (const auto x : std::get<Index>(outer._M_range).value)
	{
	  try
	    {
	      const auto f1 = outer._M_testfun.function(x);
	      const auto f2 = outer._M_basefun.function(x);

	      if (std::abs(f1) == inf || std::abs(f2) == inf)
		{
		  ++num_divergences;
		  last_divergence = std::make_tuple(f1, f2, x);
		  if (num_divergences <= 3)
		    output << "\n// Divergence at"
			   << " " << arg << "=" << x
			   << " f=" << f1
			   << " f_" << baseline << "=" << f2;
		  continue;
		}

	      if (std::isnan(std::real(f1)) || std::isnan(std::real(f2)))
		{
		  output << "\n// Failure at"
			 << " " << arg << "=" << x
			 << " f=" << f1
			 << " f_" << baseline << "=" << f2;
		  break;
		}

	      const auto diff = f1 - f2;
	      raw_stats << diff;
	      abs_stats << std::abs(diff);
	      if (std::abs(f2) > Val{10} * eps && std::abs(f1) > Val{10} * eps)
		{
		  const auto frac = diff / f2;
		  if (std::abs(frac) > max_abs_frac)
		    max_abs_frac = std::abs(frac);
		}
	      crud.emplace_back(f2, x);
	    }
	  catch (...)
	    {
	      continue;
	    }
	}
      if (num_divergences > 0)
	{
	  if (num_divergences > 4)
	    output << "\n// ...";
	  output << "\n// Divergence at"
		 << " " << arg << "=" << std::get<2>(last_divergence)
		 << " f=" << std::get<0>(last_divergence)
		 << " f_" << baseline << "=" << std::get<1>(last_divergence);
	  num_divergences = 0;
	}

      if (abs_stats.max() >= Val{0} && max_abs_frac >= Val{0})
	{
	  bool tol_ok = false;
	  const auto min_tol = Val{1.0e-3L};
	  const auto frac_toler = get_tolerance(max_abs_frac, min_tol, tol_ok);
	  std::string tname;
	  if (ret_complex)
	    tname = type_strings<Ret>::type().to_string();
	  std::ostringstream dataname;
	  dataname.fill('0');
	  dataname << "data" << std::setw(3) << outer._M_start_test;
	  dataname.fill(' ');
	  output << '\n';
	  output << "// Test data.\n";
	  output << "// max(|f - f_" << baseline << "|): " << abs_stats.max() << '\n';
	  output << "// max(|f - f_" << baseline << "| / |f_" << baseline << "|): " << max_abs_frac << '\n';
	  output << "// mean(f - f_" << baseline << "): " << raw_stats.mean() << '\n';
	  output << "// variance(f - f_" << baseline << "): " << raw_stats.variance() << '\n';
	  output << "// stddev(f - f_" << baseline << "): " << raw_stats.std_deviation() << '\n';
	  output.fill('0');
	  output << "const " << structname << '\n' << dataname.str() << '[' << crud.size() << "] =\n{\n";
	  output.fill(' ');
	  for (unsigned int i = 0; i < crud.size(); ++i)
	    {
	      output << "  { " << tname << std::get<0>(crud[i]) << type_strings<Ret>::suffix();
	      output << ", " << std::get<1>(crud[i]) << type_strings<Ret>::suffix();
	      output << " },\n";
	    }
	  output << "};\n";
	  output.fill('0');
	  output << "const " << numname << " toler" << std::setw(3) << outer._M_start_test << " = " << frac_toler << ";\n";
	  output.fill(' ');
	  ++outer._M_start_test;
	}
    }

/**
 * Build the function call string for the test suite.
 */
template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
  std::ostringstream
  testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>::
  get_funcall() const
  {
    std::ostringstream funcall;
    funcall << _M_testfun.funcname << _M_templparm << '(';
    get_funcall_help(funcall, std::index_sequence_for<Arg...>{});
    funcall << ')';
    return funcall;
  }

template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
  template<std::size_t... Index>
    void
    testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>::
    get_funcall_help(std::ostringstream & funcall, std::index_sequence<Index...>) const
    {
      using swallow = int[]; // guaranties left to right order
      (void)swallow
      {
	(void(funcall << (Index > 0 ? ", " : "")
		      << "data[i]."
		      << std::get<Index>(_M_range).name), 0)...
      };
      // This actually should work:
      //((funcall << (Index > 0 ? ", " : "") << "data[i].") << std::get<Index>(_M_range).name << ...);
    }

/**
 * Write the test case main function.
 */
template<typename MaskFun, typename TestFun, typename BaselineFun,
	 typename Ret, typename... Arg>
  void
  testcase2<MaskFun, TestFun, BaselineFun, Ret, Arg...>::
  write_main(std::ostream & output,
	     bool riemann_zeta_limits,
	     std::string_view funcall) const
  {
    //using Ret = decltype(TestFun(Arg{}.value[0]...));
    constexpr auto ret_complex = emsr::is_complex_v<Ret>;
    std::string ret_tname = "Ret";
    if (ret_complex)
      ret_tname = "std::complex<Ret>";

    output << '\n';
    output << "template<typename Ret, unsigned int Num>\n";
    output.fill('0');
    output << "  void\n";
    output << "  test(const " << _M_structname << "<Ret> (&data)[Num], Ret toler)\n";
    output.fill(' ');
    output << "  {\n";
    output << "    bool test __attribute__((unused)) = true;\n";
    output << "    const Ret eps = std::numeric_limits<Ret>::epsilon();\n";
    output << "    Ret max_abs_diff = Ret(-1);\n";
    output << "    Ret max_abs_frac = Ret(-1);\n";
    if (riemann_zeta_limits)
      output << "    unsigned int num_datum = MAX_ITERATIONS;\n";
    else
      output << "    unsigned int num_datum = Num;\n";
    output << "    for (unsigned int i = 0; i < num_datum; ++i)\n";
    output << "      {\n";
    output << "\tconst " << ret_tname << " f = " << funcall << ";\n";
    output << "\tconst " << ret_tname << " f0 = data[i].f0;\n";
    output << "\tconst " << ret_tname << " diff = f - f0;\n";
    output << "\tif (std::abs(diff) > max_abs_diff)\n";
    output << "\t  max_abs_diff = std::abs(diff);\n";
    output << "\tif (std::abs(f0) > Ret(10) * eps\n";
    output << "\t && std::abs(f) > Ret(10) * eps)\n";
    output << "\t  {\n";
    output << "\t    const " << ret_tname << " frac = diff / f0;\n";
    output << "\t    if (std::abs(frac) > max_abs_frac)\n";
    output << "\t      max_abs_frac = std::abs(frac);\n";
    output << "\t  }\n";
    output << "      }\n";
    output << "    VERIFY(max_abs_frac < toler);\n";
    output << "  }\n";
    output << '\n';

    output << "int\n";
    output << "main()\n";
    output << "{\n";
    output.fill('0');
    for (unsigned int t = 1; t < _M_start_test; ++t)
      output << "  test(data" << std::setw(3) << t << ", toler" << std::setw(3) << t << ");\n";
    output.fill(' ');
    output << "  return 0;\n";
    output << "}\n";
  }

#endif // TESTCASE2_TCC