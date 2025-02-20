#ifndef TESTCASE_TCC
#define TESTCASE_TCC 1

#include <cmath>
#include <complex>
#include <sstream>
#include <type_traits> // For *_v type traits.
#include <string_view>
#include <../laboratories/complex_tools/complex_compare.h> // For the Statistics min/max (maybe rethink that there)
#include <../statistics.h>

/// A class to abstract the scalar data type in a generic way.
template<typename Tp>
  struct num_traits
  {
    using value_type = Tp;
  };

template<typename Tp>
  struct num_traits<std::complex<Tp>>
  {
    using value_type = typename std::complex<Tp>::value_type;
  };

template<typename Tp>
  using __num_traits_t = typename num_traits<Tp>::value_type;

template<typename Tp>
  struct type_strings
  {
    static std::string_view
    type()
    { return std::string_view(""); }

    static std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<float>
  {
    static std::string_view
    type()
    { return std::string_view("float"); }

    static std::string_view
    suffix()
    { return std::string_view("F"); }
  };

template<>
  struct type_strings<double>
  {
    static std::string_view
    type()
    { return std::string_view("double"); }

    static std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<long double>
  {
    static std::string_view
    type()
    { return std::string_view("long double"); }

    static std::string_view
    suffix()
    { return std::string_view("L"); }
  };

template<>
  struct type_strings<__float128>
  {
    static std::string_view
    type()
    { return std::string_view("__float128"); }

    static std::string_view
    suffix()
    { return std::string_view("Q"); }
  };

template<>
  struct type_strings<std::complex<float>>
  {
    static std::string_view
    type()
    { return std::string_view("std::complex<float>"); }

    static std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<std::complex<double>>
  {
    static std::string_view
    type()
    { return std::string_view("std::complex<double>"); }

    static std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<std::complex<long double>>
  {
    static std::string_view
    type()
    { return std::string_view("std::complex<long double>"); }

    static std::string_view
    suffix()
    { return std::string_view(""); }
  };

template<>
  struct type_strings<std::complex<__float128>>
  {
    static std::string_view
    type()
    { return std::string_view("std::complex<__float128>"); }

    static std::string_view
    suffix()
    { return std::string_view(""); }
  };


///
///  @brief  Append two testcase vectors.
///
template<typename Tp>
  std::vector<Tp>
  fill_argument(std::vector<Tp>&& arg, std::vector<Tp>&& more)
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
  fill_argument(const std::pair<Tp,Tp> & range,
		const std::pair<bool,bool> & inclusive,
		unsigned int num_points)
  {
    std::vector<Tp> argument;

    if (num_points == 1)
      {
	argument.push_back(range.first);
	return argument;
      }

    auto delta = (range.second - range.first) / (num_points - 1);
    auto min = std::min(range.first, range.second);
    auto max = std::max(range.first, range.second);
    auto clamp = [min, max](Tp x)
		 -> Tp
		 {
		   if (x < min)
		     return min;
		   else if (x > max)
		     return max;
		   else
		     return x;
		 };
    argument.reserve(num_points);
    for (unsigned int i = 0; i < num_points; ++i)
      {
	if (i == 0 && ! inclusive.first)
	  continue;
	if (i == num_points - 1 && ! inclusive.second)
	  continue;

	argument.push_back(clamp(range.first + i * delta));
      }

    return argument;
  }


///
///  @brief  Find a nice round number limit.
///
template<typename Tp>
  Tp
  get_tolerance(Tp delta, Tp min_tol, bool& ok)
  {
    using Val = __num_traits_t<Tp>;

    const auto abs_delta = std::abs(delta);
    //  Make this some number larger because you lose some accuracy writing and reading.
    const auto eps = Tp{10} * std::numeric_limits<Val>::epsilon();
    Tp tol = min_tol;
    while (tol > std::abs(delta)
	&& tol > eps)
      {
	if (Tp(0.5Q) * tol > abs_delta
	 && Tp(0.5Q) * tol > eps)
	  {
	    if (Tp(0.2Q) * tol > abs_delta
	     && Tp(0.2Q) * tol > eps)
	      {
		if (Tp(0.1Q) * tol > abs_delta
		 && Tp(0.1Q) * tol > eps)
		  {
		    tol *= Tp(0.1Q);
		  }
		else
		  {
		    tol *= Tp(0.2Q);
		    break;
		  }
	      }
	    else
	      {
		tol *= Tp(0.5Q);
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
    return Tp(50.0Q) * tol;
  }


///
///  @brief  Difference two one-argument functions.
///
template<typename Ret, typename Arg1>
  unsigned int
  maketest(Ret function1(Arg1),
	   Ret function2(Arg1),
	   const std::string& test_struct_name,
	   const std::string& nsname,
	   const std::string& funcname,
	   const std::string& arg1, const std::vector<Arg1>& argument1,
	   const std::string& baseline,
	   std::ostream & output,
	   bool write_main = true, unsigned int start_test = 1)
  {
    using Val = __num_traits_t<Ret>;

    if (output.bad() || output.fail())
      throw std::runtime_error("maketest: Bad output stream.");

    output.precision(std::numeric_limits<Val>::max_digits10);
    output.flags(std::ios::showpoint);

    std::string templparm;
    if (!std::is_floating_point_v<Arg1>)
      templparm += "<Ret>";

    output << "\n//  " << funcname << '\n';

    if (start_test == 1)
      output << "\n#include \"verify.h\"\n";

    constexpr auto eps = std::numeric_limits<Val>::epsilon();
    constexpr auto inf = std::numeric_limits<Val>::infinity();
    constexpr auto ret_complex = emsr::is_complex_v<Ret>;

    std::string numname(type_strings<Val>::type());

    std::string structname = test_struct_name;
    structname += '<' + numname + '>';

    std::vector<std::tuple<Ret, Arg1>> crud;
    _Statistics<Ret> raw_stats;
    _Statistics<decltype(std::abs(Ret{}))> abs_stats;
    auto max_abs_frac = Val{-1};
    auto num_divergences = 0;
    std::tuple<Ret, Ret, Arg1> last_divergence;
    auto num_failures = 0;
    for (const auto x : argument1)
      {
	try
	  {
	    const auto f1 = function1(x);
	    const bool inf_f1 = std::abs(f1) == inf;
	    const bool nan_f1 = std::isnan(std::real(f1));
	    const auto f2 = function2(x);
	    const bool inf_f2 = std::abs(f2) == inf;
	    const bool nan_f2 = std::isnan(std::real(f2));

	    if (inf_f1 || inf_f2)
	      {
		++num_divergences;
		last_divergence = std::make_tuple(f1, f2, x);
		if (num_divergences <= 3)
		  output << "\n// Divergence at"
			 << " " << arg1 << "=" << x
			 << " f=" << f1
			 << " f_" << baseline << "=" << f2;
		if (inf_f2)
		  continue;
	      }

	    if (nan_f1 || nan_f2)
	      {
		++num_failures;
		if (num_failures <= 3)
		  output << "\n// Failure at"
		         << " " << arg1 << "=" << x
		         << " f=" << f1
		         << " f_" << baseline << "=" << f2;
		if (nan_f2)
		  continue;
	      }

	    if (!(inf_f1 || inf_f2 || nan_f1 || nan_f2))
	      {
	        const auto diff = f1 - f2;
	        raw_stats << diff;
	        abs_stats << std::abs(diff);
	        if (std::abs(f2) > Val{10} * eps && std::abs(f1) > Val{10} * eps)
	          {
		    const auto frac = diff / f2;
		    if (std::abs(frac) > max_abs_frac)
		      max_abs_frac = std::abs(frac);
	          }
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
	  output << "\n//  ...";
	output << "\n// Divergence at"
	       << " " << arg1 << "=" << std::get<2>(last_divergence)
	       << " f=" << std::get<0>(last_divergence)
	       << " f_" << baseline << "=" << std::get<1>(last_divergence);
	num_divergences = 0;
      }

    if (abs_stats.max() >= Val{0} && max_abs_frac >= Val{0})
      {
	bool tol_ok = false;
	const auto min_tol = Val{1.0e-3Q};
	const auto frac_toler = get_tolerance(max_abs_frac, min_tol, tol_ok);
	std::string tname;
	if (ret_complex)
	  tname = type_strings<Ret>::type();
	std::ostringstream dataname;
	dataname.fill('0');
	dataname << "data" << std::setw(3) << start_test;
	dataname.fill(' ');
	output << '\n';
	output << "// Test data.\n";
	output << "// max(|f - f_" << baseline << "|): " << abs_stats.max() << " at index " << abs_stats.max_index() << '\n';
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
	    output << ", 0.0" << type_strings<Ret>::suffix();
	    output << " },\n";
	  }
	output << "};\n";
	output.fill('0');
	output << "const " << numname << " toler" << std::setw(3) << start_test << " = " << frac_toler << ";\n";
	output.fill(' ');
	++start_test;
      }

    if (write_main)
      {
	std::string structname = test_struct_name;
	structname += "<Ret>";

	std::string tname = "Ret";
	if (ret_complex)
	  tname = "std::complex<Ret>";

	output << '\n';
	output << "template<typename Ret, unsigned int Num>\n";
	output.fill('0');
	output << "  int\n";
	output << "  test(const " << structname << " (&data)[Num], Ret toler)\n";
	output.fill(' ');
	output << "  {\n";
	output << "    const Ret eps = std::numeric_limits<Ret>::epsilon();\n";
	output << "    Ret max_abs_diff = Ret(-1);\n";
	output << "    Ret max_abs_frac = Ret(-1);\n";
	output << "    bool failure = false;\n";
	output << "    unsigned int num_datum = Num;\n";
	output << "    for (unsigned int i = 0; i < num_datum; ++i)\n";
	output << "      {\n";
	output << "\tconst " << tname << " f = " << nsname << "::" << funcname << templparm << "(data[i]." << arg1 << ");\n";
	output << "\tconst bool failure_f = std::isnan(f);\n";
	output << "\tif (!failure && failure_f)\n";
	output << "\t  failure = true;\n";
	output << "\tif (!failure_f)\n";
	output << "\t  {\n";
	output << "\t    const " << tname << " f0 = data[i].f0;\n";
	output << "\t    const " << tname << " diff = f - f0;\n";
	output << "\t    if (std::abs(diff) > max_abs_diff)\n";
	output << "\t      max_abs_diff = std::abs(diff);\n";
	output << "\t    if (std::abs(f0) > Ret(10) * eps\n";
	output << "\t     && std::abs(f) > Ret(10) * eps)\n";
	output << "\t      {\n";
	output << "\t\tconst " << tname << " frac = diff / f0;\n";
	output << "\t\tif (std::abs(frac) > max_abs_frac)\n";
	output << "\t\t  max_abs_frac = std::abs(frac);\n";
	output << "\t      }\n";
	output << "\t  }\n";
	output << "      }\n";
	output << "    int num_errors = 0;\n";
	output << "    VERIFY(!failure && max_abs_frac < toler);\n";
	output << "    return num_errors;\n";
	output << "  }\n";

	output << '\n';
	output << "int\n";
	output << "main()\n";
	output << "{\n";
	output.fill('0');
        if (start_test - 1 == 1)
          output << "  return test(data" << std::setw(3) << 1 << ", toler" << std::setw(3) << 1 << ");\n";
        else
          {
            output << "  int num_errors = 0;\n";
	    for (unsigned int t = 1; t < start_test; ++t)
	      output << "  num_errors += test(data" << std::setw(3) << t << ", toler" << std::setw(3) << t << ");\n";
            output << "  return num_errors;\n";
          }
	output.fill(' ');
	output << "}\n";
      }

    output.flush();

    return start_test;
  }


///
///  @brief  Difference two two-argument functions.
///
template<typename Ret, typename Arg1, typename Arg2>
  unsigned int
  maketest(Ret function1(Arg1,Arg2),
	   Ret function2(Arg1,Arg2),
	   const std::string& test_struct_name,
	   const std::string& nsname,
	   const std::string& funcname,
	   const std::string& arg1, const std::vector<Arg1>& argument1,
	   const std::string& arg2, const std::vector<Arg2>& argument2,
	   const std::string& baseline,
	   std::ostream & output,
	   bool write_main = true, unsigned int start_test = 1)
  {
    using Val = __num_traits_t<Ret>;

    if (output.bad() || output.fail())
      throw std::runtime_error("maketest: Bad output stream.");

    output.precision(std::numeric_limits<Val>::max_digits10);
    output.flags(std::ios::showpoint);

    std::string templparm;
    if (!std::is_floating_point_v<Arg1>
     && !std::is_floating_point_v<Arg2>)
      templparm += "<Ret>";

    output << "\n//  " << funcname << '\n';

    if (start_test == 1)
      output << "\n#include \"verify.h\"\n";

    constexpr auto eps = std::numeric_limits<Val>::epsilon();
    constexpr auto inf = std::numeric_limits<Val>::infinity();
    constexpr auto ret_complex = emsr::is_complex_v<Ret>;

    std::string numname(type_strings<Val>::type());

    std::string structname = test_struct_name;
    structname += '<' + numname + '>';

    for (const auto x : argument1)
      {
	std::vector<std::tuple<Ret, Arg1, Arg2>> crud;
	_Statistics<Ret> raw_stats;
	_Statistics<decltype(std::abs(Ret{}))> abs_stats;
	auto num_divergences = 0;
	std::tuple<Ret, Ret, Arg1, Arg2> last_divergence;
	auto num_failures = 0;

	auto max_abs_frac = Val{-1};
	for (const auto y : argument2)
	  {
	    try
	      {
		const auto f1 = function1(x, y);
		const bool inf_f1 = std::abs(f1) == inf;
		const bool nan_f1 = std::isnan(std::real(f1));
		const auto f2 = function2(x, y);
		const bool inf_f2 = std::abs(f2) == inf;
		const bool nan_f2 = std::isnan(std::real(f2));

		if (inf_f1 || inf_f2)
		  {
		    ++num_divergences;
		    last_divergence = std::make_tuple(f1, f2, x, y);
		    if (num_divergences <= 3)
		      output << "\n// Divergence at"
			     << " " << arg1 << "=" << x
			     << " " << arg2 << "=" << y
			     << " f=" << f1
			     << " f_" << baseline << "=" << f2;
		    if (inf_f2)
		      continue;
		  }

		if (nan_f1 || nan_f2)
		  {
		    ++num_failures;
		    if (num_failures <= 3)
		      output << "\n// Failure at"
			     << " " << arg1 << "=" << x
			     << " " << arg2 << "=" << y
			     << " f=" << f1
			     << " f_" << baseline << "=" << f2;
		    if (nan_f2)
		      continue;
		  }

		if (!(inf_f1 || inf_f2 || nan_f1 || nan_f2))
		  {
		    const auto diff = f1 - f2;
		    raw_stats << diff;
		    abs_stats << std::abs(diff);
		    if (std::abs(f2) > Val{10} * eps && std::abs(f1) > Val{10} * eps)
		      {
		        const auto frac = diff / f2;
		        if (std::abs(frac) > max_abs_frac)
		          max_abs_frac = std::abs(frac);
		      }
		  }

		crud.emplace_back(f2, x, y);
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
		   << " " << arg1 << "=" << std::get<2>(last_divergence)
		   << " " << arg2 << "=" << std::get<3>(last_divergence)
		   << " f=" << std::get<0>(last_divergence)
		   << " f_" << baseline << "=" << std::get<1>(last_divergence);
	    num_divergences = 0;
	  }

	if (abs_stats.max() >= Val{0} && max_abs_frac >= Val{0})
	  {
	    bool tol_ok = false;
	    const auto min_tol = Val{1.0e-3Q};
	    const auto frac_toler = get_tolerance(max_abs_frac, min_tol, tol_ok);
	    std::string tname;
	    if (ret_complex)
	      tname = type_strings<Ret>::type();
	    std::ostringstream dataname;
	    dataname.fill('0');
	    dataname << "data" << std::setw(3) << start_test;
	    dataname.fill(' ');
	    output << '\n';
	    output << "// Test data for " << arg1 << '=' << std::get<1>(crud[0]) << ".\n";
	    output << "// max(|f - f_" << baseline << "|): " << abs_stats.max() << " at index " << abs_stats.max_index() << '\n';
	    output << "// max(|f - f_" << baseline << "| / |f_" << baseline << "|): " << max_abs_frac << '\n';
	    output << "// mean(f - f_" << baseline << "): " << raw_stats.mean() << '\n';
	    output << "// variance(f - f_" << baseline << "): " << raw_stats.variance() << '\n';
	    output << "// stddev(f - f_" << baseline << "): " << raw_stats.std_deviation() << '\n';
	    output.fill('0');
	    output << "const " << structname << '\n' << dataname.str() << '[' << crud.size() << "] =\n{\n";
	    output.fill(' ');
	    for (unsigned int j = 0; j < crud.size(); ++j)
	      {
		output << "  { " << tname << std::get<0>(crud[j]) << type_strings<Ret>::suffix();
		output << ", " << std::get<1>(crud[j]) << type_strings<Ret>::suffix();
		output << ", " << std::get<2>(crud[j]) << type_strings<Ret>::suffix();
		output << ", 0.0" << type_strings<Ret>::suffix();
		output << " },\n";
	      }
	    output << "};\n";
	    output.fill('0');
	    output << "const " << numname << " toler" << std::setw(3) << start_test << " = " << frac_toler << ";\n";
	    output.fill(' ');
	    ++start_test;
	  }
      }

    if (write_main)
      {
	std::string structname = test_struct_name;
	structname += "<Ret>";

	std::string tname = "Ret";
	if (ret_complex)
	  tname = "std::complex<Ret>";

	output << '\n';
	output << "template<typename Ret, unsigned int Num>\n";
	output.fill('0');
	output << "  int\n";
	output << "  test(const " << structname << " (&data)[Num], Ret toler)\n";
	output.fill(' ');
	output << "  {\n";
	output << "    const Ret eps = std::numeric_limits<Ret>::epsilon();\n";
	output << "    Ret max_abs_diff = Ret(-1);\n";
	output << "    Ret max_abs_frac = Ret(-1);\n";
	output << "    bool failure = false;\n";
	output << "    unsigned int num_datum = Num;\n";
	output << "    for (unsigned int i = 0; i < num_datum; ++i)\n";
	output << "      {\n";
	output << "\tconst " << tname << " f = " << nsname << "::" << funcname << templparm << '('
	       << "data[i]." << arg1 << ", data[i]." << arg2 << ");\n";
	output << "\tconst bool failure_f = std::isnan(f);\n";
	output << "\tif (!failure && failure_f)\n";
	output << "\t  failure = true;\n";
	output << "\tif (!failure_f)\n";
	output << "\t  {\n";
	output << "\t    const " << tname << " f0 = data[i].f0;\n";
	output << "\t    const " << tname << " diff = f - f0;\n";
	output << "\t    if (std::abs(diff) > max_abs_diff)\n";
	output << "\t      max_abs_diff = std::abs(diff);\n";
	output << "\t    if (std::abs(f0) > Ret(10) * eps\n";
	output << "\t     && std::abs(f) > Ret(10) * eps)\n";
	output << "\t      {\n";
	output << "\t\tconst " << tname << " frac = diff / f0;\n";
	output << "\t\tif (std::abs(frac) > max_abs_frac)\n";
	output << "\t\t  max_abs_frac = std::abs(frac);\n";
	output << "\t      }\n";
	output << "\t  }\n";
	output << "      }\n";
	output << "    int num_errors = 0;\n";
	output << "    VERIFY(!failure && max_abs_frac < toler);\n";
	output << "    return num_errors;\n";
	output << "  }\n";
	output << '\n';

	output << "int\n";
	output << "main()\n";
	output << "{\n";
	output.fill('0');
        if (start_test - 1 == 1)
          output << "  return test(data" << std::setw(3) << 1 << ", toler" << std::setw(3) << 1 << ");\n";
        else
          {
            output << "  int num_errors = 0;\n";
	    for (unsigned int t = 1; t < start_test; ++t)
	      output << "  num_errors += test(data" << std::setw(3) << t << ", toler" << std::setw(3) << t << ");\n";
            output << "  return num_errors;\n";
          }
	output.fill(' ');
	output << "}\n";
      }

    output.flush();

    return start_test;
  }


///
///  @brief  Difference two three-argument functions.
///
template<typename Ret, typename Arg1, typename Arg2, typename Arg3>
  unsigned int
  maketest(Ret function1(Arg1,Arg2,Arg3),
	   Ret function2(Arg1,Arg2,Arg3),
	   const std::string& test_struct_name,
	   const std::string& nsname,
	   const std::string& funcname,
	   const std::string& arg1, const std::vector<Arg1>& argument1,
	   const std::string& arg2, const std::vector<Arg2>& argument2,
	   const std::string& arg3, const std::vector<Arg3>& argument3,
	   const std::string& baseline,
	   std::ostream & output,
	   bool write_main = true, unsigned int start_test = 1)
  {
    using Val = __num_traits_t<Ret>;

    if (output.bad() || output.fail())
      throw std::runtime_error("maketest: Bad output stream.");

    output.precision(std::numeric_limits<Val>::max_digits10);
    output.flags(std::ios::showpoint);

    bool do_assoc_laguerre = (funcname == "assoc_laguerre");

    std::string templparm;
    if (!std::is_floating_point_v<Arg1>
     && !std::is_floating_point_v<Arg2>
     && !std::is_floating_point_v<Arg3>)
      templparm += "<Ret>";

    output << "\n//  " << funcname << '\n';

    if (start_test == 1)
      output << "\n#include \"verify.h\"\n";

    constexpr auto eps = std::numeric_limits<Val>::epsilon();
    constexpr auto inf = std::numeric_limits<Val>::infinity();
    constexpr auto ret_complex = emsr::is_complex_v<Ret>;

    std::string numname(type_strings<Val>::type());

    std::string structname = test_struct_name;
    if (do_assoc_laguerre)
      structname += "<unsigned int, " + numname + '>';
    else
      structname += '<' + numname + '>';

    for (const auto x : argument1)
      {
	for (const auto y : argument2)
	  {
	    std::vector<std::tuple<Ret, Arg1, Arg2, Arg3>> crud;
	    _Statistics<Ret> raw_stats;
	    _Statistics<decltype(std::abs(Ret{}))> abs_stats;
	    auto num_divergences = 0;
	    std::tuple<Ret, Ret, Arg1, Arg2, Arg3> last_divergence;
	    auto num_failures = 0;

	    auto max_abs_frac = Val{-1};
	    for (const auto z : argument3)
	      {
		try
		  {
		    const auto f1 = function1(x, y, z);
		    const bool inf_f1 = std::abs(f1) == inf;
		    const bool nan_f1 = std::isnan(std::real(f1));
		    const auto f2 = function2(x, y, z);
		    const bool inf_f2 = std::abs(f2) == inf;
		    const bool nan_f2 = std::isnan(std::real(f2));

		    if (inf_f1 || inf_f2)
		      {
			++num_divergences;
			last_divergence = std::make_tuple(f1, f2, x, y, z);
			if (num_divergences <= 3)
			  output << "\n// Divergence at"
				 << " " << arg1 << "=" << x
				 << " " << arg2 << "=" << y
				 << " " << arg3 << "=" << z
				 << " f=" << f1
				 << " f_" << baseline << "=" << f2;
			if (inf_f2)
			  continue;
		      }

		    if (nan_f1 || nan_f2)
		      {
			++num_failures;
			if (num_failures <= 3)
			  output << "\n// Failure at"
			         << " " << arg1 << "=" << x
			         << " " << arg2 << "=" << y
			         << " " << arg3 << "=" << z
			         << " f=" << f1
			         << " f_" << baseline << "=" << f2;
			if (nan_f2)
			  continue;
		      }

		    if (!(inf_f1 || inf_f2 || nan_f1 || nan_f2))
		      {
		        const auto diff = f1 - f2;
		        raw_stats << diff;
		        abs_stats << std::abs(diff);
		        if (std::abs(f2) > Val{10} * eps && std::abs(f1) > Val{10} * eps)
		          {
			    const auto frac = diff / f2;
			    if (std::abs(frac) > max_abs_frac)
			      max_abs_frac = std::abs(frac);
		          }
		      }

		    crud.emplace_back(f2, x, y, z);
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
		       << " " << arg1 << "=" << std::get<2>(last_divergence)
		       << " " << arg2 << "=" << std::get<3>(last_divergence)
		       << " " << arg3 << "=" << std::get<4>(last_divergence)
		       << " f=" << std::get<0>(last_divergence)
		       << " f_" << baseline << "=" << std::get<1>(last_divergence);
		num_divergences = 0;
	      }

	    if (abs_stats.max() >= Val{0} && max_abs_frac >= Val{0})
	      {
		bool tol_ok = false;
		const auto min_tol = Val{1.0e-3Q};
		const auto frac_toler = get_tolerance(max_abs_frac, min_tol, tol_ok);
		std::string tname;
		if (ret_complex)
		  tname = type_strings<Ret>::type();
		std::ostringstream dataname;
		dataname.fill('0');
		dataname << "data" << std::setw(3) << start_test;
		dataname.fill(' ');
		output << '\n';
		output << "// Test data for " << arg1 << '=' << std::get<1>(crud[0]);
		output << ", " << arg2 << '=' << std::get<2>(crud[0]) << ".\n";
		output << "// max(|f - f_" << baseline << "|): " << abs_stats.max() << " at index " << abs_stats.max_index() << '\n';
		output << "// max(|f - f_" << baseline << "| / |f_" << baseline << "|): " << max_abs_frac << '\n';
		output << "// mean(f - f_" << baseline << "): " << raw_stats.mean() << '\n';
		output << "// variance(f - f_" << baseline << "): " << raw_stats.variance() << '\n';
		output << "// stddev(f - f_" << baseline << "): " << raw_stats.std_deviation() << '\n';
		output.fill('0');
		output << "const " << structname << '\n' << dataname.str() << '[' << crud.size() << "] =\n{\n";
		output.fill(' ');
		for (unsigned int k = 0; k < crud.size(); ++k)
		  {
		    output << "  { " << tname << std::get<0>(crud[k]) << type_strings<Ret>::suffix();
		    output << ", " << std::get<1>(crud[k]) << type_strings<Ret>::suffix();
		    output << ", " << std::get<2>(crud[k]) << type_strings<Ret>::suffix();
		    output << ", \n";
		    output << "\t  " << std::get<3>(crud[k]) << type_strings<Ret>::suffix();
		    output << ", 0.0" << type_strings<Ret>::suffix();
		    output << " },\n";
		  }
		output << "};\n";
		output.fill('0');
		output << "const " << numname << " toler" << std::setw(3) << start_test << " = " << frac_toler << ";\n";
		output.fill(' ');
		++start_test;
	      }
	  }
      }


    if (write_main)
      {
	std::string structname = test_struct_name;
	if (do_assoc_laguerre)
	  structname += "<unsigned int, Ret>";
	else
	  structname += "<Ret>";

	std::string tname = "Ret";
	if (ret_complex)
	  tname = "std::complex<Ret>";

	output << '\n';
	output << "template<typename Ret, unsigned int Num>\n";
	output.fill('0');
	output << "  int\n";
	output << "  test(const " << structname << " (&data)[Num], Ret toler)\n";
	output.fill(' ');
	output << "  {\n";
	output << "    const Ret eps = std::numeric_limits<Ret>::epsilon();\n";
	output << "    Ret max_abs_diff = Ret(-1);\n";
	output << "    Ret max_abs_frac = Ret(-1);\n";
	output << "    bool failure = false;\n";
	output << "    unsigned int num_datum = Num;\n";
	output << "    for (unsigned int i = 0; i < num_datum; ++i)\n";
	output << "  	 {\n";
	output << "\tconst " << tname << " f = " << nsname << "::" << funcname << templparm << '('
	       << "data[i]." << arg1 << ", data[i]." << arg2 << ",\n";
	output << "\t\t     data[i]." << arg3 << ");\n";
	output << "\tconst bool failure_f = std::isnan(f);\n";
	output << "\tif (!failure && failure_f)\n";
	output << "\t  failure = true;\n";
	output << "\tif (!failure_f)\n";
	output << "\t  {\n";
	output << "\t    const " << tname << " f0 = data[i].f0;\n";
	output << "\t    const " << tname << " diff = f - f0;\n";
	output << "\t    if (std::abs(diff) > max_abs_diff)\n";
	output << "\t      max_abs_diff = std::abs(diff);\n";
	output << "\t    if (std::abs(f0) > Ret(10) * eps\n";
	output << "\t     && std::abs(f) > Ret(10) * eps)\n";
	output << "\t      {\n";
	output << "\t\tconst " << tname << " frac = diff / f0;\n";
	output << "\t\tif (std::abs(frac) > max_abs_frac)\n";
	output << "\t\t  max_abs_frac = std::abs(frac);\n";
	output << "\t      }\n";
	output << "\t  }\n";
	output << "      }\n";
	output << "    int num_errors = 0;\n";
	output << "    VERIFY(!failure && max_abs_frac < toler);\n";
	output << "    return num_errors;\n";
	output << "  }\n";

	output << '\n';
	output << "int\n";
	output << "main()\n";
	output << "{\n";
	output.fill('0');
        if (start_test - 1 == 1)
          output << "  return test(data" << std::setw(3) << 1 << ", toler" << std::setw(3) << 1 << ");\n";
        else
          {
            output << "  int num_errors = 0;\n";
	    for (unsigned int t = 1; t < start_test; ++t)
	      output << "  num_errors += test(data" << std::setw(3) << t << ", toler" << std::setw(3) << t << ");\n";
            output << "  return num_errors;\n";
          }
	output.fill(' ');
	output << "}\n";
      }

    output.flush();

    return start_test;
  }


///
///  @brief  Difference two four-argument functions.
///
template<typename Ret, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  unsigned int
  maketest(Ret function1(Arg1,Arg2,Arg3,Arg4),
	   Ret function2(Arg1,Arg2,Arg3,Arg4),
	   const std::string& test_struct_name,
	   const std::string& nsname,
	   const std::string& funcname,
	   const std::string& arg1, const std::vector<Arg1>& argument1,
	   const std::string& arg2, const std::vector<Arg2>& argument2,
	   const std::string& arg3, const std::vector<Arg3>& argument3,
	   const std::string& arg4, const std::vector<Arg4>& argument4,
	   const std::string& baseline,
	   std::ostream & output,
	   bool write_main = true, unsigned int start_test = 1)
  {
    using Val = __num_traits_t<Ret>;

    if (output.bad() || output.fail())
      throw std::runtime_error("maketest: Bad output stream.");

    output.precision(std::numeric_limits<Val>::max_digits10);
    output.flags(std::ios::showpoint);

    std::string templparm;
    if (!std::is_floating_point_v<Arg1>
     && !std::is_floating_point_v<Arg2>
     && !std::is_floating_point_v<Arg3>
     && !std::is_floating_point_v<Arg4>)
      templparm += "<Ret>";

    output << "\n//  " << funcname << '\n';

    if (start_test == 1)
      output << "\n#include \"verify.h\"\n";

    constexpr auto eps = std::numeric_limits<Val>::epsilon();
    constexpr auto inf = std::numeric_limits<Val>::infinity();
    constexpr auto ret_complex = emsr::is_complex_v<Ret>;

    std::string numname(type_strings<Val>::type());

    std::string structname = test_struct_name;
    structname += '<' + numname + '>';

    for (const auto w : argument1)
      {
	for (const auto x : argument2)
	  {
	    for (const auto y : argument3)
	      {
		std::vector<std::tuple<Ret, Arg1, Arg2, Arg3, Arg4>> crud;
		_Statistics<Ret> raw_stats;
		_Statistics<decltype(std::abs(Ret{}))> abs_stats;
		auto num_divergences = 0;
		std::tuple<Ret, Ret, Arg1, Arg2, Arg3, Arg4> last_divergence;
		auto num_failures = 0;

		auto max_abs_frac = Val{-1};
		for (const auto z : argument4)
		  {
		    try
		      {
			const auto f1 = function1(w, x, y, z);
			const bool inf_f1 = std::abs(f1) == inf;
			const bool nan_f1 = std::isnan(std::real(f1));
			const auto f2 = function2(w, x, y, z);
			const bool inf_f2 = std::abs(f2) == inf;
			const bool nan_f2 = std::isnan(std::real(f2));
			if (inf_f1 || inf_f2)
			  {
			    ++num_divergences;
			    last_divergence = std::make_tuple(f1, f2, w, x, y, z);
			    if (num_divergences <= 3)
			      output << "\n// Divergence at"
				     << " " << arg1 << "=" << w
				     << " " << arg2 << "=" << x
				     << " " << arg3 << "=" << y
				     << " " << arg4 << "=" << z
				     << " f=" << f1
				     << " f_" << baseline << "=" << f2;
			    if (inf_f2)
			      continue;
			  }

			if (nan_f1 || nan_f2)
			  {
			    ++num_failures;
			    if (num_failures <= 3)
			      output << "\n// Failure at"
				     << " " << arg1 << "=" << w
				     << " " << arg2 << "=" << x
				     << " " << arg3 << "=" << y
				     << " " << arg4 << "=" << z
				     << " f=" << f1
				     << " f_" << baseline << "=" << f2 << '\n';
			    if (nan_f2)
			      continue;
			  }

			if (!(inf_f1 || inf_f2 || nan_f1 || nan_f2))
			  {
			    const auto diff = f1 - f2;
			    raw_stats << diff;
			    abs_stats << std::abs(diff);
			    if (std::abs(f2) > Val{10} * eps && std::abs(f1) > Val{10} * eps)
			      {
			        const auto frac = diff / f2;
			        if (std::abs(frac) > max_abs_frac)
			          max_abs_frac = std::abs(frac);
			      }
			  }

			crud.emplace_back(f2, w, x, y, z);
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
			   << " " << arg1 << "=" << std::get<2>(last_divergence)
			   << " " << arg2 << "=" << std::get<3>(last_divergence)
			   << " " << arg3 << "=" << std::get<4>(last_divergence)
			   << " " << arg4 << "=" << std::get<5>(last_divergence)
			   << " f=" << std::get<0>(last_divergence)
			   << " f_" << baseline << "=" << std::get<1>(last_divergence);
		    num_divergences = 0;
		  }

		if (abs_stats.max() >= Val{0} && max_abs_frac >= Val{0})
		 {
		    bool tol_ok = false;
		    const auto min_tol = Val{1.0e-3Q};
		    const auto frac_toler = get_tolerance(max_abs_frac, min_tol, tol_ok);
		    std::string tname;
		    if (ret_complex)
		      tname = type_strings<Ret>::type();
		    std::ostringstream dataname;
		    dataname.fill('0');
		    dataname << "data" << std::setw(3) << start_test;
		    dataname.fill(' ');
		    output << '\n';
		    output << "// Test data for " << arg1 << '=' << std::get<1>(crud[0]);
		    output << ", " << arg2 << '=' << std::get<2>(crud[0]);
		    output << ", " << arg3 << '=' << std::get<3>(crud[0]) << ".\n";
		    output << "// max(|f - f_" << baseline << "|): " << abs_stats.max() << " at index " << abs_stats.max_index() << '\n';
		    output << "// max(|f - f_" << baseline << "| / |f_" << baseline << "|): " << max_abs_frac << '\n';
		    output << "// mean(f - f_" << baseline << "): " << raw_stats.mean() << '\n';
		    output << "// variance(f - f_" << baseline << "): " << raw_stats.variance() << '\n';
		    output << "// stddev(f - f_" << baseline << "): " << raw_stats.std_deviation() << '\n';
		    output.fill('0');
		    output << "const " << structname << '\n' << dataname.str() << '[' << crud.size() << "] =\n{\n";
		    output.fill(' ');
		    for (unsigned int l = 0; l < crud.size(); ++l)
		      {
			output << "  { " << tname << std::get<0>(crud[l]) << type_strings<Ret>::suffix();
			output << ", " << std::get<1>(crud[l]) << type_strings<Ret>::suffix();
			output << ", " << std::get<2>(crud[l]) << type_strings<Ret>::suffix();
			output << ", \n";
			output << "\t  " << std::get<3>(crud[l]) << type_strings<Ret>::suffix();
			output << ", " << std::get<4>(crud[l]) << type_strings<Ret>::suffix();
			output << ", 0.0" << type_strings<Ret>::suffix();
			output << " },\n";
		      }
		    output << "};\n";
		    output.fill('0');
		    output << "const " << numname << " toler" << std::setw(3) << start_test << " = " << frac_toler << ";\n";
		    output.fill(' ');
		    ++start_test;
		  }
	      }
	  }
      }

    if (write_main)
      {
	std::string structname = test_struct_name;
	structname += "<Ret>";

	std::string tname = "Ret";
	if (ret_complex)
	  tname = "std::complex<Ret>";

	output << '\n';
	output << "template<typename Ret, unsigned int Num>\n";
	output.fill('0');
	output << "  int\n";
	output << "  test(const " << structname << " (&data)[Num], Ret toler)\n";
	output.fill(' ');
	output << "  {\n";
	output << "    const Ret eps = std::numeric_limits<Ret>::epsilon();\n";
	output << "    Ret max_abs_diff = Ret(-1);\n";
	output << "    Ret max_abs_frac = Ret(-1);\n";
	output << "    bool failure = false;\n";
	output << "    unsigned int num_datum = Num;\n";
	output << "    for (unsigned int i = 0; i < num_datum; ++i)\n";
	output << "      {\n";
	output << "\tconst " << tname << " f = " << nsname << "::" << funcname << templparm << '('
	       << "data[i]." << arg1 << ", data[i]." << arg2 << ",\n";
	output << "\t\t     data[i]." << arg3 << ", data[i]." << arg4 << ");\n";
	output << "\tconst bool failure_f = std::isnan(f);\n";
	output << "\tif (!failure && failure_f)\n";
	output << "\t  failure = true;\n";
	output << "\tif (!failure_f)\n";
	output << "\t  {\n";
	output << "\t    const " << tname << " f0 = data[i].f0;\n";
	output << "\t    const " << tname << " diff = f - f0;\n";
	output << "\t    if (std::abs(diff) > max_abs_diff)\n";
	output << "\t      max_abs_diff = std::abs(diff);\n";
	output << "\t    if (std::abs(f0) > Ret(10) * eps\n";
	output << "\t     && std::abs(f) > Ret(10) * eps)\n";
	output << "\t      {\n";
	output << "\t\tconst " << tname << " frac = diff / f0;\n";
	output << "\t\tif (std::abs(frac) > max_abs_frac)\n";
	output << "\t\t  max_abs_frac = std::abs(frac);\n";
	output << "\t      }\n";
	output << "\t  }\n";
	output << "      }\n";
	output << "    int num_errors = 0;\n";
	output << "    VERIFY(!failure && max_abs_frac < toler);\n";
	output << "    return num_errors;\n";
	output << "  }\n";

	output << '\n';
	output << "int\n";
	output << "main()\n";
	output << "{\n";
	output.fill('0');
        if (start_test - 1 == 1)
          output << "  return test(data" << std::setw(3) << 1 << ", toler" << std::setw(3) << 1 << ");\n";
        else
          {
            output << "  int num_errors = 0;\n";
	    for (unsigned int t = 1; t < start_test; ++t)
	      output << "  num_errors += test(data" << std::setw(3) << t << ", toler" << std::setw(3) << t << ");\n";
          }
	output.fill(' ');
	output << "}\n";
      }

    output.flush();

    return start_test;
  }

#endif // TESTCASE_TCC
