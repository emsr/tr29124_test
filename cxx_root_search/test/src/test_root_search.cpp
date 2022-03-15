/*
$HOME/bin/bin/g++ -g -Wall -Wextra -Wno-psabi -Iinclude -I../include -I../polynomial/include -I../cxx_fp_utils/include -o test_root_search test_root_search.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_root_search > test_root_search.txt
*/

#include <emsr/root_search.h>
#include <emsr/polynomial.h>
#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <test_function.h>

template<typename Tp>
  Tp
  signum(Tp x)
  {
    if (x > Tp{0})
      return Tp{1};
    else if (x < Tp{0})
      return Tp{-1};
    else
      return Tp{0};
  }

template<typename Tp>
  emsr::root_state<Tp>
  func1_state(Tp x)
  { return {std::pow(x, Tp{20}) - Tp{1}, Tp{20} * std::pow(x, Tp{19})}; }

template<typename Tp>
  emsr::root_state<Tp>
  func2_state(Tp x)
  {
    return {signum(x) * std::sqrt(std::abs(x)),
	    Tp{1} / std::sqrt(std::abs(x))};
  }

template<typename Tp>
  emsr::root_state<Tp>
  func3_state(Tp x)
  { return {x * x - Tp{1.0e-8}, Tp{2} * x}; }

template<typename Tp>
  emsr::root_state<Tp>
  func4_state(Tp x)
  { return {x * std::exp(-x), std::exp(-x) - x * std::exp(-x)}; }

template<typename Tp>
  emsr::root_state<Tp>
  func5_state(Tp x)
  {
    return {Tp{1} / (Tp{1} + std::exp(x)),
	    -std::exp(x) / std::pow(Tp{1} + std::exp(x), Tp{2})};
  }

template<typename Tp>
  emsr::root_state<Tp>
  func6_state(Tp x)
  {
    return {std::pow(x - Tp{1}, Tp{7}),
	    Tp{7} * std::pow(x - Tp{1}, Tp{6})};
  }

template<typename Tp>
  emsr::root_state<Tp>
  sin_state(Tp x)
  { return {std::sin(x), std::cos(x)}; }

template<typename Tp>
  emsr::root_state<Tp>
  cos_state(Tp x)
  { return {std::cos(x), -std::sin(x)}; }

template<typename Tp>
  emsr::root_state<Tp>
  func7_state(Tp x)
  { return {-M_PI * x + M_E, -M_PI}; }

template<typename Tp>
  void
  test_roots(Tp proto = Tp{})
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();
    auto eps = std::numeric_limits<Tp>::epsilon();
    eps = std::sqrt(std::sqrt(eps));

    // We need to cast to resolve the ambiguity in overloads in std::cos.
    using fun_t = Tp (*)(Tp);

    try
      {
	Tp x_lower = Tp{1}, x_upper = Tp{2};
	if (emsr::root_bracket((fun_t)&std::cos, x_lower, x_upper))
	  {
	    std::cout << "cos root bracket: " << std::setw(w) << x_lower << " <= x <= " << std::setw(w) << x_upper << '\n';

	    std::cout << "epsilon                  : e = " << std::setw(w) << eps << '\n';

	    auto x_bisect = emsr::root_bisect((fun_t)&std::cos, x_lower, x_upper, eps);
	    std::cout << "cos root (bisect)        : x = " << std::setw(w) << x_bisect << '\n';

	    auto x_secant = emsr::root_secant((fun_t)&std::cos, x_lower, x_upper, eps);
	    std::cout << "cos root (secant)        : x = " << std::setw(w) << x_secant << '\n';

	    auto x_false_position = emsr::root_false_position((fun_t)&std::cos, x_lower, x_upper, eps);
	    std::cout << "cos root (false position): x = " << std::setw(w) << x_false_position << '\n';

	    auto x_ridder = emsr::root_ridder((fun_t)&std::cos, x_lower, x_upper, eps);
	    std::cout << "cos root (ridder)        : x = " << std::setw(w) << x_ridder << '\n';

	    auto x_brent = emsr::root_brent((fun_t)&std::cos, x_lower, x_upper, eps);
	    std::cout << "cos root (brent)         : x = " << std::setw(w) << x_brent << '\n';

	    auto x_newton = emsr::root_newton(cos_state<Tp>, x_lower, x_upper, eps);
	    std::cout << "cos root (newton)        : x = " << std::setw(w) << x_newton << '\n';

	    auto x_steffensen = emsr::root_steffensen((fun_t)&std::cos, x_lower, x_upper, eps);
	    std::cout << "cos root (steffensen)    : x = " << std::setw(w) << x_steffensen << '\n';

	    auto x_safe = emsr::root_safe(cos_state<Tp>, x_lower, x_upper, eps);
	    std::cout << "cos root (safe)          : x = " << std::setw(w) << x_safe << '\n';
	  }
	else
	  std::cout << "No cos root found in " << std::setw(w) << x_lower << " <= x <= " << std::setw(w) << x_upper << '\n';
      }
    catch(std::exception& err)
      {
	std::cout << "Error: " << err.what() << '\n';
      }

    try
      {
	Tp x_lower = Tp{2}, x_upper = Tp{3};
	if (emsr::root_bracket((fun_t)&std::cos, x_lower, x_upper, 40))
	  {
	    std::cout << "cos root bracket: " << std::setw(w) << x_lower << " <= x <= " << std::setw(w) << x_upper << '\n';
	  }
	else
	  std::cout << "No cos root found in " << std::setw(w) << x_lower << " <= x <= " << std::setw(w) << x_upper << '\n';
      }
    catch(std::exception& err)
      {
	std::cout << "Error: " << err.what() << '\n';
      }

    try
      {
	Tp x_lower = Tp{0}, x_upper = Tp{30};
	auto bracket = emsr::root_brackets((fun_t)&std::cos, x_lower, x_upper, 40);
	std::cout << "cos root brackets:\n";
	for (auto& br : bracket)
	  std::cout << "  " << std::setw(w) << br.first << " <= x <= " << std::setw(w) << br.second
		    << ": x = " << std::setw(w) << emsr::root_newton(cos_state<Tp>, br.first, br.second, eps) << '\n';
	std::cout << '\n';
      }
    catch(std::exception& err)
      {
	std::cout << "Error: " << err.what() << '\n';
      }

    try
      {
	emsr::Polynomial<Tp> P({-4, 1, -2, 3});
	Tp x_lower = Tp{-10}, x_upper = Tp{10};
	auto bracket = emsr::root_brackets(P, x_lower, x_upper, 50 * P.degree());
	std::cout << "root brackets:\n";
	for (auto& br : bracket)
	  {
	    std::cout << "  " << std::setw(w) << br.first
	       << " <= x <= " << std::setw(w) << br.second;
	    auto r = emsr::root_bisect(P, br.first, br.second, eps);
	    std::cout << ": x = " << std::setw(w) << r << '\n';
	    P /= emsr::Polynomial<Tp>(-r, Tp{1});
	  }
	std::cout << '\n';
      }
    catch(std::exception& err)
      {
	std::cout << "Error: " << err.what() << '\n';
      }
  }

template<typename Tp, typename _Searcher, typename _Test>
  void
  run_test(_Searcher srchr, _Test test, const char* str)
  {
    auto w = 8 + std::cout.precision();
    try
      {
	auto x = srchr(test.func, test.lower, test.upper, test.eps);
	std::cout << str << std::setw(w) << x << '\n';
      }
    catch (std::exception& err)
      {
	std::cout << "Error: " << err.what() << '\n';
      }
  }

/*
 * Run the GSL-derived tests.
 */
template<typename Tp>
  void
  test_moar(Tp proto = Tp{})
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    auto eps = std::numeric_limits<Tp>::epsilon();
    eps = std::sqrt(eps);

    for (const auto test : test_value_funcs<Tp>)
      {
	std::cout << '\n' << test.str << '\n';
	std::cout << "epsilon       : e = " << std::setw(w) << eps << '\n';
	std::cout << "exact         : x = " << std::setw(w) << test.root << '\n';

	//FIXME: run_test<Tp>(&emsr::root_bisect, test, "bisect        : x = ");
	try
	  {
	    auto x_bisect = emsr::root_bisect(test.func, test.lower, test.upper, eps);
	    std::cout << "bisect        : x = " << std::setw(w) << x_bisect << '\n';
	  }
	catch (std::exception& err)
	  {
	    std::cout << "Error: " << err.what() << '\n';
	  }

	try
	  {
	    auto x_secant = emsr::root_secant(test.func, test.lower, test.upper, eps);
	    std::cout << "secant        : x = " << std::setw(w) << x_secant << '\n';
	  }
	catch (std::exception& err)
	  {
	    std::cout << "Error: " << err.what() << '\n';
	  }

	try
	  {
	    auto x_false_position = emsr::root_false_position(test.func, test.lower, test.upper, eps);
	    std::cout << "false position: x = " << std::setw(w) << x_false_position << '\n';
	  }
	catch (std::exception& err)
	  {
	    std::cout << "Error: " << err.what() << '\n';
	  }

	try
	  {
	    auto x_ridder = emsr::root_ridder(test.func, test.lower, test.upper, eps);
	    std::cout << "ridder        : x = " << std::setw(w) << x_ridder << '\n';
	  }
	catch (std::exception& err)
	  {
	    std::cout << "Error: " << err.what() << '\n';
	  }

	try
	  {
	    auto x_brent = emsr::root_brent(test.func, test.lower, test.upper, eps);
	    std::cout << "brent         : x = " << std::setw(w) << x_brent << '\n';
	  }
	catch (std::exception& err)
	  {
	    std::cout << "Error: " << err.what() << '\n';
	  }

	try
	  {
	    auto x_steffensen = emsr::root_steffensen(test.func, test.lower, test.upper, eps);
	    std::cout << "steffensen    : x = " << std::setw(w) << x_steffensen << '\n';
	  }
	catch (std::exception& err)
	  {
	    std::cout << "Error: " << err.what() << '\n';
	  }
      }

    for (const auto test : test_fdf_funcs<Tp>)
      {
	std::cout << '\n' << test.str << '\n';
	std::cout << "epsilon       : e = " << std::setw(w) << eps << '\n';
	std::cout << "exact         : x = " << std::setw(w) << test.root << '\n';

	try
	  {
	  auto x_newton = emsr::root_newton(test.func, test.lower, test.upper, eps);
	  std::cout << "newton        : x = " << std::setw(w) << x_newton << '\n';
	  }
	catch (std::exception& err)
	  {
	    std::cout << "Error: " << err.what() << '\n';
	  }

	try
	  {
	  auto x_safe = emsr::root_safe(test.func, test.lower, test.upper, eps);
	  std::cout << "safe          : x = " << std::setw(w) << x_safe << '\n';
	  }
	catch (std::exception& err)
	  {
	    std::cout << "Error: " << err.what() << '\n';
	  }
      }
  }

int
main()
{
  std::cout << "\n  float";
  std::cout << "\n===============\n\n";
  test_roots<float>();

  std::cout << "\n  double";
  std::cout << "\n===============\n\n";
  test_roots<double>();

  std::cout << "\n  long double";
  std::cout << "\n===============\n\n";
  test_roots<long double>();

  std::cout << "\n  moar float";
  std::cout << "\n====================\n\n";
  test_moar<float>();

  std::cout << "\n  moar double";
  std::cout << "\n====================\n\n";
  test_moar<double>();

  std::cout << "\n  moar long double";
  std::cout << "\n====================\n\n";
  test_moar<long double>();
}
