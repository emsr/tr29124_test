/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -I. -o test_factorial test_factorial.cpp -lquadmath
./test_factorial > test_factorial.txt

$HOME/bin/bin/g++ -std=gnu++17 -I. -o test_factorial test_factorial.cpp -lquadmath
./test_factorial > test_factorial.txt
*/

#include <bits/numeric_limits.h>
#include <limits>
#include <iostream>
#include <iomanip>
#include <bits/float128_io.h>

template<typename Tp>
  std::string
  suffix()
  { return ""; }

template<>
  std::string
  suffix<float>()
  { return "F"; }

template<>
  std::string
  suffix<long double>()
  { return "L"; }

template<>
  std::string
  suffix<__float128>()
  { return "L"; }

template<typename Tp>
  void
  factorial()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto prev = Tp{1};
    auto fact = Tp{1};
    std::cout << "{\n";
    for (auto i = 1ULL; i < 1000LL; ++i)
      {
	std::cout << "  {" << std::setw(4) << (i - 1)
		  << ", " << fact << suffix<Tp>()
		  << ", " << std::log(fact) << suffix<Tp>() << "},\n";

	prev = fact;
	fact *= i;
	if (fact < prev || __isinf(fact))
	  break;
      }
    std::cout << "};\n";
  }

template<typename Tp>
  void
  double_factorial()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    long long i_o = -1LL;
    auto prev_o = Tp{1};
    auto fact_o = Tp{1};
    long long i_e = 0LL;
    auto prev_e = Tp{1};
    auto fact_e = Tp{1};
    std::cout << "{\n";
    for (auto i = 1LL; i < 1000LL; ++i)
      {
	std::cout << "  {" << std::setw(4) << i_e
		  << ", " << fact_e << suffix<Tp>()
		  << ", " << std::log(fact_e) << suffix<Tp>() << "},\n";

	i_o += 2LL;
	prev_o = fact_o;
	fact_o *= i_o;
	if (fact_o < prev_o || __isinf(fact_o))
	  break;

	std::cout << "  {" << std::setw(4) << i_o
		  << ", " << fact_o << suffix<Tp>()
		  << ", " << std::log(fact_o) << suffix<Tp>() << "},\n";

	i_e += 2LL;
	prev_e = fact_e;
	fact_e *= i_e;
	if (fact_e < prev_e || __isinf(fact_e))
	  break;
      }
    std::cout << "};\n";
  }


template<typename Tp>
  void
  neg_double_factorial()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    std::cout << "{\n";
    auto i_o = -1LL;
    auto prev_o = Tp{1};
    auto fact_o = Tp{1};
    for (auto i = 1LL; i < 1000LL; ++i)
      {
	std::cout << "  {" << std::setw(4) << i_o
		  << ", " << fact_o << suffix<Tp>()
		  << ", " << std::log(std::abs(fact_o)) << suffix<Tp>() << "},\n";
	i_o -= 2LL;
	fact_o /= i_o;
	if (fact_o > prev_o || __isinf(fact_o) || std::abs(fact_o) < std::numeric_limits<Tp>::min())
	  break;
      }
    std::cout << "};\n";
  }

int
main()
{
  std::cout << "\n_S_factorial_table<float>\n";
  factorial<float>();
  std::cout << "\n_S_factorial_table<double>\n";
  factorial<double>();
  std::cout << "\n_S_factorial_table<long double>\n";
  factorial<long double>();
  std::cout << "\n_S_factorial_table<__float128>\n";
  factorial<__float128>();

  std::cout << "\n_S_double_factorial_table<float>\n";
  double_factorial<float>();
  std::cout << "\n_S_double_factorial_table<double>\n";
  double_factorial<double>();
  std::cout << "\n_S_double_factorial_table<long double>\n";
  double_factorial<long double>();
  std::cout << "\n_S_double_factorial_table<__float128>\n";
  double_factorial<__float128>();

  std::cout << "\n_S_neg_double_factorial_table<float>\n";
  neg_double_factorial<float>();
  std::cout << "\n_S_neg_double_factorial_table<double>\n";
  neg_double_factorial<double>();
  std::cout << "\n_S_neg_double_factorial_table<long double>\n";
  neg_double_factorial<long double>();
  std::cout << "\n_S_neg_double_factorial_table<__float128>\n";
  neg_double_factorial<__float128>();
}
