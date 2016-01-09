
// $HOME/bin_specfun/bin/g++ -o test_factorial test_factorial.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_factorial > test_factorial.txt

#include "numeric_limits.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include "float128.h"

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
    Tp prev = 1ULL;
    Tp fact = 1ULL;
    std::cout << '\n';
    for (auto i = 1ULL; i < 600ULL; ++i)
      {
	std::cout << "  {" << std::setw(4) << (i - 1)
		  << ", " << fact << suffix<Tp>()
		  << ", " << std::log(fact) << suffix<Tp>() << "},\n";

	prev = fact;
	fact *= i;
	if (fact < prev || __isinf(fact))
	  break;
      }
  }

template<typename Tp>
  void
  double_factorial()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    long long i_o = -1LL;
    Tp prev_o = 1LL;
    Tp fact_o = 1LL;
    long long i_e = 0LL;
    Tp prev_e = 1LL;
    Tp fact_e = 1LL;
    std::cout << '\n';
    for (auto i = 1LL; i < 600LL; ++i)
      {
	std::cout << "  {" << std::setw(4) << i_e << ", "
		  << fact_e << suffix<Tp>() << ", "
		  << std::log(fact_e) << suffix<Tp>() << "},\n";

	i_o += 2LL;
	prev_o = fact_o;
	fact_o *= i_o;
	if (fact_o < prev_o || __isinf(fact_o))
	  break;

	std::cout << "  {" << std::setw(4) << i_o
		  << ", " << fact_o << suffix<Tp>() << ", "
		  << std::log(fact_o) << suffix<Tp>() << "},\n";

	i_e += 2LL;
	prev_e = fact_e;
	fact_e *= i_e;
	if (fact_e < prev_e || __isinf(fact_e))
	  break;
      }
    std::cout << '\n';
    i_o = -1LL;
    prev_o = 1LL;
    fact_o = 1LL;
    for (auto i = 1LL; i < 300LL; ++i)
      {
	std::cout << "  {" << std::setw(4) << i_o
		  << ", " << fact_o << suffix<Tp>() << "},\n";
	i_o -= 2LL;
	fact_o /= i_o;
	if (fact_o > prev_o || __isinf(fact_o) || std::abs(fact_o) < std::numeric_limits<Tp>::min())
	  break;
      }
  }

int
main()
{
  factorial<float>();
  factorial<double>();
  factorial<long double>();
  factorial<__float128>();

  double_factorial<float>();
  double_factorial<double>();
  double_factorial<long double>();
  double_factorial<__float128>();
}
