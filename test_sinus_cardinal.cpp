/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -I. -o test_sinus_cardinal test_sinus_cardinal.cpp wrap_boost.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_sinus_cardinal > test_sinus_cardinal.txt

g++ -std=c++14 -o test_sinus_cardinal test_sinus_cardinal.cpp wrap_boost.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
./test_sinus_cardinal > test_sinus_cardinal.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "wrap_boost.h"
#include "wrap_gsl.h"


template<typename _Tp>
  void
  test_sinc()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    constexpr auto pi = __gnu_cxx::__math_constants<_Tp>::__pi;

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sinc"
	      << std::setw(width) << "sinc GSL"
	      << std::setw(width) << "sinc Boost"
	      << std::setw(width) << "delta GSL"
	      << std::setw(width) << "delta Boost"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    for (int i = -40; i <= +40; ++i)
      {
	auto x = _Tp(0.1Q * i) * pi;
	auto sinc = __gnu_cxx::sinc(x);
	auto sinc_gsl = gsl::sinc(x);
	auto sinc_boost = beast::sinc(x);
	auto delta_gsl = sinc - sinc_gsl;
	auto delta_boost = sinc - sinc_boost;
	std::cout << std::setw(width) << x
		  << std::setw(width) << sinc
		  << std::setw(width) << sinc_gsl
		  << std::setw(width) << sinc_boost
		  << std::setw(width) << delta_gsl
		  << std::setw(width) << delta_boost
		  << '\n';
      }
  }

template<typename _Tp>
  void
  test_sinc_pi()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    constexpr auto pi = __gnu_cxx::__math_constants<_Tp>::__pi;

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sinc_pi"
	      << std::setw(width) << "sinc_pi GSL"
	      << std::setw(width) << "sinc_pi Boost"
	      << std::setw(width) << "delta GSL"
	      << std::setw(width) << "delta Boost"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    for (int i = -40; i <= +40; ++i)
      {
	auto x = _Tp(0.1Q * i);
	auto sinc_pi = __gnu_cxx::sinc_pi(x);
	auto sinc_pi_gsl = gsl::sinc_pi(x);
	auto sinc_pi_boost = beast::sinc_pi(x);
	auto delta_gsl = sinc_pi - sinc_pi_gsl;
	auto delta_boost = sinc_pi - sinc_pi_boost;
	std::cout << std::setw(width) << x
		  << std::setw(width) << sinc_pi
		  << std::setw(width) << sinc_pi_gsl
		  << std::setw(width) << sinc_pi_boost
		  << std::setw(width) << delta_gsl
		  << std::setw(width) << delta_boost
		  << '\n';
      }
  }

template<typename _Tp>
  void
  test_sinhc()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sinhc"
	      << std::setw(width) << "sinhc GSL"
	      << std::setw(width) << "sinhc Boost"
	      << std::setw(width) << "delta GSL"
	      << std::setw(width) << "delta Boost"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    constexpr auto pi = __gnu_cxx::__math_constants<_Tp>::__pi;

    for (int i = -40; i <= +40; ++i)
      {
	auto x = _Tp(0.1Q * i);
	auto sinhc = __gnu_cxx::sinhc(x);
	auto sinhc_gsl = gsl::sinhc(x);
	auto sinhc_boost = beast::sinhc(x);
	auto delta_gsl = sinhc - sinhc_gsl;
	auto delta_boost = sinhc - sinhc_boost;
	std::cout << std::setw(width) << x
		  << std::setw(width) << sinhc
		  << std::setw(width) << sinhc_gsl
		  << std::setw(width) << sinhc_boost
		  << std::setw(width) << delta_gsl
		  << std::setw(width) << delta_boost
		  << '\n';
      }
  }

template<typename _Tp>
  void
  test_sinhc_pi()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    constexpr auto pi = __gnu_cxx::__math_constants<_Tp>::__pi;

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sinhc_pi"
	      << std::setw(width) << "sinhc_pi GSL"
	      << std::setw(width) << "sinhc_pi Boost"
	      << std::setw(width) << "delta GSL"
	      << std::setw(width) << "delta Boost"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    for (int i = -40; i <= +40; ++i)
      {
	auto x = _Tp(0.1Q * i);
	auto sinhc_pi = __gnu_cxx::sinhc_pi(x);
	auto sinhc_pi_gsl = gsl::sinhc_pi(x);
	auto sinhc_pi_boost = beast::sinhc_pi(x);
	auto delta_gsl = sinhc_pi - sinhc_pi_gsl;
	auto delta_boost = sinhc_pi - sinhc_pi_boost;
	std::cout << std::setw(width) << x
		  << std::setw(width) << sinhc_pi
		  << std::setw(width) << sinhc_pi_gsl
		  << std::setw(width) << sinhc_pi_boost
		  << std::setw(width) << delta_gsl
		  << std::setw(width) << delta_boost
		  << '\n';
      }
  }

int
main()
{
  test_sinc<double>();

  test_sinc_pi<double>();

  test_sinhc<double>();

  test_sinhc_pi<double>();
}