// $HOME/bin_specfun/bin/g++ -std=gnu++1z -o test_clausen test_clausen.cpp gsl_wrap.cpp gslextras/Fresnel/fresnel.c -lgsl -lgslcblas -ljacobi

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_clausen > test_clausen.txt

// g++ -std=c++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_clausen test_clausen.cpp gsl_wrap.cpp -lgsl -lgslcblas

// ./test_clausen

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include "gsl_wrap.h"


template<typename _Tp>
  void
  run_clausen()
  {
    using __cmplx = std::complex<_Tp>;

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << std::setw(width) << "t"
	      << std::setw(width) << "Re(Cl_1)"
	      << std::setw(width) << "Im(Cl_1)"
	      << std::setw(width) << "Re(Cl_2)"
	      << std::setw(width) << "Im(Cl_2)"
	      << '\n';
    std::cout << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << '\n';
    for (int i = -1000; i <= +1000; ++i)
      {
	auto w = __cmplx{_Tp(0.01Q * i)};
	auto clausen1 = __gnu_cxx::clausen(1, w);
	auto clausen2 = __gnu_cxx::clausen(2, w);
	std::cout << std::setw(width) << std::real(w)
		  << std::setw(width) << std::imag(clausen1)
		  << std::setw(width) << std::real(clausen1)
		  << std::setw(width) << std::real(clausen2)
		  << std::setw(width) << std::imag(clausen2)
		  << std::setw(width) << gsl::clausen_2(std::real(w))
		  << std::setw(width) << __gnu_cxx::clausen(2, std::real(w))
		  << std::setw(width) << __gnu_cxx::clausen_c(2, std::real(w))
		  << '\n';
      }
      std::cout << std::endl;
  }


int
main()
{
  std::cout << "\nfloat\n=====\n\n";
  run_clausen<float>();

  std::cout << "\ndouble\n======\n";
  run_clausen<double>();

  std::cout << "\nlong double\n===========\n";
  run_clausen<long double>();

  std::cout << "\n__float128\n==========\n";
  //run_clausen<__float128>();
}