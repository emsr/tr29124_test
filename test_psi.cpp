/*
$HOME/bin_tr29124/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_psi test_psi.cpp wrap_gsl.cpp gslextras/Fresnel/fresnel.c -L/usr/local/lib -lgsl -lgslcblas
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:/usr/local/lib:$LD_LIBRARY_PATH ./test_psi
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include "wrap_gsl.h"

int
main()
{

  std::cout.precision(16);
  std::cout.flags(std::ios::showpoint);

  double x_start = -9.9375;
  double x_stop = +10.0625;
  const unsigned int max = 801;
  double delta = (x_stop - x_start) / (max - 1);
  std::cout << std::endl;
  for (unsigned int i = 0; i < max; ++i)
    {
      double x = x_start + i * delta;
      double y_gcc = __gnu_cxx::psi(x);
      double y_gsl = gsl::psi(x);
      std::cout << "psi(" << std::setw(5) << x << ") ="
                << ' ' << std::setw(20) << y_gcc
                << ' ' << std::setw(20) << y_gsl
                << ' ' << std::setw(20) << y_gcc - y_gsl
                << std::endl;
    }

  std::cout << std::endl;
  for (unsigned int i = 1; i < 100; ++i)
    {
      double x = i * 0.5;
      double y_gcc = __gnu_cxx::psi(x);
      double y_gsl = gsl::psi(x);
      std::cout << "psi(" << std::setw(5) << x << ") ="
                << ' ' << std::setw(20) << y_gcc
                << ' ' << std::setw(20) << y_gsl
                << ' ' << std::setw(20) << y_gcc - y_gsl
                << std::endl;
    }

  return 0;
}


