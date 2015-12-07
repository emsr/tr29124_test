// $HOME/bin_specfun/bin/g++ -std=gnu++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_chebyshev test_chebyshev.cpp gsl_wrap.cpp -lgsl -lgslcblas -lquadmath 2> err.txt

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_chebyshev

// $HOME/bin/bin/g++ -std=gnu++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_chebyshev test_chebyshev.cpp gsl_wrap.cpp -lgsl -lgslcblas -lquadmath 2> err.txt

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_chebyshev

#include <iostream>
#include <iomanip>

#include "chebyshev.h"
#include <cmath>

int
main()
{
  _Chebyshev<double> c1{-1.0, 1.0,
   {-1.142022680371168e0, 6.5165112670737e-3,
     3.087090173086e-4, -3.4706269649e-6,
     6.9437664e-9, 3.67765e-11, -1.356e-13}};
  _Chebyshev<double> c2{-1.0, 1.0,
   { 1.843740587300905e0, -7.86528408447867e-2,
     1.2719271366546e-3, -4.9717367042e-6,
    -3.31261198e-8, 2.423096e-10, -1.702e-13, -1.49e-15}};

  _Chebyshev<double> cdilog(-4.0, 1.0, 40, __gnu_cxx::dilog);
  for (int i = 0; i <= 500; ++i)
    {
      auto x = -4.0 + i * 0.01;
      std::cout << "x = " << x
		<< "  dilog = " << cdilog(x) << '\n';

    }
}
