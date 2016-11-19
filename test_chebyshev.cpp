/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_chebyshev test_chebyshev.cpp -lquadmath
./test_chebyshev

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_chebyshev test_chebyshev.cpp -lquadmath
./test_chebyshev

g++ -std=gnu++14 -o test_chebyshev test_chebyshev.cpp -lquadmath
./test_chebyshev
*/

#include <iostream>
#include <iomanip>

#include <bits/float128_io.h>
#include "chebyshev.h"
#include <cmath>

int
main()
{
  __gnu_cxx::_Chebyshev<double> c1{-1.0, 1.0,
   {-1.142022680371168e0, 6.5165112670737e-3,
     3.087090173086e-4, -3.4706269649e-6,
     6.9437664e-9, 3.67765e-11, -1.356e-13}};
  __gnu_cxx::_Chebyshev<double> c2{-1.0, 1.0,
   { 1.843740587300905e0, -7.86528408447867e-2,
     1.2719271366546e-3, -4.9717367042e-6,
    -3.31261198e-8, 2.423096e-10, -1.702e-13, -1.49e-15}};

  //
  //__gnu_cxx::dilog
  __gnu_cxx::_Chebyshev<double> cdilog(-4.0, 1.0, 40, [](double x){return x * x;});
  std::cout << cdilog << '\n';
  for (int i = 0; i <= 500; ++i)
    {
      auto x = -4.0 + i * 0.01;
      std::cout << "x = " << x
		<< "  dilog(x) = " << cdilog(x) << '\n';
    }
  cdilog.truncate(1.0e-12);
  std::cout << cdilog << '\n';

  //
  __gnu_cxx::_Chebyshev<__float128> cquad(-4.0Q, +7.0Q, 40, [](__float128 x){return 3.5Q + x * (13.25Q + x * 6.375Q);});
  std::cout << cquad << '\n';
  for (int i = 0; i <= 500; ++i)
    {
      auto x = -4.0Q + i * 0.01Q;
      std::cout << "x = " << x
		<< "  quad(x) = " << cquad(x) << '\n';
    }
  cquad.truncate<float>();
  std::cout << cquad << '\n';
}
