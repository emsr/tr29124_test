
// $HOME/bin_specfun/bin/g++ -std=gnu++1z -Wall -Wextra -o test_hankel test_hankel.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_hankel > test_hankel.txt

#include <iostream>
#include <iomanip>

#include <ext/cmath>
#include "complex_util.h"
//#include "hankel.h"
#include "float128.h"


  template<typename _Tp>
    inline std::complex<_Tp>
    operator/(const std::complex<_Tp>& __x, const __complex__ & __y)
    {
      std::complex<_Tp> __r = __x;
      __r /= std::complex<_Tp>(__y);
      return __r;
    }









int
main()
{
  using namespace std::literals::complex_literals;

  std::complex<double> h1, h2, h1p, h2p;
  std::complex<double> z, nu;

  nu = 5.0;
  z = 1.0 - 3.0i;

  std::__detail::__hankel_uniform(nu, z, h1, h2, h1p, h2p);

  std::cout.precision(std::numeric_limits<double>::max_digits10);
  auto width = 6 + std::cout.precision();

  std::cout << '\n';
  std::cout << " nu     = " << std::setw(width) << nu << '\n';
  std::cout << " z      = " << std::setw(width) << z << '\n';
  std::cout << " H1(z)  = " << std::setw(width) << h1 << '\n';
  std::cout << " H1'(z) = " << std::setw(width) << h1p << '\n';
  std::cout << " H2(z)  = " << std::setw(width) << h2 << '\n';
  std::cout << " H2'(z) = " << std::setw(width) << h2p << '\n';
  std::cout << " J(z)   = " << std::setw(width) << (h1 + h2) / 2.0 << '\n';
  std::cout << " J'(z)  = " << std::setw(width) << (h1p + h2p) / 2.0 << '\n';
  std::cout << " Y(z)   = " << std::setw(width) << (h1 - h2) / 2.0i << '\n';
  std::cout << " Y'(z)  = " << std::setw(width) << (h1p - h2p) / 2.0i << '\n';

  std::complex<double> ai, aip, bi, bip;
  std::__detail::__airy(z, 1.0e-16, ai, aip, bi, bip);
  std::cout << '\n';
  std::cout << " z      = " << std::setw(width) << z << '\n';
  std::cout << " Ai(z)  = " << std::setw(width) << ai << '\n';
  std::cout << " Ai'(z) = " << std::setw(width) << aip << '\n';
  std::cout << " Bi(z)  = " << std::setw(width) << bi << '\n';
  std::cout << " Bi'(z) = " << std::setw(width) << bip << '\n';
}
