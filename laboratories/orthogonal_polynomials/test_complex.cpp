/*
$HOME/bin/bin/g++ -std=c++2a -g -Wall -Wextra -Wno-psabi -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -I../../quadrature/include -I../../cxx_summation/include -o test_complex test_complex.cpp -lquadmath
./test_complex > test_complex.txt
*/

#include <cmath>
#include <complex>

int
main()
{
  std::legendre(1, std::complex<double>{});
  __gnu_cxx::legendre_q(1, std::complex<double>{});
  std::assoc_legendre(1, 0, std::complex<double>{});
  __gnu_cxx::assoc_legendre_q(1, 0, std::complex<double>{});
}
