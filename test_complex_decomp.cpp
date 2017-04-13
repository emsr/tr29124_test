/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_complex_decomp test_complex_decomp.cpp -lquadmath
./test_complex_decomp > test_complex_decomp.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_complex_decomp test_complex_decomp.cpp -lquadmath
./test_complex_decomp
*/

#include <complex>
#include "tuple_complex.h"
#include <iostream>
#include <iomanip>

template<typename Tp>
  void
  test_complex_decomp()
  {
    std::complex<Tp> z(1, 2);

    auto [x, y] = z;
    std::cout << z << '\n';
    std::cout << ' ' << x << ' ' << y << '\n';

    z.real(3);
    z.imag(4);
    std::cout << z << '\n';
    std::cout << ' ' << x << ' ' << y << '\n';

    std::complex<Tp> w(1, 2);
    auto& [p, q] = w;

    w.real(3);
    w.imag(4);
    std::cout << w << '\n';
    std::cout << ' ' << p << ' ' << q << '\n';

    return;
  }

int
main()
{
  test_complex_decomp<double>();
  return 0;
}
