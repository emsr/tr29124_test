/**
 *
 */

#include <complex>
#include <iostream>
#include <iomanip>

#include <tuple_complex.h>

template<typename Tp>
  void
  test_complex_decomp()
  {
    std::complex<Tp> z(1, 2);

    auto [x, y] = z;
    std::cout << "z = " << z << '\n';
    std::cout << "z parts:" << ' ' << x << ' ' << y << '\n';

    z.real(3);
    z.imag(4);
    std::cout << "z = " << z << '\n';
    std::cout << "z parts:" << ' ' << x << ' ' << y << '\n';

    std::complex<Tp> w(1, 2);
    std::cout << "w = " << w << '\n';
    auto& [p, q] = w;

    w.real(3);
    w.imag(4);
    std::cout << "w = " << w << '\n';
    std::cout << "w parts:" << ' ' << p << ' ' << q << '\n';

    p = -5;
    q = -6;
    std::cout << "w = " << w << '\n';
    std::cout << "w parts:" << ' ' << p << ' ' << q << '\n';

    return;
  }

int
main()
{
  test_complex_decomp<double>();
  return 0;
}
