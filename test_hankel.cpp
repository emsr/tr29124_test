// $HOME/bin/bin/g++ -std=c++14 -Wall -Wextra -o test_hankel test_hankel.cpp > /dev/null

#include <iomanip>

#include "hankel.h"

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
  std::cout << '\n';
  std::cout << " nu     = " << nu << '\n';
  std::cout << " z      = " << z << '\n';
  std::cout << " H1(z)  = " << h1 << '\n';
  std::cout << " H1'(z) = " << h1p << '\n';
  std::cout << " H2(z)  = " << h2 << '\n';
  std::cout << " H2'(z) = " << h2p << '\n';
  std::cout << " J(z)   = " << (h1 + h2) / 2.0 << '\n';
  std::cout << " J'(z)  = " << (h1p + h2p) / 2.0 << '\n';
  std::cout << " Y(z)   = " << (h1 - h2) / 2.0i << '\n';
  std::cout << " Y'(z)  = " << (h1p - h2p) / 2.0i << '\n';

  std::complex<double> ai, aip;
  int ier;
  std::__detail::__airy(z, 1.0e-16, ai, aip, ier);
  std::cout << '\n';
  std::cout << " z      = " << z << '\n';
  std::cout << " Ai(z)  = " << ai << '\n';
  std::cout << " Ai'(z) = " << aip << '\n';
}
