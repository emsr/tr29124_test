/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_complex128 test_complex128.cpp -lquadmath
./test_complex128
*/

#include <limits>
#include <iostream>

#include <complex>

int
main()
{
  __float128 x{2.5Q};
  std::cout.precision(std::numeric_limits<__float128>::max_digits10);
  std::complex<__float128> w{0.125Q, -0.5Q}, z{3.5Q, 1.5Q};

  std::cout << "abs(z)    = " << std::abs(z) << '\n';
  std::cout << "arg(z)    = " << std::arg(z) << '\n';
  std::cout << "imag(z)   = " << std::imag(z) << '\n';
  std::cout << "real(z)   = " << std::real(z) << '\n';
  std::cout << "acos(z)   = " << std::acos(z) << '\n';
  std::cout << "acosh(z)  = " << std::acosh(z) << '\n';
  std::cout << "asin(z)   = " << std::asin(z) << '\n';
  std::cout << "asinh(z)  = " << std::asinh(z) << '\n';
  std::cout << "atan(z)   = " << std::atan(z) << '\n';
  std::cout << "atanh(z)  = " << std::atanh(z) << '\n';
  std::cout << "cos(z)    = " << std::cos(z) << '\n';
  std::cout << "cosh(z)   = " << std::cosh(z) << '\n';
  std::cout << "exp(z)    = " << std::exp(z) << '\n';
  std::cout << "expi(x)   = " << std::expi(x) << '\n';
  std::cout << "log(z)    = " << std::log(z) << '\n';
  std::cout << "log10(z)  = " << std::log10(z) << '\n';
  std::cout << "conj(z)   = " << std::conj(z) << '\n';
  std::cout << "pow(z, w) = " << std::pow(z, w) << '\n';
  std::cout << "proj(z)   = " << std::proj(z) << '\n';
  std::cout << "sin(z)    = " << std::sin(z) << '\n';
  std::cout << "sinh(z)   = " << std::sinh(z) << '\n';
  std::cout << "sqrt(z)   = " << std::sqrt(z) << '\n';
  std::cout << "tan(z)    = " << std::tan(z) << '\n';
  std::cout << "tanh(z)   = " << std::tanh(z) << '\n';

  std::cout << "cabsq(z)    = " << cabsq(z) << '\n';
  std::cout << "cargq(z)    = " << cargq(z) << '\n';
  std::cout << "cimagq(z)   = " << cimagq(z) << '\n';
  std::cout << "crealq(z)   = " << crealq(z) << '\n';
  std::cout << "cacosq(z)   = " << cacosq(z) << '\n';
  std::cout << "cacoshq(z)  = " << cacoshq(z) << '\n';
  std::cout << "casinq(z)   = " << casinq(z) << '\n';
  std::cout << "casinhq(z)  = " << casinhq(z) << '\n';
  std::cout << "catanq(z)   = " << catanq(z) << '\n';
  std::cout << "catanhq(z)  = " << catanhq(z) << '\n';
  std::cout << "ccosq(z)    = " << ccosq(z) << '\n';
  std::cout << "ccoshq(z)   = " << ccoshq(z) << '\n';
  std::cout << "cexpq(z)    = " << cexpq(z) << '\n';
  std::cout << "cexpiq(x)   = " << cexpiq(x) << '\n';
  std::cout << "clogq(z)    = " << clogq(z) << '\n';
  std::cout << "clog10q(z)  = " << clog10q(z) << '\n';
  std::cout << "conjq(z)    = " << conjq(z) << '\n';
  std::cout << "cpowq(z, w) = " << cpowq(z, w) << '\n';
  std::cout << "cprojq(z)   = " << cprojq(z) << '\n';
  std::cout << "csinq(z)    = " << csinq(z) << '\n';
  std::cout << "csinhq(z)   = " << csinhq(z) << '\n';
  std::cout << "csqrtq(z)   = " << csqrtq(z) << '\n';
  std::cout << "ctanq(z)    = " << ctanq(z) << '\n';
  std::cout << "ctanhq(z)   = " << ctanhq(z) << '\n';
}
