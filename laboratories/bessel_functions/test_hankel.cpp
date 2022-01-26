/**
 *
 */

#include <iostream>
#include <iomanip>

#include <cmath>
#include <emsr/complex_util.h>
#include <emsr/complex_safe_math.h>
//#include <hankel.h>
#include <emsr/float128_io.h>
#include <emsr/special_functions.h>

/*
  template<typename _Tp>
    inline std::complex<_Tp>
    operator/(const std::complex<_Tp>& x, const complex__ & y)
    {
      std::complex<_Tp> r = x;
      r /= std::complex<_Tp>(y);
      return r;
    }
*/
int
main()
{
  using namespace std::literals::complex_literals;
  std::complex<double> z, nu;
  emsr::detail::Airy<std::complex<double>> airy_thing;

  std::cout.precision(std::numeric_limits<double>::max_digits10);
  auto width = 6 + std::cout.precision();

  nu = 5.0;
  z = 1.0 - 3.0i;

  try
  {
    // Cool but we nead an eater.
    auto [zx, nux, h1, h2, h1p, h2p] = emsr::detail::hankel_uniform(nu, z);

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
  }
  catch (std::exception & err)
  {
    std::cerr << "Error in test_hankel (hankel_uniform): " << err.what() << '\n';
  }

  try
  {
    auto Airy = airy_thing(z);
    std::cout << '\n';
    std::cout << " z      = " << std::setw(width) << Airy.z << '\n';
    std::cout << " Ai(z)  = " << std::setw(width) << Airy.Ai_value << '\n';
    std::cout << " Ai'(z) = " << std::setw(width) << Airy.Ai_deriv << '\n';
    std::cout << " Bi(z)  = " << std::setw(width) << Airy.Bi_value << '\n';
    std::cout << " Bi'(z) = " << std::setw(width) << Airy.Bi_deriv << '\n';
  }
  catch (std::exception & err)
  {
    std::cerr << "Error in test_hankel (airy): " << err.what() << '\n';
  }
}
