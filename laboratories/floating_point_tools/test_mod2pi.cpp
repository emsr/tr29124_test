/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/math_constants.h>
#include <emsr/numeric_limits.h>

template<typename Tp>
  Tp
  mod2pi_cheap(Tp x)
  {
    const auto s_2pi = emsr::tau_v<Tp>;
    return x - s_2pi * std::floor(x / s_2pi);
  }

template<typename Tp>
  Tp
  mod2pi_cephes_wtf(Tp x)
  {
    const auto s_2pi = emsr::tau_v<Tp>;
    const auto n = std::floor(x / s_2pi);
    auto a = x - ldexp(n, 2);  /* 4n */
    a -= ldexp( n, 1);    /* 2n */
    a -= ldexp( n, -2 );  /* n/4 */
    a -= ldexp( n, -5 );  /* n/32 */
    a -= ldexp( n, -9 );  /* n/512 */
    a += ldexp( n, -15 ); /* add n/32768 */
    a -= ldexp( n, -17 ); /* n/131072 */
    a -= ldexp( n, -18 );
    a -= ldexp( n, -20 );
    a -= ldexp( n, -22 );
    a -= ldexp( n, -24 );
    a -= ldexp( n, -28 );
    a -= ldexp( n, -32 );
    a -= ldexp( n, -37 );
    a -= ldexp( n, -39 );
    a -= ldexp( n, -40 );
    a -= ldexp( n, -42 );
    a -= ldexp( n, -46 );
    a -= ldexp( n, -47 );
    return a;
  }

template<typename Tp>
  void
  test_mod2pi(Tp proto = Tp{})
  {
    std::cout.precision(emsr::max_digits10(proto));
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto x = Tp{12.34567895432Q};
    std::cout << '\n';
    std::cout << "x         = " << std::setw(w) << x << '\n';
    std::cout << "mod2pi(x) = " << std::setw(w) << mod2pi_cheap(x) << '\n';
    std::cout << "mod2pi(x) = " << std::setw(w) << mod2pi_cephes_wtf(x) << '\n';

    auto y = Tp{123456.78954326521Q};
    std::cout << '\n';
    std::cout << "y         = " << std::setw(w) << y << '\n';
    std::cout << "mod2pi(y) = " << std::setw(w) << mod2pi_cheap(y) << '\n';
    std::cout << "mod2pi(y) = " << std::setw(w) << mod2pi_cephes_wtf(y) << '\n';
  }

int
main()
{
  test_mod2pi<float>();
  test_mod2pi<double>();
  test_mod2pi<long double>();
  //test_mod2pi<__float128>();
  //test_mod2pi<mpfr::mpreal>(mpfr::mpreal(1,  256));
}
