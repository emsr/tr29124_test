/**
 *
 */

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <bits/float128_io.h>
#include <mpreal.h>

template<std::size_t _NumBits>
  struct int_t;

template<>
  struct int_t<16>
  {
    using type = short;
  };

template<>
  struct int_t<32>
  {
    using type = int;
  };

template<>
  struct int_t<64>
  {
    using type = long;
  };

template<>
  struct int_t<128>
  {
    using type = __int128;
  };

template<std::size_t _NumBits>
  struct cordic
  {
    mpfr::mpreal
    k_circular()
    {
      mpfr::mpreal one(1, 4 * _NumBits);
      mpfr::mpreal two(2, 4 * _NumBits);
      mpfr::mpreal K(1, 4 * _NumBits);
      for (auto k = 0u; k < _NumBits; ++k)
	K *= mpfr::sqrt(one + mpfr::pow(two, -2 * k));
      return K;
    }

    mpfr::mpreal
    k_hyperbolic()
    {
      mpfr::mpreal one(1, 4 * _NumBits);
      mpfr::mpreal two(2, 4 * _NumBits);
      mpfr::mpreal K(1, 4 * _NumBits);
      for (auto k = 1u; k < _NumBits; ++k)
	K *= mpfr::sqrt(one - mpfr::pow(two, -2 * k));
      for (auto k = 1u, i = 3*k + 1u; i < _NumBits; ++k, i = 3*i + 1)
	K *= mpfr::sqrt(one - mpfr::pow(two, -2 * i));
      return K;
    }

    void
    atan_table()
    {
      auto one = std::numeric_limits<typename int_t<_NumBits>::type>::max() >> 2;// + 1;
      std::cout.precision(std::numeric_limits<double>::digits10);
      auto w = std::cout.precision() + 6;

      mpq_t p;
      mpq_init(p);
      unsigned long int i = 1;
      for (auto k = 0u; k < _NumBits; ++k)
	{
	  mpq_set_ui(p, 1ul, i);
	  mpfr::mpreal x(p, 2 * _NumBits);
	  auto a = mpfr::atan(x);
	  std::cout << ' ' << std::setw(w) << x
		    << ' ' << std::hex << std::setw(w) << a
		    << ' ' << std::hex << std::setw(w) << one * a
		    << '\n';
	  i <<= 1;
	}
    }

    void
    atanh_table()
    {
      auto one = std::numeric_limits<typename int_t<_NumBits>::type>::max() >> 2;// + 1;
      std::cout.precision(std::numeric_limits<double>::digits10);
      auto w = std::cout.precision() + 6;

      mpq_t p;
      mpq_init(p);
      unsigned long int i = 1;
      for (auto k = 0u; k < _NumBits; ++k)
	{
	  mpq_set_ui(p, 1ul, i);
	  mpfr::mpreal x(p, 2 * _NumBits);
	  auto a = mpfr::atanh(x);
	  std::cout << ' ' << std::setw(w) << x
		    << ' ' << std::hex << std::setw(w) << a
		    << ' ' << std::hex << std::setw(w) << one * a
		    << '\n';
	  i <<= 1;
	}
    }
  };

int
main()
{
  std::cout << "\n\nCORDIC Algorithm\n\n";

  std::cout << "\ncordic<16>\n";
  cordic<16> c16;
  std::cout << "K = " << c16.k_circular() << '\n';
  c16.atan_table();
  std::cout << "K' = " << c16.k_hyperbolic() << '\n';
  c16.atanh_table();

  std::cout << "\ncordic<32>\n";
  cordic<32> c32;
  std::cout << "K = " << c32.k_circular() << '\n';
  c32.atan_table();
  std::cout << "K' = " << c32.k_hyperbolic() << '\n';
  c32.atanh_table();

  std::cout << "\ncordic<64>\n";
  cordic<64> c64;
  std::cout << "K = " << c64.k_circular() << '\n';
  c64.atan_table();
  std::cout << "K' = " << c64.k_hyperbolic() << '\n';
  c64.atanh_table();

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //std::cout << "\ncordic<128>\n";
  //cordic<128> c128;
  //std::cout << "K = " << c128.k_circular() << '\n';
  //c128.atan_table();
  //std::cout << "K' = " << c128.k_hyperbolic() << '\n';
  //c128.atanh_table();
#endif
}
