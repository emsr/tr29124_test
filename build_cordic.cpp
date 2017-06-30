/*
$HOME/bin_tr29124/bin/g++ -I. -o build_cordic build_cordic.cpp -lquadmath -lmpfr -lgmp
./build_cordic > build_cordic.txt

$HOME/bin/bin/g++ -std=gnu++14 -I. -o build_cordic build_cordic.cpp -lquadmath -lmpfr -lgmp
./build_cordic > build_cordic.txt
*/

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include <bits/float128_io.h>
#include <mpreal.h>

template<std::size_t _NumBits>
  struct cordic
  {
    void operator()()
    {
      mpq_t p;
      mpq_init(p);
      unsigned long int i = 1;
      for (int k = 0; k < _NumBits; ++k)
	{
	  mpq_set_ui(p, 1ul, i);
	  mpfr::mpreal x(p, _NumBits);
	  auto a = mpfr::atan(x);
	  std::cout << ' ' << x << ' ' << std::hex << a << '\n';
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
  c16();
  std::cout << "\ncordic<32>\n";
  cordic<32> c32;
  c32();
  std::cout << "\ncordic<64>\n";
  cordic<64> c64;
  c64();
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\ncordic<128>\n";
  cordic<64>();
#endif
}
