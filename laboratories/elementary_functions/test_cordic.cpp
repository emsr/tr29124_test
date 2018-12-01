/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_cordic test_cordic.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_cordic > test_cordic.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_cordic test_cordic.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
PATH=wrappers/debug:$PATH ./test_cordic > test_cordic.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "cordic16.h"
#include "cordic32.h"

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  auto w = std::cout.precision() + 6;

  std::cout << '\n';
  for (int i = -0x6487ED51; i < 0x6487ED51; i += 16384)
    {
      auto sc = sincos_cordic(32, i);
      std::cout << ' ' << std::setw(2) << i
		<< ' ' << std::setw(9) << sc.first
		<< ' ' << std::setw(9) << sc.second
		<< ' ' << std::setw(w) << i / double(cordic32_one)
		<< ' ' << std::setw(w) << sc.first / double(cordic32_one)
		<< ' ' << std::setw(w) << sc.second / double(cordic32_one)
		<< '\n';
    }

  std::cout << '\n';
  short x0 = 1;
  for (short i = -0x6487; i < 0x6487; i += 128)
    {
      short y0 = 0x0001 << i;
      auto at2 = atan2_cordic(16, y0, x0);
      std::cout << ' ' << std::setw(2) << i
		<< ' ' << std::setw(5) << at2
		<< ' ' << std::setw(w) << at2 / double(cordic16_one)
		<< '\n';
    }

  std::cout << '\n';
  for (short i = -0x6487; i < 0x6487; i += 128)
    {
      auto sc = sincos_cordic(16, i);
      std::cout << ' ' << std::setw(2) << i
		<< ' ' << std::setw(5) << sc.first
		<< ' ' << std::setw(5) << sc.second
		<< ' ' << std::setw(w) << i / double(cordic16_one)
		<< ' ' << std::setw(w) << sc.first / double(cordic16_one)
		<< ' ' << std::setw(w) << sc.second / double(cordic16_one)
		<< '\n';
    }
}
