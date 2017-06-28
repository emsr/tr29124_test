/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_cordic test_cordic.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_cordic > test_cordic.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_cordic test_cordic.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
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
  std::cout << '\n';
  for (int i = -0x6487ED51; i < 0x6487ED51; i += 16384)
    {
      auto sc = cordic(i, 32);
      std::cout << ' ' << i << ' ' << sc.first << ' ' << sc.second << '\n';
    }

  std::cout << '\n';
  for (short i = -0x6487; i < 0x6487; i += 128)
    {
      auto sc = cordic(i, 32);
      std::cout << ' ' << i << ' ' << sc.first << ' ' << sc.second << '\n';
    }
}
