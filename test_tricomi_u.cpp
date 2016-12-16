/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_tricomi_u test_tricomi_u.cpp -lquadmath -L. -lwgsl -lburkhardt
./test_tricomi_u > test_tricomi_u.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_tricomi_u test_tricomi_u.cpp -lquadmath -L. -lwgsl -lburkhardt
./test_tricomi_u > test_tricomi_u.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  void
  test_tricomi_u(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto a = _Tp{1.2Q};
    auto c = _Tp{0.2Q};
    for (int i = -200; i < +200; ++i)
    {
      auto z = _Tp{0.1Q} * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << __gnu_cxx::tricomi_u(a, c, z)
		<< '\n';
    }
  }

int
main()
{
}
