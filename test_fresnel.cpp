/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_fresnel test_fresnel.cpp -lquadmath -Lwrappers/debug -lwrap_boost
./test_fresnel > test_fresnel.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_fresnel test_fresnel.cpp -lquadmath -Lwrappers/debug -lwrap_boost
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_fresnel > test_fresnel.txt
*/

#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <bits/numeric_limits.h>

#include "fresnel.tcc"

template<typename _Tp>
  void
  test_fresnel(_Tp proto = _Tp{})
  {
    //using _Val = _Tp;
    //using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "  " << std::setw(width) << "x";
    std::cout << "  " << std::setw(width) << "C(x)";
    std::cout << "  " << std::setw(width) << "S(x)";
    std::cout << '\n';
    const auto del = _Tp{1} / _Tp{100};
    for (int i = 0; i <= 1000; ++i)
      {
	auto x = i * del;
	auto frnl = __fresnel(x);
	std::cout << "  " << std::setw(width) << x;
	std::cout << "  " << std::setw(width) << frnl.first;
	std::cout << "  " << std::setw(width) << frnl.second;
	std::cout << '\n';
      }
  }

int
main()
{
  test_fresnel(1.0);

  return 0;
}

