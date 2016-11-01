/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_fresnel test_fresnel.cpp wrap_boost.cpp -lquadmath
./test_fresnel > test_fresnel.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_fresnel test_fresnel.cpp wrap_boost.cpp -lquadmath
./test_fresnel > test_fresnel.txt

g++ -std=gnu++17 -g -Wall -Wextra -DNO_LOGBQ -I. -o test_fresnel test_fresnel.cpp wrap_boost.cpp -lquadmath
./test_fresnel > test_fresnel.txt
*/

#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "fresnel.tcc"

template<typename _Tp>
  void
  test_fresnel(_Tp proto = _Tp{})
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "  " << std::setw(width) << "x";
    std::cout << "  " << std::setw(width) << "C(x)";
    std::cout << "  " << std::setw(width) << "S(x)";
    std::cout << '\n';
    for (int i = 0; i <= 1000; ++i)
      {
	auto x = i * _Tp{0.01Q};
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

