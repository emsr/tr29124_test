/**
 *
 */

#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <emsr/numeric_limits.h>

#include <fresnel.tcc>

template<typename Tp>
  void
  test_fresnel(Tp proto = Tp{})
  {
    //using _Val = Tp;
    //using _Real = emsr::num_traits_t<_Val>;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "  " << std::setw(width) << "x";
    std::cout << "  " << std::setw(width) << "C(x)";
    std::cout << "  " << std::setw(width) << "S(x)";
    std::cout << '\n';
    const auto del = Tp{1} / Tp{100};
    for (int i = 0; i <= 1000; ++i)
      {
	auto x = i * del;
	auto frnl = fresnel(x);
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

