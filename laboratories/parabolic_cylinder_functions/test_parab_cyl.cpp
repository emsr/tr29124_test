/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <sf_parab_cyl.tcc>

template<typename Tp>
  void
  test_parab_cyl(Tp proto = Tp{})
  {
    std::cout.precision(emsr::max_digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    // Tp a = 0.5; // blows up a factor!
    Tp a = 1.25;
    for (int i = 0; i <= +200; ++i)
      {
	auto z = Tp{0.1L} * i;
	auto [U, V] = emsr::detail::parabolic_cylinder(a, z);
	std::cout << ' ' << std::setw(width) << z
		  << ' ' << std::setw(width) << U
		  << ' ' << std::setw(width) << V;
	std::cout << '\n';
      }
  }

int
main()
{
  test_parab_cyl(1.0);
}
