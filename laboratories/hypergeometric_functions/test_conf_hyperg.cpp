/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include <emsr/numeric_limits.h>
#include <emsr/sf_hyperg.h>

template<typename Tp>
  void
  test_conf_hyperg(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto a = Tp{6} / Tp{5};
    auto c = Tp{1} / Tp{5};
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i < +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << emsr::conf_hyperg(a, c, z)
		<< '\n';
    }
  }

int
main()
{
}
