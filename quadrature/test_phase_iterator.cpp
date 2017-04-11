
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "fourier_transform.h"

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  auto w = 4 + std::cout.precision();

  __gnu_cxx::__phase_iterator<double> iter(+0.001);
  for (int k = 0; k < 1000; ++k, ++iter)
    {
      auto phi = k * 0.001;
      std::cout << ' ' << std::setw(w) << phi
		<< ' ' << std::setw(w) << std::sin(phi)
		<< ' ' << std::setw(w) << std::cos(phi)
		<< ' ' << std::setw(w) << iter.sin()
		<< ' ' << std::setw(w) << iter.cos()
		<< ' ' << std::setw(w) << iter.sin() - std::sin(phi)
		<< ' ' << std::setw(w) << iter.cos() - std::cos(phi)
		<< '\n';
    }
}
