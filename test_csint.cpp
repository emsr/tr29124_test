#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include <cmath>

int
main()
{
  std::cout.precision(8);
  std::cout.flags(std::ios::showpoint);

  std::cout << "  " << std::setw(16) << "x";
  std::cout << "  " << std::setw(16) << "Ci(x)";
  std::cout << "  " << std::setw(16) << "Si(x)";
  std::cout << std::endl;
  for (int i = 0; i <= 1000; ++i)
    {
      double x = i * 0.01;
      std::pair<double, double> cisi = __csint(x);
      std::cout << "  " << std::setw(16) << x;
      std::cout << "  " << std::setw(16) << cisi.first;
      std::cout << "  " << std::setw(16) << cisi.second;
      std::cout << std::endl;
    }

  return 0;
}

