#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "fresnel.tcc"

int main(int, char **)
{

  std::cout.precision(8);
  std::cout.flags(std::ios::showpoint);

  std::cout << "  " << std::setw(16) << "x";
  std::cout << "  " << std::setw(16) << "C(x)";
  std::cout << "  " << std::setw(16) << "S(x)";
  std::cout << std::endl;
  for (int i = 0; i <= 1000; ++i)
    {
      double x = i * 0.01;
      std::pair<double,double> frnl = __fresnel(x);
      std::cout << "  " << std::setw(16) << x;
      std::cout << "  " << std::setw(16) << frnl.first;
      std::cout << "  " << std::setw(16) << frnl.second;
      std::cout << std::endl;
    }

  return 0;
}

