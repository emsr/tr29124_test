#include <iostream>
#include <iomanip>
#include <cmath>

#include "new_hermite.tcc"

int
main()
{

  std::cout.precision(8);
  std::cout.flags(std::ios::showpoint);

  double power = 1.0;
  for (int n = 0; n <= 50; ++n)
    {
      // sqrt(factrl(n) * 2**n sqrt(pi))
      double factor = std::exp(0.5 * ::lgamma(n + 1)) * std::sqrt(power * std::sqrt(M_PI));
      power *= 2.0;
      std::cout << "  " << std::setw(16) << "factor = " << factor << std::endl;
      std::cout << "  " << std::setw(16) << "x";
      std::cout << "  " << std::setw(16) << "H_" << n << "(x)";
      std::cout << "  " << std::setw(16) << "H~_" << n << "(x)";
      std::cout << "  " << std::setw(16) << "ratio";
      std::cout << std::endl;
      for (int i = 0; i <= 100; ++i)
        {
          double x = i * 0.1;
          double h = __hermite(n, x);
          double ht = __hermite_norm(n, x);
          std::cout << "  " << std::setw(16) << x;
          std::cout << "  " << std::setw(16) << h;
          std::cout << "  " << std::setw(16) << ht;
          std::cout << "  " << std::setw(16) << h / ht;
          std::cout << std::endl;
        }
    }
  return 0;
}
