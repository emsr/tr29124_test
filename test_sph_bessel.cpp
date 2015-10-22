#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>

//#include "new_bessel.tcc.bad"
#include "new_bessel.tcc"

int main(int, char **)
{

  std::cout.precision(8);
  std::cout.flags(std::ios::showpoint);

  for (int n = 0; n <= 50; ++n)
    {
      std::cout << "  " << std::setw(16) << "x";
      std::cout << "  " << std::setw(16) << "j_" << n << "(x)";
      std::cout << "  " << std::setw(16) << "n_" << n << "(x)";
      std::cout << "  " << std::setw(16) << "i_" << n << "(x)";
      std::cout << "  " << std::setw(16) << "k_" << n << "(x)";
      std::cout << std::endl;
      for (int i = 0; i <= 100; ++i)
        {
          double x = i * 0.1;
          double j_n, n_n, jp_n, np_n;
          try
            {
              __sph_bessel_jn(n, x, j_n, n_n, jp_n, np_n);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          double i_n, k_n, ip_n, kp_n;
          try
            {
              __sph_bessel_ik(n, x, i_n, k_n, ip_n, kp_n);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          std::cout << "  " << std::setw(16) << x;
          std::cout << "  " << std::setw(16) << j_n;
          std::cout << "  " << std::setw(16) << n_n;
          std::cout << "  " << std::setw(16) << i_n;
          std::cout << "  " << std::setw(16) << k_n;
          std::cout << std::endl;
        }
    }
  return 0;
}

