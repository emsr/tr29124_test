/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_legendre test_legendre.cpp -lquadmath
./test_legendre > test_legendre.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_legendre test_legendre.cpp -lquadmath
./test_legendre > test_legendre.txt
*/

#include <iostream>
#include <iomanip>
#include <tr1/cmath>
#include <cmath>
#include <stdexcept>
#include <limits>

#include "legendre.tcc"
#include "bits/sf_legendre.tcc"

int main(int, char **)
{

  std::cout.precision(8);
  std::cout.flags(std::ios::showpoint);

  for (int l = 0; l <= 50; ++l)
    {
      std::cout << "  " << std::setw(16) << "x";
      std::cout << "  " << std::setw(16) << "P_" << l << "(x)";
      std::cout << std::endl;
      for (int i = -100; i <= 100; ++i)
        {
          double x = i * 0.01;
          double P, P_l, P_l0;
          try
            {
              P = __legendre_p(l, x);
              P_l = std::tr1::__detail::__poly_legendre_p(l, x);
              P_l0 = std::tr1::__detail::__assoc_legendre_p(l, 0, x);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          std::cout << "  " << std::setw(16) << x;
          std::cout << "  " << std::setw(16) << P;
          std::cout << "  " << std::setw(16) << P_l;
          std::cout << "  " << std::setw(16) << P_l0;
          std::cout << "  " << std::setw(16) << P_l - P;
          std::cout << std::endl;
        }
    }

  for (unsigned int l = 0; l <= 3; ++l)
    {
      for (unsigned int m = l; m <= l; --m)
        {
          std::cout << "  " << std::setw(16) << "x";
          std::cout << "  " << std::setw(16) << "P_" << l << "_" << m << "(x)";
          std::cout << std::endl;
          for (int i = -100; i <= 100; ++i)
            {
              double x = i * 0.01;
              double P_lm;
              try
                {
                  P_lm = std::tr1::__detail::__assoc_legendre_p(l, m, x);

                }
              catch (std::exception & err)
                {
                  std::cerr << err.what() << std::endl;
                }
              std::cout << "  " << std::setw(16) << x;
              std::cout << "  " << std::setw(16) << P_lm;
              std::cout << std::endl;
            }
        }
    }

  return 0;
}

