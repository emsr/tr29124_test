#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>

#include "nric_bessel.tcc"

int main(int, char **)
{

  std::cout.precision(8);
  std::cout.flags(std::ios::showpoint);

  for (int n = 0; n <= 50; ++n)
    {
      double nu = n * 1.0;
      std::cout << "  " << std::setw(16) << "x";
      std::cout << "  " << std::setw(16) << "J_" << n << "(x)";
      std::cout << "  " << std::setw(16) << "N_" << n << "(x)";
      std::cout << "  " << std::setw(16) << "I_" << n << "(x)";
      std::cout << "  " << std::setw(16) << "K_" << n << "(x)";
      std::cout << std::endl;
      for (int i = 0; i <= 1000; ++i)
        {
          double x = i * 0.1;
          double J_nu, N_nu, Jp_nu, Np_nu;
          try
            {
              __bessel_jn(nu, x, J_nu, N_nu, Jp_nu, Np_nu);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          double I_nu, K_nu, Ip_nu, Kp_nu;
          try
            {
              __bessel_ik(nu, x, I_nu, K_nu, Ip_nu, Kp_nu);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          if (i == 0)
            {
              std::cout << "  Jp_nu(0) = " << std::setw(16) << Jp_nu;
              std::cout << "  Np_nu(0) = " << std::setw(16) << Np_nu;
              std::cout << "  Ip_nu(0) = " << std::setw(16) << Ip_nu;
              std::cout << "  Kp_nu(0) = " << std::setw(16) << Kp_nu;
              std::cout << std::endl;
            }
          std::cout << "  " << std::setw(16) << x;
          std::cout << "  " << std::setw(16) << J_nu;
          std::cout << "  " << std::setw(16) << N_nu;
          std::cout << "  " << std::setw(16) << I_nu;
          std::cout << "  " << std::setw(16) << K_nu;
          std::cout << std::endl;
        }
    }
  return 0;
}

