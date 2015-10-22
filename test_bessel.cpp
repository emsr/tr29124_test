#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

//#include "new_bessel.tcc.bad"
#include "new_bessel.tcc"

template <typename Type>
void
test_bessel( void )
{

  std::cout.precision(8);
  std::cout.flags(std::ios::showpoint);

  std::vector<Type> nu;
  nu.push_back(Type(0));
  nu.push_back(Type(1)/Type(3));
  nu.push_back(Type(1)/Type(2));
  nu.push_back(Type(2)/Type(3));
  nu.push_back(Type(1));
  nu.push_back(Type(2));
  nu.push_back(Type(5));
  nu.push_back(Type(10));
  nu.push_back(Type(20));
  nu.push_back(Type(50));
  nu.push_back(Type(100));

  for (std::size_t n = 0; n <= nu.size(); ++n)
    {
      std::cout << "  nu = " << nu[n] << std::endl;
      std::cout << "  " << std::setw(16) << "x";
      std::cout << "  " << std::setw(16) << "J_nu(x)";
      std::cout << "  " << std::setw(16) << "N_nu(x)";
      std::cout << "  " << std::setw(16) << "I_nu(x)";
      std::cout << "  " << std::setw(16) << "K_nu(x)";
      std::cout << std::endl;
      for (unsigned int i = 0; i <= 100; ++i)
        {
          Type x = i / Type(10);
          Type J_nu, N_nu, Jp_nu, Np_nu;
          try
            {
              __bessel_jn(nu[n], x, J_nu, N_nu, Jp_nu, Np_nu);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          Type I_nu, K_nu, Ip_nu, Kp_nu;
          try
            {
              __bessel_ik(nu[n], x, I_nu, K_nu, Ip_nu, Kp_nu);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          std::cout << "  " << std::setw(16) << x;
          std::cout << "  " << std::setw(16) << J_nu;
          std::cout << "  " << std::setw(16) << N_nu;
          std::cout << "  " << std::setw(16) << I_nu;
          std::cout << "  " << std::setw(16) << K_nu;
          std::cout << std::endl;
        }
    }

  return;
}


///
///
///
int main(int, char **)
{

  test_bessel<double>();

  return 0;
}

