/*
$HOME/bin_tr29124/bin/g++ -std=gnu++1z -Wall -Wextra -o test_bessel test_bessel.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_bessel > test_bessel.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

#include <cmath>
//#include "new_bessel.tcc.bad"
//#include "new_bessel.tcc"

template <typename Type>
void
test_bessel( void )
{
  std::cout.precision(std::numeric_limits<Type>::digits10);
  std::cout.flags(std::ios::showpoint);
  auto width = 8 + std::cout.precision();
  constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<Type>::__pi_half;

  std::cyl_neumann(1.0, 0.01);

  std::vector<Type> nu
  {
    Type{0},
    Type{1} / Type{3},
    Type{1} / Type{2},
    Type{2} / Type{3},
    Type{1},
    Type{2},
    Type{5},
    Type{10},
    Type{20},
    Type{50},
    Type{100},
  };

  for (std::size_t n = 0; n <= nu.size(); ++n)
    {
      std::cout << "\n  nu = " << nu[n] << std::endl;
      std::cout << "  " << std::setw(width) << "x";
      std::cout << "  " << std::setw(width) << "J_nu(x)";
      std::cout << "  " << std::setw(width) << "N_nu(x)";
      std::cout << "  " << std::setw(width) << "pi x W{J,N} / 2";
      std::cout << "  " << std::setw(width) << "I_nu(x)";
      std::cout << "  " << std::setw(width) << "K_nu(x)";
      std::cout << "  " << std::setw(width) << "-x W{I,K}";
      std::cout << std::endl;
      for (unsigned int i = 0; i <= 1000; ++i)
        {
          auto x = Type(i) / Type{100};
          Type J_nu, N_nu, Jp_nu, Np_nu;
          try
            {
              std::__detail::__cyl_bessel_jn(nu[n], x, J_nu, N_nu, Jp_nu, Np_nu);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          Type I_nu, K_nu, Ip_nu, Kp_nu;
          try
            {
              std::__detail::__cyl_bessel_ik(nu[n], x, I_nu, K_nu, Ip_nu, Kp_nu);
            }
          catch (std::exception & err)
            {
              std::cerr << err.what() << std::endl;
            }
          auto W_JN = _S_pi_2 * x * (J_nu * Np_nu - N_nu * Jp_nu);
          auto W_IK = -x * (I_nu * Kp_nu - K_nu * Ip_nu);
          std::cout << "  " << std::setw(width) << x;
          std::cout << "  " << std::setw(width) << J_nu;
          std::cout << "  " << std::setw(width) << N_nu;
          std::cout << "  " << std::setw(width) << W_JN;
          std::cout << "  " << std::setw(width) << I_nu;
          std::cout << "  " << std::setw(width) << K_nu;
          std::cout << "  " << std::setw(width) << W_IK;
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

