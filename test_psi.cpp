// $HOME/bin/bin/g++ -o test_psi test_psi.cpp

#include <cmath>
#include <iostream>
#include <iomanip>
#include "psi.tcc"
#include "sf_gamma.tcc"

int
main()
{

  std::cout.precision(16);
  std::cout.flags(std::ios::showpoint);

//  for (unsigned int i = 1; i <= 100; ++i)
//    {
//      double x = i * 1.0;
//      std::cout << "psi(0," << std::setw(5) << x << ") = " << std::__detail::__psi(0,x);
//      std::cout << "psi(1," << std::setw(5) << x << ") = " << std::__detail::__psi(1,x);
//      std::cout << std::endl;
//    }

  const unsigned int max = 4001;
  for (unsigned int i = 1; i < max; ++i)
    {
      double x = -200.0 + i * 0.1;
      //double x = i * 0.1;
      double y_tr1 = std::__detail::__psi(x);
      double y_ser = std::__detail::__psi_series(x);
      //double y_asym = std::__detail::__psi_asymp(x);
      std::cout << "psi(" << std::setw(5) << x << ") = "
                << std::setw(20) << y_tr1
                << std::setw(20) << y_ser
                //<< std::setw(20) << y_asym
                << std::endl;
    }

  return 0;
}


