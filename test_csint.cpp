/*
$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_csint test_csint.cpp -L$HOME/bin/lib64 -lquadmath
./test_csint > test_csint.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -DNO_LOGBQ -I. -o test_csint test_csint.cpp -lquadmath
./test_csint > test_csint.txt
*/

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <complex>

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
      auto [Ci, Si] = std::__detail::__sincosint(x);
      std::cout << "  " << std::setw(16) << x;
      std::cout << "  " << std::setw(16) << Ci;
      std::cout << "  " << std::setw(16) << Si;
      std::cout << std::endl;
    }

  return 0;
}

