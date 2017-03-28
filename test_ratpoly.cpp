/*
$HOME/bin_tr29124/bin/g++ -g -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_ratpoly test_ratpoly.cpp -lquadmath
./test_ratpoly > test_ratpoly.txt

$HOME/bin/bin/g++ -g -std=gnu++17 -g -Wall -Wextra -I. -o test_ratpoly test_ratpoly.cpp -lquadmath
./test_ratpoly > test_ratpoly.txt
*/

#include "ratpoly.h"
#include <iostream>
#include <complex>
#include <sstream>

int
main()
{
  __gnu_cxx::_Polynomial<double> P({0.0, 1.0, 2.0, 3.0});
  __gnu_cxx::_Polynomial<double> Q({2.0, 1.0});
  __gnu_cxx::_RationalPolynomial<double> R(P, Q);

  std::cout << "R = " << R << '\n';

  std::cout << "P = " << R.numer() << '\n';
  std::cout << "+P = " << +R.numer() << '\n';
  std::cout << "-P = " << -R.numer() << '\n';
  std::cout << "P = " << R.numer() << '\n';
  std::cout << "degree(P) = " << R.numer().degree() << '\n';

  std::cout << "Q = " << R.denom() << '\n';
  std::cout << "+Q = " << +R.denom() << '\n';
  std::cout << "-Q = " << -R.denom() << '\n';
  std::cout << "Q = " << R.denom() << '\n';
  std::cout << "degree(Q) = " << R.denom().degree() << '\n';

}

