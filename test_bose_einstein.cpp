/*
$HOME/bin_tr29124/bin/g++ -std=gnu++1z -o test_bose_einstein test_bose_einstein.cpp wrap_gsl.cpp -lgsl -lgslcblas
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_bose_einstein
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>


template<typename _Tp>
  void
  run_bose_einstein()
  {
  }


int
main()
{
  std::cout << "\nfloat\n=====\n\n";
  run_bose_einstein<float>();

  std::cout << "\ndouble\n======\n";
  run_bose_einstein<double>();

  std::cout << "\nlong double\n===========\n";
  run_bose_einstein<long double>();

  std::cout << "\n__float128\n==========\n";
  //run_bose_einstein<__float128>();
}
