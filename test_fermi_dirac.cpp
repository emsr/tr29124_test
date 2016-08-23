/*
$HOME/bin_tr29124/bin/g++ -std=gnu++1z -o test_fermi_dirac test_fermi_dirac.cpp wrap_gsl.cpp -lgsl -lgslcblas
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_fermi_dirac
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>


template<typename _Tp>
  void
  run_fermi_dirac()
  {
  }


int
main()
{
  std::cout << "\nfloat\n=====\n\n";
  run_fermi_dirac<float>();

  std::cout << "\ndouble\n======\n";
  run_fermi_dirac<double>();

  std::cout << "\nlong double\n===========\n";
  run_fermi_dirac<long double>();

  std::cout << "\n__float128\n==========\n";
  //run_fermi_dirac<__float128>();
}
