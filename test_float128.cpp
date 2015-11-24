
// $HOME/bin_specfun/bin/g++ -o test_float128 test_float128.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_float128

#include "float128.h"

int
main()
{
  std::cout.precision(std::numeric_limits<__float128>::max_digits10);

  auto x = 0.00004472229441850588228136889483397204368247Q;

  std::cout << x << '\n';
}
