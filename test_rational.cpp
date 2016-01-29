// $HOME/bin/bin/g++ -o test_rational test_rational.cpp

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_rational

#include "rational.h"
#include "polynomial.h"

int
main()
{
  __gnu_cxx::_Polynomial<__gnu_cxx::_Rational<int>> rpoly{{1,2}, {3, 4}, {5, 6}, {7, 8}};
  std::cout << rpoly << '\n';
}
