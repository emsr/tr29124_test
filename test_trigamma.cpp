/*
$HOME/bin_tr29124/bin/g++ -std=gnu++14 -g -Wall -Wextra -o test_trigamma test_trigamma.cpp -L$HOME/bin/lib64 -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_trigamma > test_trigamma.txt

g++ -std=gnu++14 -Wall -Wextra -DNO_LOGBQ -I. -o test_trigamma test_trigamma.cpp -lquadmath
./test_trigamma > test_trigamma.txt
*/

#include <limits>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  _Tp
  __trigamma(_Tp __x)
  {
  }

template<typename _Tp>
  _Tp
  __tetragamma(_Tp __x)
  {
  }

int
main()
{
}
