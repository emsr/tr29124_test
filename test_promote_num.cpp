// $HOME/bin/bin/g++ -o test_promote_num test_promote_num.cpp

#include "promote_num.h"

int
main()
{
  int i;
  short j;
  float f;
  double d;
  long double ld;

  __promote_num_t<int, short> x;
  __promote_num_t<int, long double> y;
}
