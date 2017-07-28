/*
$HOME/bin/bin/g++ -std=gnu++17 -g -I. -o test_solvers test_solvers.cpp -lquadmath
./test_solvers > test_solvers.txt
*/

#include "solver.h"

int
main()
{
  auto quad = quadratic(std::experimental::make_array(-2.0, 1.0, 1.0));

}
