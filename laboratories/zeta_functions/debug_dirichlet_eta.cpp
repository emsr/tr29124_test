/*
$HOME/bin_specfun/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -Izeta -o debug_dirichlet_eta debug_dirichlet_eta.cpp -lquadmath -Lwrappers/debug -lwrap_cephes
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./debug_dirichlet_eta > debug_dirichlet_eta.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -Izeta -o debug_dirichlet_eta debug_dirichlet_eta.cpp -lquadmath -Lwrappers/debug -lwrap_cephes
PATH=wrappers/debug:$PATH ./debug_dirichlet_eta > debug_dirichlet_eta.txt

PATH=wrappers/debug:$PATH gdb ./debug_dirichlet_eta
br 'std::__detail::__dirichlet_eta<double>(double)' 
*/

#include <ext/cmath>
#include <iostream>

int
main()
{
  std::cout.precision(__gnu_cxx::__digits10(1.0));
  std::cout << std::scientific;

  double s = -1.0;
  std::cout << __gnu_cxx::dirichlet_eta(s) << '\n';
}
