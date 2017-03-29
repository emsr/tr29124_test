
#include <cmath>

int
main()
{
  // Answer from Wolfram Alpha.
  long double ans = 2.8001067714281046381729339089728201718367278248825236L;
  auto Pi_comp = std::comp_ellint_3(0.75L, -0.5L);
  auto diff = Pi_comp - ans;
}

/*
Thread 1 "a" hit Breakpoint 1, main () at pr66689.cc:7
7	  long double ans = 2.8001067714281046381729339089728201718367278248825236L;
(gdb) n
8	  auto Pi_comp = std::comp_ellint_3(0.75L, -0.5L);
(gdb) 
9	  auto diff = Pi_comp - ans;
(gdb) p Pi_comp 
$3 = 2.8001067714281046381469275630138327
(gdb) p ans 
$4 = 2.8001067714281046381469275630138327
(gdb) n
10	}
(gdb) p diff
$5 = 0

*/
