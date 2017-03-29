
#include <cmath>

int
main()
{
  // Answers from Wolfram Alpha.
  long double ans_ok = -0.10001943365331651406888645149537315243646135979573L;
  long double ans_bomb = -0.10777727809650077516264612749163100483995270163783L;
  auto Ei_ok = std::expint(-1.500001L);
  auto diff_ok = Ei_ok - ans_ok;
  auto Ei_bomb = std::expint(-1.450001L);
  auto diff_bomb = Ei_bomb - ans_bomb;
}

/*
Thread 1 "pr68397" hit Breakpoint 1, main () at pr68397.cc:8
8	  long double ans_ok = -0.10001943365331651406888645149537315243646135979573L;
(gdb) n
9	  long double ans_bomb = -0.10777727809650077516264612749163100483995270163783L;
(gdb) 
10	  auto Ei_ok = std::expint(-1.500001L);
(gdb) 
11	  auto diff_ok = Ei_ok - ans_ok;
(gdb) 
12	  auto Ei_bomb = std::expint(-1.450001L);
(gdb) 
13	  auto diff_bomb = Ei_bomb - ans_bomb;
(gdb) 
14	}
(gdb) p diff_bomb 
$1 = -1.7618285302889447052621108014136553e-19
(gdb) p diff_ok 
$2 = -9.48676900924816379756521200761199e-20
*/
