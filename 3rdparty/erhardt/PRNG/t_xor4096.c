/* t_xor4096.c   we Apr.2004
   test program for xor4096i in xorgens.c V3.04 by R. P. Brent, 20060628
   calculate 1st and 10000th random longint from seed 1
   
   gcc -O3 -o t_xor4096 -lm t_xor4096.c
*/

#include <stdio.h>
#include <math.h>

#include "xorgens.h"
#include "xorgens.c"

const int trials=10000;
const int seed0=1;

int main() {
  UINT i,f1,f2;
  f1 = xor4096i(seed0);
  f2 = f1;
  printf("R1 = %d\n", f1);
  for (i=1; i<trials; i++) f2 = xor4096i(0);
  printf("R%d = %d\n", trials, f2);
  return 0;
}
