/* Program to generate selftest values for taus113 prng from
**
**  [1] P. L'Ecuyer, "Tables of Maximally-Equidistributed Combined LFSR
**      Generators", Mathematics of Computation, 68, 225 (1999), 261-269.
**  [2] http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps
**      Online version of [1]
**
** Copied from [2] Fig.1, modified for unsigned long, uses minimum seeds
** Compile with: gcc -Wall t_taus113.c -o t_taus113.exe
** W.Ehrhardt Aug.2005
**
*/

#include <stdio.h>

static unsigned long z1=2, z2=8, z3=16, z4=128;


unsigned long lsfr113() {
   unsigned long b;
   b  = (((z1 <<  6) ^ z1) >> 13);
   z1 = (((z1 & 4294967294U) << 18) ^ b);
   b  = (((z2 <<  2) ^ z2) >> 27);
   z2 = (((z2 & 4294967288U) << 2) ^ b);
   b  = (((z3 << 13) ^ z3) >> 21);
   z3 = (((z3 & 4294967280U) << 7) ^ b);
   b  = (((z4 <<  3) ^ z4) >> 12);
   z4 = (((z4 & 4294967168U) << 13) ^ b);
   return (z1 ^ z2 ^ z3 ^ z4);
}

int main() {
  int j=0;
  printf("First: %08lx\n", lsfr113());
  for (++j;   j<499; j++) {lsfr113();}  printf("  500: %08lx\n", lsfr113());
  for (++j;   j<999; j++) {lsfr113();}  printf(" 1000: %08lx\n", lsfr113());
  for (++j;  j<2499; j++) {lsfr113();}  printf(" 2500: %08lx\n", lsfr113());
  for (++j;  j<9999; j++) {lsfr113();}  printf("10000: %08lx\n", lsfr113());
  for (++j; j<19999; j++) {lsfr113();}  printf("20000: %08lx\n", lsfr113());
  for (++j; j<29999; j++) {lsfr113();}  printf("30000: %08lx\n", lsfr113());
  return(0);
}

