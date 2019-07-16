/* Program to generate selftest values for taus88 prng from
**
** [1] P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe Generators",
**     Mathematics of Computation 65, 213 (1996), 203-213
** [2] http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps
**     Online version of [1] with corrections
**
** Copied from [2] Fig.1, modified for unsigned long, uses minimum seeds
** Compile with: gcc -Wall t_taus88.c -o t_taus88.exe
** W.Ehrhardt May 2005
**
*/

#include <stdio.h>

static unsigned long s1=2, s2=8, s3=16;

unsigned long taus88() {
  unsigned long b;
  b = (((s1 << 13) ^ s1) >> 19); s1 = (((s1 & 4294967294U) << 12) ^ b);
  b = (((s2 <<  2) ^ s2) >> 25); s2 = (((s2 & 4294967288U) <<  4) ^ b);
  b = (((s3 <<  3) ^ s3) >> 11); s3 = (((s3 & 4294967280U) << 17) ^ b);
  return (s1 ^ s2 ^ s3) ;
}

int main() {
  int j;
  for (j=0;   j<499; j++) {taus88();}  printf("  500: %08lx\n", taus88());
  for (++j;   j<999; j++) {taus88();}  printf(" 1000: %08lx\n", taus88());
  for (++j;  j<2499; j++) {taus88();}  printf(" 2500: %08lx\n", taus88());
  for (++j;  j<9999; j++) {taus88();}  printf("10000: %08lx\n", taus88());
  for (++j; j<19999; j++) {taus88();}  printf("20000: %08lx\n", taus88());
  for (++j; j<29999; j++) {taus88();}  printf("30000: %08lx\n", taus88());
  return(0);
}

