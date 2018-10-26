/*
  g++ -o mpfr_bernoulli mpfr_bernoulli.c -lmpfr -lgmp
  ./mpfr_bernoulli
*/

#include "mpfr.h"
#include "mpfr.h"

// \zeta(-s) = -B_{s+1} / (s + 1)
// s+1 == 2n
// -2n \zeta(1 - 2n) = B_{2n}
int
main()
{
  mpz_t B;  /* variable B declared as initialized */
  &B = mpfr_bernoulli_internal((mpz_t *)0, 0);
  for (int n = 1; n < 100; ++n)
    {
      mpfr_zeta();
      &B = mpfr_bernoulli_internal(&B, n);
      mpfr_printf("%f", B);
    }
}
