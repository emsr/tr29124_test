#include <quadmath.h>
#include <assert.h>

int
main(void)
{
  __float128 result = tgammaq(-1.5Q);
  assert(result > 0.0Q);

  result = tgammaq(-2.5Q);
  assert(result < 0.0Q);

  return 0;
}
