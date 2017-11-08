#include <quadmath.h>
#include <assert.h>

void
test_alt_signs()
{
  __float128 result = tgammaq(-1.5Q);
  assert(result > 0.0Q);

  result = tgammaq(-2.5Q);
  assert(result < 0.0Q);

  result = tgammaq(-3.5Q);
  assert(result > 0.0Q);

  result = tgammaq(-4.5Q);
  assert(result < 0.0Q);
}

/*
 * Return |\Gamma(x) \Gamma(1 - x) - \frac{\pi}{sin(\pi x)}|
 */
__float128
abs_delta(__float128 x)
{
  return fabsq(tgammaq(x) * tgammaq(1.0Q - x) - M_PIq / sinq(M_PIq * x));
}

/*
 * Test reflection: \Gamma(x) \Gamma(1 - x) = \frac{\pi}{sin(\pi x)}
 */
void
test_reflection()
{
  assert(abs_delta(+1.5Q) < 100 * FLT128_EPSILON);
  assert(abs_delta(+0.5Q) < 100 * FLT128_EPSILON);
  assert(abs_delta(-0.5Q) < 100 * FLT128_EPSILON);
  assert(abs_delta(-1.5Q) < 100 * FLT128_EPSILON);
  assert(abs_delta(-2.5Q) < 100 * FLT128_EPSILON);
}

int
main()
{
  test_alt_signs();

  test_reflection();

  return 0;
}
