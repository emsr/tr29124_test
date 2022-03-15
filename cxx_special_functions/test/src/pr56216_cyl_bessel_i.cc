
// PR libstdc++/56216 - Crash of Bessel functions at x==0!

#include <emsr/special_functions.h>

int
test01()
{
  double j0 = emsr::cyl_bessel_j(0.0, 0.0);
  double i0 = emsr::cyl_bessel_i(0.0, 0.0);
  double j1 = emsr::cyl_bessel_j(1.0, 0.0);
  double i1 = emsr::cyl_bessel_i(1.0, 0.0);

  int num_errors = 0;
  if (j0 != 1.0)
    ++num_errors;
  if (i0 != 1.0)
    ++num_errors;
  if (j1 != 0.0)
    ++num_errors;
  if (i1 != 0.0)
    ++num_errors;
  return num_errors;
}

int
main()
{
  return test01();
}
