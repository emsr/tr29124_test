
#include <limits>
#include <emsr/sf_bessel.h>

int
test01()
{
  const double inf = std::numeric_limits<double>::infinity();
  double jm1o2 = std::cyl_bessel_j(-0.5, 0.0);
  double jm1 = std::cyl_bessel_j(-1.0, 0.0);
  double jm3o2 = std::cyl_bessel_j(-1.5, 0.0);
  double jm2 = std::cyl_bessel_j(-2.0, 0.0);

  int num_errors = 0;
  if (jm1o2 != inf)
    ++num_errors;
  if (jm1 != 0.0)
    ++num_errors;
  if (jm3o2 != -inf)
    ++num_errors;
  if (jm2 != 0.0)
    ++num_errors;
  return num_errors;
}

int
main()
{
  return test01();
}
