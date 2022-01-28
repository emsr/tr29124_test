
#include <limits>
#include <emsr/sf_bessel.h>

int
test01()
{
  const double inf = std::numeric_limits<double>::infinity();
  double nm1o4 = emsr::cyl_neumann(-0.25, 0.0);
  double nm1o2 = emsr::cyl_neumann(-0.5, 0.0);
  double nm1 = emsr::cyl_neumann(-1.0, 0.0);
  double nm3o2 = emsr::cyl_neumann(-1.5, 0.0);
  double nm2 = emsr::cyl_neumann(-2.0, 0.0);

  int num_errors = 0;
  if (nm1o4 != -inf)
    ++num_errors;
  if (nm1o2 != 0.0)
    ++num_errors;
  if (nm1 != inf)
    ++num_errors;
  if (nm3o2 != 0.0)
    ++num_errors;
  if (nm2 != -inf)
    ++num_errors;
  return num_errors;
}

int
main()
{
  return test01();
}
