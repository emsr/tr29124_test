
#include <limits>
#include "wrap_lerch.h"
#include "lerchphi/Source/lerchphi.h"

namespace lurch
{

double
lerch(double z, double s, double a)
{
  const double acc = 10.0 * std::numeric_limits<double>::epsilon();
  double phi = 0.0;
  int iter = 0;
  auto ok = lerchphi(&z, &s, &a, &acc, &phi, &iter);

  return phi;
}

} // namespace lurch
