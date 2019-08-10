#include <iostream>
#include <iomanip>
#include <iterator>
#include <deque>
#include <algorithm>
#include <e_float/e_float.h>
#include <functions/functions.h>

int main(void)
{
  // Show e_float STL usage.
  std::deque<e_float> z;

  // Compute three zeros of a Bessel function with very high order
  // and store these in a container.
  ef::cyl_bessel_j_zero(ef::million() + ef::pi(), 3u, z);

  std::cout.precision(std::numeric_limits<e_float>::digits10);

  // Copy the zeros to the output using an algorithm.
  std::copy(z.begin(), z.end(), std::ostream_iterator<e_float>(std::cout, "\n"));
}
