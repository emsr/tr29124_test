#include <iomanip>
#include <iterator>
#include <fstream>
#include <vector>
#include <e_float/e_float.h>
#include <functions/functions.h>
#include <examples/examples.h>
#include <test/spot/test_spot.h>
#include <utility/util_timer.h>

void test::spot::test_spot(void)
{
  // Calculate some real function values...
  static const e_float x = ef::pi() / 7;
  static const e_float s = ef::sin(x);
  static const e_float g = ef::gamma(x);
  static const e_float b = ef::cyl_bessel_j(ef::third(), x);

  // ...set the output precision...
  std::cout.precision(std::numeric_limits<e_float>::digits10);

  // ...and print the values.
  std::cout << s << std::endl << g << std::endl << b << std::endl;
}
