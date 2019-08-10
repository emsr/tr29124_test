#include <deque>
#include <iostream>
#include <algorithm>
#include <iterator>

#include <e_float/e_float.h>
#include <examples/examples.h>
#include <functions/functions.h>
#include <utility/util_timer.h>

void examples::nr_001::basic_usage_real(void)
{
  // Print 21 values of the function legendre_p[(222/10) + k, (111/10) + k, euler_gamma]
  // for 0 <= k < 21 to the standard output using e_float. Also compute the computation
  // time for the calculation using e_float in [ms].
  // A comparable Mathematica code is:
  // Timing[Table[N[LegendreP[(222/10) + k, (111/10) + k, EulerGamma], 100],{k, 0, 20, 1}]].

  const e_float v(22.2);
  const e_float u(11.1);

  std::vector<e_float> values(21u);

  const Util::timer tm;
  for(INT32 k = 0; k < static_cast<INT32>(values.size()); k++)
  {
    values[static_cast<std::size_t>(k)] = ef::legendre_p(v + k, u + k, ef::euler_gamma());
  }
  const double elapsed = tm.elapsed();

  std::cout << "Elapsed time: " << elapsed << "\n";

  const std::streamsize         original_prec = std::cout.precision(std::numeric_limits<e_float>::digits10);
  const std::ios_base::fmtflags original_flag = std::cout.setf     (std::ios_base::showpos | std::ios_base::scientific);

  std::copy(values.begin(), values.end(), std::ostream_iterator<e_float>(std::cout, "\n"));

  std::cout.precision(original_prec);
  std::cout.unsetf(std::ios_base::showpos | std::ios_base::scientific);
  std::cout.setf(original_flag);
}
