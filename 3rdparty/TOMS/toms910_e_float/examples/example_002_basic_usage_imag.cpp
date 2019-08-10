#include <deque>
#include <fstream>
#include <algorithm>
#include <iterator>

#include <e_float/e_float.h>
#include <examples/examples.h>
#include <functions/functions.h>
#include <utility/util_timer.h>

void examples::nr_002::basic_usage_imag(void)
{
  // Print 21 values of the function riemann_zeta[(1/2) + (((1234/10) + k) I)]
  // for 0 <= k < 21 to the standard output using e_float. Also compute the
  // computation time for the calculation using e_float in [ms].
  // A comparable Mathematica code is:
  // Timing[Table[N[Zeta[(1/2) + (((1234/10) + k) I)], 100],{k, 0, 20, 1}]].

  static const e_float y(123.4);

  std::vector<ef_complex> values(21u);

  const Util::timer tm;
  for(INT32 k = 0; k < static_cast<INT32>(values.size()); k++)
  {
    const ef_complex z(ef::half(), y + k);

    values[static_cast<std::size_t>(k)] = efz::riemann_zeta(z);
  }
  const double elapsed = tm.elapsed();

  std::cout << "Elapsed time: " << elapsed << "\n";

  const std::streamsize         original_prec = std::cout.precision(std::numeric_limits<e_float>::digits10);
  const std::ios_base::fmtflags original_flag = std::cout.setf     (std::ios_base::showpos | std::ios_base::scientific);

  std::copy(values.begin(), values.end(), std::ostream_iterator<ef_complex>(std::cout, "\n"));

  std::cout.precision(original_prec);
  std::cout.unsetf(std::ios_base::showpos | std::ios_base::scientific);
  std::cout.setf(original_flag);
}
