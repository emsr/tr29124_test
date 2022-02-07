
//  cyl_bessel_i_scaled

#include "verify.h"

// Test data for nu=100.00000000000000.
// max(|f - f_GSL|): 2.3104348295666099e-16 at index 10
// max(|f - f_GSL| / |f_GSL|): 1.2043137544512246e-12
// mean(f - f_GSL): -3.7063697676001627e-17
// variance(f - f_GSL): 4.1390972807283784e-33
// stddev(f - f_GSL): 6.4335816468965233e-17
const testcase_cyl_bessel_i_scaled<double>
data001[11] =
{
  { 8.5155875815486127e-05, 100.00000000000000, 1000.0000000000000, 0.0 },
  { 0.00012783770668059490, 100.00000000000000, 1100.0000000000000, 0.0 },
  { 0.00017868879940604636, 100.00000000000000, 1200.0000000000000, 0.0 },
  { 0.00023648188537317519, 100.00000000000000, 1300.0000000000000, 0.0 },
  { 0.00029987385710643594, 100.00000000000000, 1400.0000000000000, 0.0 },
  { 0.00036754116875294015, 100.00000000000000, 1500.0000000000000, 0.0 },
  { 0.00043825961543828511, 100.00000000000000, 1600.0000000000000, 0.0 },
  { 0.00051094422898306715, 100.00000000000000, 1700.0000000000000, 0.0 },
  { 0.00058466302033197875, 100.00000000000000, 1800.0000000000000, 0.0 },
  { 0.00065863487030955832, 100.00000000000000, 1900.0000000000000, 0.0 },
  { 0.00073221866792298469, 100.00000000000000, 2000.0000000000000, 0.0 },
};
const double toler001 = 1.0000000000000006e-10;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_cyl_bessel_i_scaled<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::cyl_bessel_i_scaled(data[i].nu, data[i].x);
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const Ret f0 = data[i].f0;
	    const Ret diff = f - f0;
	    if (std::abs(diff) > max_abs_diff)
	      max_abs_diff = std::abs(diff);
	    if (std::abs(f0) > Ret(10) * eps
	     && std::abs(f) > Ret(10) * eps)
	      {
		const Ret frac = diff / f0;
		if (std::abs(frac) > max_abs_frac)
		  max_abs_frac = std::abs(frac);
	      }
	  }
      }
    int num_errors = 0;
    VERIFY(!failure && max_abs_frac < toler);
    return num_errors;
  }

int
main()
{
  return test(data001, toler001);
}
