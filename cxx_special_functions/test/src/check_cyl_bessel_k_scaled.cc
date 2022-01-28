
//  cyl_bessel_k_scaled

#include "verify.h"

// Test data for nu=100.00000000000000.
// max(|f - f_GSL|): 6.6613381477509392e-16 at index 4
// max(|f - f_GSL| / |f_GSL|): 5.7960745000079599e-16
// mean(f - f_GSL): 8.0743492700011387e-17
// variance(f - f_GSL): 7.1714627747364710e-34
// stddev(f - f_GSL): 2.6779586954873801e-17
const testcase_cyl_bessel_k_scaled<double>
data001[11] =
{
  { 5.8424465521265159, 100.00000000000000, 1000.0000000000000, 0.0 },
  { 3.5410426710741936, 100.00000000000000, 1100.0000000000000, 0.0 },
  { 2.3237462838842249, 100.00000000000000, 1200.0000000000000, 0.0 },
  { 1.6216147846822155, 100.00000000000000, 1300.0000000000000, 0.0 },
  { 1.1879504112001300, 100.00000000000000, 1400.0000000000000, 0.0 },
  { 0.90491922756991028, 100.00000000000000, 1500.0000000000000, 0.0 },
  { 0.71165910479993821, 100.00000000000000, 1600.0000000000000, 0.0 },
  { 0.57464221239235203, 100.00000000000000, 1700.0000000000000, 0.0 },
  { 0.47437600624377624, 100.00000000000000, 1800.0000000000000, 0.0 },
  { 0.39899827109763819, 100.00000000000000, 1900.0000000000000, 0.0 },
  { 0.34100208493029188, 100.00000000000000, 2000.0000000000000, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_cyl_bessel_k_scaled<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::cyl_bessel_k_scaled(data[i].nu, data[i].x);
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
    VERIFY(max_abs_frac < toler);
    return num_errors;
  }

int
main()
{
  return test(data001, toler001);
}
