
//  comp_ellint_1

#include "verify.h"

// Test data.
// max(|f - f_Boost|): 2.2204460492503131e-16 at index 1
// max(|f - f_Boost| / |f_Boost|): 1.3992633682840101e-16
// mean(f - f_Boost): -7.0119348923694093e-17
// variance(f - f_Boost): 2.8832635424744582e-34
// stddev(f - f_Boost): 1.6980175330291670e-17
const testcase_comp_ellint_1<double>
data001[19] =
{
  { 2.2805491384227703, -0.90000000000000002, 0.0 },
  { 1.9953027776647294, -0.80000000000000004, 0.0 },
  { 1.8456939983747234, -0.69999999999999996, 0.0 },
  { 1.7507538029157526, -0.59999999999999998, 0.0 },
  { 1.6857503548125961, -0.50000000000000000, 0.0 },
  { 1.6399998658645112, -0.39999999999999991, 0.0 },
  { 1.6080486199305128, -0.29999999999999993, 0.0 },
  { 1.5868678474541662, -0.19999999999999996, 0.0 },
  { 1.5747455615173560, -0.099999999999999978, 0.0 },
  { 1.5707963267948966, 0.0000000000000000, 0.0 },
  { 1.5747455615173560, 0.10000000000000009, 0.0 },
  { 1.5868678474541662, 0.20000000000000018, 0.0 },
  { 1.6080486199305128, 0.30000000000000004, 0.0 },
  { 1.6399998658645112, 0.40000000000000013, 0.0 },
  { 1.6857503548125961, 0.50000000000000000, 0.0 },
  { 1.7507538029157526, 0.60000000000000009, 0.0 },
  { 1.8456939983747238, 0.70000000000000018, 0.0 },
  { 1.9953027776647294, 0.80000000000000004, 0.0 },
  { 2.2805491384227707, 0.90000000000000013, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_comp_ellint_1<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::comp_ellint_1(data[i].k);
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
