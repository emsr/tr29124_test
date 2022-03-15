
//  airy_ai_scaled

#include "verify.h"

// Test data.
// max(|f - f_GSL|): 1.3877787807814457e-15 at index 4
// max(|f - f_GSL| / |f_GSL|): 6.0294974095362732e-15
// mean(f - f_GSL): -5.8154539385127250e-17
// variance(f - f_GSL): 1.7755239868255660e-34
// stddev(f - f_GSL): 1.3324878936881814e-17
const testcase_airy_ai_scaled<double>
data001[21] =
{
  { 0.35502805388781722, 0.0000000000000000, 0.0 },
  { 0.29327715912994734, 0.50000000000000000, 0.0 },
  { 0.26351364474914007, 1.0000000000000000, 0.0 },
  { 0.24418489767140847, 1.5000000000000000, 0.0 },
  { 0.23016491865251160, 2.0000000000000000, 0.0 },
  { 0.21932220512871203, 2.5000000000000000, 0.0 },
  { 0.21057204278597699, 3.0000000000000000, 0.0 },
  { 0.20329208081635175, 3.5000000000000000, 0.0 },
  { 0.19709480264306647, 4.0000000000000000, 0.0 },
  { 0.19172396872398537, 4.5000000000000000, 0.0 },
  { 0.18700211893594343, 5.0000000000000000, 0.0 },
  { 0.18280173946240377, 5.5000000000000000, 0.0 },
  { 0.17902840741321008, 6.0000000000000000, 0.0 },
  { 0.17561043019266195, 6.5000000000000000, 0.0 },
  { 0.17249220797740278, 7.0000000000000000, 0.0 },
  { 0.16962983096364936, 7.5000000000000000, 0.0 },
  { 0.16698807106393279, 8.0000000000000000, 0.0 },
  { 0.16453827306790608, 8.5000000000000000, 0.0 },
  { 0.16225684290423317, 9.0000000000000000, 0.0 },
  { 0.16012414238108222, 9.5000000000000000, 0.0 },
  { 0.15812366685434615, 10.000000000000000, 0.0 },
};
const double toler001 = 5.0000000000000039e-13;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_airy_ai_scaled<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::airy_ai_scaled(data[i].x);
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
