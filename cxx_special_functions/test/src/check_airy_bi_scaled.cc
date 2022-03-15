
//  airy_bi_scaled

#include "verify.h"

// Test data.
// max(|f - f_GSL|): 0.38947881234488979 at index 4
// max(|f - f_GSL| / |f_GSL|): 0.86033214250877665
// mean(f - f_GSL): 0.30116939386853875
// variance(f - f_GSL): 3.9103589272068649e-05
// stddev(f - f_GSL): 0.0062532862777957519
const testcase_airy_bi_scaled<double>
data001[21] =
{
  { 0.61492662744600068, 0.0000000000000000, 0.0 },
  { 0.67489241111563025, 0.50000000000000000, 0.0 },
  { 0.61991194357267854, 1.0000000000000000, 0.0 },
  { 0.55209437228578417, 1.5000000000000000, 0.0 },
  { 0.50043725430409502, 2.0000000000000000, 0.0 },
  { 0.46475048019609250, 2.5000000000000000, 0.0 },
  { 0.43938402355009643, 3.0000000000000000, 0.0 },
  { 0.42017718823530570, 3.5000000000000000, 0.0 },
  { 0.40480946788929806, 4.0000000000000000, 0.0 },
  { 0.39202734094459063, 4.5000000000000000, 0.0 },
  { 0.38110853108887738, 5.0000000000000000, 0.0 },
  { 0.37160000660099485, 5.5000000000000000, 0.0 },
  { 0.36319693054542684, 6.0000000000000000, 0.0 },
  { 0.35568337591227883, 6.5000000000000000, 0.0 },
  { 0.34890049029582110, 7.0000000000000000, 0.0 },
  { 0.34272793654369132, 7.5000000000000000, 0.0 },
  { 0.33707237582041633, 8.0000000000000000, 0.0 },
  { 0.33185997650946220, 8.5000000000000000, 0.0 },
  { 0.32703135827743030, 9.0000000000000000, 0.0 },
  { 0.32253807502213006, 9.5000000000000000, 0.0 },
  { 0.31834010533673446, 10.000000000000000, 0.0 },
};
const double toler001 = 0.050000000000000003;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_airy_bi_scaled<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::airy_bi_scaled(data[i].x);
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
