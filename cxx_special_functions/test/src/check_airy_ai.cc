
//  airy_ai

#include "verify.h"

// Test data.
// max(|f - f_GSL|): 1.5126788710517758e-15 at index 4
// max(|f - f_GSL| / |f_GSL|): 5.7527019513802774e-14
// mean(f - f_GSL): -5.3665833949484320e-17
// variance(f - f_GSL): 7.3800556920774120e-35
// stddev(f - f_GSL): 8.5907250520997424e-18
const testcase_airy_ai<double>
data001[41] =
{
  { 0.040241238486444071, -10.000000000000000, 0.0 },
  { 0.31910324771912801, -9.5000000000000000, 0.0 },
  { -0.022133721547341240, -9.0000000000000000, 0.0 },
  { -0.33029023763020882, -8.5000000000000000, 0.0 },
  { -0.052705050356385910, -8.0000000000000000, 0.0 },
  { 0.32177571638064789, -7.5000000000000000, 0.0 },
  { 0.18428083525050609, -7.0000000000000000, 0.0 },
  { -0.23802030199711663, -6.5000000000000000, 0.0 },
  { -0.32914517362982321, -6.0000000000000000, 0.0 },
  { 0.017781541276574383, -5.5000000000000000, 0.0 },
  { 0.35076100902411411, -5.0000000000000000, 0.0 },
  { 0.29215278105595921, -4.5000000000000000, 0.0 },
  { -0.070265532949289680, -4.0000000000000000, 0.0 },
  { -0.37553382314043182, -3.5000000000000000, 0.0 },
  { -0.37881429367765823, -3.0000000000000000, 0.0 },
  { -0.11232506769296607, -2.5000000000000000, 0.0 },
  { 0.22740742820168561, -2.0000000000000000, 0.0 },
  { 0.46425657774886947, -1.5000000000000000, 0.0 },
  { 0.53556088329235207, -1.0000000000000000, 0.0 },
  { 0.47572809161053958, -0.50000000000000000, 0.0 },
  { 0.35502805388781722, 0.0000000000000000, 0.0 },
  { 0.23169360648083348, 0.50000000000000000, 0.0 },
  { 0.13529241631288141, 1.0000000000000000, 0.0 },
  { 0.071749497008105428, 1.5000000000000000, 0.0 },
  { 0.034924130423274372, 2.0000000000000000, 0.0 },
  { 0.015725923380470481, 2.5000000000000000, 0.0 },
  { 0.0065911393574607175, 3.0000000000000000, 0.0 },
  { 0.0025840987869896349, 3.5000000000000000, 0.0 },
  { 0.00095156385120480195, 4.0000000000000000, 0.0 },
  { 0.00033025032351430934, 4.5000000000000000, 0.0 },
  { 0.00010834442813607434, 5.0000000000000000, 0.0 },
  { 3.3685311908599812e-05, 5.5000000000000000, 0.0 },
  { 9.9476943602528973e-06, 6.0000000000000000, 0.0 },
  { 2.7958823432049148e-06, 6.5000000000000000, 0.0 },
  { 7.4921288639971570e-07, 7.0000000000000000, 0.0 },
  { 1.9172560675134295e-07, 7.5000000000000000, 0.0 },
  { 4.6922076160992236e-08, 8.0000000000000000, 0.0 },
  { 1.0997009755195515e-08, 8.5000000000000000, 0.0 },
  { 2.4711684308724904e-09, 9.0000000000000000, 0.0 },
  { 5.3302637046174900e-10, 9.5000000000000000, 0.0 },
  { 1.1047532552898652e-10, 10.000000000000000, 0.0 },
};
const double toler001 = 5.0000000000000029e-12;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_airy_ai<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::airy_ai(data[i].x);
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