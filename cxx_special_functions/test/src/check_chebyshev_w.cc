
//  chebyshev_w

#include "verify.h"

// Test data for n=0.
// max(|f - f_GSL|): 0.0000000000000000 at index 0
// max(|f - f_GSL| / |f_GSL|): 0.0000000000000000
// mean(f - f_GSL): 0.0000000000000000
// variance(f - f_GSL): 0.0000000000000000
// stddev(f - f_GSL): 0.0000000000000000
const testcase_chebyshev_w<double>
data001[21] =
{
  { 1.0000000000000000, 0, -1.0000000000000000, 0.0 },
  { 1.0000000000000000, 0, -0.90000000000000002, 0.0 },
  { 1.0000000000000000, 0, -0.80000000000000004, 0.0 },
  { 1.0000000000000000, 0, -0.69999999999999996, 0.0 },
  { 1.0000000000000000, 0, -0.59999999999999998, 0.0 },
  { 1.0000000000000000, 0, -0.50000000000000000, 0.0 },
  { 1.0000000000000000, 0, -0.39999999999999991, 0.0 },
  { 1.0000000000000000, 0, -0.29999999999999993, 0.0 },
  { 1.0000000000000000, 0, -0.19999999999999996, 0.0 },
  { 1.0000000000000000, 0, -0.099999999999999978, 0.0 },
  { 1.0000000000000000, 0, 0.0000000000000000, 0.0 },
  { 1.0000000000000000, 0, 0.10000000000000009, 0.0 },
  { 1.0000000000000000, 0, 0.20000000000000018, 0.0 },
  { 1.0000000000000000, 0, 0.30000000000000004, 0.0 },
  { 1.0000000000000000, 0, 0.40000000000000013, 0.0 },
  { 1.0000000000000000, 0, 0.50000000000000000, 0.0 },
  { 1.0000000000000000, 0, 0.60000000000000009, 0.0 },
  { 1.0000000000000000, 0, 0.70000000000000018, 0.0 },
  { 1.0000000000000000, 0, 0.80000000000000004, 0.0 },
  { 1.0000000000000000, 0, 0.90000000000000013, 0.0 },
  { 1.0000000000000000, 0, 1.0000000000000000, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

// Test data for n=1.
// max(|f - f_GSL|): 3.7138001787229613e-16 at index 5
// max(|f - f_GSL| / |f_GSL|): 1.3877787807814465e-15
// mean(f - f_GSL): 2.6936621294366595e-17
// variance(f - f_GSL): 3.8093032254696509e-35
// stddev(f - f_GSL): 6.1719553023897140e-18
const testcase_chebyshev_w<double>
data002[21] =
{
  { -1.0000000000000000, 1, -1.0000000000000000, 0.0 },
  { -0.80000000000000016, 1, -0.90000000000000002, 0.0 },
  { -0.60000000000000009, 1, -0.80000000000000004, 0.0 },
  { -0.39999999999999986, 1, -0.69999999999999996, 0.0 },
  { -0.19999999999999971, 1, -0.59999999999999998, 0.0 },
  { -3.7138001787229613e-16, 1, -0.50000000000000000, 0.0 },
  { 0.19999999999999990, 1, -0.39999999999999991, 0.0 },
  { 0.40000000000000002, 1, -0.29999999999999993, 0.0 },
  { 0.59999999999999998, 1, -0.19999999999999996, 0.0 },
  { 0.80000000000000016, 1, -0.099999999999999978, 0.0 },
  { 1.0000000000000002, 1, 0.0000000000000000, 0.0 },
  { 1.2000000000000000, 1, 0.10000000000000009, 0.0 },
  { 1.4000000000000004, 1, 0.20000000000000018, 0.0 },
  { 1.6000000000000001, 1, 0.30000000000000004, 0.0 },
  { 1.8000000000000003, 1, 0.40000000000000013, 0.0 },
  { 2.0000000000000000, 1, 0.50000000000000000, 0.0 },
  { 2.2000000000000002, 1, 0.60000000000000009, 0.0 },
  { 2.4000000000000004, 1, 0.70000000000000018, 0.0 },
  { 2.6000000000000001, 1, 0.80000000000000004, 0.0 },
  { 2.8000000000000003, 1, 0.90000000000000013, 0.0 },
  { 3.0000000000000000, 1, 1.0000000000000000, 0.0 },
};
const double toler002 = 2.5000000000000020e-13;

// Test data for n=5.
// max(|f - f_GSL|): 2.2204460492503131e-15 at index 19
// max(|f - f_GSL| / |f_GSL|): 4.4796955301271334e-15
// mean(f - f_GSL): 4.4937598615780146e-17
// variance(f - f_GSL): 1.0601785789103069e-34
// stddev(f - f_GSL): 1.0296497360317764e-17
const testcase_chebyshev_w<double>
data003[21] =
{
  { -1.0000000000000000, 5, -1.0000000000000000, 0.0 },
  { 0.80992000000000008, 5, -0.90000000000000002, 0.0 },
  { 0.97184000000000004, 5, -0.80000000000000004, 0.0 },
  { 0.35935999999999785, 5, -0.69999999999999996, 0.0 },
  { -0.42272000000000015, 5, -0.59999999999999998, 0.0 },
  { -0.99999999999999944, 5, -0.50000000000000000, 0.0 },
  { -1.1900799999999998, 5, -0.39999999999999991, 0.0 },
  { -0.96415999999999946, 5, -0.29999999999999993, 0.0 },
  { -0.40864000000000017, 5, -0.19999999999999996, 0.0 },
  { 0.31328000000000084, 5, -0.099999999999999978, 0.0 },
  { 1.0000000000000013, 5, 0.0000000000000000, 0.0 },
  { 1.4499199999999999, 5, 0.10000000000000009, 0.0 },
  { 1.4998399999999998, 5, 0.20000000000000018, 0.0 },
  { 1.0633599999999990, 5, 0.30000000000000004, 0.0 },
  { 0.16927999999999876, 5, 0.40000000000000013, 0.0 },
  { -0.99999999999999933, 5, 0.50000000000000000, 0.0 },
  { -2.0700800000000008, 5, 0.60000000000000009, 0.0 },
  { -2.4361599999999992, 5, 0.70000000000000018, 0.0 },
  { -1.2246399999999982, 5, 0.80000000000000004, 0.0 },
  { 2.7452800000000068, 5, 0.90000000000000013, 0.0 },
  { 11.000000000000000, 5, 1.0000000000000000, 0.0 },
};
const double toler003 = 2.5000000000000020e-13;

// Test data for n=8.
// max(|f - f_GSL|): 2.5604518505417673e-15 at index 4
// max(|f - f_GSL| / |f_GSL|): 8.1720443183085951e-14
// mean(f - f_GSL): 3.0233752009881494e-16
// variance(f - f_GSL): 4.7989187431238202e-33
// stddev(f - f_GSL): 6.9274228563902614e-17
const testcase_chebyshev_w<double>
data004[21] =
{
  { 1.0000000000000000, 8, -1.0000000000000000, 0.0 },
  { -0.78988543999999905, 8, -0.90000000000000002, 0.0 },
  { 0.72417535999999993, 8, -0.80000000000000004, 0.0 },
  { 0.96322815999999889, 8, -0.69999999999999996, 0.0 },
  { -0.031331840000002782, 8, -0.59999999999999998, 0.0 },
  { -0.99999999999999967, 8, -0.50000000000000000, 0.0 },
  { -1.0868710400000003, 8, -0.39999999999999991, 0.0 },
  { -0.28722944000000117, 8, -0.29999999999999993, 0.0 },
  { 0.77578495999999975, 8, -0.19999999999999996, 0.0 },
  { 1.3454617600000003, 8, -0.099999999999999978, 0.0 },
  { 1.0000000000000004, 8, 0.0000000000000000, 0.0 },
  { -0.098352640000000116, 8, 0.10000000000000009, 0.0 },
  { -1.2638182400000013, 8, 0.20000000000000018, 0.0 },
  { -1.6443622400000000, 8, 0.30000000000000004, 0.0 },
  { -0.75960063999999949, 8, 0.40000000000000013, 0.0 },
  { 0.99999999999999956, 8, 0.50000000000000000, 0.0 },
  { 2.2351897600000004, 8, 0.60000000000000009, 0.0 },
  { 1.1870489599999958, 8, 0.70000000000000018, 0.0 },
  { -2.2978534400000030, 8, 0.80000000000000004, 0.0 },
  { -2.8540390399999929, 8, 0.90000000000000013, 0.0 },
  { 17.000000000000000, 8, 1.0000000000000000, 0.0 },
};
const double toler004 = 5.0000000000000029e-12;

// Test data for n=10.
// max(|f - f_GSL|): 3.1124501748083564e-15 at index 5
// max(|f - f_GSL| / |f_GSL|): 6.0032321926169395e-15
// mean(f - f_GSL): 1.6192448913430982e-16
// variance(f - f_GSL): 1.3765258595238791e-33
// stddev(f - f_GSL): 3.7101561416251460e-17
const testcase_chebyshev_w<double>
data005[21] =
{
  { 1.0000000000000000, 10, -1.0000000000000000, 0.0 },
  { 0.023998054399998959, 10, -0.90000000000000002, 0.0 },
  { 0.93808220160000000, 10, -0.80000000000000004, 0.0 },
  { -0.51782512640000378, 10, -0.69999999999999996, 0.0 },
  { -1.0641181695999991, 10, -0.59999999999999998, 0.0 },
  { -3.1124501748083564e-15, 10, -0.50000000000000000, 0.0 },
  { 1.1036806144000004, 10, -0.39999999999999991, 0.0 },
  { 0.92616028160000097, 10, -0.29999999999999993, 0.0 },
  { -0.30930032639999966, 10, -0.19999999999999996, 0.0 },
  { -1.3008490496000005, 10, -0.099999999999999978, 0.0 },
  { -1.0000000000000002, 10, 0.0000000000000000, 0.0 },
  { 0.39238717440000015, 10, 0.10000000000000009, 0.0 },
  { 1.5350895616000009, 10, 0.20000000000000018, 0.0 },
  { 1.1243380735999995, 10, 0.30000000000000004, 0.0 },
  { -0.70076712960000320, 10, 0.40000000000000013, 0.0 },
  { -2.0000000000000000, 10, 0.50000000000000000, 0.0 },
  { -0.68601026559999589, 10, 0.60000000000000009, 0.0 },
  { 2.2687420416000044, 10, 0.70000000000000018, 0.0 },
  { 1.4422260735999965, 10, 0.80000000000000004, 0.0 },
  { -4.4709124096000030, 10, 0.90000000000000013, 0.0 },
  { 21.000000000000000, 10, 1.0000000000000000, 0.0 },
};
const double toler005 = 5.0000000000000039e-13;

// Test data for n=20.
// max(|f - f_GSL|): 7.0499162063697440e-15 at index 9
// max(|f - f_GSL| / |f_GSL|): 1.8173872761738678e-14
// mean(f - f_GSL): 2.2204460492503131e-16
// variance(f - f_GSL): 2.5884498452564449e-33
// stddev(f - f_GSL): 5.0876810486276013e-17
const testcase_chebyshev_w<double>
data006[21] =
{
  { 1.0000000000000000, 20, -1.0000000000000000, 0.0 },
  { -1.0096350973538502, 20, -0.90000000000000002, 0.0 },
  { 0.85458211259118788, 20, -0.80000000000000004, 0.0 },
  { -0.89660062573515276, 20, -0.69999999999999996, 0.0 },
  { 1.1037543614593988, 20, -0.59999999999999998, 0.0 },
  { -0.99999999999999611, 20, -0.50000000000000000, 0.0 },
  { 0.24130140035435260, 20, -0.39999999999999991, 0.0 },
  { 0.84402681417869085, 20, -0.29999999999999993, 0.0 },
  { -1.2650429316715133, 20, -0.19999999999999996, 0.0 },
  { 0.40203624022563073, 20, -0.099999999999999978, 0.0 },
  { 0.99999999999999623, 20, 0.0000000000000000, 0.0 },
  { -1.4229092060125628, 20, 0.10000000000000009, 0.0 },
  { 0.31543552675929620, 20, 0.20000000000000018, 0.0 },
  { 1.2386077195392600, 20, 0.30000000000000004, 0.0 },
  { -1.7881475926508603, 20, 0.40000000000000013, 0.0 },
  { 0.99999999999999423, 20, 0.50000000000000000, 0.0 },
  { 0.35623761485475713, 20, 0.60000000000000009, 0.0 },
  { -1.4530226431858577, 20, 0.70000000000000018, 0.0 },
  { 1.8512711080640289, 20, 0.80000000000000004, 0.0 },
  { 0.79504869890744312, 20, 0.90000000000000013, 0.0 },
  { 41.000000000000000, 20, 1.0000000000000000, 0.0 },
};
const double toler006 = 1.0000000000000008e-12;

// Test data for n=40.
// max(|f - f_GSL|): 9.9364960703951510e-15 at index 16
// max(|f - f_GSL| / |f_GSL|): 3.1039930526159287e-14
// mean(f - f_GSL): 1.5646802725698727e-15
// variance(f - f_GSL): 1.2853177865688990e-31
// stddev(f - f_GSL): 3.5851328937277889e-16
const testcase_chebyshev_w<double>
data007[21] =
{
  { 1.0000000000000000, 40, -1.0000000000000000, 0.0 },
  { 0.85651884958432190, 40, -0.90000000000000002, 0.0 },
  { 0.63097169179115631, 40, -0.80000000000000004, 0.0 },
  { 0.75745186284307398, 40, -0.69999999999999996, 0.0 },
  { 1.1065174331497729, 40, -0.59999999999999998, 0.0 },
  { -8.4888903334451430e-15, 40, -0.50000000000000000, 0.0 },
  { -1.1773725987499921, 40, -0.39999999999999991, 0.0 },
  { 0.65788833124633084, 40, -0.29999999999999993, 0.0 },
  { 0.60116875594653152, 40, -0.19999999999999996, 0.0 },
  { -1.3370585083307165, 40, -0.099999999999999978, 0.0 },
  { 0.99999999999999500, 40, 0.0000000000000000, 0.0 },
  { 0.19293637359528409, 40, 0.10000000000000009, 0.0 },
  { -1.3992477229963809, 40, 0.20000000000000018, 0.0 },
  { 1.4329479238333860, 40, 0.30000000000000004, 0.0 },
  { 0.31440756245615642, 40, 0.40000000000000013, 0.0 },
  { -2.0000000000000000, 40, 0.50000000000000000, 0.0 },
  { -0.32011979092611192, 40, 0.60000000000000009, 0.0 },
  { 1.8481101593324996, 40, 0.70000000000000018, 0.0 },
  { 2.5331546572256292, 40, 0.80000000000000004, 0.0 },
  { -2.4619369906292432, 40, 0.90000000000000013, 0.0 },
  { 81.000000000000000, 40, 1.0000000000000000, 0.0 },
};
const double toler007 = 2.5000000000000015e-12;

// Test data for n=100.
// max(|f - f_GSL|): 3.4999780851308060e-14 at index 3
// max(|f - f_GSL| / |f_GSL|): 2.1815036587121280e-12
// mean(f - f_GSL): 1.9429940336578212e-16
// variance(f - f_GSL): 1.9819935527856925e-33
// stddev(f - f_GSL): 4.4519586170422708e-17
const testcase_chebyshev_w<double>
data008[21] =
{
  { 1.0000000000000000, 100, -1.0000000000000000, 0.0 },
  { 0.22880117622182833, 100, -0.90000000000000002, 0.0 },
  { -0.28035903752475666, 100, -0.80000000000000004, 0.0 },
  { -0.18674254467748605, 100, -0.69999999999999996, 0.0 },
  { 0.55182443700502926, 100, -0.59999999999999998, 0.0 },
  { -1.9241770650718719e-14, 100, -0.50000000000000000, 0.0 },
  { -1.1523528299091899, 100, -0.39999999999999991, 0.0 },
  { -0.011059589108426603, 100, -0.29999999999999993, 0.0 },
  { 1.0643966872442276, 100, -0.19999999999999996, 0.0 },
  { -1.3345690929659346, 100, -0.099999999999999978, 0.0 },
  { 0.99999999999999845, 100, 0.0000000000000000, 0.0 },
  { -0.21296288072517755, 100, 0.10000000000000009, 0.0 },
  { -0.89475669222160348, 100, 0.20000000000000018, 0.0 },
  { 1.6902413146121975, 100, 0.30000000000000004, 0.0 },
  { -0.48462380179611242, 100, 0.40000000000000013, 0.0 },
  { -2.0000000000000000, 100, 0.50000000000000000, 0.0 },
  { -1.9447259865843485, 100, 0.60000000000000009, 0.0 },
  { -2.5434335568364355, 100, 0.70000000000000018, 0.0 },
  { 3.0483748605943708, 100, 0.80000000000000004, 0.0 },
  { 4.3595126348515461, 100, 0.90000000000000013, 0.0 },
  { 201.00000000000000, 100, 1.0000000000000000, 0.0 },
};
const double toler008 = 2.5000000000000017e-10;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_chebyshev_w<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::chebyshev_w(data[i].n, data[i].x);
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
  int num_errors = 0;
  num_errors += test(data001, toler001);
  num_errors += test(data002, toler002);
  num_errors += test(data003, toler003);
  num_errors += test(data004, toler004);
  num_errors += test(data005, toler005);
  num_errors += test(data006, toler006);
  num_errors += test(data007, toler007);
  num_errors += test(data008, toler008);
  return num_errors;
}