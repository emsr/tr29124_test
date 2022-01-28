
//  eulerian_1

#include "verify.h"

// Test data for n=0.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data001[11] =
{
  { 1.0000000000000000, 0, 0, 0.0 },
  { 0.0000000000000000, 0, 1, 0.0 },
  { 0.0000000000000000, 0, 2, 0.0 },
  { 0.0000000000000000, 0, 3, 0.0 },
  { 0.0000000000000000, 0, 4, 0.0 },
  { 0.0000000000000000, 0, 5, 0.0 },
  { 0.0000000000000000, 0, 6, 0.0 },
  { 0.0000000000000000, 0, 7, 0.0 },
  { 0.0000000000000000, 0, 8, 0.0 },
  { 0.0000000000000000, 0, 9, 0.0 },
  { 0.0000000000000000, 0, 10, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

// Test data for n=1.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data002[11] =
{
  { 1.0000000000000000, 1, 0, 0.0 },
  { 0.0000000000000000, 1, 1, 0.0 },
  { 0.0000000000000000, 1, 2, 0.0 },
  { 0.0000000000000000, 1, 3, 0.0 },
  { 0.0000000000000000, 1, 4, 0.0 },
  { 0.0000000000000000, 1, 5, 0.0 },
  { 0.0000000000000000, 1, 6, 0.0 },
  { 0.0000000000000000, 1, 7, 0.0 },
  { 0.0000000000000000, 1, 8, 0.0 },
  { 0.0000000000000000, 1, 9, 0.0 },
  { 0.0000000000000000, 1, 10, 0.0 },
};
const double toler002 = 2.5000000000000020e-13;

// Test data for n=2.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data003[11] =
{
  { 1.0000000000000000, 2, 0, 0.0 },
  { 1.0000000000000000, 2, 1, 0.0 },
  { 0.0000000000000000, 2, 2, 0.0 },
  { 0.0000000000000000, 2, 3, 0.0 },
  { 0.0000000000000000, 2, 4, 0.0 },
  { 0.0000000000000000, 2, 5, 0.0 },
  { 0.0000000000000000, 2, 6, 0.0 },
  { 0.0000000000000000, 2, 7, 0.0 },
  { 0.0000000000000000, 2, 8, 0.0 },
  { 0.0000000000000000, 2, 9, 0.0 },
  { 0.0000000000000000, 2, 10, 0.0 },
};
const double toler003 = 2.5000000000000020e-13;

// Test data for n=3.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data004[11] =
{
  { 1.0000000000000000, 3, 0, 0.0 },
  { 4.0000000000000000, 3, 1, 0.0 },
  { 1.0000000000000000, 3, 2, 0.0 },
  { 0.0000000000000000, 3, 3, 0.0 },
  { 0.0000000000000000, 3, 4, 0.0 },
  { 0.0000000000000000, 3, 5, 0.0 },
  { 0.0000000000000000, 3, 6, 0.0 },
  { 0.0000000000000000, 3, 7, 0.0 },
  { 0.0000000000000000, 3, 8, 0.0 },
  { 0.0000000000000000, 3, 9, 0.0 },
  { 0.0000000000000000, 3, 10, 0.0 },
};
const double toler004 = 2.5000000000000020e-13;

// Test data for n=4.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data005[11] =
{
  { 1.0000000000000000, 4, 0, 0.0 },
  { 11.000000000000000, 4, 1, 0.0 },
  { 11.000000000000000, 4, 2, 0.0 },
  { 1.0000000000000000, 4, 3, 0.0 },
  { 0.0000000000000000, 4, 4, 0.0 },
  { 0.0000000000000000, 4, 5, 0.0 },
  { 0.0000000000000000, 4, 6, 0.0 },
  { 0.0000000000000000, 4, 7, 0.0 },
  { 0.0000000000000000, 4, 8, 0.0 },
  { 0.0000000000000000, 4, 9, 0.0 },
  { 0.0000000000000000, 4, 10, 0.0 },
};
const double toler005 = 2.5000000000000020e-13;

// Test data for n=5.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data006[11] =
{
  { 1.0000000000000000, 5, 0, 0.0 },
  { 26.000000000000000, 5, 1, 0.0 },
  { 66.000000000000000, 5, 2, 0.0 },
  { 26.000000000000000, 5, 3, 0.0 },
  { 1.0000000000000000, 5, 4, 0.0 },
  { 0.0000000000000000, 5, 5, 0.0 },
  { 0.0000000000000000, 5, 6, 0.0 },
  { 0.0000000000000000, 5, 7, 0.0 },
  { 0.0000000000000000, 5, 8, 0.0 },
  { 0.0000000000000000, 5, 9, 0.0 },
  { 0.0000000000000000, 5, 10, 0.0 },
};
const double toler006 = 2.5000000000000020e-13;

// Test data for n=6.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data007[11] =
{
  { 1.0000000000000000, 6, 0, 0.0 },
  { 57.000000000000000, 6, 1, 0.0 },
  { 302.00000000000000, 6, 2, 0.0 },
  { 302.00000000000000, 6, 3, 0.0 },
  { 57.000000000000000, 6, 4, 0.0 },
  { 1.0000000000000000, 6, 5, 0.0 },
  { 0.0000000000000000, 6, 6, 0.0 },
  { 0.0000000000000000, 6, 7, 0.0 },
  { 0.0000000000000000, 6, 8, 0.0 },
  { 0.0000000000000000, 6, 9, 0.0 },
  { 0.0000000000000000, 6, 10, 0.0 },
};
const double toler007 = 2.5000000000000020e-13;

// Test data for n=7.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data008[11] =
{
  { 1.0000000000000000, 7, 0, 0.0 },
  { 120.00000000000000, 7, 1, 0.0 },
  { 1191.0000000000000, 7, 2, 0.0 },
  { 2416.0000000000000, 7, 3, 0.0 },
  { 1191.0000000000000, 7, 4, 0.0 },
  { 120.00000000000000, 7, 5, 0.0 },
  { 1.0000000000000000, 7, 6, 0.0 },
  { 0.0000000000000000, 7, 7, 0.0 },
  { 0.0000000000000000, 7, 8, 0.0 },
  { 0.0000000000000000, 7, 9, 0.0 },
  { 0.0000000000000000, 7, 10, 0.0 },
};
const double toler008 = 2.5000000000000020e-13;

// Test data for n=8.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data009[11] =
{
  { 1.0000000000000000, 8, 0, 0.0 },
  { 247.00000000000000, 8, 1, 0.0 },
  { 4293.0000000000000, 8, 2, 0.0 },
  { 15619.000000000000, 8, 3, 0.0 },
  { 15619.000000000000, 8, 4, 0.0 },
  { 4293.0000000000000, 8, 5, 0.0 },
  { 247.00000000000000, 8, 6, 0.0 },
  { 1.0000000000000000, 8, 7, 0.0 },
  { 0.0000000000000000, 8, 8, 0.0 },
  { 0.0000000000000000, 8, 9, 0.0 },
  { 0.0000000000000000, 8, 10, 0.0 },
};
const double toler009 = 2.5000000000000020e-13;

// Test data for n=9.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data010[11] =
{
  { 1.0000000000000000, 9, 0, 0.0 },
  { 502.00000000000000, 9, 1, 0.0 },
  { 14608.000000000000, 9, 2, 0.0 },
  { 88234.000000000000, 9, 3, 0.0 },
  { 156190.00000000000, 9, 4, 0.0 },
  { 88234.000000000000, 9, 5, 0.0 },
  { 14608.000000000000, 9, 6, 0.0 },
  { 502.00000000000000, 9, 7, 0.0 },
  { 1.0000000000000000, 9, 8, 0.0 },
  { 0.0000000000000000, 9, 9, 0.0 },
  { 0.0000000000000000, 9, 10, 0.0 },
};
const double toler010 = 2.5000000000000020e-13;

// Test data for n=10.
// max(|f - f_Burkhardt|): 0.0000000000000000 at index 0
// max(|f - f_Burkhardt| / |f_Burkhardt|): 0.0000000000000000
// mean(f - f_Burkhardt): 0.0000000000000000
// variance(f - f_Burkhardt): 0.0000000000000000
// stddev(f - f_Burkhardt): 0.0000000000000000
const testcase_eulerian_1<double>
data011[11] =
{
  { 1.0000000000000000, 10, 0, 0.0 },
  { 1013.0000000000000, 10, 1, 0.0 },
  { 47840.000000000000, 10, 2, 0.0 },
  { 455192.00000000000, 10, 3, 0.0 },
  { 1310354.0000000000, 10, 4, 0.0 },
  { 1310354.0000000000, 10, 5, 0.0 },
  { 455192.00000000000, 10, 6, 0.0 },
  { 47840.000000000000, 10, 7, 0.0 },
  { 1013.0000000000000, 10, 8, 0.0 },
  { 1.0000000000000000, 10, 9, 0.0 },
  { 0.0000000000000000, 10, 10, 0.0 },
};
const double toler011 = 2.5000000000000020e-13;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_eulerian_1<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::eulerian_1<Ret>(data[i].n, data[i].m);
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
  num_errors += test(data009, toler009);
  num_errors += test(data010, toler010);
  num_errors += test(data011, toler011);
  return num_errors;
}
