
//  bernoulli

#include "verify.h"

// Test data.
// max(|f - f_Boost|): 1.1929908040578744e+64 at index 100
// max(|f - f_Boost| / |f_Boost|): 4.2734227959298992e-15
// mean(f - f_Boost): -1.1764257338735019e+62
// variance(f - f_Boost): inf
// stddev(f - f_Boost): inf
const testcase_bernoulli<double>
data001[101] =
{
  { 1.0000000000000000, 0, 0.0 },
  { -0.50000000000000000, 1, 0.0 },
  { 0.16666666666666666, 2, 0.0 },
  { 0.0000000000000000, 3, 0.0 },
  { -0.033333333333333333, 4, 0.0 },
  { 0.0000000000000000, 5, 0.0 },
  { 0.023809523809523808, 6, 0.0 },
  { 0.0000000000000000, 7, 0.0 },
  { -0.033333333333333333, 8, 0.0 },
  { 0.0000000000000000, 9, 0.0 },
  { 0.075757575757575760, 10, 0.0 },
  { 0.0000000000000000, 11, 0.0 },
  { -0.25311355311355310, 12, 0.0 },
  { 0.0000000000000000, 13, 0.0 },
  { 1.1666666666666667, 14, 0.0 },
  { 0.0000000000000000, 15, 0.0 },
  { -7.0921568627450977, 16, 0.0 },
  { 0.0000000000000000, 17, 0.0 },
  { 54.971177944862156, 18, 0.0 },
  { 0.0000000000000000, 19, 0.0 },
  { -529.12424242424242, 20, 0.0 },
  { 0.0000000000000000, 21, 0.0 },
  { 6192.1231884057970, 22, 0.0 },
  { 0.0000000000000000, 23, 0.0 },
  { -86580.253113553117, 24, 0.0 },
  { 0.0000000000000000, 25, 0.0 },
  { 1425517.1666666667, 26, 0.0 },
  { 0.0000000000000000, 27, 0.0 },
  { -27298231.067816094, 28, 0.0 },
  { 0.0000000000000000, 29, 0.0 },
  { 601580873.90064240, 30, 0.0 },
  { 0.0000000000000000, 31, 0.0 },
  { -15116315767.092157, 32, 0.0 },
  { 0.0000000000000000, 33, 0.0 },
  { 429614643061.16669, 34, 0.0 },
  { 0.0000000000000000, 35, 0.0 },
  { -13711655205088.332, 36, 0.0 },
  { 0.0000000000000000, 37, 0.0 },
  { 488332318973593.19, 38, 0.0 },
  { 0.0000000000000000, 39, 0.0 },
  { -19296579341940068., 40, 0.0 },
  { 0.0000000000000000, 41, 0.0 },
  { 8.4169304757368256e+17, 42, 0.0 },
  { 0.0000000000000000, 43, 0.0 },
  { -4.0338071854059454e+19, 44, 0.0 },
  { 0.0000000000000000, 45, 0.0 },
  { 2.1150748638081993e+21, 46, 0.0 },
  { 0.0000000000000000, 47, 0.0 },
  { -1.2086626522296526e+23, 48, 0.0 },
  { 0.0000000000000000, 49, 0.0 },
  { 7.5008667460769642e+24, 50, 0.0 },
  { 0.0000000000000000, 51, 0.0 },
  { -5.0387781014810688e+26, 52, 0.0 },
  { 0.0000000000000000, 53, 0.0 },
  { 3.6528776484818122e+28, 54, 0.0 },
  { 0.0000000000000000, 55, 0.0 },
  { -2.8498769302450882e+30, 56, 0.0 },
  { 0.0000000000000000, 57, 0.0 },
  { 2.3865427499683627e+32, 58, 0.0 },
  { 0.0000000000000000, 59, 0.0 },
  { -2.1399949257225335e+34, 60, 0.0 },
  { 0.0000000000000000, 61, 0.0 },
  { 2.0500975723478097e+36, 62, 0.0 },
  { 0.0000000000000000, 63, 0.0 },
  { -2.0938005911346379e+38, 64, 0.0 },
  { 0.0000000000000000, 65, 0.0 },
  { 2.2752696488463515e+40, 66, 0.0 },
  { 0.0000000000000000, 67, 0.0 },
  { -2.6257710286239577e+42, 68, 0.0 },
  { 0.0000000000000000, 69, 0.0 },
  { 3.2125082102718032e+44, 70, 0.0 },
  { 0.0000000000000000, 71, 0.0 },
  { -4.1598278166794712e+46, 72, 0.0 },
  { 0.0000000000000000, 73, 0.0 },
  { 5.6920695482035283e+48, 74, 0.0 },
  { 0.0000000000000000, 75, 0.0 },
  { -8.2183629419784578e+50, 76, 0.0 },
  { 0.0000000000000000, 77, 0.0 },
  { 1.2502904327166994e+53, 78, 0.0 },
  { 0.0000000000000000, 79, 0.0 },
  { -2.0015583233248370e+55, 80, 0.0 },
  { 0.0000000000000000, 81, 0.0 },
  { 3.3674982915364376e+57, 82, 0.0 },
  { 0.0000000000000000, 83, 0.0 },
  { -5.9470970503135450e+59, 84, 0.0 },
  { 0.0000000000000000, 85, 0.0 },
  { 1.1011910323627977e+62, 86, 0.0 },
  { 0.0000000000000000, 87, 0.0 },
  { -2.1355259545253502e+64, 88, 0.0 },
  { 0.0000000000000000, 89, 0.0 },
  { 4.3328896986641194e+66, 90, 0.0 },
  { 0.0000000000000000, 91, 0.0 },
  { -9.1885528241669332e+68, 92, 0.0 },
  { 0.0000000000000000, 93, 0.0 },
  { 2.0346896776329074e+71, 94, 0.0 },
  { 0.0000000000000000, 95, 0.0 },
  { -4.7003833958035730e+73, 96, 0.0 },
  { 0.0000000000000000, 97, 0.0 },
  { 1.1318043445484249e+76, 98, 0.0 },
  { 0.0000000000000000, 99, 0.0 },
  { -2.8382249570693707e+78, 100, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_bernoulli<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::bernoulli<Ret>(data[i].n);
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