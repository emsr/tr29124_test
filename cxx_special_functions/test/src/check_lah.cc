
//  lah

//  Compare against values from teh interwebs.
//  https://en.wikipedia.org/wiki/Lah_number
//  I start with degree 0 lah numbers and include order 0.
//  Their L_11^8 is wrong! 11880 should be 118800!

#include "verify.h"

long long
lah_number_0[]
{
  1ll
};

long long
lah_number_1[]
{
  0ll,
  1ll
};

long long
lah_number_2[]
{
  0ll,
  2ll,
  1ll
};

long long
lah_number_3[]
{
  0ll,
  6ll,
  6ll,
  1ll
};

long long
lah_number_4[]
{
  0ll,
  24ll,
  36ll,
  12ll,
  1ll
};

long long
lah_number_5[]
{
  0ll,
  120ll,
  240ll,
  120ll,
  20ll,
  1ll
};

long long
lah_number_6[]
{
  0ll,
  720ll,
  1800ll,
  1200ll,
  300ll,
  30ll,
  1ll
};

long long
lah_number_7[]
{
  0ll,
  5040ll,
  15120ll,
  12600ll,
  4200ll,
  630ll,
  42ll,
  1ll
};

long long
lah_number_8[]
{
  0ll,
  40320ll,
  141120ll,
  141120ll,
  58800ll,
  11760ll,
  1176ll,
  56ll,
  1ll
};

long long
lah_number_9[]
{
  0ll,
  362880ll,
  1451520ll,
  1693440ll,
  846720ll,
  211680ll,
  28224ll,
  2016ll,
  72ll,
  1ll
};

long long
lah_number_10[]
{
  0ll,
  3628800ll,
  16329600ll,
  21772800ll,
  12700800ll,
  3810240ll,
  635040ll,
  60480ll,
  3240ll,
  90ll,
  1ll
};

long long
lah_number_11[]
{
  0ll,
  39916800ll,
  199584000ll,
  299376000ll,
  199584000ll,
  69854400ll,
  13970880ll,
  1663200ll,
  118800ll,
  4950ll,
  110ll,
  1ll
};

long long
lah_number_12[]
{
  0ll,
  479001600ll,
  2634508800ll,
  4390848000ll,
  3293136000ll,
  1317254400ll,
  307359360ll,
  43908480ll,
  3920400ll,
  217800ll,
  7260ll,
  132ll,
  1ll
};

template<typename Ret, unsigned int Num>
  int
  test(const long long (&lah_numbers)[Num])
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    const Ret toler = 100 * eps;
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    const auto lah = emsr::lah<Ret>(Num - 1);
    for (unsigned int i = 0; i < Num; ++i)
      {
	const Ret f0 = lah_numbers[i];
	const Ret f = lah[i];
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
    int num_errors = 0;
    VERIFY(max_abs_frac < toler);
    return num_errors;
  }

int
main()
{
  int num_errors = 0;
  num_errors = test<long double>(lah_number_0);
  num_errors = test<long double>(lah_number_1);
  num_errors = test<long double>(lah_number_2);
  num_errors = test<long double>(lah_number_3);
  num_errors = test<long double>(lah_number_4);
  num_errors = test<long double>(lah_number_5);
  num_errors = test<long double>(lah_number_6);
  num_errors = test<long double>(lah_number_7);
  num_errors = test<long double>(lah_number_8);
  num_errors = test<long double>(lah_number_9);
  num_errors = test<long double>(lah_number_10);
  num_errors = test<long double>(lah_number_11);
  num_errors = test<long double>(lah_number_12);
  return num_errors;
}
