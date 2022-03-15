
//  bell

//  Compare against values from The On-Line Encyclopedia of Integer Sequences:
//  https://oeis.org/A000110

#include "verify.h"

template<typename Ret>
  int
  test()
  {
    constexpr unsigned int num_bells = 20;
    constexpr uint64_t
    bell[num_bells]
    {
      1ull,
      1ull,
      2ull,
      5ull,
      15ull,
      52ull,
      203ull,
      877ull,
      4140ull,
      21147ull,
      115975ull,
      678570ull,
      4213597ull,
      27644437ull,
      190899322ull,
      1382958545ull,
      10480142147ull,
      82864869804ull,
      682076806159ull,
      5832742205057ull
    };

    const Ret eps = std::numeric_limits<Ret>::epsilon();
    const Ret toler = 100 * eps;
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    const auto bells = emsr::bell<Ret>(num_bells);
    for (unsigned int i = 0; i < num_bells; ++i)
      {
	const Ret f0 = bell[i];
	const Ret f = bells[i];
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
  num_errors += test<long double>();
  return num_errors;
}
