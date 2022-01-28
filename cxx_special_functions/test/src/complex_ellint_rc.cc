
// ellint_rc

#include <complex>
#include <vector>

#include "verify.h"

#include <emsr/sf_ellint.h>

int
test01()
{
  using cmplx = std::complex<double>;

  std::vector<cmplx> rcvals
  {
    emsr::ellint_rc(0.0, 0.25),
    emsr::ellint_rc(2.25, 2.0),  //  ln 2 = 0.69314718055995
    emsr::ellint_rc(cmplx(0.0, 0.0), cmplx(0.0, 1.0)),
    //emsr::ellint_rc(cmplx(0.0, -1.0), cmplx(0.0, 1.0)),
    emsr::ellint_rc(0.25, -2.0),
    emsr::ellint_rc(cmplx(0.0, 1.0), cmplx(-1.0, 0.0))
  };

  std::vector<cmplx> rctests
  {
    cmplx(3.1415926535898),
    cmplx(0.69314718055995),
    cmplx(1.1107207345396, -1.1107207345396),
    //cmplx(1.2260849569072, -0.34471136988768),
    cmplx(0.23104906018665),
    cmplx(0.77778596920447, 0.19832484993429)
  };

  double max_abs_diff = double(-1);
  double max_abs_frac = double(-1);
  double eps = 1.0e-12;
  bool failure = false;
  for (auto i = 0ul; i < rcvals.size(); ++i)
    {
        const auto f = rcvals[i];
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const auto f0 = rctests[i];
	    const auto diff = f - f0;
	    if (std::abs(diff) > max_abs_diff)
	      max_abs_diff = std::abs(diff);
	    if (std::abs(f0) > double(10) * eps
	     && std::abs(f) > double(10) * eps)
	      {
		const auto frac = diff / f0;
		if (std::abs(frac) > max_abs_frac)
		  max_abs_frac = std::abs(frac);
	      }
	  }
    }
  int num_errors = 0;
  VERIFY(max_abs_frac < eps);
  return num_errors;
}

int
main()
{
  return test01();
}
