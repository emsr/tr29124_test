
// ellint_rf

#include <complex>
#include <vector>

#include "verify.h"

#include <emsr/sf_ellint.h>

int
test01()
{
  using cmplx = std::complex<double>;

  std::vector<cmplx> rfvals
  {
    emsr::ellint_rf(1.0, 2.0, 0.0),
    emsr::ellint_rf(cmplx(0.0, 1.0), cmplx(0.0, -1.0), cmplx(0.0, 0.0)),
    emsr::ellint_rf(cmplx(-1.0, 1.0), cmplx(0.0, 1.0), cmplx(0.0, 0.0)),
    emsr::ellint_rf(0.5, 1.0, 0.0),
    emsr::ellint_rf(2.0, 3.0, 4.0),
    emsr::ellint_rf(cmplx(0.0, 1.0), cmplx(0.0, -1.0), cmplx(2.0, 0.0)),
    emsr::ellint_rf(cmplx(-1.0, 1.0), cmplx(1.0, -1.0), cmplx(0.0, 1.0))
  };

  std::vector<cmplx> rftests
  {
    cmplx(1.3110287771461),
    cmplx(1.8540746773014),
    cmplx(0.79612586584234, -1.2138566698365),
    cmplx(1.8540746773014),
    cmplx(0.58408284167715),
    cmplx(1.0441445654064),
    cmplx(0.93912050218619, -0.53296252018635),
  };

  double max_abs_diff = double(-1);
  double max_abs_frac = double(-1);
  double eps = 1.0e-12;
  bool failure = false;
  for (auto i = 0ul; i < rfvals.size(); ++i)
    {
        const auto f = rfvals[i];
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const auto f0 = rftests[i];
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
