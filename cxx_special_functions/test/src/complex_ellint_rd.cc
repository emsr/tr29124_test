
// ellint_rd

#include <complex>
#include <vector>

#include "verify.h"

#include <emsr/sf_ellint.h>

int
test01()
{
  using cmplx = std::complex<double>;

  std::vector<cmplx> rdvals
  {
    emsr::ellint_rd(0.0, 2.0, 1.0),
    emsr::ellint_rd(2.0, 3.0, 4.0),
    emsr::ellint_rd(cmplx(0.0, 1.0), cmplx(0.0, -1.0), cmplx(2.0, 0.0)),
    emsr::ellint_rd(cmplx(0.0, 0.0), cmplx(0.0, 1.0), cmplx(0.0, -1.0)),
    emsr::ellint_rd(cmplx(0.0, 0.0), cmplx(-1.0, 1.0), cmplx(0.0, 1.0)),
    emsr::ellint_rd(cmplx(-2.0, -1.0), cmplx(0.0, -1.0), cmplx(-1.0, 1.0))
  };

  std::vector<cmplx> rdtests
  {
    cmplx(1.7972103521034),
    cmplx(0.16510527294261),
    cmplx(0.65933854154220),
    cmplx(1.2708196271910, +2.7811120159521),
    cmplx(-1.8577235439239, -0.96193450888839),
    cmplx(1.8249027393704, -1.2218475784827),
  };

  double max_abs_diff = double(-1);
  double max_abs_frac = double(-1);
  double eps = 1.0e-12;
  bool failure = false;
  for (auto i = 0ul; i < rdvals.size(); ++i)
    {
        const auto f = rdvals[i];
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const auto f0 = rdtests[i];
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
