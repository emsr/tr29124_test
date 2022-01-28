
// ellint_rg

#include <complex>
#include <vector>

#include "verify.h"

#include <emsr/sf_ellint.h>

int
test01()
{
  using cmplx = std::complex<double>;

  std::vector<cmplx> rgvals
  {
    emsr::ellint_rg(0.0, 16.0, 16.0),
    emsr::ellint_rg(2, 3, 4),
    emsr::ellint_rg(cmplx(0.0, 0.0), cmplx(0.0, 1.0), cmplx(0.0, -1.0)),
    emsr::ellint_rg(cmplx(-1.0, 1.0), cmplx(0.0, 1.0), cmplx(0.0, 0.0)),
    emsr::ellint_rg(cmplx(0.0, -1.0), cmplx(-1.0, 1.0), cmplx(0.0, 1.0)),
    emsr::ellint_rg(0, 0.0796, 4)
  };

  std::vector<cmplx> rgtests
  {
    cmplx(3.1415926535898),
    cmplx(1.7255030280692),
    cmplx(0.42360654239699),
    cmplx(0.44660591677018, +0.70768352357515),
    cmplx(0.36023392184473, +0.40348623401722),
    cmplx(1.0284758090288)
  };


  double max_abs_diff = double(-1);
  double max_abs_frac = double(-1);
  double eps = 1.0e-12;
  bool failure = false;
  for (auto i = 0ul; i < rgvals.size(); ++i)
    {
        const auto f = rgvals[i];
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const auto f0 = rgtests[i];
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
