
// ellint_rj

#include <complex>
#include <vector>

#include "verify.h"

#include <emsr/sf_ellint.h>

int
test01()
{
  using cmplx = std::complex<double>;

  std::vector<cmplx> rjvals
  {
     emsr::ellint_rj(0.0, 1.0, 2.0, 3.0),
     emsr::ellint_rj(2.0, 3.0, 4.0, 5.0),
     emsr::ellint_rj(cmplx(2.0,0.0), cmplx(3.0,0.0),
                     cmplx(4.0,0.0), cmplx(-1.0,1.0)),
     emsr::ellint_rj(cmplx(0.0,1.0), cmplx(0.0,-1.0),
                     cmplx(0.0,0.0), cmplx(2.0,0.0)),
     emsr::ellint_rj(cmplx(-1.0,1.0), cmplx(-1.0,-1.0),
                     cmplx(1.0,0.0), cmplx(2.0,0.0)),
     emsr::ellint_rj(cmplx(0.0,1.0), cmplx(0.0,-1.0),
                     cmplx(0.0,0.0), cmplx(1.0,-1.0)),
     emsr::ellint_rj(cmplx(-1.0,1.0), cmplx(-1.0,-1.0),
                     cmplx(1.0,0.0), cmplx(-3.0,1.0)),
     emsr::ellint_rj(cmplx(-1.0,1.0), cmplx(-2.0,-1.0),
                     cmplx(0.0,-1.0), cmplx(-1.0,1.0))
  };

  std::vector<cmplx> rjtests
  {
    cmplx(0.77688623778582),
    cmplx(0.14297579667157),
    cmplx(0.13613945827771, -0.38207561624427),
    cmplx(1.6490011662711),
    cmplx(0.94148358841220),
    cmplx(1.8260115229009, +1.2290661908643),
    cmplx(-0.61127970812028, -1.0684038390007),
    cmplx(1.8249027393704, -1.2218475784827)
  };

  double max_abs_diff = double(-1);
  double max_abs_frac = double(-1);
  double eps = 1.0e-12;
  bool failure = false;
  for (auto i = 0ul; i < rjvals.size(); ++i)
    {
        const auto f = rjvals[i];
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const auto f0 = rjtests[i];
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
