
// airy_ai

#include <complex>
#include <vector>

#include "verify.h"

#include <emsr/sf_airy.h>

int
test01()
{
  using cmplx = std::complex<double>;

  std::vector<cmplx> aivals
  {
    emsr::airy_ai(cmplx( 1.0000000000000000,  0.0000000000000000)),
    emsr::airy_ai(cmplx( 0.8090169943749474,  0.5877852522924731)),
    emsr::airy_ai(cmplx( 0.3090169943749474,  0.9510565162951536)),
    emsr::airy_ai(cmplx(-0.3090169943749474,  0.9510565162951536)),
    emsr::airy_ai(cmplx(-0.8090169943749474,  0.5877852522924731)),
    emsr::airy_ai(cmplx(-1.0000000000000000,  0.0000000000000000)),
    emsr::airy_ai(cmplx(-0.8090169943749474, -0.5877852522924731)),
    emsr::airy_ai(cmplx(-0.3090169943749474, -0.9510565162951536)),
    emsr::airy_ai(cmplx( 0.3090169943749474, -0.9510565162951536)),
    emsr::airy_ai(cmplx( 0.8090169943749474, -0.5877852522924731))
  };

  std::vector<cmplx> aitests
  {
    cmplx( 0.1352924163128814,  0.0000000000000000),
    cmplx( 0.1433824486882056, -0.1092193342707378),
    cmplx( 0.2215404472324631, -0.2588711788891803),
    cmplx( 0.4763929771766866, -0.3036484220291284),
    cmplx( 0.5983692170633874, -0.08154602160771214),
    cmplx( 0.5355608832923521,  0.00000000000000000),
    cmplx( 0.5983692170633874,  0.08154602160771214),
    cmplx( 0.4763929771766866,  0.3036484220291284),
    cmplx( 0.2215404472324631,  0.2588711788891803),
    cmplx( 0.1433824486882056,  0.1092193342707378)
  };

  double max_abs_diff = double(-1);
  double max_abs_frac = double(-1);
  double eps = 1.0e-12;
  bool failure = false;
  for (auto i = 0ul; i < aivals.size(); ++i)
    {
        const auto f = aivals[i];
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const auto f0 = aitests[i];
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
