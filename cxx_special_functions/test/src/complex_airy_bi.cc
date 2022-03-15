
// airy_bi

#include <complex>
#include <vector>

#include "verify.h"

#include <emsr/sf_airy.h>

int
test01()
{
  using cmplx = std::complex<double>;

  std::vector<cmplx> bivals
  {
    emsr::airy_bi(cmplx( 1.0000000000000000,  0.0000000000000000)),
    emsr::airy_bi(cmplx( 0.8090169943749474,  0.5877852522924731)),
    emsr::airy_bi(cmplx( 0.3090169943749474,  0.9510565162951536)),
    emsr::airy_bi(cmplx(-0.3090169943749474,  0.9510565162951536)),
    emsr::airy_bi(cmplx(-0.8090169943749474,  0.5877852522924731)),
    emsr::airy_bi(cmplx(-1.0000000000000000,  0.0000000000000000)),
    emsr::airy_bi(cmplx(-0.8090169943749474, -0.5877852522924731)),
    emsr::airy_bi(cmplx(-0.3090169943749474, -0.9510565162951536)),
    emsr::airy_bi(cmplx( 0.3090169943749474, -0.9510565162951536)),
    emsr::airy_bi(cmplx( 0.8090169943749474, -0.5877852522924731))
  };

  std::vector<cmplx> bitests
  {
    cmplx( 1.207423594952871,   0.0000000000000000),
    cmplx( 0.9127160108293936,  0.3800456133135556),
    cmplx( 0.6824453575635721,  0.3343047153635002),
    cmplx( 0.5726265660086474,  0.3988641086982559),
    cmplx( 0.2511841251049547,  0.3401447690712719),
    cmplx( 0.1039973894969446,  0.0000000000000000),
    cmplx( 0.2511841251049547, -0.3401447690712719),
    cmplx( 0.5726265660086474, -0.3988641086982559),
    cmplx( 0.6824453575635721, -0.3343047153635002),
    cmplx( 0.9127160108293936, -0.3800456133135556)
  };

  double max_abs_diff = double(-1);
  double max_abs_frac = double(-1);
  double eps = 1.0e-12;
  bool failure = false;
  for (auto i = 0ul; i < bivals.size(); ++i)
    {
        const auto f = bivals[i];
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const auto f0 = bitests[i];
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
