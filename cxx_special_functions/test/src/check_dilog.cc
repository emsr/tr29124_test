
//  dilog

#include "verify.h"

// Test data.
// max(|f - f_GSL|): 2.2204460492503131e-15 at index 4
// max(|f - f_GSL| / |f_GSL|): 1.9068116597202735e-15
// mean(f - f_GSL): -9.8069700508555491e-17
// variance(f - f_GSL): 2.2355112453651842e-34
// stddev(f - f_GSL): 1.4951626150239258e-17
const testcase_dilog<double>
data001[45] =
{
  { -4.1982778868581043, -10.000000000000000, 0.0 },
  { -4.1378595698085165, -9.7500000000000000, 0.0 },
  { -4.0764759702627735, -9.5000000000000000, 0.0 },
  { -4.0140903713193570, -9.2500000000000000, 0.0 },
  { -3.9506637782441563, -9.0000000000000000, 0.0 },
  { -3.8861547188822385, -8.7500000000000000, 0.0 },
  { -3.8205190211678008, -8.5000000000000000, 0.0 },
  { -3.7537095644470124, -8.2500000000000000, 0.0 },
  { -3.6856760007574056, -8.0000000000000000, 0.0 },
  { -3.6163644415187579, -7.7500000000000000, 0.0 },
  { -3.5457171042558460, -7.5000000000000000, 0.0 },
  { -3.4736719129571010, -7.2500000000000000, 0.0 },
  { -3.4001620444283858, -7.0000000000000000, 0.0 },
  { -3.3251154114686625, -6.7500000000000000, 0.0 },
  { -3.2484540717954049, -6.5000000000000000, 0.0 },
  { -3.1700935492807281, -6.2500000000000000, 0.0 },
  { -3.0899420510880309, -6.0000000000000000, 0.0 },
  { -3.0078995605434260, -5.7500000000000000, 0.0 },
  { -2.9238567807919029, -5.5000000000000000, 0.0 },
  { -2.8376938981442601, -5.2500000000000000, 0.0 },
  { -2.7492791260608072, -5.0000000000000000, 0.0 },
  { -2.6584669803084551, -4.7500000000000000, 0.0 },
  { -2.5650962220765550, -4.5000000000000000, 0.0 },
  { -2.4689873874722430, -4.2500000000000000, 0.0 },
  { -2.3699397969983651, -4.0000000000000000, 0.0 },
  { -2.2677279046461170, -3.7500000000000000, 0.0 },
  { -2.1620967990779754, -3.5000000000000000, 0.0 },
  { -2.0527566029048092, -3.2500000000000000, 0.0 },
  { -1.9393754207667084, -3.0000000000000000, 0.0 },
  { -1.8215703477370619, -2.7500000000000000, 0.0 },
  { -1.6988958419950144, -2.5000000000000000, 0.0 },
  { -1.5708284488702069, -2.2500000000000000, 0.0 },
  { -1.4367463668836815, -2.0000000000000000, 0.0 },
  { -1.2959015448891098, -1.7500000000000000, 0.0 },
  { -1.1473806603755707, -1.5000000000000000, 0.0 },
  { -0.99004900126010376, -1.2500000000000000, 0.0 },
  { -0.82246703342411320, -1.0000000000000000, 0.0 },
  { -0.64276126883997886, -0.75000000000000000, 0.0 },
  { -0.44841420692364631, -0.50000000000000000, 0.0 },
  { -0.23590029768626336, -0.25000000000000000, 0.0 },
  { 0.0000000000000000, 0.0000000000000000, 0.0 },
  { 0.26765263908273251, 0.25000000000000000, 0.0 },
  { 0.58224052646501256, 0.50000000000000000, 0.0 },
  { 0.97846939293030599, 0.75000000000000000, 0.0 },
  { 1.6449340668482264, 1.0000000000000000, 0.0 },
};
const double toler001 = 2.5000000000000020e-13;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_dilog<Ret> (&data)[Num], Ret toler)
  {
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = emsr::dilog(data[i].x);
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