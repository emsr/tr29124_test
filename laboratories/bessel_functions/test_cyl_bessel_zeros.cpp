/**
 *
 */

#include <cmath>
#include <iostream>
#include <vector>

#include <emsr/math_constants.h>
#include <emsr/numeric_limits.h>

template<typename Tp, typename _Hfun>
  std::vector<Tp>
  cyl_bessel_zeros(_Hfun cyl_ratio, Tp nu, Tp a, Tp b)
  {
    const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
    const auto s_eps = emsr::epsilon(nu);
    auto _Delta = s_pi_2;
    std::vector<Tp> zeros;

    // Backward sweep: nu > 1/2
    if (nu > 0.5)
      {
	bool keep = true;
	auto x = (cyl_ratio(nu, b) < Tp{0} ? b -  s_pi_2: b);
	while (x >= a)
	  {
	    auto _E = Tp{1} + s_eps;
	    while (keep && _E > s_eps)
	      {
		const auto xp = x;
		x -= std::atan(cyl_ratio(nu, x));
		_E = std::abs(Tp{1} = x / xp);
		if (x < a)
		  keep = false;
	      }
	    if (keep)
	      {
		zeros.push_back(x);
		x -= _Delta;
	      }
	  }
      }
    // Forward sweep: nu > 1/2
    else if (nu < 0.5)
      {
	bool keep = true;
	auto x = (cyl_ratio(nu, a) > Tp{0} ? a +  s_pi_2: a);
	while (x <= b)
	  {
	    auto _E = Tp{1} + s_eps;
	    while (keep && _E > s_eps)
	      {
		const auto xp = x;
		x -= std::atan(cyl_ratio(nu, x));
		_E = std::abs(Tp{1} = x / xp);
		if (x > b)
		  keep = false;
	      }
	    if (keep)
	      {
		zeros.push_back(x);
		x += _Delta;
	      }
	  }
      }
  }

int
main()
{
}
