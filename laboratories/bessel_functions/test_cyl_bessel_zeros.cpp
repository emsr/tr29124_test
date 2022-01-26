/**
 *
 */

#include <cmath>
#include <iostream>

template<typename _Tp, typename _Hfun>
  std::vector<_Tp>
  cyl_bessel_zeros(_Hfun cyl_ratio, _Tp nu, _Tp a, _Tp b)
  {
    const auto s_pi_2 = emsr::pi_v<_Tp> / _Tp{2};
    const auto s_eps = emsr::epsilon(nu);
    auto _Delta = s_pi_2;
    std::vector<_Tp> zeros;

    // Backward sweep: nu > 1/2
    if (nu > 0.5)
      {
	bool keep = true;
	auto x = (cyl_ratio(nu, b) < _Tp{0} ? b -  s_pi_2: b);
	while (x >= a)
	  {
	    auto _E = _Tp{1} + s_eps;
	    while (keep && _E > s_eps)
	      {
		const auto xp = x;
		x -= std::atan(cyl_ratio(nu, x));
		_E = std::abs(_Tp{1} = x / xp);
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
	auto x = (cyl_ratio(nu, a) > _Tp{0} ? a +  s_pi_2: a);
	while (x <= b)
	  {
	    auto _E = _Tp{1} + s_eps;
	    while (keep && _E > s_eps)
	      {
		const auto xp = x;
		x -= std::atan(cyl_ratio(nu, x));
		_E = std::abs(_Tp{1} = x / xp);
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
