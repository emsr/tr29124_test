/**
 *
 */

#include <iostream>
#include <iomanip>

#include <emsr/special_functions.h>
#include <emsr/root_search.h>

template<typename Tp>
  Tp
  ibeta_inv(Tp a, Tp b, Tp _Ibeta)
  {
    constexpr auto _S_eps = std::numeric_limits<Tp>::epsilon();
    constexpr int _S_max_inner_iter = 12;
    constexpr int _S_max_outer_iter = 35;
    if (a < Tp{0} || b < Tp{0})
      throw std::domain_error("ibeta_inv: parameters a and b must be positive");
    else if (_Ibeta < Tp{0} || _Ibeta > Tp{1})
      throw std::domain_error("ibeta_inv: incomplete beta function out of range");
    const auto _Beta = emsr::beta(a, b);

    auto ibeta = [a, b](Tp x)
		   -> Tp
		   { return emsr::ibeta(a, b, x); };
    auto xl = Tp{0};
    auto xr = Tp{1};
    emsr::root_bracket(ibeta, xl, xr); // This API blows.
/*
    auto _Ix = _Ibeta;
    auto _Iy = Tp{1} - _Ix;
    auto _Ixy_prev = std::min(_Ix, _Iy);
    Tp x_prev, y_prev;
    if (a <= b)
      {
	x_prev = std::pow(a * _Ix * _Beta, Tp{1} / a);
	y_prev = Tp{1} - x_prev;
      }
    else
      {
	y_prev = std::pow(b * _Iy * _Beta, Tp{1} / b);
	x_prev = Tp{1} - y_prev;
      }
    auto xy_prev = std::min(x_prev, y_prev);
    auto xy_lower = Tp{0};
    auto xy_upper = Tp{1};
    int outer_iter = 0;
    while (outer_iter++ < _S_max_outer_iter)
      {
	int iter = 0;
	while (iter++ < _S_max_inner_iter && _Ixy_prev < _Ixy)
	  {
	    xy = Tp{2} * xy_prev;
	    xy_lower = xy;
	  }
	iter = 0;
	while (iter++ < _S_max_inner_iter && _Ixy_prev >= _Ixy)
	  {
	    xy = xy_prev / Tp{2};
	    xy_upper = xy;
	  }
	xy = (xy_lower + xy_upper) / Tp{2};
	if (std::abs(xy_upper - xy_lower) < _S_eps * xy)
	  break;
      }
*/
    auto thing = [a, b](Tp x){ return emsr::ibeta(a, b, x); };
    auto xy_lower = Tp{0};
    auto xy_upper = Tp{1};
    if (emsr::root_bracket(thing, xy_lower, xy_upper))
      {
	return emsr::root_brent(thing, xy_lower, xy_upper, Tp{10} * _S_eps);
      }

    return Tp{0};
  }

int
main()
{
  double a = 2.4, b = 5.6;
  double x = 0.8;
  double ibet = emsr::ibeta(a, b, x);
  auto y = ibeta_inv(a, b, ibet);
}
