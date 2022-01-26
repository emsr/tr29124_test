/**
 *
 */

#include <iostream>
#include <iomanip>

#include <emsr/specfun.h>
#include <emsr/root_search.h>

template<typename _Tp>
  _Tp
  ibeta_inv(_Tp a, _Tp b, _Tp _Ibeta)
  {
    constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    constexpr int _S_max_inner_iter = 12;
    constexpr int _S_max_outer_iter = 35;
    if (a < _Tp{0} || b < _Tp{0})
      throw std::domain_error("ibeta_inv: parameters a and b must be positive");
    else if (_Ibeta < _Tp{0} || _Ibeta > _Tp{1})
      throw std::domain_error("ibeta_inv: incomplete beta function out of range");
    const auto _Beta = emsr::beta(a, b);

    auto ibeta = [a, b](_Tp x)
		   -> _Tp
		   { return emsr::ibeta(a, b, x); };
    auto xl = _Tp{0};
    auto xr = _Tp{1};
    emsr::root_bracket(ibeta, xl, xr); // This API blows.
/*
    auto _Ix = _Ibeta;
    auto _Iy = _Tp{1} - _Ix;
    auto _Ixy_prev = std::min(_Ix, _Iy);
    _Tp x_prev, y_prev;
    if (a <= b)
      {
	x_prev = std::pow(a * _Ix * _Beta, _Tp{1} / a);
	y_prev = _Tp{1} - x_prev;
      }
    else
      {
	y_prev = std::pow(b * _Iy * _Beta, _Tp{1} / b);
	x_prev = _Tp{1} - y_prev;
      }
    auto xy_prev = std::min(x_prev, y_prev);
    auto xy_lower = _Tp{0};
    auto xy_upper = _Tp{1};
    int outer_iter = 0;
    while (outer_iter++ < _S_max_outer_iter)
      {
	int iter = 0;
	while (iter++ < _S_max_inner_iter && _Ixy_prev < _Ixy)
	  {
	    xy = _Tp{2} * xy_prev;
	    xy_lower = xy;
	  }
	iter = 0;
	while (iter++ < _S_max_inner_iter && _Ixy_prev >= _Ixy)
	  {
	    xy = xy_prev / _Tp{2};
	    xy_upper = xy;
	  }
	xy = (xy_lower + xy_upper) / _Tp{2};
	if (std::abs(xy_upper - xy_lower) < _S_eps * xy)
	  break;
      }
*/
    auto thing = [a, b](_Tp x){ return emsr::ibeta(a, b, x); };
    auto xy_lower = _Tp{0};
    auto xy_upper = _Tp{1};
    if (emsr::root_bracket(thing, xy_lower, xy_upper))
      {
	return emsr::root_brent(thing, xy_lower, xy_upper, _Tp{10} * _S_eps);
      }

    return _Tp{0};
  }

int
main()
{
  double a = 2.4, b = 5.6;
  double x = 0.8;
  double ibet = emsr::ibeta(a, b, x);
  auto y = ibeta_inv(a, b, ibet);
}
