/**
 *
 */

#include <iostream>
#include <iomanip>

#include <emsr/special_functions.h>
#include <emsr/root_search.h>

template<typename Tp>
  Tp
  ibeta_inv(Tp a, Tp b, Tp Ibeta)
  {
    constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
    constexpr int s_max_inner_iter = 12;
    constexpr int s_max_outer_iter = 35;
    if (a < Tp{0} || b < Tp{0})
      throw std::domain_error("ibeta_inv: parameters a and b must be positive");
    else if (Ibeta < Tp{0} || Ibeta > Tp{1})
      throw std::domain_error("ibeta_inv: incomplete beta function out of range");
    const auto _Beta = emsr::beta(a, b);

    auto ibeta = [a, b](Tp x)
		   -> Tp
		   { return emsr::ibeta(a, b, x); };
    auto xl = Tp{0};
    auto xr = Tp{1};
    emsr::root_bracket(ibeta, xl, xr); // This API blows.
/*
    auto Ix = Ibeta;
    auto Iy = Tp{1} - Ix;
    auto Ixy_prev = std::min(Ix, Iy);
    Tp x_prev, y_prev;
    if (a <= b)
      {
	x_prev = std::pow(a * Ix * _Beta, Tp{1} / a);
	y_prev = Tp{1} - x_prev;
      }
    else
      {
	y_prev = std::pow(b * Iy * _Beta, Tp{1} / b);
	x_prev = Tp{1} - y_prev;
      }
    auto xy_prev = std::min(x_prev, y_prev);
    auto xy_lower = Tp{0};
    auto xy_upper = Tp{1};
    int outer_iter = 0;
    while (outer_iter++ < s_max_outer_iter)
      {
	int iter = 0;
	while (iter++ < s_max_inner_iter && Ixy_prev < Ixy)
	  {
	    xy = Tp{2} * xy_prev;
	    xy_lower = xy;
	  }
	iter = 0;
	while (iter++ < s_max_inner_iter && Ixy_prev >= Ixy)
	  {
	    xy = xy_prev / Tp{2};
	    xy_upper = xy;
	  }
	xy = (xy_lower + xy_upper) / Tp{2};
	if (std::abs(xy_upper - xy_lower) < s_eps * xy)
	  break;
      }
*/
    auto thing = [a, b](Tp x){ return emsr::ibeta(a, b, x); };
    auto xy_lower = Tp{0};
    auto xy_upper = Tp{1};
    if (emsr::root_bracket(thing, xy_lower, xy_upper))
      {
	return emsr::root_brent(thing, xy_lower, xy_upper, Tp{10} * s_eps);
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
