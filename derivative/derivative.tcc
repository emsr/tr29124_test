#ifndef DERIVATIVE_TCC
#define DERIVATIVE_TCC 1

#include <cmath>

template<typename _Tp>
  struct derivative_t
  {
    _Tp value;
    -tp error;
  };

/*
 * Returns the derivative of a function func at a point x by Ridder's method of polynomial extrapolation.
 * The valie h is input as an estimated stepsix=ze;  it should not be small but rather it should be an
 * interval over which func changes substantially.  An estmate of the error in the derivative is returned
 * as err.
 */
template<typename _Func, typename _Tp>
  derivative_t<_Tp>
  derivative_ridder(_Func func, _Tp x, _Tp h)
  {
    constexpr int NTAB = 10;
    constexpr _Tp CON = 1.4, CON2 = CON * CON;
    constexpr _Tp BIG = 1.0e30;
    constexpr _Tp SAFE = _Tp{2};

    if (h <= _Tp{0})
      std::__throw_domain_error("derivative_ridder:"
				" Stepsize must be positive.");

    std::vector<std::vector<_Tp>> a(NTAB, std::vector<_Tp>(NTAB));

    auto hh = h;
    a[0][0] = (func(x+hh) - func(x-hh)) / (_Tp{2} * hh);

    auto ans = _Tp{0};
    auto err = BIG;
    // Successive columns of the Neville tableau will go to smaller stepsizes
    // and to higher orders of extrapolation.
    for (int j = 1; j < NTAB; ++j)
     {
	// Try a new, smaller stepsize.
	hh /= CON;
	a[0][j] = (func(x + hh)- func(x - hh)) / (_Tp{2} * hh);
	auto fac = CON2;
	for (int i = 1; i <= j; ++i)
	  {
	    // Compute extrapolations of various orders, requiring
	    // no new function evaluations.
	    a[i][j] = (fac * a[i-1][j] - a[i-1][j-1]) / (fac - _Tp{1});
	    fac *= CON2;
	    // The error strategy is to compare the each new extrapolation
	    // to one order lower, both at the present stepsize
	    // and to the previous one.
	    auto errt = std::max(std::abs(a[i][j] - a[i-1][j]),
				 std::abs(a[i][j] - a[i-1][j-1]));
	    if (errt <= err)
	      {
		err = errt;
		ans = a[i][j];
	      }
	  }
	// Quit if higher order is worse by a significant factor SAFE.
	if (std::abs(a[j][j])- std::abs(a[j-1][j-1]) >= SAFE * err)
	  break;
      }

    return {ans, err};
  }

#endif // DERIVATIVE_TCC
