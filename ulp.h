
#include <limits>
#include <cmath>
#include <algorithm>

/**
 * Returns the "Unit in the Last Place" or ulp.
 * This is radix inedpendent.
 *
 * For radix @f$ \beta @f$, exponent e, and precision p value
 * @f[
 *    x \in \left[ \beta^{e}, \beta^{e+1} \right)
 * @f]
 * @f[
 *    ulp(x) = \beta^{max(e,e_{min}) - p + 1}
 * @f[
 */
template<typename _Tp>
  constexpr _Tp
  ulp(_Tp __x)
  {
    if (std::isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (std::isinf(__x))
      return std::numeric_limits<_Tp>::infinity();
    else
      {
	const auto __exp = std::ilogb(__t);
	const auto __exp_min = std::numeric_limits<_Tp>::min_exponent;
	const auto _prec = std::numeric_limits<_Tp>::digits;

	return std::scalbn(_Tp{1}, std::max(__exp, __exp_min) - _prec + 1);
      }
  }
