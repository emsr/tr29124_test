
#ifndef ULP_H
#define ULP_H 1

#include <limits>
#include <cmath>
#include <algorithm>

namespace emsr
{

  /**
   * Returns the "Unit in the Last Place" or ulp.
   * This is radix independent.
   *
   * For radix @f$ \beta @f$, exponent e, and precision p value
   * @f[
   *    x \in \left[ \beta^{e}, \beta^{e+1} \right)
   * @f]
   * @f[
   *    ulp(x) = \beta^{max(e,e_{min}) - p + 1}
   * @f[
   *
   * @see Handbook of Floating-Point Arithmetic, Muller, et. al.
   *      Birkhauser, 2010, Chapter 2.6, pp. 32-39.
   */
  template<typename Tp>
    constexpr Tp
    ulp(Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (std::isinf(x))
	return std::numeric_limits<Tp>::infinity();
      else
	{
	  const auto exp = std::ilogb(std::abs(x));
	  const auto exp_min = std::numeric_limits<Tp>::min_exponent;
	  const auto prec = std::numeric_limits<Tp>::digits;

	  return std::scalbn(Tp{1}, std::max(exp, exp_min) - prec + 1);
	}
    }

} // namespace emsr

#endif // ULP_H
