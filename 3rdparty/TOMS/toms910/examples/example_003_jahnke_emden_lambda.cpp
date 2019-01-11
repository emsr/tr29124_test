#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

#include <e_float/e_float.h>
#include <functions/functions.h>
#include <examples/examples.h>
#include <utility/util_coefficient_expansion.h>

namespace examples
{
  namespace nr_003
  {
    struct divides_n : std::binary_function<const e_float&, INT32, e_float> 
    {
      result_type operator() (first_argument_type a, second_argument_type b) { return a / b; }
    };

    namespace JahnkeEmden_Series
    {
      static e_float AtZero(const e_float& v, const e_float& x)
      {
        static const std::tr1::array<INT32, 8u> A047053 = {{ 1, -4, 32, -384, 6144, -122880, 2949120, -82575360 }};

        const std::tr1::array<e_float,8u> coef_nv =
        {{
          ef::one(),
          ef::one()   / (1 + v),
          coef_nv[1u] / (2 + v),
          coef_nv[2u] / (3 + v),
          coef_nv[3u] / (4 + v),
          coef_nv[4u] / (5 + v),
          coef_nv[5u] / (6 + v),
          coef_nv[6u] / (7 + v)
        }};

        std::vector<e_float> coef(coef_nv.size());

        std::transform(coef_nv.begin(),
                       coef_nv.end(),
                       A047053.begin(),
                       coef.begin(),
                       divides_n());

        return std::accumulate(coef.begin(),
                               coef.end(),
                               ef::zero(),
                               Util::coefficient_expansion<e_float>(x * x));
      }
    }
  }
}

e_float examples::nr_003::jahnke_emden_lambda(const e_float& v, const e_float& x)
{
  if(!ef::small_arg(x))
  {
    const e_float gv = ef::gamma(ef::one() + v);
    const e_float pv = ef::pow  (ef::half(), v);

    return (gv * ef::cyl_bessel_j(v, x)) / pv;
  }
  else
  {
    return JahnkeEmden_Series::AtZero(v, x);
  }
}
