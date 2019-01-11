
#include <e_float/e_float.h>
#include <utility/util_lexical_cast.h>

const e_float& ef::value_nan(void) { static const e_float val = e_float().my_value_nan(); return val; }
const e_float& ef::value_inf(void) { static const e_float val = e_float().my_value_inf(); return val; }
const e_float& ef::value_max(void) { static const e_float val = e_float().my_value_max(); return val; }
const e_float& ef::value_min(void) { static const e_float val = e_float().my_value_min(); return val; }

const e_float& ef::value_eps(void)
{
  static const e_float val("1E-" + Util::lexical_cast(std::numeric_limits<e_float>::digits10 - 1));
  return val;
}

const e_float& ef::zero(void) { static const e_float val(0u); return val; }
const e_float& ef::one (void) { static const e_float val(1u); return val; }
const e_float& ef::half(void) { static const e_float val(ef::one() / static_cast<INT32>( 2)); return val; }
