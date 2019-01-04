#include <sstream>
#include <iomanip>

#include <functions/complex/e_float_complex.h>
#include <functions/elementary/elementary_trans.h>

// NOCOVER_BLK_BEG
e_float ef_complex::get_abs(void) const
{
  return ef::sqrt(norm());
}

std::string ef_complex::get_str(void) const
{
  const std::streamsize my_prec = static_cast<std::streamsize>(std::numeric_limits<e_float>::digits10);

  std::stringstream ss;

  static_cast<void>(ss.setf(std::ios_base::showpos | std::ios_base::scientific));
  static_cast<void>(ss.precision(my_prec));

  ss << *this;

  return ss.str();  
}
// NOCOVER_BLK_END
