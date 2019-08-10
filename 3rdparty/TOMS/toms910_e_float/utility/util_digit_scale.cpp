#include <vector>

#include <e_float/e_float.h>

#include <utility/util_digit_scale.h>
#include <utility/util_interpolate.h>

const double& Util::DigitScale(void)
{
  static const std::tr1::array<Util::point<double>, static_cast<std::size_t>(4u)> scale_data =
  {{
    Util::point<double>(static_cast<double>( 50.0), static_cast<double>(1.0 / 6.0)),
    Util::point<double>(static_cast<double>(100.0), static_cast<double>(1.0 / 3.0)),
    Util::point<double>(static_cast<double>(200.0), static_cast<double>(1.0 / 2.0)),
    Util::point<double>(static_cast<double>(300.0), static_cast<double>(1.0)),
  }};

  static const std::vector<Util::point<double> > scale(scale_data.begin(), scale_data.end());

  static const double the_scale = Util::linear_interpolate<double>::interpolate(static_cast<double>(std::numeric_limits<e_float>::digits10), scale);

  return the_scale;
}
