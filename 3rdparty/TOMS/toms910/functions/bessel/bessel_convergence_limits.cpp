
#include <fstream>
#include <vector>

#include <e_float/e_float.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <utility/util_point.h>
#include <utility/util_interpolate.h>

double BesselConvergenceLimits::UpperLimitSmallX_Iv(const double v)
{
  return UpperLimitSmallX_Jv(v);
}

double BesselConvergenceLimits::LowerLimitLargeX_Iv(const double v)
{
  static const std::tr1::array<Util::point<double>, 19u> limit_data =
  {{
    Util::point<double>(static_cast<double>(    2.0), static_cast<double>(     400 * 1.3)),
    Util::point<double>(static_cast<double>(    4.0), static_cast<double>(     400 * 1.3)),
    Util::point<double>(static_cast<double>(    8.0), static_cast<double>(     400 * 1.3)),
    Util::point<double>(static_cast<double>(   10.0), static_cast<double>(     400 * 1.3)),
    Util::point<double>(static_cast<double>(   20.0), static_cast<double>(     400 * 1.3)),
    Util::point<double>(static_cast<double>(   50.0), static_cast<double>(     400 * 1.3)),
    Util::point<double>(static_cast<double>(  100.0), static_cast<double>(     500 * 1.3)),
    Util::point<double>(static_cast<double>(  200.0), static_cast<double>(     900 * 1.3)),
    Util::point<double>(static_cast<double>(  300.0), static_cast<double>(    1900 * 1.3)),
    Util::point<double>(static_cast<double>(  400.0), static_cast<double>(    3000 * 1.3)),
    Util::point<double>(static_cast<double>(  500.0), static_cast<double>(    5000 * 1.3)),
    Util::point<double>(static_cast<double>( 1000.0), static_cast<double>(   22000 * 1.3)),
    Util::point<double>(static_cast<double>( 2000.0), static_cast<double>(   85000 * 1.3)),
    Util::point<double>(static_cast<double>( 3000.0), static_cast<double>(  180000 * 1.3)),
    Util::point<double>(static_cast<double>( 4000.0), static_cast<double>(  340000 * 1.3)),
    Util::point<double>(static_cast<double>( 5000.0), static_cast<double>(  520000 * 1.3)),
    Util::point<double>(static_cast<double>(10000.0), static_cast<double>( 2000000 * 1.3)),
    Util::point<double>(static_cast<double>(20000.0), static_cast<double>( 8500000 * 1.3)),
    Util::point<double>(static_cast<double>(30000.0), static_cast<double>(19000000 * 1.3)),
  }};

  static const std::vector<Util::point<double> > limit_values(limit_data.begin(), limit_data.end());

  return Util::linear_interpolate<double>::interpolate(v, limit_values);
}

double BesselConvergenceLimits::UpperLimitSmallX_Jv(const double v)
{
  static const std::tr1::array<Util::point<double>, 19u> limit_data =
  {{
    Util::point<double>(static_cast<double>(    2.0), static_cast<double>(  24 * 0.8)),
    Util::point<double>(static_cast<double>(    4.0), static_cast<double>(  24 * 0.8)),
    Util::point<double>(static_cast<double>(    8.0), static_cast<double>(  24 * 0.8)),
    Util::point<double>(static_cast<double>(   10.0), static_cast<double>(  24 * 0.8)),
    Util::point<double>(static_cast<double>(   20.0), static_cast<double>(  32 * 0.8)),
    Util::point<double>(static_cast<double>(   50.0), static_cast<double>(  50 * 0.8)),
    Util::point<double>(static_cast<double>(  100.0), static_cast<double>(  70 * 0.8)),
    Util::point<double>(static_cast<double>(  200.0), static_cast<double>( 100 * 0.8)),
    Util::point<double>(static_cast<double>(  300.0), static_cast<double>( 125 * 0.8)),
    Util::point<double>(static_cast<double>(  400.0), static_cast<double>( 145 * 0.8)),
    Util::point<double>(static_cast<double>(  500.0), static_cast<double>( 170 * 0.8)),
    Util::point<double>(static_cast<double>( 1000.0), static_cast<double>( 230 * 0.8)),
    Util::point<double>(static_cast<double>( 2000.0), static_cast<double>( 300 * 0.8)),
    Util::point<double>(static_cast<double>( 3000.0), static_cast<double>( 390 * 0.8)),
    Util::point<double>(static_cast<double>( 4000.0), static_cast<double>( 460 * 0.8)),
    Util::point<double>(static_cast<double>( 5000.0), static_cast<double>( 510 * 0.8)),
    Util::point<double>(static_cast<double>(10000.0), static_cast<double>( 700 * 0.8)),
    Util::point<double>(static_cast<double>(20000.0), static_cast<double>( 900 * 0.8)),
    Util::point<double>(static_cast<double>(30000.0), static_cast<double>(1200 * 0.8)),
  }};

  static const std::vector<Util::point<double> > limit_values(limit_data.begin(), limit_data.end());

  return Util::linear_interpolate<double>::interpolate(v, limit_values);
}

double BesselConvergenceLimits::LowerLimitLargeX_Jv(const double v)
{
  static const std::tr1::array<Util::point<double>, 19u> limit_data =
  {{
    Util::point<double>(static_cast<double>(    2.0), static_cast<double>(    400 * 1.3)),
    Util::point<double>(static_cast<double>(    4.0), static_cast<double>(    400 * 1.3)),
    Util::point<double>(static_cast<double>(    8.0), static_cast<double>(    400 * 1.3)),
    Util::point<double>(static_cast<double>(   10.0), static_cast<double>(    400 * 1.3)),
    Util::point<double>(static_cast<double>(   20.0), static_cast<double>(    400 * 1.3)),
    Util::point<double>(static_cast<double>(   50.0), static_cast<double>(    400 * 1.3)),
    Util::point<double>(static_cast<double>(  100.0), static_cast<double>(    500 * 1.3)),
    Util::point<double>(static_cast<double>(  200.0), static_cast<double>(    600 * 1.3)),
    Util::point<double>(static_cast<double>(  300.0), static_cast<double>(   1000 * 1.3)),
    Util::point<double>(static_cast<double>(  400.0), static_cast<double>(   1500 * 1.3)),
    Util::point<double>(static_cast<double>(  500.0), static_cast<double>(   2400 * 1.3)),
    Util::point<double>(static_cast<double>( 1000.0), static_cast<double>(  20000 * 1.3)),
    Util::point<double>(static_cast<double>( 2000.0), static_cast<double>(  40000 * 1.3)),
    Util::point<double>(static_cast<double>( 3000.0), static_cast<double>(  90000 * 1.3)),
    Util::point<double>(static_cast<double>( 4000.0), static_cast<double>( 150000 * 1.3)),
    Util::point<double>(static_cast<double>( 5000.0), static_cast<double>( 230000 * 1.3)),
    Util::point<double>(static_cast<double>(10000.0), static_cast<double>(1000000 * 1.3)),
    Util::point<double>(static_cast<double>(20000.0), static_cast<double>(3500000 * 1.3)),
    Util::point<double>(static_cast<double>(30000.0), static_cast<double>(9000000 * 1.3)),
  }};

  static const std::vector<Util::point<double> > limit_values(limit_data.begin(), limit_data.end());

  return Util::linear_interpolate<double>::interpolate(v, limit_values);
}

double BesselConvergenceLimits::UpperLimitSmallX_Kv(const double v)
{
  static const std::tr1::array<Util::point<double>, 19u> limit_data =
  {{
    Util::point<double>(static_cast<double>(    2.0), static_cast<double>(  12 * 0.8)),
    Util::point<double>(static_cast<double>(    4.0), static_cast<double>(  12 * 0.8)),
    Util::point<double>(static_cast<double>(    8.0), static_cast<double>(  15 * 0.8)),
    Util::point<double>(static_cast<double>(   10.0), static_cast<double>(  15 * 0.8)),
    Util::point<double>(static_cast<double>(   20.0), static_cast<double>(  22 * 0.8)),
    Util::point<double>(static_cast<double>(   50.0), static_cast<double>(  40 * 0.8)),
    Util::point<double>(static_cast<double>(  100.0), static_cast<double>(  73 * 0.8)),
    Util::point<double>(static_cast<double>(  200.0), static_cast<double>( 140 * 0.8)),
    Util::point<double>(static_cast<double>(  300.0), static_cast<double>( 170 * 0.8)),
    Util::point<double>(static_cast<double>(  400.0), static_cast<double>( 200 * 0.8)),
    Util::point<double>(static_cast<double>(  500.0), static_cast<double>( 228 * 0.8)),
    Util::point<double>(static_cast<double>( 1000.0), static_cast<double>( 312 * 0.8)),
    Util::point<double>(static_cast<double>( 2000.0), static_cast<double>( 448 * 0.8)),
    Util::point<double>(static_cast<double>( 3000.0), static_cast<double>( 536 * 0.8)),
    Util::point<double>(static_cast<double>( 4000.0), static_cast<double>( 660 * 0.8)),
    Util::point<double>(static_cast<double>( 5000.0), static_cast<double>( 702 * 0.8)),
    Util::point<double>(static_cast<double>(10000.0), static_cast<double>( 996 * 0.8)),
    Util::point<double>(static_cast<double>(20000.0), static_cast<double>(1370 * 0.8)),
    Util::point<double>(static_cast<double>(30000.0), static_cast<double>(1700 * 0.8)),
  }};

  static const std::vector<Util::point<double> > limit_values(limit_data.begin(), limit_data.end());

  return Util::linear_interpolate<double>::interpolate(v, limit_values);
}

double BesselConvergenceLimits::LowerLimitLargeX_Kv(const double v)
{
  return LowerLimitLargeX_Jv(v);
}

double BesselConvergenceLimits::UpperLimitSmallX_Yv(const double v)
{
  return UpperLimitSmallX_Jv(v);
}

double BesselConvergenceLimits::LowerLimitLargeX_Yv(const double v)
{
  return LowerLimitLargeX_Jv(v);
}
