
#include <vector>
#include <algorithm>
#include <numeric>

#include <e_float/e_float.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <utility/util_find_root_bisect.h>
#include <utility/util_interpolate.h>
#include <utility/util_point.h>

const INT32& HypergeometricUtil::AsympConvergeMaxOrder(void)
{
  static const double d_max = static_cast<double>(static_cast<double>(std::numeric_limits<e_float>::digits10) * static_cast<double>(0.02));
  static const INT32  n_max = static_cast<INT32>(static_cast<INT64>(d_max));
  static const INT32  n_two = static_cast<INT32>(2);

  return std::max(n_two, n_max);
}

namespace HypergeometricUtil
{
  static double approx_log_n_fact(const INT32 n)
  {
    if(n < static_cast<INT32>(0))
    {
      return static_cast<double>(0.0);
    }
    else
    {
      switch(n)
      {
        case 0:
        case 1:
          return static_cast<double>(0.0);
        case 2:
          { static const double lg_2_fact = ::log(static_cast<double>(2.0)); return lg_2_fact; }
        case 3:
          { static const double lg_3_fact = ::log(static_cast<double>(6.0)); return lg_3_fact; }
        case 4:
          { static const double lg_4_fact = ::log(static_cast<double>(24.0)); return lg_4_fact; }
        case 5:
          { static const double lg_5_fact = ::log(static_cast<double>(120.0)); return lg_5_fact; }
        case 6:
          { static const double lg_6_fact = ::log(static_cast<double>(720.0)); return lg_6_fact; }
        case 7:
          { static const double lg_7_fact = ::log(static_cast<double>(5040.0)); return lg_7_fact; }
        case 8:
          { static const double lg_8_fact = ::log(static_cast<double>(40320.0)); return lg_8_fact; }
        default:
          {
            static const double half_lg_two_pi = ::log(3.14159265358979323 * 2.0) / static_cast<double>(2.0);

            const double x = static_cast<double>(n + static_cast<INT32>(1));

            return (((((x - static_cast<double>(0.5)) * ::log(x)) - x) + half_lg_two_pi) + (1.0 / (12.0 * x))) - (1.0 / (360.0 * ((x * x) * x)));
          }
      }
    }
  }

  struct log_value
  {
    log_value() { }

    double operator()(double& x) const { return ::log(::fabs(x)); }
  };
}

bool HypergeometricUtil::AsympConverge(const std::deque<double>& a,
                                       const std::deque<double>& b,
                                       const double x,
                                       const INT32  prec,
                                       const INT32  maxo)
{
  std::vector<double> ap(a.begin(), a.end());
  std::vector<double> bp(b.begin(), b.end());

  std::vector<double> ap_log(a.size());
  std::vector<double> bp_log(b.size());
  
  const double lgx = ::log(::fabs(x));

  bool bo_converge = true;

  // Take the log of each of the pochhammer elements in a and b
  // and store these in the logarithmic containers.
  std::transform(ap.begin(), ap.end(), ap_log.begin(), HypergeometricUtil::log_value());
  std::transform(bp.begin(), bp.end(), bp_log.begin(), HypergeometricUtil::log_value());

  double ap_log_sum = std::accumulate(ap_log.begin(), ap_log.end(), static_cast<double>(0.0));
  double bp_log_sum = std::accumulate(bp_log.begin(), bp_log.end(), static_cast<double>(0.0));

  static const double one_d = static_cast<double>(1.0);

  for(INT32 n = static_cast<INT32>(2); ; n++)
  {
    // Increment each of the pochhammer elements in a and b.
    std::transform(ap.begin(), ap.end(), ap.begin(), std::bind1st(std::plus<double>(), one_d));
    std::transform(bp.begin(), bp.end(), bp.begin(), std::bind1st(std::plus<double>(), one_d));

    // Take the log of each of the pochhammer elements in a and b
    // and store these in the logarithmic containers.
    std::transform(ap.begin(), ap.end(), ap_log.begin(), HypergeometricUtil::log_value());
    std::transform(bp.begin(), bp.end(), bp_log.begin(), HypergeometricUtil::log_value());

    // Increment the logarithmic sums, which is equivalent to the
    // multiplication of the pochhammer elements.
    ap_log_sum = std::accumulate(ap_log.begin(), ap_log.end(), ap_log_sum);
    bp_log_sum = std::accumulate(bp_log.begin(), bp_log.end(), bp_log_sum);

    const double lg_top  = ap_log_sum + (static_cast<double>(n) * lgx);
    const double lg_bot  = bp_log_sum + HypergeometricUtil::approx_log_n_fact(n);
    const double lg_term = lg_top - lg_bot;

    // Check if the term exceeds the maximum allowed order.
    // Check if the maximum allowed number of iterations have been reached.

    static const INT32 max_iter = static_cast<INT32>(300);
    static const double lg10    = ::log(static_cast<double>(10.0));

    if((lg_term > static_cast<double>(maxo) * lg10) || (n >= max_iter))
    {
      bo_converge = false;
      break;
    }
    
    if(lg_term < -static_cast<double>(prec) * lg10)
    {
      break;
    }
  }

  return bo_converge;
}
