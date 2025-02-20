/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include <emsr/math_constants.h>


constexpr size_t s_num = 11;

template<typename Tp>
  constexpr Tp
  zeta_even[s_num]
  {
                                       -0.5L,
    1.6449340668482264364724151666460251892L,
    1.0823232337111381915160036965411679028L,
    1.0173430619844491397145179297909205279L,
    1.0040773561979443393786852385086524653L,
     1.000994575127818085337145958900319017L,
     1.000246086553308048298637998047739671L,
    1.0000612481350587048292585451051353337L,
     1.000015282259408651871732571487636722L,
    1.0000038172932649998398564616446219397L,
    1.0000009539620338727961131520386834493L,
  };

template<typename Tp>
  constexpr Tp
  eta_even[s_num]
  {
                                	 0.5L,
    0.82246703342411321823620758332301259461L,
    0.94703282949724591757650323447352191493L,
     0.9855510912974351040984392444849542614L,
    0.99623300185264789922728926008280361787L,
    0.99903950759827156563922184569934183142L,
    0.99975768514385819085317967871275542308L,
    0.99993917034597971817095419225539105453L,
    0.99998476421490610644168277496140724092L,
    0.99999618786961011347968922641160768328L,
    0.99999904661158152211505084255772634432L,
  };

template<typename Tp>
  constexpr Tp
  lambda_even[s_num]
  {
                                        0.0L,
    1.2337005501361698273543113749845188919L,
    1.0146780316041920545462534655073449089L,
    1.0014470766409421219064785871379373947L,
    1.0001551790252961193029872492957280416L,
    1.0000170413630448254881839022998304242L,
     1.000001885848583119575908838380247547L,
    1.0000002092405192115001063686802631941L,
    1.0000000232371573791567076732245219815L,
    1.0000000025814375566597728440281148115L,
    1.0000000002868076974555819972982048968L,
  };

template<typename Tp>
  constexpr Tp
  beta_odd[s_num]
  {
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
                                	0.0L,
  };

/**
 * Return the reperiodized cotangent for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The cotangent function is given by:
 * @f[
 *    cot(\pi x) = \frac{1}{\pi x}
 *      - \frac{2}{\pi x}\sum_{k=1}^{\infty}\zeta(2k)x^{2k}
 * @f]
 */
template<typename Tp>
  Tp
  cot_pi_zeta(Tp x)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto pi_x = s_pi * x;
    const auto xx = x * x;
    auto xxk = xx;
    unsigned int k = 2;
    auto res = Tp{0};
    for (unsigned int k = 0; k < s_num; ++k, xxk *= xx)
      res += zeta_even<Tp>[k] * xx;
    return (Tp{1} - Tp{2} * res) / pi_x;
  }

/**
 * Return the reperiodized cosecant for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The cosecant function is given by:
 * @f[
 *    csc(\pi x) = \frac{1}{\pi x}
 *      + \frac{2}{\pi x}\sum_{k=1}^{\infty}\eta(2k)x^{2k}
 * @f]
 */
template<typename Tp>
  Tp
  csc_pi_zeta(Tp x)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto pi_x = s_pi * x;
    const auto xx = x * x;
    auto xxk = xx;
    auto res = Tp{0};
    for (unsigned int k = 0; k < s_num; ++k, xxk *= xx)
      res += eta_even<Tp>[k] * xx;
    return (Tp{1} + Tp{2} * res) / pi_x;
  }

/**
 * Return the reperiodized tangent for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The tangent function is given by:
 * @f[
 *    tan\left(\frac{\pi x}{2}\right) =
 *      -\frac{4}{\pi x}\sum_{k=1}^{\infty}(-1)^k\lambda(2k)x^{2k}
 * @f]
 */
template<typename Tp>
  Tp
  tan_pi_zeta(Tp x)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto pi_x = s_pi * x;
    const auto xx = x * x;
    auto xxk = xx;
    auto res = Tp{0};
    for (unsigned int k = 0; k < s_num; ++k, xxk *= xx)
      res += lambda_even<Tp>[k] * xx;
    return Tp{4} * res / pi_x;
  }

/**
 * Return the reperiodized secant for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The secant function is given by:
 * @f[
 *    sec\left(\frac{\pi x}{2}\right) =
 *      \frac{4}{\pi x}\sum_{k=1}^{\infty}\beta(2k-1)x^{2k-1}
 * @f]
 */
template<typename Tp>
  Tp
  sec_pi_zeta(Tp x)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto pi_x = s_pi * x;
    const auto xx = x * x;
    auto xxk = x;
    unsigned int k = 1;
    auto res = Tp{0};
    for (unsigned int k = 0; k < s_num; ++k, xxk *= xx)
      res += beta_odd<Tp>[k] * xx;
    return Tp{4} * res / pi_x;
  }

/**
 * Return the reperiodized hyperbolic cotangent for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The hyperbolic cotangent function is given by:
 * @f[
 *    coth(\pi x) = \frac{1}{\pi x}
 *      - \frac{2}{\pi x}\sum_{k=1}^{\infty}(-1)^k\zeta(2k)x^{2k}
 * @f]
 */
template<typename Tp>
  Tp
  coth_pi_zeta(Tp x)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto pi_x = s_pi * x;
    const auto xx = -x * x;
    auto xxk = xx;
    unsigned int k = 2;
    auto res = Tp{0};
    for (unsigned int k = 0; k < s_num; ++k, xxk *= xx)
      res += zeta_even<Tp>[k] * xx;
    return (Tp{1} - Tp{2} * res) / pi_x;
  }

/**
 * Return the reperiodized csch for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The csch function is given by:
 * @f[
 *    csch(\pi x) = \frac{1}{\pi x}
 *      + \frac{2}{\pi x}\sum_{k=1}^{\infty}(-1)^k\eta(2k)x^{2k}
 * @f]
 */
template<typename Tp>
  Tp
  csch_pi_zeta(Tp x)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto pi_x = s_pi * x;
    const auto xx = -x * x;
    auto xxk = xx;
    auto res = Tp{0};
    for (unsigned int k = 0; k < s_num; ++k, xxk *= xx)
      res += eta_even<Tp>[k] * xx;
    return (Tp{1} + Tp{2} * res) / pi_x;
  }

/**
 * Return the reperiodized hyperbolic tangent for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The hyperbolic tangent function is given by:
 * @f[
 *    tan\left(\frac{\pi x}{2}\right) =
 *      -\frac{4}{\pi x}\sum_{k=1}^{\infty}(-1)^k\lambda(2k)x^{2k}
 * @f]
 */
template<typename Tp>
  Tp
  tanh_pi_zeta(Tp x)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto pi_x = s_pi * x;
    const auto xx = -x * x;
    auto xxk = xx;
    auto res = Tp{0};
    for (unsigned int k = 0; k < s_num; ++k, xxk *= xx)
      res += lambda_even<Tp>[k] * xx;
    return -Tp{4} * res / pi_x;
  }

/**
 * Return the reperiodized hyperbolic secant function for real argument @c x.
 * This function uses series involving integer argument Riemann zeta function
 * and friends.
 *
 * The hyperbolic secant function is given by:
 * @f[
 *    sech\left(\frac{\pi x}{2}\right) =
 *      -\frac{4}{\pi x}\sum_{k=1}^{\infty}(-1)^k\beta(2k-1)x^{2k-1}
 * @f]
 */
template<typename Tp>
  Tp
  sech_pi_zeta(Tp x)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto pi_x = s_pi * x;
    const auto xx = -x * x;
    auto xxk = x;
    auto res = Tp{0};
    for (unsigned int k = 0; k < s_num; ++k, xxk *= xx)
      res += beta_odd<Tp>[k] * xx;
    return -Tp{4} * res / pi_x;
  }

/**
 * 
 */
template<typename Tp>
  Tp
  cot_pi()
  {
  }

/**
 * 
 */
template<typename Tp>
  Tp
  csc_pi()
  {
  }

/**
 * 
 */
template<typename Tp>
  Tp
  tan_pi()
  {
  }

/**
 * 
 */
template<typename Tp>
  Tp
  sec_pi()
  {
  }

/**
 * 
 */
template<typename Tp>
  Tp
  coth_pi()
  {
  }

/**
 * 
 */
template<typename Tp>
  Tp
  csch_pi()
  {
  }

/**
 * 
 */
template<typename Tp>
  Tp
  tanh_pi()
  {
  }

/**
 * 
 */
template<typename Tp>
  Tp
  sech_pi()
  {
  }

int
main()
{
}
