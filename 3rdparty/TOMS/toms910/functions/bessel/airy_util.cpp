
#include <vector>

#include <functions/bessel/bessel_recursion_order.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/bessel/airy_util.h>

e_float AiryUtil::make_zeta(const e_float& sqrt_x)
{
  return ef::two_third() * ef::pown(sqrt_x, static_cast<INT64>(3));
}

const std::tr1::array<e_float, 2u>& AiryUtil::Three_Pow_1_3and2_3(void)
{
  static const e_float p3_13 = ef::pow(ef::three(), ef::third());
  
  // Create a vector containing 3^(1/3) and 3^(2/3)
  static const std::tr1::array<e_float, 2u> p3_13_23 =
  {{
    p3_13,
    p3_13 * p3_13
  }};

  return p3_13_23;
}

const std::tr1::array<e_float, 2u>& AiryUtil::Ai_and_AiPrime_OfZero(void)
{
  static const std::tr1::array<e_float, 2u> ai_aip =
  {{
    ef::one()       / (AiryUtil::Three_Pow_1_3and2_3()[1] * AiryUtil::Gamma_Of_1_3and2_3()[1]),
    ef::one_minus() / (AiryUtil::Three_Pow_1_3and2_3()[0] * AiryUtil::Gamma_Of_1_3and2_3()[0])
  }};
  
  return ai_aip;
}

const std::tr1::array<e_float, 2u>& AiryUtil::Bi_and_BiPrime_OfZero(void)
{
  static const std::tr1::array<e_float, 2u> bi_bip =
  {{
     ef::sqrt3() * Ai_and_AiPrime_OfZero().front(),
    -ef::sqrt3() * Ai_and_AiPrime_OfZero().back()
  }};

  return bi_bip;
}

const std::tr1::array<e_float, 2u>& AiryUtil::Gamma_Of_1_3and2_3(void)
{
  // Create a vector containing gamma(1/3) and gamma(2/3)

  static std::vector<e_float> g_data;
  
  if(g_data.empty())
  {
    g_data.resize(static_cast<std::size_t>(2u));

    GammaUtil::GammaOfPlusXMinusX(ef::third(), g_data[0], g_data[1]);

    g_data[1] *= -ef::third();
  }

  static const std::tr1::array<e_float, 2u> g_values =
  {{
    g_data.front(),
    g_data.back()
  }};

  return g_values;
}

void AiryUtil::BesselJ_Of_1_3and2_3(const e_float& x, e_float& J1_3, e_float& Jm1_3, e_float& J2_3, e_float& Jm2_3)
{
  // Implement recursion algorithms for Jp1_3, Jm1_3, Jp2_3, Jm2_3 which are used
  // for the expansion of the Airy function in the negative transition region.
  // Use downward recursion in combination with Neumann sums.

  const double xd = ef::to_double(x);

  // Find the recursion order which is necessary for downward recursion and
  // ensure that the value of the order is even numbered.
  const INT32 N0 = BesselRecursionOrder::RecursionStartOrderJ0(   xd);
  const INT32 Nm = static_cast<INT32>(static_cast<INT32>(2) * (N0 / static_cast<INT32>(2)));

  INT32   k                = static_cast<INT32>(Nm / static_cast<INT32>(2));
  e_float one_over_k_fact  = ef::one() / ef::factorial(static_cast<UINT32>(k));
  e_float k_plus_1_3       = ef::third() + k;
  e_float k_plus_2_3       = k_plus_1_3 + ef::third();
  e_float gamma_k_plus_1_3 = ef::gamma(k_plus_1_3);
  e_float gamma_k_plus_2_3 = ef::gamma(k_plus_2_3);

  e_float J1_3_p2 = ef::zero();
  e_float J1_3_p1 = ef::one();

  e_float J2_3_p2 = ef::zero();
  e_float J2_3_p1 = ef::one();

  const e_float third_plus_two_k_over_k_fact  = (ef::third() + static_cast<INT32>(k * static_cast<INT32>(2))) * one_over_k_fact;

  e_float sum_1_3 = third_plus_two_k_over_k_fact * gamma_k_plus_1_3;
  e_float sum_2_3 = third_plus_two_k_over_k_fact * gamma_k_plus_2_3;

  static const e_float eight_thirds = ef::two()   + ef::two_third();
  static const e_float ten_thirds   = ef::three() + ef::third();
  
  const e_float one_over_x = ef::one() / x;

  // Do the downward recursion of Jv:
  //
  //                  Jv+1
  //   Jv = [ 2 (v+1) ---- ] - Jv+2
  //                   x 

  const e_float two_over_x = ef::two() / x;

  for(INT32 m = static_cast<INT32>(Nm - static_cast<INT32>(1)); m >= static_cast<INT32>(0); m--)
  {
    const e_float n_plus_1_3_plus_one = ef::third() + static_cast<INT32>(m + static_cast<INT32>(1));

    // Downward recursion for J_1/3.
    J1_3    = ((n_plus_1_3_plus_one * J1_3_p1) * two_over_x) - J1_3_p2;
    J1_3_p2 = J1_3_p1;
    J1_3_p1 = J1_3;

    const e_float n_plus_2_3_plus_one = ef::two_third() + static_cast<INT32>(m + static_cast<INT32>(1));

    // Downward recursion for J_2/3.
    J2_3    = ((n_plus_2_3_plus_one * J2_3_p1) * two_over_x) - J2_3_p2;
    J2_3_p2 = J2_3_p1;
    J2_3_p1 = J2_3;

    // Do the normalization using a Neumann sum which is:
    // (x/2)^v = Sum_k { [((v + 2k) gamma(v + k)) / k!] * J_v+2k
    if((m % static_cast<INT32>(2)) == static_cast<INT32>(0))
    {
      one_over_k_fact  *= k--;
      gamma_k_plus_1_3 /= --k_plus_1_3;
      gamma_k_plus_2_3 /= --k_plus_2_3;

      // Increment the normalization sums for J_1/3 and J_2/3.
      sum_1_3 += (((ef::third()     + m) * gamma_k_plus_1_3) * J1_3) * one_over_k_fact;
      sum_2_3 += (((ef::two_third() + m) * gamma_k_plus_2_3) * J2_3) * one_over_k_fact;
    }
  }

  const e_float x_half_pow_1_3 = ef::pow(x / static_cast<INT32>(2), ef::third());

  // Compute the normalizations using the Neumann sums in combination with (x/2)^v.
  const e_float norm_1_3 =  x_half_pow_1_3 / sum_1_3;
  const e_float norm_2_3 = (x_half_pow_1_3 * x_half_pow_1_3) / sum_2_3;

  // Normalize J_1/3 and J_2/3.
  J1_3 *= norm_1_3;
  J2_3 *= norm_2_3;

  const e_float three_x = x * static_cast<INT32>(3);

  // Recur down one single time in order to obtain J-1/3 from J2/3 and J5/3,
  // as well as to obtain J-2/3 from J1/3 and J4/3.
  Jm1_3 = ((J2_3 * static_cast<INT32>(4)) / three_x) - (J2_3_p2 * norm_2_3);
  Jm2_3 = ((J1_3 * static_cast<INT32>(2)) / three_x) - (J1_3_p2 * norm_1_3);
}
