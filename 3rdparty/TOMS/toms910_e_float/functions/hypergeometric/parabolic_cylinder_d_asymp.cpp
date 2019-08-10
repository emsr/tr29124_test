
#include <deque>
#include <numeric>

#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/parabolic_cylinder_d_asymp.h>
#include <functions/tables/tables.h>

e_float ParabolicCylinder_Asymp::WeberU_TemmeSiamBook_Eq_12_119_Asymp(const e_float& a, const e_float& x)
{
  const e_float x_half                     = x / static_cast<INT32>(2);
  const e_float sqrt_x_half_squared_plus_a = ef::sqrt((x_half * x_half) + a);
  const e_float tau                        = ((x_half / sqrt_x_half_squared_plus_a) - ef::one()) / static_cast<INT32>(2);

  std::deque<e_float> tau_powers(static_cast<std::size_t>(1u), ef::one());

  const e_float one_over_minus_two_a       = ef::one() / (a * static_cast<INT32>(-2));
        e_float one_over_minus_two_a_pow_s = ef::one();

  e_float Fa = ef::one();

  for(std::size_t s = static_cast<std::size_t>(1u); s < Tables::A001164().size(); s++)
  {
    while(tau_powers.size() < Tables::A158503()[s]().size())
    {
      tau_powers.push_back(tau_powers.back() * tau);
    }

    one_over_minus_two_a_pow_s *= one_over_minus_two_a;

    const e_float phi_s = (tau_powers[s] * std::inner_product(tau_powers.begin(),
                                                              tau_powers.end(),
                                                              Tables::A158503()[s]().begin(),
                                                              ef::zero())) / Tables::A001164()[s]();

    const e_float term = phi_s * one_over_minus_two_a_pow_s;

    if(term.order() < -ef::tol())
    {
      break;
    }

    Fa += term;
  }

  const e_float Uax_tilde = Fa / ef::sqrt(sqrt_x_half_squared_plus_a * static_cast<INT32>(2));

  const e_float Fax = ef::pow(x_half + sqrt_x_half_squared_plus_a, a) * ef::exp((x_half * sqrt_x_half_squared_plus_a) - (a / static_cast<INT32>(2)));

  return Uax_tilde / Fax;
}
