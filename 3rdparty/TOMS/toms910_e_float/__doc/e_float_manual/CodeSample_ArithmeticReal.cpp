#include <iostream>
#include <iomanip>
#include <e_float/e_float.h>

// Compute the sin(x) for x near zero.
e_float sin_small_arg(const e_float& x)
{
  const e_float x_squared            = x * x;
        e_float sum                  = x;
        e_float k_fact               = ef::one();
        e_float x_pow_two_k_plus_one = x;
   
  bool b_neg_term = true;

  // Do the Taylor series expansion.
  for(INT32 k = 3; k < ef::max_iteration(); k += 2)
  {
    k_fact *= static_cast<INT32>((k - 1) * k);
    x_pow_two_k_plus_one *= x_squared;

    const e_float term = x_pow_two_k_plus_one / k_fact;

    if((term.order() - sum.order()) < -ef::tol())
    {
      // The required precision has been reached.
      break;
    }

    !b_neg_term ? sum += term : sum -= term;
    b_neg_term = !b_neg_term;
  }

  return sum;
}

int main(void)
{
  // Calculate the value of sin(1/10).
  static const e_float x = ef::tenth();
  static const e_float y = sin_small_arg(x);

  // Set the output stream precision and print the result.
  std::cout << std::setprecision(std::numeric_limits<e_float>::digits10)
            << y << std::endl;
}
