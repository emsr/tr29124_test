
// expint fails in expint_En_cont_frac
// for some long double arguments due to low max_iter value

#include <emsr/special_functions.h>

int
test01()
{
  int num_errors = 0;

  // Answers from Wolfram Alpha.
  long double ans_ok = -0.10001943365331651406888645149537315243646135979573L;
  long double ans_bomb = -0.10777727809650077516264612749163100483995270163783L;

  auto Ei_ok = emsr::expint(-1.500001L);
  auto diff_ok = Ei_ok - ans_ok;
  if (std::abs(diff_ok) > 1.0e-15) ++num_errors;

  auto Ei_bomb = emsr::expint(-1.450001L);
  auto diff_bomb = Ei_bomb - ans_bomb;
  if (std::abs(diff_bomb) > 1.0e-15) ++num_errors;

  return num_errors;
}

int
main()
{
  return test01();
}
