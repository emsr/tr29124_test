
#include <emsr/special_functions.h>

template<typename Tp>
  int
  test_m_gt_l()
  {
    int num_errors = 0;
    for (auto l : {0u, 1u, 2u, 5u})
      for (auto m : {l + 1u, l + 2u})
	for (auto i : {-2, -1, 0, 1, 2})
	  {
	    auto x = Tp(i * 0.5L);
	    if (emsr::assoc_legendre(l, m, x) != Tp(0)) ++num_errors;
	  }
    return num_errors;
  }

int
main()
{
  int num_errors = 0;
  test_m_gt_l<float>();
  test_m_gt_l<double>();
  test_m_gt_l<long double>();
  return num_errors;
}