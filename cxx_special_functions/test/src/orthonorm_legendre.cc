
#include <iostream>

#include <emsr/integration.h>
#include <emsr/sf_legendre.h>

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename Tp>
  Tp
  norm_legendre(int l1, int l2, Tp x)
  {
    return (Tp(l1 + l2 + 1) / Tp{2})
	 * emsr::legendre(l1,x)
	 * emsr::legendre(l2,x);
  }

template<typename Tp>
  Tp
  delta(int l1, int l2)
  { return l1 == l2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  int
  test_legendre()
  {
    const auto eps_factor = 1 << (std::numeric_limits<Tp>::digits / 3);
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = Tp{10} * rel_prec;

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int num_errors = 0;

    for (const auto l1 : degree)
      {
	for (const auto l2 : degree)
	  {
	    auto func = [l1, l2](Tp x)
			-> Tp
			{ return norm_legendre(l1, l2, x); };

	    auto [result, error]
		= emsr::integrate_tanh_sinh(func, Tp{-1}, Tp{1}, abs_prec, rel_prec, 6);

            auto del = delta<Tp>(l1, l2) - result;
	    if (std::abs(del) > cmp_prec)
              {
		++num_errors;
        	std::cout << "l1 = " << l1 << "; l2 = " << l2
                          << "; res = " << result << "; del = " << del << '\n';
              }
	  }
      }
    return num_errors;
  }

int
main()
{
  int num_errors = 0;
  num_errors += test_legendre<float>();
  num_errors += test_legendre<double>();
  num_errors += test_legendre<long double>();
  return num_errors;
}
