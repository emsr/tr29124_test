
#include <emsr/integration.h>
#include <emsr/sf_gegenbauer.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_gegenbauer(int n1, int n2, _Tp lambda, _Tp x)
  {
    const auto
    _S_pi = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    auto gama = std::tgamma(lambda);
    auto gamn2a = std::tgamma(n1 + _Tp{2} * lambda);
    auto norm = _S_pi * std::pow(_Tp{2}, _Tp{1} - _Tp{2} * lambda) * gamn2a
	      / emsr::factorial<_Tp>(n1) / (_Tp(n1) + lambda) / gama / gama;
    return std::pow(_Tp{1} - x * x, lambda - _Tp{0.5})
	 * emsr::gegenbauer(n1, lambda, x)
	 * emsr::gegenbauer(n2, lambda, x) / norm;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  int
  test_gegenbauer(_Tp lambda)
  {
    const auto eps_factor = 1 << (std::numeric_limits<_Tp>::digits / 3);
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto abs_prec = eps_factor * eps;
    const auto rel_prec = eps_factor * eps;
    const auto cmp_prec = _Tp{10} * rel_prec;

    const bool singular = (lambda < _Tp{0.5});

    const std::array<int, 10> degree{{0, 1, 2, 4, 8, 16, 32, 64, 81, 128}};
    int num_errors = 0;

    for (const auto n1 : degree)
      {
	for (const auto n2 : degree)
	  {
	    auto func = [n1, n2, lambda](_Tp x)
			-> _Tp
			{ return norm_gegenbauer<_Tp>(n1, n2, lambda, x); };

	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func,
					_Tp{-1}, _Tp{1},
					lambda - _Tp{0.5}, lambda - _Tp{0.5},
					0, 0,
					abs_prec, rel_prec)
		: emsr::integrate_tanh_sinh(func, _Tp{-1}, _Tp{1}, abs_prec, rel_prec);

	    if (std::abs(delta<_Tp>(n1, n2) - result) > cmp_prec)
	      ++num_errors;
	  }
      }
    return num_errors;
  }

int
main()
{
  int num_errors = 0;
  num_errors += test_gegenbauer<float>(0.5F);
  num_errors += test_gegenbauer<double>(0.5);
  num_errors += test_gegenbauer<long double>(0.5L);
  return num_errors;
}
