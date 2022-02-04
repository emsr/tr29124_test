/**
 *
 */

#include <cmath>
#include <algorithm> // For clamp.
#include <limits>
#include <stdexcept>
#include <iostream>
#include <iomanip>

#include <emsr/numeric_limits.h>

  //  Evaluates the continued fraction for the incomplete beta function
  //  by the modified Lentz's method
  template<typename Tp>
    Tp
    ibeta_cont_frac(Tp a, Tp b, Tp x)
    {
      constexpr auto _S_itmax = 100;
      const auto _S_fpmin = 1000 * std::numeric_limits<Tp>::min();
      const auto _S_eps = std::numeric_limits<Tp>::epsilon();
      auto apb = a + b;
      auto ap1 = a + Tp{1};
      auto am1 = a - Tp{1};
      auto c = Tp{1};
      //  First step of Lentz's method.
      auto d = Tp{1} - apb * x / ap1;
      if (std::abs(d) < _S_fpmin)
	d = _S_fpmin;
      d = Tp{1} / d;
      auto h = d;
      for (int m = 1; m <= _S_itmax; ++m)
	{
	  auto m2 = 2 * m;

	  //  Even step of the recurrence.
	  auto aa = Tp(m) * (b - Tp(m)) * x
		     / ((am1 + Tp(m2)) * (a + Tp(m2)));
	  d = Tp{1} + aa * d;
	  if (std::abs(d) < _S_fpmin)
	    d = _S_fpmin;
	  c = Tp{1} + aa / c;
	  if (std::abs(c) < _S_fpmin)
	    c = _S_fpmin;
	  d = Tp{1} / d;
	  h *= d * c;

	  //  Odd step of the recurrence.
	  aa = -(a + Tp(m)) * (apb + Tp(m)) * x
		/ ((a + Tp(m2)) * (ap1 + Tp(m2)));
	  d = Tp{1} + aa * d;
	  if (std::abs(d) < _S_fpmin)
	    d = _S_fpmin;
	  c = Tp{1} + aa / c;
	  if (std::abs(c) < _S_fpmin)
	    c = _S_fpmin;
	  d = Tp{1} / d;
	  auto del = d * c;
	  h *= del;

	  if (std::abs(del - Tp{1}) < _S_eps)
	    return h;
	}
      throw std::runtime_error("ibeta_cont_frac: a or b too big, or _S_itmax too small");
    }

  ///  Returns the incomplete beta function I_x(a;b).
  template<typename Tp>
    Tp
    ibeta(Tp a, Tp b, Tp x)
    {
      const auto _S_NaN = emsr::quiet_NaN(x);
      if (x < Tp{0} || x > Tp{1})
	throw std::domain_error("ibeta: argument out of range");
      else if (std::isnan(x) || std::isnan(a) || std::isnan(b))
	return _S_NaN;
      else if (a == Tp{0} && b == Tp{0})
	return _S_NaN;
      else if (a == Tp{0})
	{
	  if (x > Tp{0})
	    return Tp{1};
	  else
	    return Tp{0};
	}
      else if (b == Tp{0})
	{
	  if (x < Tp{1})
	    return Tp{0};
	  else
	    return Tp{1};
	}
      else
	{
	  Tp fact;
	  if (x == Tp{0} || x == Tp{1})
	    fact = Tp{0};
	  else
	    fact = std::exp(std::lgamma(a + b)
			    - std::lgamma(a) - std::lgamma(b)
			    + a * std::log(x) + b * std::log(Tp{1} - x));
	  if (x < (a + Tp{1}) / (a + b + Tp{2}))
	    return fact * ibeta_cont_frac(a, b, x) / a;
	  else
	    return Tp{1}
		 - fact * ibeta_cont_frac(b, a, Tp{1} - x) / b;
	}
    }

template<typename Tp>
  void
  test_ibeta(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 6;

    for (int ia = 0; ia <= 10; ++ia)
      {
	auto a = Tp(ia);
	for (int ib = 10; ib >= 0; --ib)
	  {
	    auto b = Tp(ib);
	    std::cout << "a = " << std::setw(6) << a << '\n';
	    std::cout << "b = " << std::setw(6) << b << '\n';
	    for (int ix = 0; ix <= 100; ++ix)
	      {
		//auto x = ix * Tp{0.01L};
		auto x = std::clamp(ix * Tp{0.01L}, Tp{0}, Tp{1});
		std::cout << ' ' << std::setw(6) << x
		          << ' ' << std::setw(width) << ibeta(a, b, x) << '\n';
	      }
	  }
      }
  }

template<typename Tp>
  void
  stress_test_ibeta(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;

    std::cout << "a = " << std::setw(6) << Tp{0.001L} << '\n';
    std::cout << "b = " << std::setw(6) << Tp{20L} << '\n';
    for (int ix = 0; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.00000000000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.000000000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.0000000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.00000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.0000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.00000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.0001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix <= 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.01L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{0.001L}, Tp{20L}, x) << '\n';
      }

    std::cout << "a = " << std::setw(6) << Tp{20L} << '\n';
    std::cout << "b = " << std::setw(6) << Tp{0.001L} << '\n';
    for (int ix = 0; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * Tp{0.01L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(Tp{0.99L} + ix * Tp{0.0001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(Tp{0.9999L} + ix * Tp{0.000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(Tp{0.999999L} + ix * Tp{0.00000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(Tp{0.99999999L} + ix * Tp{0.0000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(Tp{0.9999999999L} + ix * Tp{0.000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(Tp{0.999999999999L} + ix * Tp{0.00000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(Tp{0.99999999999999L} + ix * Tp{0.0000000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix <= 100; ++ix)
      {
	auto x = std::clamp(Tp{0.9999999999999999L} + ix * Tp{0.000000000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix <= 100; ++ix)
      {
	auto x = std::clamp(Tp{0.999999999999999999L} + ix * Tp{0.00000000000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix <= 100; ++ix)
      {
	auto x = std::clamp(Tp{0.99999999999999999999L} + ix * Tp{0.0000000000000000000001L}, Tp{0}, Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << ibeta(Tp{20L}, Tp{0.001L}, x) << '\n';
      }
  }

int
main()
{
  test_ibeta<long double>();
  stress_test_ibeta<long double>();
}
