/**
 *
 */

#include <cmath>
#include <limits>

#undef STANDALONE
#include "../error_functions/test_inv_erf.cpp"

#include <emsr/math_constants.h>
#include <emsr/sf_gamma.h>

// A fake Binet function.
template<typename Tp>
  Tp
  lgamma_scaled(Tp a)
  {
    constexpr auto s_ln2pi
      = emsr::ln2_v<Tp>
      + emsr::lnpi_v<Tp>;
    constexpr auto half = Tp{1} / Tp{2};
    return std::lgamma(a)
	 - (a - half) * std::log(a) + a - half * s_ln2pi;
  }

// A fake Binet function.
template<typename Tp>
  Tp
  tgamma_scaled(Tp a)
  {
    constexpr auto s_sqrt2pi
      = emsr::sqrt2_v<Tp>
      + emsr::sqrtpi_v<Tp>;
    return std::tgamma(a) * std::pow(a, -a + Tp{0.5}) * std::exp(a)
	 / s_sqrt2pi;
  }

/**
 * Get lambda from eta by Newton's method.
 */
template<typename Tp>
  Tp
  lambda(Tp eta)
  {
    auto func = [eta](Tp lambdac)
		-> Tp
		{
		  auto sum = Tp{1 / Tp{10}};
		  for (int n = 9; n >= 2; --n)
		    sum += Tp((n & 1) ? -1 : +1) / Tp(n) + lambdac * sum;
		  return lambdac * lambdac * sum - eta * eta / Tp{2};
		};
    auto funcp = [eta](Tp lambdac)
		-> Tp
		{
		  //auto sum = Tp{1};
		  //for (int n = 9; n >= 2; --n)
		  //  sum += Tp((n & 1) ? -1 : +1) + lambdac * sum;
		  //return lambdac * sum;
		  return lambdac / (Tp{1} + lambdac);
		};

    auto lambdac = eta;
    for (int i = 0; i < 50; ++i)
      lambdac -= func(eta) / funcp(eta);

    return Tp{1} + lambdac;
  }


template<typename Tp>
  Tp
  epsilon1(Tp eta)
  {
    return Tp{-1.0L} / Tp{3.0L}
	 + eta * (Tp{1.0L} / Tp{36.0L}
 	 + eta * (Tp{1.0L} / Tp{1'620.0L}
	 + eta * (Tp{-7.0L} / Tp{6'480.0L}
	 + eta * (Tp{5.0L} / Tp{18'144.0L}
	 + eta * Tp{-11.0L} / Tp{382'725.0L}))));
  }

template<typename Tp>
  Tp
  epsilon2(Tp eta)
  {
    return Tp{-7.0L} / Tp{405.0L}
	 + eta * (Tp{-7.0L} / Tp{2'592.0L}
	 + eta * (Tp{533.0L} / Tp{204'120.0L}
	 + eta * (Tp{-1'579.0L} / Tp{2'099'520.0L}
	 + eta * (Tp{109.0L} / Tp{1'749'600.0L}))));
  }

template<typename Tp>
  Tp
  epsilon3(Tp eta)
  {
    return Tp{449.0L} / Tp{102'060.0L}
	 + eta * (Tp{-63'149.0L} / Tp{20'995'200.0L}
	 + eta * (Tp{29'233.0L} /Tp{36'741'600.0L} 
	 + eta * (Tp{346'793.0L} / Tp{5'290'790'400.0L})));
  }

template<typename Tp>
  Tp
  epsilon4(Tp eta)
  {
    return Tp{319.0L} / Tp{183'708.0L}
	 + eta * (Tp{-269'383.0L} / Tp{4'232'632'320.0L}
	 + eta * (Tp{-449'882'243.0L} / Tp{982'102'968'000.0L}));
  }

/**
 * Return the inverse of the incomplete gamma functions equations
 *  P(a,xr) = p and Q(a,xr) = q
 * for positive a and the two incomplete gamma functions.
 * In most cases, we invert the equation with min(p,q)
 *
 * This material is described in Numerical Methods for special functions
 * Amparo Gil, Javier Segura, Nico M. Temme, Chapter 10.
 *
 * @param[in] a Argument of the functions.
 * @param[in] p The function value P(a,x).
 * @param[in] q The function value Q(a,x).
 * @return The argument x.
 */
template<typename Tp>
  Tp
  inv_gamma(Tp a, Tp p, Tp q)
  {
    constexpr auto s_giant = std::numeric_limits<Tp>::max() / Tp{10};
    constexpr auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
    constexpr auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;
    constexpr auto s_sqrt_2pi = s_sqrt_2 * s_sqrt_pi;
    bool pcase;
    Tp porq, s;
    if (p < 0.5)
      {
	pcase = true;
	porq = p;
	s = -1;
      }
    else
      {
	pcase = false;
	porq = q;
	s = 1;
      }

    const auto ra = Tp{1} / a;
    const auto ap1 = a + Tp{1};
    Tp ck[5];
    Tp eta{};
    Tp x0{};
    int m{};
    const auto logr = ra * (std::log(p) + std::lgamma(ap1));
    if (porq == Tp{0})
      return Tp(s) * std::numeric_limits<Tp>::infinity();
    else if (logr < std::log(0.2 * ap1))
      {
	const auto r = std::exp(logr);
	m = 0;
	const auto a2 = a * a;
	const auto a3 = a * a2;
	const auto a4 = a * a3;
	const auto ap12 = ap1 * ap1;
	const auto ap13 = ap1 * ap12;
	const auto ap14 = ap12 * ap12;
	const auto ap2 = a + Tp{2};
	const auto ap22 = ap2 * ap2;
	ck[0] = Tp{1};
	ck[1] = Tp{1} / ap1;
	ck[2] = Tp{0.5L} * (3 * a + 5) / (ap12 * (a + 2));
	ck[3] = (Tp{1} / Tp{3})
	      * (31 + 8 * a2 + 33 * a)
	      / (ap13 * ap2 * (a + 3));
	ck[4] = (Tp{1} / Tp{24})
	      * (2888 + 1179 * a3 + 125 * a4 + 3971 * a2 + 5661 * a)
	      / (ap14 * ap22 * (a + 3) * (a + 4));
	x0 = r * (1 + r * (ck[1] + r * (ck[2] + r * (ck[3] + r * ck[4]))));
      }
    else if (q < std::min(Tp{0.02L}, std::exp(-1.5 * a) / std::tgamma(a))
	     && (a < 10))
      {
	m = 0;
	const auto b = Tp{1} - a;
	const auto b2 = b * b;
	const auto b3 = b2 * b;
	eta = std::sqrt(-2 / a * std::log(q * tgamma_scaled(a) * s_sqrt_2pi / std::sqrt(a)));
	x0 = a * lambda(eta);
	const auto L = std::log(x0);
	if ((a > Tp{0.12L}) || (x0 > Tp{5}))
	  {
	    const auto L2 = L * L;
	    const auto L3 = L * L2;
	    const auto L4 = L * L3;
	    const auto r = Tp{1} / x0;
	    ck[0] = L - Tp{1};
	    ck[1] = (3 * b - 2 * b * L + L2 - 2 * L + 2) / Tp{2};
	    ck[2] = (24 * b * L - 11 * b2 - 24 * b - 6 * L2 + 12 * L
		     - 12 - 9 * b * L2 + 6 * b2 * L + 2 * L3)
		  / Tp{6};
	    ck[3] = (-12 * b3 * L + 84 * b * L2
		     - 114 * b2 * L + 72 + 36 * L2 + 3 * L4
		     - 72 * L + 162 * b - 168 * b * L - 12 * L3 + 25 * b3
		     - 22 * b * L3 + 36 * b2 * L2 + 120 * b2) / Tp{12};
	    x0 += -L + b * r * (ck[0] + r * (ck[1] + r * (ck[2] + r * ck[3])));
	  }
	else
	  {
	    const auto rx0 = Tp{1} / x0;
	    ck[0] = L - Tp{1};
	    x0 += -L + b * rx0 * ck[0];
	  }
      }
    else if (std::abs(porq - 0.5) < Tp{1.0e-5L})
      {
	m = 0;
	x0 = a - Tp{1} / Tp{3}
	   + (Tp{8} / Tp{405} + Tp{184} / Tp{25515} * ra) * ra;
      }
    else if (std::abs(a - Tp{1}) < Tp{1.0e-4L})
      {
	m = 0;
	if (pcase)
	  x0 = -std::log(Tp{1} - p);
	else
	  x0 = -std::log(q);
      }
    else if (a < Tp{1})
      {
	m = 0;
	if (pcase)
	  x0 = std::exp(ra * (std::log(porq) + std::lgamma(ap1)));
	else
	  x0 = std::exp(ra * (std::log(Tp{1} - porq) + std::lgamma(ap1)));
      }
    else
      {
	m = 1;
	const auto r = erfc_inv(Tp{2} * porq);
	eta = s * r / std::sqrt(a * Tp{0.5L});
	eta += ra * (epsilon1(eta)
	     + ra * (epsilon2(eta)
	     + ra * (epsilon3(eta)
	     + ra * epsilon4(eta))));
	x0 = a * lambda(eta);
      }

    auto delta = Tp{1};
    auto x = x0;
    int n = 0;
    const auto a2 = a * a;
    const int max_num_iter = 50;
    const auto max_delta = 50 * std::numeric_limits<Tp>::epsilon();
    // Implementation of the high order Newton-like method
    while (delta > max_delta && n < max_num_iter)
      {
	x = x0;
	const auto x2 = x * x;
	if (m == 0)
	  {
	    const auto dlnr = (Tp{1} - a) * std::log(x) + x + std::lgamma(a);
	    if (dlnr > std::log(s_giant))
	      throw std::runtime_error("inv_gamma: Overflow in computation");
	    else
	      {
		auto r = std::exp(dlnr);
		const auto [px, qx] = emsr::detail::gamma(a, x);
		if (pcase)
		  ck[0] = -r * (px - p);
		else
		  ck[0] = r * (qx - q);
		ck[1] = (x - a + Tp{1}) / (Tp{2} * x);
		ck[2] = (2 * x2 - 4 * x * a + 4 * x + 2 * a2 - 3 * a + 1)
		      / (6 * x2);

		r = ck[0];
		if (a > Tp{0.1L})
		  x0 = x + r * (1 + r * (ck[1] + r * ck[2]));
		else
		  {
		    if (a > Tp{0.05L})
		      x0 = x + r * (1 + r * (ck[1]));
		    else
		      x0 = x + r;
		  }
	      }
	  }
	else
	  {
	    const auto y = eta;
	    const auto fp = -std::sqrt(a) / s_sqrt_2pi
			  * std::exp(-0.5 * a * y * y)
			  / tgamma_scaled(a);
	    auto r = -x / fp;
	    const auto [px, qx] = emsr::detail::gamma(a, x);
	    if (pcase)
	      ck[0] = -r * (px - p);
	    else
	      ck[0] = r * (qx - q);
	    ck[1] = (x - a + Tp{1}) / (Tp{2} * x);
	    ck[2] = (2 * x2 - 4 * x * a + 4 * x + 2 * a2 - 3 * a + 1)
		  / (6 * x2);

	    r = ck[0];
	    if (a > Tp{0.1L})
	      x0 = x + r * (1 + r * (ck[1] + r * ck[2]));
	    else
	      {
		if (a > Tp{0.05L})
		  x0 = x + r * (1 + r * ck[1]);
		else
		  x0 = x + r;
	      }
	  }
	delta = std::abs(x / x0 - Tp{1});
	x = x0;
	++n;
      }
    if (n == max_num_iter)
      throw std::runtime_error("inv_gamma: Too many iterations in Newton root search.");

    return x;
  }

int
main()
{
  for (int ia = 1; ia <= 100; ++ia)
    {
      auto a = ia * 0.1;
      std::cout << "\na = " << a << '\n';
      for (int ix = 0; ix <= 100; ++ix)
	{
	  auto x = ix * 0.1;
	  const auto [p, q] = emsr::detail::gamma(a, x);
	  const auto xr = inv_gamma(a, p, q);
	  std::cout << ' ' << x
		    << ' ' << p
		    << ' ' << q
		    << ' ' << xr
		    << ' ' << xr - x
		    << '\n';
	}
    }
}
