/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -Wno-psabi -I. -o test_inv_gamma test_inv_gamma.cpp
*/

#include <ext/cmath>
#include <limits>

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
template<typename RealTp>
  RealTp
  inv_gamma(RealTp a, RealTp p, RealTp q)
  {
    constexpr auto _S_giant = std::numeric_limits<RealTp>::max() / RealTp{10};
    constexpr auto _S_sqrt_2 = __gnu_cxx::__const_root_2<RealTp>();
    constexpr auto _S_sqrt_pi = __gnu_cxx::__const_root_pi<RealTp>();
    constexpr auto _S_sqrt_2pi = _S_sqrt_2 * _S_sqrt_pi;
    bool pcase;
    RealTp porq, s;
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

    const auto ra = RealTp{1} / a;
    const auto ap1 = a + RealTp{1};
    RealTp ck[5];
    RealTp eta{};
    RealTp x0{};
    int m{};
    const auto logr = ra * (std::log(p) + std::lgamma(ap1));
    if (logr < std::log(0.2 * ap1))
      {
	const auto r = std::exp(logr);
	m = 0;
	const auto a2 = a * a;
	const auto a3 = a * a2;
	const auto a4 = a * a3;
	const auto ap12 = ap1 * ap1;
	const auto ap13 = ap1 * ap12;
	const auto ap14 = ap12 * ap12;
	const auto ap2 = a + RealTp{2};
	const auto ap22 = ap2 * ap2;
	ck[0] = RealTp{1};
	ck[1] = RealTp{1} / ap1;
	ck[2] = RealTp{0.5L} * (3 * a + 5) / (ap12 * (a + 2));
	ck[3] = (RealTp{1} / RealTp{3})
	      * (31 + 8 * a2 + 33 * a)
	      / (ap13 * ap2 * (a + 3));
	ck[4] = (RealTp{1} / RealTp{24})
	      * (2888 + 1179 * a3 + 125 * a4 + 3971 * a2 + 5661 * a)
	      / (ap14 * ap22 * (a + 3) * (a + 4));
	x0 = r * (1 + r * (ck[1] + r * (ck[2] + r * (ck[3] + r * ck[4]))));
      }
    else if (q < std::min(RealTp{0.02L}, std::exp(-1.5 * a) / std::tgamma(a))
	     && (a < 10))
      {
	m = 0;
	const auto b = RealTp{1} - a;
	const auto b2 = b * b;
	const auto b3 = b2 * b;
	eta = std::sqrt(-2 / a * std::log(q * gamma_scaled(a) * _S_sqrt_2pi / std::sqrt(a)));
	x0 = a * lambdaeta(eta);
	const auto L = std::log(x0);
	if ((a > RealTp{0.12L}) || (x0 > RealTp{5}))
	  {
	    const auto L2 = L * L;
	    const auto L3 = L * L2;
	    const auto L4 = L * L3;
	    const auto r = RealTp{1} / x0;
	    ck[0] = L - 1;
	    ck[1] = (3 * b - 2 * b * L + L2 - 2 * L + 2) / RealTp{2};
	    ck[2] = (24 * b * L - 11 * b2 - 24 * b - 6 * L2 + 12 * L
		     - 12 - 9 * b * L2 + 6 * b2 * L + 2 * L3)
		  / RealTp{6};
	    ck[3] = (-12 * b3 * L + 84 * b * L2
		     - 114 * b2 * L + 72 + 36 * L2 + 3 * L4
		     - 72 * L + 162 * b - 168 * b * L -12 * L3 + 25 * b3
		     - 22 * b * L3 + 36 * b2 * L2 + 120 * b2) / RealTp{12};
	    x0 = x0 - L + b * r * (ck[0] + r * (ck[1] + r * (ck[2] + r * ck[3])));
	  }
	else
	  {
	    const auto rx0 = RealTp{1} / x0;
	    ck[0] = L - RealTp{1};
	    x0 = x0 - L + b * rx0 * ck[0];
	  }
      }
    else if (std::abs(porq - 0.5) < RealTp{1.0e-5L})
      {
	m = 0;
	x0 = a - RealTp{1} / RealTp{3}
	   + (RealTp{8} / RealTp{405} + RealTp{184} / RealTp{25515} * ra) * ra;
      }
    else if (std::abs(a - RealTp{1}) < RealTp{1.0e-4L})
      {
	m = 0;
	if (pcase)
	  x0 = -std::log(RealTp{1} - p);
	else
	  x0 = -std::log(q);
      }
    else if (a < RealTp{1})
      {
	m = 0;
	if (pcase)
	  x0 = std::exp(ra * (std::log(porq) + std::lgamma(ap1)));
	else
	  x0 = std::exp(ra * (std::log(RealTp{1} - porq) + std::lgamma(ap1)));
      }
    else
      {
	m = 1;
	const auto r = inv_erfc(2 * porq);
	eta = s * r / std::sqrt(a * RealTp{0.5L});
	eta += (eps1(eta) + (eps2(eta) + eps3(eta) * ra) * ra) * ra;
	x0 = a * lambdaeta(eta);
      }

    auto delta = RealTp{1};
    auto x = x0;
    int n = 0;
    const auto a2 = a * a;
    const auto a3 = a * a2;
    // Implementation of the high order Newton-like method
    while (delta > RealTp{1.0e-15L} && n < 15)
      {
	x = x0;
	const auto x2 = x * x;
	if (m == 0)
	{
	  const auto dlnr = (RealTp{1} - a) * std::log(x) + x + std::lgamma(a);
	  if (dlnr > std::log(_S_giant))
	    std::__throw_runtime_error("inv_gamma: Overflow in computation");
	  else
	    {
	      auto r = std::exp(dlnr);
	      const auto [px, qx] = std::__detail::__gamma(a, x);
	      if (pcase)
		ck[0] = -r * (px - p);
	      else
		ck[0] = r * (qx - q);
	      ck[1] = (x - a + RealTp{1}) / (RealTp{2} * x);
	      ck[2] = (2 * x2 - 4 * x * a + 4 * x + 2 * a2 - 3 * a + 1)
		    / (6 * x2);

	      r = ck[0];
	      if (a > RealTp{0.1L})
		x0 = x + r * (1 + r * (ck[1] + r * ck[2]));
	      else
		{
		  if (a > RealTp{0.05L})
		    x0 = x + r * (1 + r * (ck[1]));
		  else
		    x0 = x + r;
		}
	    }
	}
	else
	{
	  const auto y = eta;
	  const auto fp = -std::sqrt(a) / _S_sqrt_2pi
			* std::exp(-0.5 * a * y * y)
			/ gamma_scaled(a);
	  auto r = -x / fp;
	  const auto [px, qx] = std::__detail::__gamma(a, x);
	  if (pcase)
	    ck[0] = -r * (px - p);
	  else
	    ck[0] = r * (qx - q);
	  ck[1] = (x - a + RealTp{1}) / (RealTp{2} * x);
	  ck[2] = (2 * x2 - 4 * x * a + 4 * x + 2 * a2 - 3 * a + 1)
		/ (6 * x2);

	  r = ck[0];
	  if (a > RealTp{0.1L})
	    x0 = x + r * (1 + r * (ck[1] + r * ck[2]));
	  else
	    {
	      if (a > RealTp{0.05L})
		x0 = x + r * (1 + r * ck[1]);
	      else
		x0 = x + r;
	    }
	}
	delta = std::abs(x / x0 - RealTp{1});
	x = x0;
	++n;
      }
    if (n == 15)
      std::__throw_runtime_error("inv_gamma:"
				" Too many iterations in Newton root search.");

    return x;
  }

int
main()
{
  for (int ia = 1; ia <= 100; ++ia)
    {
      auto a = ia * 0.1;
      for (int ix = 0; ix <= 100; ++ix)
	{
	  auto x = ix * 0.1;
	  const auto [p, q] = std::__detail::__gamma(a, x);
	  const auto xr = inv_gamma(a, p, q);
	}
    }
}
