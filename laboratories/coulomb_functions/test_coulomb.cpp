/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <complex>

enum Func
{
  Coulomb,
  SphBessel,
  CylBessel
};

template<typename Tp>
  void
  coulomb_jwkb(Tp lambda, Tp eta, Tp rho, Tp& fjwkb, Tp& gjwkb, int& iexp);

/**
 * Revised Coulomb wavefunction program using Steed's method
 *
 * A. R. Barnett           Manchester  March   1981
 *
 * Return F, G, F', G', for real rho>0, real eta1 (including 0),
 * and real lambda > -1 for integer-spaced lambda values
 * thus giving positive-energy solutions to the Coulomb Schrodinger
 * equation, to the Klein-Gordon equation and to suitable forms of
 * the Dirac equation, and also spherical + cylindrical Bessel equations.
 *
 * If func = 0 Real Coulomb functions are returned
 *         = 1 Spherical Bessel
 *         = 2 Cylindrical Bessel
 */
template<typename Tp>
  void
  coulomb_steed(Func func, Tp lambda, Tp eta, Tp rho,
		  Tp& fc, Tp& gc, Tp& fcp, Tp& gcp, int& ifail)
  {
    constexpr auto abort = Tp{20000};
    //constexpr auto s_tiny = Tp{1.0e-30};
    // This constant is  sqrt(Tp{2}/pi)
    constexpr auto rt2epi = Tp{0.79788'45608'02865e0};
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    const auto s_i = std::complex<Tp>{0, 1};
    ifail = 0;
    int iexp = 1;
    if (func != Coulomb)
      eta = Tp{0};
    const auto acc4 = std::pow(s_eps, 0.75);// * Tp{10000};
    const auto acch = std::sqrt(s_eps);
    // Test range of rho, exit if <= sqrt(epsilon) or if negative
    if (rho <= acch)
      {
	ifail = -1;
	return;
      }
    if (func == 2)
      lambda -= Tp{0.5};
    if (lambda <= -Tp{1})
      {
	ifail = -2;
	return;
      }
    const auto e2mm1 = eta * eta + lambda * (lambda + Tp{1});
    bool turn = rho * (rho - Tp{2} * eta) < lambda * (lambda + Tp{1});
    // lxtra is number of additional lambda values to be computed
    // lamax is max lambda value, or 0.5 smaller for J, Y Bessels
    int lxtra = std::numeric_limits<Tp>::digits10;
    const auto lamax = lambda + Tp(lxtra);

    //  Evaluate CF1 = f = F'(lambda,eta,x) / F(lambda,eta,x)
    const auto xi = Tp{1} / rho;
    auto f = Tp{0};
    auto fcl = Tp{1};
    auto ek = Tp{0};
    auto pk = lamax + Tp{1};
    auto pk1 = Tp{0};
    auto px = pk + abort;
    while (true)
      {
	ek = eta / pk;
	f = (ek + pk * xi) * fcl + (fcl - Tp{1}) * xi;
	pk1 = pk + Tp{1};
	const auto ek1 = eta / pk1;
	// Test ensures b1 != Tp{0} for negative eta; fixup is exact.
	if (std::abs(eta * rho + pk * pk1) > s_eps)
	  break;
	fcl = (Tp{1} + ek * ek) / (Tp{1} + ek1 * ek1);
	pk +=  Tp{2};
      }
    auto d = Tp{1} / ((pk + pk1) * (xi + ek / pk1));
    auto df = -fcl * (Tp{1} + ek * ek) * d;
    if (fcl != Tp{1} )
      fcl = -Tp{1};
    if (d < Tp{0})
      fcl = -fcl;
    f += df;

    // Begin CF1 loop on pk = k = lambda + 1
    auto test = Tp{1};
    auto tk = Tp{0};
    do
      {
	pk = pk1;
	pk1 += Tp{1};
	ek = eta / pk;
	tk = (pk + pk1) * (xi + ek / pk1);
	d = tk - d * (Tp{1} + ek * ek);
	if (std::abs(d) < acch)
	  {
	    test += Tp{1};
	    if (test > Tp{2})
	      {
		ifail = 1;
		return;
	      }
	  }
	d = Tp{1} / d;
	if (d < Tp{0})
	  fcl = -fcl;
	df *= (d * tk - Tp{1});
	f += df;
	if (pk > px)
	  {
	    ifail = 1;
	    return;
	  }
      }
    while (std::abs(df) > s_eps * std::abs(f));
    //if (lxtra == 0)
    //  goto 7;
/*
    // Downward recurrence to lambda = lambda.
    fcl *= s_tiny;
    auto fcpl = fcl * f;
    fcp = fcpl;
    fc = fcl;
    auto xl = lamax;
    auto rl = Tp{1};
    auto el = Tp{0};
    for (auto lp = 1; lp <= lxtra; ++lp)
      {
	if (eta != Tp{0})
	  {
            el = eta / xl;
            rl = std::sqrt(Tp{1} + el * el);
	  }
	auto sl = el + xl * xi;
	auto fcl1 = (fcl * sl + fcpl) / rl;
	fcpl = fcl1 * sl - fcl * rl;
	fcl = fcl1;
	fc = fcl;
	fcp = fcpl;
	//if (eta != Tp{0})
          //gc = rl;
	xl -= Tp{1};
      }
    if (std::abs(fcl) < s_eps)
      fcl = s_eps;
    f = fcpl / fcl;
*/
  //7:

    // Now we have reached lambda
    // Evaluate CF2 = P + I.Q  again using Steed's algorithm
    // see text for compact complex code for sp cdc or non-ansi ibm 
    auto fjwkb = Tp{0};
    auto gjwkb = Tp{0};
    Tp gam = Tp{0}, w = Tp{0};
    Tp p, q;
    if (turn)
      coulomb_jwkb(std::max(lambda, Tp{0}), eta, rho, fjwkb, gjwkb, iexp);
    if (iexp > 1 && gjwkb > Tp{1} / (Tp{100} * acch))
      {
	// Arrive here if G > 10^6
	w = fjwkb;
	gam = gjwkb * w;
	p = f;
	q = Tp{1};
      }
    else
      {
	turn = false;
	const auto tab = Tp{2} * abort;
	pk = Tp{0};
	const auto wi = Tp{2} * eta;
	auto pq = std::complex<Tp>(Tp{0}, Tp{1} - eta * xi);
	auto a = std::complex<Tp>(-e2mm1, eta);
	auto b = Tp{2} * std::complex<Tp>(rho - eta, Tp{1});
	auto d = std::conj(b) / std::norm(b);
	auto dpq = s_i * xi * a * d;
	do
	  {
	    pq += dpq;
	    pk += Tp{2};
	    a += std::complex<Tp>(pk, wi);
	    b += Tp{2} * s_i;
	    d = b + a * d;
	    d = std::conj(d) / std::norm(d);
	    auto ab = b * d - Tp{1};
	    dpq *= ab;
	    if (pk > tab)
	      {
		ifail = 2;
		return;
	      }
	  }
	while (std::abs(dpq) > s_eps * std::abs(pq));

	// Solve for fcm = f at lambda, then find norm factor w = w / fcm.
        auto p = std::real(pq);
        auto q = std::imag(pq);
	gam = (f - p) / q;
	if (q <= acc4 * std::abs(p))
	  {
	    ifail = 3;
	    return;
	  }
	w = Tp{1} / std::sqrt((f - p) * gam + q);
      }

    // Normalise for spherical or cylindrical Bessel functions
    auto alpha = Tp{0};
    auto beta = Tp{1};
    if (func == SphBessel)
      {
	alpha = xi;
	beta = xi;
      }
    else if (func == CylBessel)
      {
	alpha = Tp{0.5} * xi;
	beta = std::sqrt(xi) * rt2epi;
      }
    const auto fcm = std::copysign(w, fcl) * beta;
    fc = fcm;

    auto gcl = turn ? gjwkb * beta : fcm * gam;
    if (func != Coulomb)
      gcl = -gcl;
    gc = gcl;
    auto gpl = gcl * (p - q / gam) - alpha * gcl;
    gcp = gpl;
    fcp = fcm * (f - alpha);
/*
    // Upward recurrence from GC, G' stored value is rl
    // Renormalise FC, F' at each lambda and correct regular derivative
    //    xl = lambda here and rl = Tp{1}, el = Tp{0} for Bessels
    w *= beta / std::abs(fcl);
    xl += Tp{1};
    if (eta != Tp{0})
      {
	el = eta / xl;
	rl = gc;
      }
    const auto sl = el + xl * xi;
    const auto gcl1 = ((sl - alpha) * gcl - gpl) / rl;
    gpl = rl * gcl - (sl + alpha) * gcl1;
    gcl = gcl1;
    gc  = gcl1;
    gcp = gpl;
    fcp = w * (fcp - alpha * fc);
    fc *= w;
*/
    return;
  }


/**
 * Computes JWKB approximations to Coulomb functions for lambda >= 0
 * as modified by Biedenharn et al. Phys Rev 97 (1955) 542-554
 */
template<typename Tp>
  void
  coulomb_jwkb(Tp lambda, Tp eta, Tp rho, Tp& fjwkb, Tp& gjwkb, int& iexp)
  {
    const auto gh2 = rho * (Tp{2} * eta - rho);
    const auto lamax = std::max(lambda * (lambda + Tp{1}), Tp{0});
    if (gh2 + lamax <= Tp{0})
      return;
    const auto hll = lamax + Tp{6} / Tp{35};
    const auto hl = std::sqrt(hll);
    const auto sl = eta / hl + hl / rho;
    const auto rl2 = Tp{1} + eta * eta / hll;
    const auto gh = std::sqrt(gh2 + hll) / rho;
    auto phi = rho * gh
	     - Tp{0.5} * (hl * std::log((gh + sl) * (gh + sl) / rl2)
			  - std::log(gh));
    phi -= eta * std::atan2(rho * gh, rho - eta);
    const auto phi10 = -phi;
    iexp = int(phi10);
    if (iexp > 250)
      gjwkb = std::exp(phi10 - Tp(iexp));
    else
      {
	gjwkb = std::exp(-phi);
	iexp = 0;
      }
    fjwkb = Tp{0.5} / (gh * gjwkb);
    return;
  }

template<typename Tp>
  void
  test_coulomb()
  {
    Func func = Coulomb;
    for (auto lambda : {Tp{0}, Tp{0.5L}, Tp{1}})
      {
	for (auto eta : {Tp{-2}, Tp{0}, Tp{2}, Tp{10}})
	  {
	    std::cout << "\n\nlambda = " << lambda << "; eta = " << eta << '\n';
	    for (int irho = 1; irho <= 200; ++irho)
	      {
		auto rho = irho * Tp{0.1};
		Tp fc = 0, gc = 0, fcp = 0, gcp = 0;
		int ifail = 0;
		coulomb_steed(func, lambda, eta, rho, fc, gc, fcp, gcp, ifail);
		std::cout << ' ' << std::setw(16) << rho
			  << ' ' << std::setw(16) << fc
			  << ' ' << std::setw(16) << gc
			  << ' ' << std::setw(16) << fcp
			  << ' ' << std::setw(16) << gcp
			  << ' ' << std::setw(4) << ifail
			  << '\n';
	      }
	  }
      }
  }

int
main()
{
  using Tp = double;
  test_coulomb<Tp>();
}
