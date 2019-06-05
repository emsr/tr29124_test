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

template<typename _Tp>
  void
  __coulomb_jwkb(_Tp lambda, _Tp eta, _Tp rho, _Tp& fjwkb, _Tp& gjwkb, int& iexp);

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
template<typename _Tp>
  void
  __coulomb_steed(Func func, _Tp lambda, _Tp eta, _Tp rho,
		  _Tp& fc, _Tp& gc, _Tp& fcp, _Tp& gcp, int& ifail)
  {
    constexpr auto abort = _Tp{20000};
    //constexpr auto _S_tiny = _Tp{1.0e-30};
    // This constant is  sqrt(_Tp{2}/pi)
    constexpr auto rt2epi = _Tp{0.79788'45608'02865e0};
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_i = std::complex<_Tp>{0, 1};
    ifail = 0;
    int iexp = 1;
    if (func != Coulomb)
      eta = _Tp{0};
    const auto acc4 = std::pow(_S_eps, 0.75);// * _Tp{10000};
    const auto acch = std::sqrt(_S_eps);
    // Test range of rho, exit if <= sqrt(epsilon) or if negative
    if (rho <= acch)
      {
	ifail = -1;
	return;
      }
    if (func == 2)
      lambda -= _Tp{0.5};
    if (lambda <= -_Tp{1})
      {
	ifail = -2;
	return;
      }
    const auto e2mm1 = eta * eta + lambda * (lambda + _Tp{1});
    bool turn = rho * (rho - _Tp{2} * eta) < lambda * (lambda + _Tp{1});
    // lxtra is number of additional lambda values to be computed
    // lamax is max lambda value, or 0.5 smaller for J, Y Bessels
    int lxtra = std::numeric_limits<_Tp>::digits10;
    const auto lamax = lambda + _Tp(lxtra);

    //  Evaluate CF1 = f = F'(lambda,eta,x) / F(lambda,eta,x)
    const auto xi = _Tp{1} / rho;
    auto f = _Tp{0};
    auto fcl = _Tp{1};
    auto ek = _Tp{0};
    auto pk = lamax + _Tp{1};
    auto pk1 = _Tp{0};
    auto px = pk + abort;
    while (true)
      {
	ek = eta / pk;
	f = (ek + pk * xi) * fcl + (fcl - _Tp{1}) * xi;
	pk1 = pk + _Tp{1};
	const auto ek1 = eta / pk1;
	// Test ensures b1 != _Tp{0} for negative eta; fixup is exact.
	if (std::abs(eta * rho + pk * pk1) > _S_eps)
	  break;
	fcl = (_Tp{1} + ek * ek) / (_Tp{1} + ek1 * ek1);
	pk +=  _Tp{2};
      }
    auto d = _Tp{1} / ((pk + pk1) * (xi + ek / pk1));
    auto df = -fcl * (_Tp{1} + ek * ek) * d;
    if (fcl != _Tp{1} )
      fcl = -_Tp{1};
    if (d < _Tp{0})
      fcl = -fcl;
    f += df;

    // Begin CF1 loop on pk = k = lambda + 1
    auto test = _Tp{1};
    auto tk = _Tp{0};
    do
      {
	pk = pk1;
	pk1 += _Tp{1};
	ek = eta / pk;
	tk = (pk + pk1) * (xi + ek / pk1);
	d = tk - d * (_Tp{1} + ek * ek);
	if (std::abs(d) < acch)
	  {
	    test += _Tp{1};
	    if (test > _Tp{2})
	      {
		ifail = 1;
		return;
	      }
	  }
	d = _Tp{1} / d;
	if (d < _Tp{0})
	  fcl = -fcl;
	df *= (d * tk - _Tp{1});
	f += df;
	if (pk > px)
	  {
	    ifail = 1;
	    return;
	  }
      }
    while (std::abs(df) > _S_eps * std::abs(f));
    //if (lxtra == 0)
    //  goto 7;
/*
    // Downward recurrence to lambda = lambda.
    fcl *= _S_tiny;
    auto fcpl = fcl * f;
    fcp = fcpl;
    fc = fcl;
    auto xl = lamax;
    auto rl = _Tp{1};
    auto el = _Tp{0};
    for (auto lp = 1; lp <= lxtra; ++lp)
      {
	if (eta != _Tp{0})
	  {
            el = eta / xl;
            rl = std::sqrt(_Tp{1} + el * el);
	  }
	auto sl = el + xl * xi;
	auto fcl1 = (fcl * sl + fcpl) / rl;
	fcpl = fcl1 * sl - fcl * rl;
	fcl = fcl1;
	fc = fcl;
	fcp = fcpl;
	//if (eta != _Tp{0})
          //gc = rl;
	xl -= _Tp{1};
      }
    if (std::abs(fcl) < _S_eps)
      fcl = _S_eps;
    f = fcpl / fcl;
*/
  //7:

    // Now we have reached lambda
    // Evaluate CF2 = P + I.Q  again using Steed's algorithm
    // see text for compact complex code for sp cdc or non-ansi ibm 
    auto fjwkb = _Tp{0};
    auto gjwkb = _Tp{0};
    _Tp gam = _Tp{0}, w = _Tp{0};
    _Tp p, q;
    if (turn)
      __coulomb_jwkb(std::max(lambda, _Tp{0}), eta, rho, fjwkb, gjwkb, iexp);
    if (iexp > 1 && gjwkb > _Tp{1} / (_Tp{100} * acch))
      {
	// Arrive here if G > 10^6
	w = fjwkb;
	gam = gjwkb * w;
	p = f;
	q = _Tp{1};
      }
    else
      {
	turn = false;
	const auto tab = _Tp{2} * abort;
	pk = _Tp{0};
	const auto wi = _Tp{2} * eta;
	auto pq = std::complex<_Tp>(_Tp{0}, _Tp{1} - eta * xi);
	auto a = std::complex<_Tp>(-e2mm1, eta);
	auto b = _Tp{2} * std::complex<_Tp>(rho - eta, _Tp{1});
	auto d = std::conj(b) / std::norm(b);
	auto dpq = _S_i * xi * a * d;
	do
	  {
	    pq += dpq;
	    pk += _Tp{2};
	    a += std::complex<_Tp>(pk, wi);
	    b += _Tp{2} * _S_i;
	    d = b + a * d;
	    d = std::conj(d) / std::norm(d);
	    auto ab = b * d - _Tp{1};
	    dpq *= ab;
	    if (pk > tab)
	      {
		ifail = 2;
		return;
	      }
	  }
	while (std::abs(dpq) > _S_eps * std::abs(pq));

	// Solve for fcm = f at lambda, then find norm factor w = w / fcm.
        auto p = std::real(pq);
        auto q = std::imag(pq);
	gam = (f - p) / q;
	if (q <= acc4 * std::abs(p))
	  {
	    ifail = 3;
	    return;
	  }
	w = _Tp{1} / std::sqrt((f - p) * gam + q);
      }

    // Normalise for spherical or cylindrical Bessel functions
    auto alpha = _Tp{0};
    auto beta = _Tp{1};
    if (func == SphBessel)
      {
	alpha = xi;
	beta = xi;
      }
    else if (func == CylBessel)
      {
	alpha = _Tp{0.5} * xi;
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
    //    xl = lambda here and rl = _Tp{1}, el = _Tp{0} for Bessels
    w *= beta / std::abs(fcl);
    xl += _Tp{1};
    if (eta != _Tp{0})
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
template<typename _Tp>
  void
  __coulomb_jwkb(_Tp lambda, _Tp eta, _Tp rho, _Tp& fjwkb, _Tp& gjwkb, int& iexp)
  {
    const auto gh2 = rho * (_Tp{2} * eta - rho);
    const auto lamax = std::max(lambda * (lambda + _Tp{1}), _Tp{0});
    if (gh2 + lamax <= _Tp{0})
      return;
    const auto hll = lamax + _Tp{6} / _Tp{35};
    const auto hl = std::sqrt(hll);
    const auto sl = eta / hl + hl / rho;
    const auto rl2 = _Tp{1} + eta * eta / hll;
    const auto gh = std::sqrt(gh2 + hll) / rho;
    auto phi = rho * gh
	     - _Tp{0.5} * (hl * std::log((gh + sl) * (gh + sl) / rl2)
			  - std::log(gh));
    phi -= eta * std::atan2(rho * gh, rho - eta);
    const auto phi10 = -phi;
    iexp = int(phi10);
    if (iexp > 250)
      gjwkb = std::exp(phi10 - _Tp(iexp));
    else
      {
	gjwkb = std::exp(-phi);
	iexp = 0;
      }
    fjwkb = _Tp{0.5} / (gh * gjwkb);
    return;
  }

template<typename _Tp>
  void
  test_coulomb()
  {
    Func func = Coulomb;
    _Tp lambda = 0;
    for (auto eta : {_Tp{-2}, _Tp{0}, _Tp{2}, _Tp{10}})
      {
	std::cout << "\n\neta = " << eta << '\n';
	for (int irho = 1; irho <= 200; ++irho)
	  {
	    auto rho = irho * _Tp{0.1};
	    _Tp fc = 0, gc = 0, fcp = 0, gcp = 0;
	    int ifail = 0;
	    __coulomb_steed(func, lambda, eta, rho, fc, gc, fcp, gcp, ifail);
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

int
main()
{
  using _Tp = double;
  test_coulomb<_Tp>();
}
