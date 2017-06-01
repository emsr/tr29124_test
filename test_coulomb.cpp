/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_coulomb test_coulomb.cpp -lquadmath -Lwrappers/debug -lwrap_boost
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_coulomb > test_coulomb.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_coulomb test_coulomb.cpp -lquadmath -Lwrappers/debug -lwrap_boost
./test_coulomb > test_coulomb.txt
*/

#include <cmath>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  void
  __coulomb_jwkb(_Tp lamb1, _Tp eta1, _Tp rho1, _Tp& fjwkb, _Tp& gjwkb, int& iexp);

/**
 * Revised Coulomb wavefunction program using Steed's method
 *
 * A. R. Barnett           Manchester  March   1981
 *
 * Return F, G, F', G', for real rho1>0, real eta1 (including 0),
 * and real lambda > -1 for integer-spaced lambda values
 * thus giving positive-energy solutions to the Coulomb Schrodinger
 * equation, to the Klein-Gordon equation and to suitable forms of
 * the Dirac equation, and also spherical + cylindrical Bessel equations.
 *
 * If kfn  = 0 Real Coulomb functions are returned
 *         = 1 Spherical Bessel
 *         = 2 Cylindrical Bessel
 */
template<typename _Tp>
  void
  __coulomb_steed(_Tp lambda1, _Tp eta1, _Tp rho1,
		  _Tp& fc, _Tp& gc, _Tp& fcp, _Tp& gcp, int kfn, int& ifail)
  {
    constexpr auto zero = _Tp{0};
    constexpr auto one = _Tp{1};
    constexpr auto two = _Tp{2};
    constexpr auto abort = _Tp{20000};
    constexpr auto tm30 = _Tp{1.0e-30};
    // This constant is  sqrt(two/pi)
    constexpr auto rt2epi = _Tp{0.79788'45608'02865e0};
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    ifail = 0;
    int iexp = 1;
    auto eta = eta1;
    if (kfn != 0)
      eta = zero;
    bool etane0 = (eta != zero);
    const auto acc4 = _S_eps * _Tp{10000};
    const auto acch = std::sqrt(_S_eps);
    // Test range of rho1, exit if <= sqrt(epsilon) or if negative
    if (rho1 <= acch)
      {
	ifail = -1;
	return;
      }
    const auto rho = rho1;
    auto lambda = lambda1;
    if (kfn == 2)
      lambda -= _Tp{0.5};
    if (lambda <= -one)
      {
	ifail = -2;
	return;
      }
    const auto e2mm1 = eta * eta + lambda * (lambda + one);
    bool xlturn = rho * (rho - two * eta) < lambda * (lambda + one);
    // lxtra is number of additional lambda values to be computed
    // lamax is max lambda value, or 0.5 smaller for J, Y Bessels
    int lxtra = std::numeric_limits<_Tp>::digits10;
    const auto lamax = lambda + _Tp(lxtra);

    //  Evaluate CF1 = f = F'(lambda,eta,x) / F(lambda,eta,x)
    const auto xi = one / rho;
    auto f = zero;
    auto fcl = one;
    auto ek = zero;
    auto pk = lamax + one;
    auto pk1 = zero;
    auto px = abort;
    while (true)
      {
	ek = eta / pk;
	f = (ek + pk * xi) * fcl + (fcl - one) * xi;
	pk1 = pk + one;
	const auto ek1 = eta / pk1;
	// Test ensures b1 != zero for negative eta; fixup is exact.
	if (std::abs(eta * rho + pk * pk1) > _S_eps)
	  break;
	fcl = (one + ek * ek) / (one + ek1 * ek1);
	pk +=  two;
      }
    auto d = one / ((pk + pk1) * (xi + ek / pk1));
    auto df = -fcl * (one + ek * ek) * d;
    if (fcl != one )
      fcl = -one;
    if (d < zero)
      fcl = -fcl;
    f += df;

    // Begin CF1 loop on pk = k = lambda + 1
    auto test = one;
    auto tk = zero;
    while (true)
      {
	pk = pk1;
	pk1 += one;
	ek = eta / pk;
	tk = (pk + pk1) * (xi + ek / pk1);
	d = tk - d * (one + ek * ek);
	if (std::abs(d) > acch)
	  break;
	test += one;
	if (test > two)
	  {
	    ifail = 1;
	    return;
	  }
	d = one / d;
	if (d < zero)
	  fcl = -fcl;
      }
    df *= (d * tk - one);
    f += df;
    if (pk > px)
      {
	ifail = 1;
	return;
      }
    //if (lxtra == 0)
    //  goto 7;

    // Downward recurrence to lambda = lambda.
    fcl *= tm30;
    auto fcpl = fcl * f;
    fcp = fcpl;
    fc = fcl;
    auto xl = lamax;
    auto rl = one;
    auto el = zero;
    for (auto lp = 1; lp <= lxtra; ++lp)
      {
	if (etane0)
	  {
            el = eta / xl;
            rl = std::sqrt(one + el * el);
	  }
	auto sl = el + xl * xi;
	auto fcl1 = (fcl * sl + fcpl) / rl;
	fcpl = fcl1 * sl - fcl * rl;
	fcl = fcl1;
	fc = fcl;
	fcp = fcpl;
	if (etane0)
          gc = rl;
	xl -= one;
      }
    if (std::abs(fcl) < _S_eps)
      fcl = _S_eps;
    f = fcpl / fcl;
  //7:

    // Now we have reached lambda
    // Evaluate CF2 = P + I.Q  again using Steed's algorithm
    // see text for compact complex code for sp cdc or non-ansi ibm 
    auto fjwkb = zero;
    auto gjwkb = zero;
    _Tp gam = zero, p = zero, q = zero, w = zero;
    if (xlturn)
      __coulomb_jwkb(std::max(lambda, zero), eta, rho, fjwkb, gjwkb, iexp);
    if (iexp <= 1 && gjwkb <= one / (_Tp{100} * acch))
      {
	// Arrive here if G > 10^6
	w = fjwkb;
	gam = gjwkb * w;
	p = f;
	q = one;
      }
    else
      {
	xlturn = false;
	const auto tab = two * abort;
	pk = zero;
	const auto wi = eta + eta;
	p = zero;
	q = one - eta * xi;
	auto ar = -e2mm1;
	auto ai = eta;
	auto br = two * (rho - eta);
	auto bi = two;
	auto dr = br / (br * br + bi * bi);
	auto di = -bi / (br * br + bi * bi);
	auto dp = -xi * (ar * di + ai * dr);
	auto dq = xi * (ar * dr - ai * di);
	do
	  {
	    p += dp;
	    q += dq;
	    pk += two;
	    ar += pk;
	    ai += wi;
	    bi += two;
	    d = ar * dr - ai * di + br;
	    di = ai * dr + ar * di + bi;
	    auto c = one / (d + d + di * di);
	    dr = c * d;
	    di = -c * di;
	    auto a = br * dr - bi * di - one;
	    auto b = bi * dr + br * di;
	    c = dp * a - dq * b;
	    dq = dp * b + dq * a;
	    dp = c;
	    if (pk > tab)
	      {
		ifail = 2;
		return;
	      }
	  }
	while (std::abs(dp) + std::abs(dq) > (std::abs(p) + std::abs(q)) * _S_eps);

	// Solve for fcm = f at lambda,then find norm factor w = w / fcm 
	gam = (f - p) / q;
	if (q <= acc4 * std::abs(p))
	  {
	    ifail = 3;
	    return;
	  }
	w = one / std::sqrt((f - p) * gam + q);
      }

    // Normalise for spherical or cylindrical Bessel functions
    auto alpha = zero;
    auto beta = one;
    if (kfn == 1)
      {
	alpha = xi;
	beta = xi;
      }
    else if (kfn == 2)
      {
	alpha = _Tp{0.5} * xi;
	beta = std::sqrt(xi) * rt2epi;
      }
    const auto fcm = std::copysign(w, fcl) * beta;
    fc = fcm;

    _Tp gcl;
    if (xlturn)
      gcl = gjwkb * beta;
    else
      gcl = fcm * gam;

    if (kfn != 0)
      gcl = -gcl;
    gc = gcl;
    const auto gpl = gcl * (p - q / gam) - alpha * gcl;
    gcp = gpl;
    fcp = fcm * (f - alpha);
  //  if (lxtra == 0)
  //    return;
/*
    // Upward recurrence from GC, G' stored value is rl
    // Renormalise FC, F' at each lambda and correct regular derivative
    //    xl = lambda here and rl = one, el = zero for Bessels
    w *= beta / std::abs(fcl);
    xl += one;
    if (etane0)
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
  __coulomb_jwkb(_Tp lamb1, _Tp eta1, _Tp rho1, _Tp& fjwkb, _Tp& gjwkb, int& iexp)
  {
    constexpr auto zero = _Tp{0};
    constexpr auto one = _Tp{1};
    constexpr auto loge = _Tp{0.43429'45e0};
    const auto rho = rho1;
    const auto eta = eta1;
    const auto gh2 = rho * (eta + eta - rho);
    const auto lamax = std::max(lamb1 * (lamb1 + one), zero);
    if (gh2 + lamax <= zero)
      return;
    const auto hll = lamax + _Tp{6} / _Tp{35};
    const auto hl = std::sqrt(hll);
    const auto sl = eta / hl + hl / rho;
    const auto rl2 = one + eta * eta / hll;
    const auto gh = std::sqrt(gh2 + hll) / rho;
    auto phi = rho * gh
	   - _Tp{0.5} * (hl * std::log((gh + sl) * (gh + sl) / rl2) - std::log(gh));
    if (eta != zero)
      phi -= eta * std::atan2(rho * gh, rho - eta);
    const auto phi10 = -phi * loge;
    iexp = int(phi10);
    if (iexp > 250)
      gjwkb = std::pow(_Tp{10}, phi10 - _Tp(iexp));
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
    int kfn = 0;
    _Tp lambda = 0;
    for (auto eta : {_Tp{-2}, _Tp{0}, _Tp{2}, _Tp{10}})
      {
	std::cout << '\n' << '\n';
	for (int irho = 1; irho <= 200; ++irho)
	  {
	    auto rho = irho * _Tp{0.1};
	    _Tp fc, gc, fcp, gcp;
	    int ifail = 0;
	    __coulomb_steed(lambda, eta, rho, fc, gc, fcp, gcp, kfn, ifail);
	    std::cout << ' ' << std::setw(12) << rho
		      << ' ' << std::setw(12) << fc
		      << ' ' << std::setw(12) << gc
		      << ' ' << std::setw(12) << fcp
		      << ' ' << std::setw(12) << gcp
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


