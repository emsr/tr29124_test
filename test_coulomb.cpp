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
  jwkb(_Tp xl, _Tp eta1, _Tp rho1, _Tp& fjwkb, _Tp& gjwkb);

/**
 * Revised Coulomb wavefunction program using Steed's method
 *
 * A. R. Barnett           Manchester  March   1981
 *
 * Return F, G, F', G', for real rho1>0, real eta1 (including 0),
 * and real lambda(xlmin) > -1 for integer-spaced lambda values
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
  __coulomb_steed(_Tp ll, _Tp eta1, _Tp rho1,// _Tp xlmin, _Tp xlmax,
		  _Tp& fc, _Tp& gc, _Tp& fcp, _Tp& gcp, int kfn, int& ifail)
  {
    constexpr auto zero = _Tp{0};
    constexpr auto one = _Tp{1};
    constexpr auto two = _Tp{2};
    constexpr auto ten2 = _Tp{100};
    constexpr auto abort = _Tp{20000};
    constexpr auto half = _Tp{0.5};
    //constexpr auto tm30 = _Tp{1.0e-30};
    // This constant is  dsqrt(two/pi)
    constexpr auto rt2epi = _Tp{0.79788'45608'02865e0};
    // Change accur to suit machine and precision required
    const auto accur = std::numeric_limits<_Tp>::epsilon();
    ifail = 0;
    auto eta = eta1;
    if (kfn != 0)
      eta = zero;
    //bool etane0 = (eta != zero);
    const auto acc = accur;
    const auto acc4 = acc * ten2 * ten2;
    const auto acch = std::sqrt(acc);
    // Test range of rho1, exit if<= sqrt(accur) or if negative
    if (rho1 <= acch)
      {
	ifail = -1;
	return;
      }
    const auto rho = rho1;
    auto xlm = ll;
    if (kfn == 2)
      xlm -= half;
    if (xlm <= -one/* || xlmax < xlmin*/)
      {
	ifail = -2;
	return;
      }
    const auto e2mm1 = eta * eta + xlm * xlm + xlm;
    bool xlturn = rho * (rho - two * eta) < xlm * xlm + xlm;
    //dell  = xlmax - xlmin + acc;
    //if (std::abs(std::fmod(dell, one)) > acc)
    //  write(6,2040) xlmax, xlmin, dell
    //lxtra = int(dell);
    const auto xll = xlm;// + _Tp(lxtra);
    //  lxtra is number of additional lambda values to be computed
    //  xll is max lambda value, or 0.5 smaller for J, Y Bessels
    //  Determine starting array element from xlm

    //  Evaluate CF1 = f = fprime(xl,eta,x) / f(xl,eta,x)
    const auto xi = one / rho;
    auto f = zero;
    auto fcl = one;
    auto ek = zero;
    auto pk = xll + one;
    auto pk1 = zero;
    auto px = abort;
    while (true)
      {
	ek = eta / pk;
	f = (ek + pk * xi) * fcl + (fcl - one) * xi;
	pk1 = pk + one;
	// Test ensures b1 != zero for negative eta; fixup is exact.
	if (std::abs(eta * rho + pk * pk1) > acc)
	  break;
	fcl = (one + ek * ek) / (one + (eta / pk1) * (eta / pk1));
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
	//write (6,1000) d, df, acch, pk, ek, eta, x
	test += one;
	if (test > two)
	  {
	    ifail = 1;
	    return;
	  }
      }
    d = one / d;
    if (d < zero)
      fcl = -fcl;
    df *= (d * tk - one);
    f += df;
    if (pk > px)
      {
	ifail = 1;
	return;
      }
  /*
    if (lxtra == 0)
      goto 7;

    // Downward recurrence to lambda = xlm. Array gc, if present, stores rl.
    fcl *= tm30
    fpl = fcl * f
    fcp = fpl
    fc = fcl
    xl = xll;
    rl = one;
    el = zero;
    for (auto lp = 1; lp <= lxtra; ++lxtra)
      {
	if (etane0)
          el = eta / xl;
	if (etane0)
          rl = std::sqrt(one + el * el);
	sl = el + xl * xi;
	l = lp;
	fcl1 = (fcl * sl + fpl) / rl;
	fpl = fcl1 * sl - fcl * rl;
	fcl = fcl1;
	fc = fcl;
	fcp  = fpl;
	if (etane0)
          gc = rl;
      }
    xl -= one;
    if (fcl == zero)
      fcl = acc;
    f  = fpl / fcl;
    // Now we have reached lambda = xlmin = xlm 
    // Evaluate CF2 = P + I.Q  again using Steed's algorithm
    // see text for compact complex code for sp cdc or non-ansi ibm 

  7:
  */
    auto fjwkb = zero;
    auto gjwkb = zero;
    _Tp gam = zero, p = zero, q = zero, w = zero;
    if (xlturn)
      jwkb(std::max(xlm, zero), eta, rho, fjwkb, gjwkb);
    if (gjwkb <= one / (acch * ten2))
      {
	// Arrive here if g > 10^6
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
	while (std::abs(dp) + std::abs(dq) >= (std::abs(p) + std::abs(q)) * acc);

	// Solve for fcm = f at lambda = xlm,then find norm factor w=w/fcm 
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
	alpha = xi * half;
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
  /*
    if (lxtra == 0)
      return;
    // Upward recurrence from GC, G' stored value is rl
    // Renormalise FC, F' at each lambda and correct regular derivative
    //    xl = xlm here and rl = one, el = zero for Bessels
    w *= beta / std::abs(fcl);
    xl += one;
    if (etane0)
      {
	el = eta / xl;
	rl = gc;
      }
    sl = el + xl * xi;
    gcl1 = ((sl - alpha) * gcl - gpl) / rl;
    gpl = rl * gcl - (sl + alpha) * gcl1;
    gcl = gcl1;
    gc  = gcl1;
    gcp = gpl;
    fcp = w * (fcp - alpha * fc);
    fc *= w;
  */
    return;

  // 1000 FORMAT(/' CF1 ACCURACY LOSS& D,DF,ACCH,K,ETA/K,ETA,X = ',1P7E9.2/)

  //  ERROR MESSAGES
  /*
    100 IFAIL = -1;
	WRITE(6,2000) rho1,ACCH;
   2000 FORMAT(' FOR rho1 = ',1PE12.3,' TRY SMALL-X  SOLUTIONS',
       *' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',E12.3/)
	RETURN;
    105 IFAIL = -2;
	WRITE (6,2005) XLMAX,XLMIN,XLM
   2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES&XLMAX,XLMIN,XLM = ',
       *1P3E15.6/)
	RETURN;
    110 IFAIL =  1
	WRITE (6,2010) ABORT,F ,DF,PK,PX,ACC
   2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
       *' F,DF,PK,PX,ACCUR =  ',1P5E12.3//)
	RETURN;
    120 IFAIL =  2
	WRITE (6,2020) ABORT,P,Q,DP,DQ,ACC
   2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
       *' P,Q,DP,DQ,ACCUR =  ',1P4E17.7,E12.3//)
	RETURN;
    130 IFAIL =  3
	WRITE (6,2030) P,Q,ACC,DELL,LXTRA
   2030 FORMAT(' FINAL Q<= ABS(P)*ACC*10**4 , P,Q,ACC = ',1P3E12.3,4X,
       *' DELL,LXTRA = ',E12.3 /)
	RETURN;
   2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P3E20.10/)
  */
  }


// Computes JWKB approximations to COULOMB functions for xl >= 0
// as modified by Biedenharn et al. Phys Rev 97 (1955) 542-554
template<typename _Tp>
  void
  jwkb(_Tp xl, _Tp eta1, _Tp rho1, _Tp& fjwkb, _Tp& gjwkb)
  {
    constexpr auto zero = _Tp{0};
    constexpr auto half = _Tp{0.5};
    constexpr auto one = _Tp{1};
    constexpr auto six = _Tp{6};
    //constexpr auto ten = _Tp{10};
    constexpr auto rl35 = _Tp{35};
    //constexpr auto aloge = _Tp{0.43429'45e0};
    const auto rho = rho1;
    const auto eta = eta1;
    const auto gh2 = rho * (eta + eta - rho);
    const auto xll1 = std::max(xl * xl + xl, zero);
    if (gh2 + xll1 <= zero)
      return;
    const auto hll = xll1 + six / rl35;
    const auto hl = std::sqrt(hll);
    const auto sl = eta / hl + hl / rho;
    const auto rl2 = one + eta * eta / hll;
    const auto gh = std::sqrt(gh2 + hll) / rho;
    auto phi = rho * gh - half * (hl * std::log((gh + sl) * (gh + sl) / rl2) - std::log(gh));
    if (eta != zero)
      phi -= eta * std::atan2(rho * gh, rho - eta);
    //const auto phi10 = -phi * aloge;
    gjwkb = std::exp(-phi);
    fjwkb = half / (gh * gjwkb);
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
	for (int irho = 0; irho <= 200; ++irho)
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


