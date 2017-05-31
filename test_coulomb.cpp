/*
  Revised Coulomb wavefunction program using Steed's method

  A. R. Barnett           Manchester  March   1981

  Return F, G, F', G', for real xx>0, real eta1 (including 0),
  and real lambda(xlmin) > -1 for integer-spaced lambda values
  thus giving positive-energy solutions to the Coulomb Schrodinger
  equation, to the Klein-Gordon equation and to suitable forms of
  the Dirac equation, and also spherical + cylindrical Bessel equations.

  If kfn  = 0 Real Coulomb functions are returned
          = 1 Spherical Bessel
          = 2 Cylindrical Bessel
*/
template<typename _Tp>
void
coulomb(_Tp xx, _Tp eta1, _Tp ll,// _Tp xlmin, _Tp xlmax,
	_Tp& fc, _Tp& gc, _Tp& fcp, _Tp& gcp, int kfn, int& ifail)
{
  constexpr auto zero = 0.0e0;
  constexpr auto one = 1.0e0;
  constexpr auto two = 2.0e0;
  constexpr auto ten2 = 1.0e2;
  constexpr auto abort = 2.0e4;
  constexpr auto half = 0.5e0;
  constexpr auto tm30 = 1.0e-30;
  // This constant is  dsqrt(two/pi)
  constexpr auto rt2epi = 0.79788'45608'02865e0;
  // Change accur to suit machine and precision required
  accur = r1mach(4);
  ifail = 0;
  auto eta = eta1;
  auto gjwkb = zero;
  if (kfn != 0)
    eta = zero;
  bool etane0 = (eta != zero);
  const auto acc = accur;
  const auto acc4 = acc * ten2 * ten2;
  const auto acch = std::sqrt(acc);
  // Test range of xx, exit if<= sqrt(accur) or if negative
  if (xx <= acch)
    {
      ifail = -1;
      return;
    }
  auto x = xx;
  auto xlm = ll;
  if (kfn == 2)
    xlm -= half;
  if (xlm <= -one/* || xlmax < xlmin*/)
    {
      ifail = -2;
      return;
    }
  auto e2mm1 = eta * eta + xlm * xlm + xlm;
  bool xlturn = x * (x - two * eta) < xlm * xlm + xlm;
  //dell  = xlmax - xlmin + acc;
  //if (std::abs(std::fmod(dell, one)) > acc)
  //  write(6,2040) xlmax, xlmin, dell
  //lxtra = int(dell);
  //xll = xlm + _Tp(lxtra);
  //  lxtra is number of additional lambda values to be computed
  //  XLL IS MAX LAMBDA VALUE, or 0.5 smaller for J, Y BESSELS
  //  Determine starting array element from xlmin

  //  Evaluate cf1  =  f =  fprime(xl,eta,x) / f(xl,eta,x)
  auto xi = one / x;
  auto fcl = one;
  auto pk = xll + one;
  px += abort;
  while (true)
    {
      ek = eta / pk;
      f = (ek + pk * xi) * fcl + (fcl - one) * xi;
      pk1 = pk + one;
      // Test ensures b1 != zero for negative eta; fixup is exact.
      if (std::abs(eta * x + pk * pk1) > acc)
	break;
      fcl = (one + ek * ek) / (one + (eta / pk1) * (eta / pk1));
      pk +=  two;
    }
  d = one / ((pk + pk1) * (xi + ek / pk1));
  df = -fcl * (one + ek * ek) * d;
  if (fcl != one )
    fcl = -one;
  if (d < zero)
    fcl = -fcl;
  f += df;

  // Begin CF1 loop on pk = k = lambda + 1
  p = one;
  while (true)
  {
    pk = pk1;
    pk1 = pk1 + one;
    ek = eta / pk;
    tk = (pk + pk1) * (xi + ek / pk1);
    d = tk - d * (one + ek * ek);
    if (std::abs(d) > acch)
      break;
    //write (6,1000) d, df, acch, pk, ek, eta, x
    p += one;
    if (p > two)
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

  // Downward recurrence to lambda = xlm. array gc,if present,stores rl
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
  if (xlturn)
    jwkb(x, eta, std::max(xlm, zero), fjwkb, gjwkb);
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
      ta = two * abort;
      pk = zero;
      wi = eta + eta;
      p = zero;
      q = one - eta * xi;
      ar = -e2mm1;
      ai = eta;
      br = two * (x - eta);
      bi = two;
      dr = br / (br * br + bi * bi);
      di = -bi / (br * br + bi * bi);
      dp = -xi * (ar * di + ai * dr);
      dq =  xi * (ar * dr - ai * di);
      do
	{
	  p += dp;
	  q  += dq;
	  pk += two;
	  ar += pk;
	  ai += wi;
	  bi += two;
	  d  = ar * dr - ai * di + br;
	  di = ai * dr + ar * di + bi;
	  c  = one /(d + d + di * di);
	  dr = c * d;
	  di = -c * di;
	  a = br * dr - bi * di - one;
	  b = bi * dr + br * di;
	  c = dp * a - dq * b;
	  dq = dp * b + dq * a;
	  dp = c;
	  if (pk > ta)
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
  alpha = zero;
  beta = one;
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
  fcm = sign(w, fcl) * beta;
  fc = fcm;

  if (xlturn)
    gcl = gjwkb * beta;
  else
    gcl = fcm * gam;

  if (kfn != 0)
    gcl = -gcl;
  gc = gcl;
  gpl = gcl * (p - q / gam) - alpha * gcl;
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
      WRITE(6,2000) XX,ACCH;
 2000 FORMAT(' FOR XX = ',1PE12.3,' TRY SMALL-X  SOLUTIONS',
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
  jwkb(_Tp xx, _Tp eta1, _Tp xl, _Tp fjwkb, _Tp gjwkb)
  {
    constexpr auto zero = 0.0e0;
    constexpr auto half = 0.5e0;
    constexpr auto one = 1.0e0;
    constexpr auto six = 6.0e0;
    constexpr auto ten = 10.0e0;
    constexpr auto rl35 = 35.0e0;
    constexpr auto aloge = 0.43429'45e0;
    auto x = xx;
    auto eta = eta1;
    auto gh2 = x * (eta + eta - x);
    auto xll1 = std::max(xl * xl + xl, zero);
    if (gh2 + xll1 <= zero)
      return;
    hll = xll1 + six / rl35;
    hl = std::sqrt(hll);
    sl = eta / hl + hl / x;
    rl2 = one + eta * eta / hll;
    gh = std::sqrt(gh2 + hll) / x;
    phi = x * gh - half * (hl * std::log((gh + sl) * (gh + sl) / rl2) - std::log(gh));
    if (eta != zero)
      phi -= eta * std::atan2(x * gh, x - eta);
    phi10 = -phi * aloge;
    gjwkb = std::exp(-phi);
    fjwkb = half / (gh * gjwkb);
    return;
  }
