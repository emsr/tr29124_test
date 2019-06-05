/**
 *
 */

#include <limits>
#include <cmath>

void
jwkb(double xx, double eta1, double xl,
     double& fjwkb, double& gjwkb, int& iexp);

struct steed_info
{
  double paccq = 0.0;
  int nfp = 0;
  int npq = 0;
  int iexp = 0;
  int m1 = 0;
};

/**
 *  Compute the Coulomb wavefunction program using Steed's method.
 *
 *  A. R. Barnett           Manchester  March   1981
 *
 *  COULFG RETURNS F,G,F',G', FOR REAL XX>0,REAL ETA1 (INCLUDING 0),
 *   AND REAL LAMBDA(XLMIN) > -1 FOR INTEGER-SPACED LAMBDA VALUES
 *   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
 *   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
 *   THE DIRAC EQUATION ,ALSO SPHERICAL + CYLINDRICAL BESSEL EQUATIONS
 *
 *  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,
 *  STARTING ARRAY ELEMENT IS M1 = max(  int(XLMIN+ACCUR),0) + 1
 *      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES
 *
 *  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES
 *            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN
 *            = 3      F               CALL TO AT LEAST LENGTH (1)
 *  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED
 *            = 1 SPHERICAL   BESSEL
 *            = 2 CYLINDRICAL BESSEL
 *  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT
 *
 *  PRECISION&  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'
 *   IN OSCILLATING REGION X >= ETA1 + SQRT(ETA1**2 + XLM(XLM+1))
 *   COULFG IS CODED FOR double ON IBM OR EQUIVALENT  ACCUR = 10**-16
 *   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33
 *   FOR MANTISSAS OF 56 + 112 BITS. FOR SINGLE PRECISION CDC (48 BITS)
 *   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION
 */
void
coulfg(double xx,double eta1,double xlmin,double xlmax,
       double (&fc)[100], double (&gc)[100], double (&fcp)[100], double (&gcp)[100],
       int mode1, int kfn, int& ifail, steed_info& steed)
{
      bool etane0, xlturn;

      constexpr double zero = 0.0;
      constexpr double one = 1.0;
      constexpr double two = 2.0;
      constexpr double ten2 = 1.0e2;
      constexpr double abort = 2.0e4;
      constexpr double half = 0.5;
      constexpr double tm30 = 1.0e-30;

      // This constant is sqrt(two/pi).
      constexpr double rt2epi = 0.79788'45608'02865;

      double acch, acc4, acc;
      double x, eta, xlm, fjwkb, gjwkb, e2mm1, dell, xll;
      double xl, xi, wi, w, tk, ta, sl, rl, p, q, px, pk1, pk;
      double gcl1, gpl, gam, fpl, fcm, fcl, fcl1, f, el, ek;
      double dr, dp, dq, di, df, d, a, ar, ai, b, br, bi, c;
      double alpha, beta, gcl;
      int maxl, l, l1, lxtra, lp;

      // change accur to suit machine and precision required
      double accur = std::numeric_limits<double>::epsilon();
      int mode = 1;
      if (mode1 == 2 || mode1 == 3)
        mode = mode1;
      ifail = 0;
      steed.iexp = 1;
      steed.npq = 0;
      eta = eta1;
      gjwkb = zero;
      steed.paccq = one;
      if (kfn != 0)
        eta = zero;
      etane0 = eta != zero;
      acc = accur;
      acc4 = acc * ten2 * ten2;
      acch = std::sqrt(acc);
      // Test range of xx, exit if <= std::sqrt(accur) or if negative
      if (xx <= acch)
        goto _100;
      x = xx;
      xlm = xlmin;
      if (kfn == 2)
        xlm = xlm - half;
      if (xlm <= -one || xlmax < xlmin)
        goto _105;
      e2mm1 = eta * eta + xlm * xlm + xlm;
      xlturn = x * (x - two * eta) < xlm * xlm + xlm;
      dell = xlmax - xlmin + acc;
      if (std::abs(std::fmod(dell, one)) > acc)
        {
	  ;//write(6,2040)xlmax,xlmin,dell;
	}
      // lxtra is number of additional lambda values to be computed.
      lxtra = int(dell);
      // xll is max lambda value, or 0.5 smaller for j, y bessels.
      xll = xlm + double(lxtra);
      // determine starting array element [m1] from xlmin
      steed.m1 = std::max(int(xlmin + acc), 0) + 1;
      l1 = steed.m1 + lxtra;

      // Evaluate CF1 = f = F'(XL,ETA,X)/F(XL,ETA,X).
      xi = one / x;
      fcl = one;
      pk = xll + one;
      px = pk + abort;
      while (true)
      {
        ek = eta / pk;
        f = (ek + pk * xi) * fcl + (fcl - one) * xi;
        pk1 = pk + one;
        // Test ensures B1 != ZERO for negative eta; Fixup is exact.
        if (std::abs(eta * x + pk * pk1) > acc)
          break;
        fcl = (one + ek * ek) / (one + (eta / pk1)*(eta / pk1));
        pk = two + pk;
      }
      d = one / ((pk + pk1) * (xi + ek / pk1));
      df = -fcl * (one + ek * ek) * d;
      if (fcl != one)
        fcl = -one;
      if (d < zero)
        fcl = -fcl;
      f += df;

      // Begin CF1 loop on pk = k = lambda + 1
      p = one;
      while (std::abs(df) >= acc * std::abs(f))
      {
        pk = pk1;
        pk1 = pk1 + one;
        ek = eta / pk;
        tk = (pk + pk1) * (xi + ek / pk1);
        d = tk - d * (one + ek * ek);
        if (std::abs(d) <= acch)
        {
          //write (6,1000) d, df, acch, pk, ek, eta, x;
          p = p + one;
          if (p > two)
            goto _110;
        }
        d = one / d;
        if (d < zero)
          fcl = -fcl;
        df = df * (d * tk - one);
        f = f + df;
        if(pk > px)
          goto _110;
      }

      steed.nfp = pk - xll - 1;
      if (lxtra > 0)
      {
        // Downward recurrence to lambda = xlm. array gc,if present,stores rl
        fcl *= tm30;
        fpl = fcl * f;
        if (mode == 1)
          fcp[l1] = fpl;
        fc[l1] = fcl;
        xl = xll;
        rl = one;
        el = zero;
        for (lp = 1; lp <= lxtra; ++lp)
        {
          if (etane0)
            el = eta / xl;
          if (etane0)
            rl = std::sqrt(one + el * el);
          sl = el + xl * xi;
          l = l1 - lp;
          fcl1 = (fcl * sl + fpl) / rl;
          fpl = fcl1 * sl - fcl * rl;
          fcl = fcl1;
          fc[l] = fcl;
          if (mode == 1)
            fcp[l]  = fpl;
          if (mode != 3 && etane0)
            gc[l + 1] = rl;
          xl -= one;
        }
        if (fcl == zero)
          fcl = acc;
        f = fpl / fcl;
      }
      // Now we have reached lambda = xlmin = xlm
      // Evaluate CF2 = p + i.q  again using Steed's algorithm
      // See text for compact complex code for sp cdc or non-ANSI IBM
      if (xlturn)
        jwkb(x, eta, std::max(xlm, zero), fjwkb, gjwkb, steed.iexp);
      if (steed.iexp > 1 || gjwkb > one / (acch * ten2))
      {
        // Arrive here if G(XLM) > 10**6 or IEXP > 250 + XLTURN = true
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
        dq = xi * (ar * dr - ai * di);
        while (std::abs(dp) + std::abs(dq) >= acc * (std::abs(p) + std::abs(q)))
        {
          p += dp;
          q += dq;
          pk += two;
          ar += pk;
          ai += wi;
          bi += two;
          d = ar * dr - ai * di + br;
          di = ai * dr + ar * di + bi;
          c = one / (d * d + di * di);
          dr = c * d;
          di = -c * di;
          a = br * dr - bi * di - one;
          b = bi * dr + br * di;
          c = dp * a - dq * b;
          dq = dp * b + dq * a;
          dp = c;
          if (pk > ta)
            goto _120;
        }
        
        steed.npq = pk / two;
        steed.paccq = half * acc / std::min(std::abs(q), one);
        if (std::abs(p) > std::abs(q))
          steed.paccq = steed.paccq * std::abs(p);

        // Solve for fcm = f at lambda = xlm, then find norm factor w = w / fcm.
        gam = (f - p) / q;
        if (q <= acc4 * std::abs(p))
          goto _130;
        w = one / std::sqrt((f - p) * gam + q);
      }
      // Normalise for spherical or cylindrical Bessel functions
      alpha = zero;
      if (kfn == 1)
        alpha = xi;
      if (kfn == 2)
        alpha = xi * half;

      beta  = one;
      if (kfn == 1)
        beta  = xi;
      if (kfn == 2)
        beta = std::sqrt(xi) * rt2epi;

      fcm = std::copysign(w, fcl) * beta;
      fc[steed.m1] = fcm;
      if (mode == 3)
        goto _11;
      if (!xlturn)
        gcl = fcm * gam;
      if (xlturn)
        gcl = gjwkb * beta;
      if (kfn != 0)
        gcl = -gcl;
      gc[steed.m1] = gcl;
      gpl = gcl * (p - q / gam) - alpha * gcl;
      if (mode == 2)
        goto _11;
      gcp[steed.m1] = gpl;
      fcp[steed.m1] = fcm * (f - alpha);
 _11: if (lxtra == 0)
        return;

      // Upward recurrence from GC[m1],GCP[m1] stored value is rl
      // Renormalise FC,FCP at each lambda and correct regular derivative
      // XL = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
      w = beta * w / std::abs(fcl);
      maxl = l1 - 1;
      for (l = steed.m1; l <= maxl; ++l)
      {
        if (mode == 3)
          goto _12;
        xl += one;
        if (etane0)
          el = eta / xl;
        if (etane0)
          rl = gc[l + 1];
        sl = el + xl * xi;
        gcl1 = ((sl - alpha) * gcl - gpl) / rl;
        gpl = rl * gcl - (sl + alpha) * gcl1;
        gcl = gcl1;
        gc[l + 1] = gcl1;
        if (mode == 2)
          goto _12;
        gcp[l + 1] = gpl;
        fcp[l + 1] = w * (fcp[l + 1] - alpha * fc[l + 1]);
 _12:   fc[l + 1] = w * fc[l + 1];
      }
      return;

//_1000 format(/' cf1 accuracy loss& d,df,acch,k,eta/k,eta,x = ',
//      1p,7e9.2/)

      // Error messages

 _100: ifail = -1;
//      write(6,2000) xx, acch;
//_2000 format(' for xx = ',1pe12.3,' try small-x  solutions',
//      ' or x negative'/ ,' square root accuracy parameter =  ',e12.3/)
      return;

 _105: ifail = -2;
//      write (6,2005) xlmax, xlmin, xlm;
//_2005 format(/' problem with input order values&xlmax,xlmin,xlm = ',
//      1p,3e15.6/)
      return;

 _110: ifail =  1;
//      write (6,2010) abort, f,  df, pk, px, acc;
//_2010 format(' cf1 has failed to converge after ',f10.0,' iterations',/
//      ' f,df,pk,px,accur =  ',1p,5e12.3//)
      return;

 _120: ifail =  2;
//      write (6,2020) abort, p, q, dp, dq, acc;
//_2020 format(' cf2 has failed to converge after ',f7.0,' iterations',/
//      ' p,q,dp,dq,accur =  ',1p,4e17.7,e12.3//)
      return;

 _130: ifail =  3;
//      write (6,2030) p, q, acc, dell, lxtra, steed.m1;
//_2030 format(' final q<= std::abs(p)*acc*10**4 , p,q,acc = ',1p,3e12.3,4x,
//      ' dell,lxtra,steed.m1 = ',e12.3,2i5 /)
//_2040 format(' xlmax - xlmin = dell not an integer ',1p,3e20.10/)

      return;
}


// Computes JWKB approximations to Coulomb functions for xl >= 0
// as modified by Biedenharn et al. Phys Rev 97 (1955) 542-554
void
jwkb(double xx, double eta1, double xl,
     double& fjwkb, double& gjwkb, int& iexp)
{
  constexpr double zero = 0.0;
  constexpr double half = 0.5;
  constexpr double one = 1.0;
  constexpr double six = 6.0;
  constexpr double ten = 10.0;
  constexpr double rl35 = 35.0;
  constexpr double aloge = 0.43429'45;

  auto x = xx;
  auto eta = eta1;
  auto gh2 = x * (eta + eta - x);
  auto xll1 = std::max(xl * xl + xl, zero);
  if (gh2 + xll1 <= zero)
    return;
  auto hll = xll1 + six / rl35;
  auto hl = std::sqrt(hll);
  auto sl = eta / hl + hl / x;
  auto rl2 = one + eta * eta / hll;
  auto gh = std::sqrt(gh2 + hll) / x;
  auto phi = x * gh - half * (hl * std::log((gh + sl)*(gh + sl) / rl2) - std::log(gh));
  if (eta != zero)
    phi -= eta * std::atan2(x * gh, x - eta);
  auto phi10 = -phi * aloge;
  iexp = int(phi10);
  if (iexp > 250)
    gjwkb = std::pow(ten, phi10 - double(iexp));
  if (iexp <= 250)
    gjwkb = std::exp(-phi);
  if (iexp <= 250)
    iexp = 0;
  fjwkb = half / (gh * gjwkb);

  return;
}
