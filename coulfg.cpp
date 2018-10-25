/*
g++ -std=c++17 coulfg.cpp 
*/

#include <limits>
#include <cmath>

void
JWKB(double XX, double ETA1, double XL,
     double& FJWKB, double& GJWKB, int& IEXP);

/**
 *  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD
 *
 *  A. R. BARNETT           MANCHESTER  MARCH   1981
 *
 *  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395
 *                 + 'RCWFF'      IN    CPC 11 (1976) 141-142
 *  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
 *  THIS VERSION WRITTEN UP       IN    CPC XX (1982) YYY-ZZZ
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
coulfg(double XX,double ETA1,double XLMIN,double XLMAX,
       double (&FC)[100], double (&GC)[100], double (&FCP)[100], double (&GCP)[100],
       int MODE1, int KFN, int& IFAIL)
{
      bool ETANE0, XLTURN;

      // Common block is for information only.  Not required in code
      int NFP, NPQ, IEXP, M1;
      double PACCQ;
      //COMMON /STEED/ PACCQ, NFP, NPQ, IEXP, M1

      constexpr double ZERO = 0.0e0;
      constexpr double ONE = 1.0e0;
      constexpr double TWO = 2.0e0;
      constexpr double TEN2 = 1.0e2;
      constexpr double ABORT = 2.0e4;
      constexpr double HALF = 0.5e0;
      constexpr double TM30 = 1.0e-30;

      constexpr double RT2EPI = 0.79788'45608'02865e0;
      // THIS CONSTANT IS  DSQRT(TWO/PI)&  USE Q0 FOR IBM REAL*16& D0 FOR
      //  double + CDC DOUBLE P&  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.

      double ACCH, ACC4, ACC;
      double X, ETA, XLM, FJWKB, GJWKB, E2MM1, DELL, XLL;
      double XL, XI, WI, W, TK, TA, SL, RL, P, Q, PX, PK1, PK;
      double GCL1, GPL, GAM, FPL, FCM, FCL, FCL1, F, EL, EK;
      double DR, DP, DQ, DI, DF, D, A, AR, AI, B, BR, BI, C;
      double ALPHA, BETA, GCL;
      int MAXL, L, L1, LXTRA, LP;

      // CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
      double ACCUR = std::numeric_limits<double>::epsilon();
      int MODE = 1;
      if (MODE1 == 2 || MODE1 == 3)
        MODE = MODE1;
      IFAIL = 0;
      IEXP = 1;
      NPQ = 0;
      ETA = ETA1;
      GJWKB = ZERO;
      PACCQ = ONE;
      if (KFN != 0)
        ETA = ZERO;
      ETANE0 = ETA != ZERO;
      ACC = ACCUR;
      ACC4 = ACC * TEN2 * TEN2;
      ACCH = std::sqrt(ACC);
      // Test range of xx, exit if <= std::sqrt(accur) or if negative
      if (XX <= ACCH)
        goto _100;
      X = XX;
      XLM = XLMIN;
      if (KFN == 2)
        XLM = XLM - HALF;
      if (XLM <= -ONE || XLMAX < XLMIN)
        goto _105;
      E2MM1 = ETA * ETA + XLM * XLM + XLM;
      XLTURN = X * (X - TWO * ETA) < XLM * XLM + XLM;
      DELL = XLMAX - XLMIN + ACC;
      if (std::abs(std::fmod(DELL, ONE)) > ACC)
        ;//WRITE(6,2040)XLMAX,XLMIN,DELL;
      // LXTRA is number of additional lambda values to be computed.
      LXTRA = int(DELL);
      // XLL is max lambda value, or 0.5 smaller for J, Y Bessels.
      XLL = XLM + double(LXTRA);
      // Determine starting array element (m1) from xlmin
      M1 = std::max(int(XLMIN + ACC), 0) + 1;
      L1 = M1 + LXTRA;

      // Evaluate CF1 = f = F'(XL,ETA,X)/F(XL,ETA,X).
      XI = ONE / X;
      FCL = ONE;
      PK = XLL + ONE;
      PX = PK + ABORT;
      while (true)
      {
        EK = ETA / PK;
        F = (EK + PK * XI) * FCL + (FCL - ONE) * XI;
        PK1 = PK + ONE;
        // Test ensures B1 != ZERO for negative eta; Fixup is exact.
        if (std::abs(ETA * X + PK * PK1) > ACC)
          break;
        FCL = (ONE + EK * EK) / (ONE + (ETA / PK1)*(ETA / PK1));
        PK = TWO + PK;
      }
      D = ONE / ((PK + PK1) * (XI + EK / PK1));
      DF = -FCL * (ONE + EK * EK) * D;
      if (FCL != ONE)
        FCL = -ONE;
      if (D < ZERO)
        FCL = -FCL;
      F = F + DF;

      // Begin CF1 loop on pk = k = lambda + 1
      P = ONE;
      while (std::abs(DF) >= ACC * std::abs(F))
      {
        PK = PK1;
        PK1 = PK1 + ONE;
        EK = ETA / PK;
        TK = (PK + PK1) * (XI + EK / PK1);
        D = TK - D * (ONE + EK * EK);
        if (std::abs(D) <= ACCH)
        {
          //WRITE (6,1000) D, DF, ACCH, PK, EK, ETA, X;
          P = P + ONE;
          if (P > TWO)
            goto _110;
        }
        D = ONE / D;
        if (D < ZERO)
          FCL = -FCL;
        DF = DF * (D * TK - ONE);
        F = F + DF;
        if(PK > PX)
          goto _110;
      }

      NFP = PK - XLL - 1;
      if (LXTRA > 0)
      {
        // Downward recurrence to lambda = xlm. array gc,if present,stores rl
        FCL *= TM30;
        FPL = FCL * F;
        if (MODE == 1)
          FCP[L1] = FPL;
        FC[L1] = FCL;
        XL = XLL;
        RL = ONE;
        EL = ZERO;
        for (LP = 1; LP <= LXTRA; ++LP)
        {
          if (ETANE0)
            EL = ETA / XL;
          if (ETANE0)
            RL = std::sqrt(ONE + EL * EL);
          SL = EL + XL * XI;
          L = L1 - LP;
          FCL1 = (FCL * SL + FPL) / RL;
          FPL = FCL1 * SL - FCL * RL;
          FCL = FCL1;
          FC[L] = FCL;
          if (MODE == 1)
            FCP[L]  = FPL;
          if (MODE != 3 && ETANE0)
            GC[L + 1] = RL;
          XL -= ONE;
        }
        if (FCL == ZERO)
          FCL = ACC;
        F = FPL / FCL;
      }
      // Now we have reached lambda = xlmin = xlm
      // Evaluate CF2 = p + i.q  again using Steed's algorithm
      // See text for compact complex code for sp cdc or non-ANSI IBM
      if (XLTURN)
        JWKB(X, ETA, std::max(XLM, ZERO), FJWKB, GJWKB, IEXP);
      if (IEXP > 1 || GJWKB > ONE / (ACCH * TEN2))
      {
        // Arrive here if G(XLM) > 10**6 or IEXP > 250 + XLTURN = true
        W = FJWKB;
        GAM = GJWKB * W;
        P = F;
        Q = ONE;
      }
      else
      {
        XLTURN = false;
        TA = TWO * ABORT;
        PK = ZERO;
        WI = ETA + ETA;
        P = ZERO;
        Q = ONE - ETA * XI;
        AR = -E2MM1;
        AI = ETA;
        BR = TWO * (X - ETA);
        BI = TWO;
        DR = BR / (BR * BR + BI * BI);
        DI = -BI / (BR * BR + BI * BI);
        DP = -XI * (AR * DI + AI * DR);
        DQ = XI * (AR * DR - AI * DI);
        while (std::abs(DP) + std::abs(DQ) >= ACC * (std::abs(P) + std::abs(Q)))
        {
          P = P + DP;
          Q = Q + DQ;
          PK = PK + TWO;
          AR = AR + PK;
          AI = AI + WI;
          BI = BI + TWO;
          D = AR * DR - AI * DI + BR;
          DI = AI * DR + AR * DI + BI;
          C = ONE / (D * D + DI * DI);
          DR = C * D;
          DI = -C * DI;
          A = BR * DR - BI * DI - ONE;
          B = BI * DR + BR * DI;
          C = DP * A - DQ * B;
          DQ = DP * B + DQ * A;
          DP = C;
          if (PK > TA)
            goto _120;
        }
        
        NPQ = PK / TWO;
        PACCQ = HALF * ACC / std::min(std::abs(Q), ONE);
        if (std::abs(P) > std::abs(Q))
          PACCQ = PACCQ * std::abs(P);

        // Solve for fcm = f at lambda = xlm, then find norm factor w = w / fcm.
        GAM = (F - P) / Q;
        if (Q <= ACC4 * std::abs(P))
          goto _130;
        W = ONE / std::sqrt((F - P) * GAM + Q);
      }
      // Normalise for spherical or cylindrical Bessel functions
      ALPHA = ZERO;
      if (KFN == 1)
        ALPHA = XI;
      if (KFN == 2)
        ALPHA = XI * HALF;

      BETA  = ONE;
      if (KFN == 1)
        BETA  = XI;
      if (KFN == 2)
        BETA = std::sqrt(XI) * RT2EPI;

      FCM = std::copysign(W, FCL) * BETA;
      FC[M1] = FCM;
      if (MODE == 3)
        goto _11;
      if (!XLTURN)
        GCL = FCM * GAM;
      if (XLTURN)
        GCL = GJWKB * BETA;
      if (KFN != 0)
        GCL = -GCL;
      GC[M1] = GCL;
      GPL = GCL * (P - Q / GAM) - ALPHA * GCL;
      if (MODE == 2)
        goto _11;
      GCP[M1] = GPL;
      FCP[M1] = FCM * (F - ALPHA);
 _11: if (LXTRA == 0)
        return;

      // Upward recurrence from GC[m1],GCP[m1] stored value is rl
      // Renormalise FC,FCP at each lambda and correct regular derivative
      // XL = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
      W = BETA * W / std::abs(FCL);
      MAXL = L1 - 1;
      for (L = M1; L <= MAXL; ++L)
      {
        if (MODE == 3)
          goto _12;
        XL += ONE;
        if (ETANE0)
          EL = ETA / XL;
        if (ETANE0)
          RL = GC[L + 1];
        SL = EL + XL * XI;
        GCL1 = ((SL - ALPHA) * GCL - GPL) / RL;
        GPL = RL * GCL - (SL + ALPHA) * GCL1;
        GCL = GCL1;
        GC[L + 1] = GCL1;
        if (MODE == 2)
          goto _12;
        GCP[L + 1] = GPL;
        FCP[L + 1] = W * (FCP[L + 1] - ALPHA * FC[L + 1]);
 _12:   FC[L + 1] = W * FC[L + 1];
      }
      return;

//_1000 FORMAT(/' CF1 ACCURACY LOSS& D,DF,ACCH,K,ETA/K,ETA,X = ',
//      1P,7E9.2/)

      // Error messages

 _100: IFAIL = -1;
//      WRITE(6,2000) XX,ACCH;
//_2000 FORMAT(' FOR XX = ',1PE12.3,' TRY SMALL-X  SOLUTIONS',
//      ' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',E12.3/)
      return;

 _105: IFAIL = -2;
//      WRITE (6,2005) XLMAX,XLMIN,XLM;
//_2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES&XLMAX,XLMIN,XLM = ',
//      1P,3E15.6/)
      return;

 _110: IFAIL =  1;
//      WRITE (6,2010) ABORT, F,  DF, PK, PX, ACC;
//_2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
//      ' F,DF,PK,PX,ACCUR =  ',1P,5E12.3//)
      return;

 _120: IFAIL =  2;
//      WRITE (6,2020) ABORT,P,Q,DP,DQ,ACC;
//_2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
//      ' P,Q,DP,DQ,ACCUR =  ',1P,4E17.7,E12.3//)
      return;

 _130: IFAIL =  3;
//      WRITE (6,2030) P,Q,ACC,DELL,LXTRA,M1;
//_2030 FORMAT(' FINAL Q<= std::abs(P)*ACC*10**4 , P,Q,ACC = ',1P,3E12.3,4X,
//      ' DELL,LXTRA,M1 = ',E12.3,2I5 /)
//_2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P,3E20.10/)

      return;
}


// COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS FOR XL >= 0
// AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
void
JWKB(double XX, double ETA1, double XL,
     double& FJWKB, double& GJWKB, int& IEXP)
{
  constexpr double ZERO = 0.0;
  constexpr double HALF = 0.5;
  constexpr double ONE = 1.0;
  constexpr double SIX = 6.0;
  constexpr double TEN = 10.0;
  constexpr double RL35 = 35.0;
  constexpr double ALOGE = 0.43429'45;

  auto X = XX;
  auto ETA = ETA1;
  auto GH2 = X * (ETA + ETA - X);
  auto XLL1 = std::max(XL * XL + XL, ZERO);
  if (GH2 + XLL1 <= ZERO)
    return;
  auto HLL = XLL1 + SIX / RL35;
  auto HL = std::sqrt(HLL);
  auto SL = ETA / HL + HL / X;
  auto RL2 = ONE + ETA * ETA / HLL;
  auto GH = std::sqrt(GH2 + HLL) / X;
  auto PHI = X * GH - HALF * (HL * std::log((GH + SL)*(GH + SL) / RL2) - std::log(GH));
  if (ETA != ZERO)
    PHI -= ETA * std::atan2(X * GH, X - ETA);
  auto PHI10 = -PHI * ALOGE;
  IEXP = int(PHI10);
  if (IEXP > 250)
    GJWKB = std::pow(TEN, PHI10 - double(IEXP));
  if (IEXP <= 250)
    GJWKB = std::exp(-PHI);
  if (IEXP <= 250)
    IEXP = 0;
  FJWKB = HALF / (GH * GJWKB);

  return;
}
