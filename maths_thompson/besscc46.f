      SUBROUTINE BESSCC(ZZ,XNU,NL, FI,FK,FIP,FKP, MODE1,ACC,IFAIL)
C BESSCC46
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  COMPLEX I,K BESSEL FUNCTIONS PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  I. J. Thompson          Bristol     JULY    1986                    C
C                                                                      C
C  original program  RCWFN       in    CPC  8 (1974) 377-395           C
C                 +  RCWFF       in    CPC 11 (1976) 141-142           C
C                 +  COULFG      in    CPC 27 (1982) 147-166           C
C                 +  COULCC      in    CPC 36 (1985) 363-372           C
C  description of real algorithm in    CPC 21 (1981) 297-314           C
C  description of complex algorithm    JCP 64 (1986) 490-509           C
C  this version written up       in    CPC 47 (1987) 245-257           C
C                                                                      C
C  BESSCC returns I,K,I',K' for complex Z, real XNU, & integer NL > 0  C
C   for integer-spaced orders  XNU  to  XMU = XNU + NL - 1             C
C  The first order XNU must be > -0.5, and > 2/ln(|Z|) if |Z| << 1     C
C              (i.e. not too negative if |Z| << 1)                     C
C                                                                      C
C  if |MODE1|= 1  get I,K,I',K'   for integer-spaced lambda values     C
C            = 2      I,K      unused arrays must be dimensioned in    C
C            = 3      I,  I'          call to at least length (1)      C
C            = 4      I                                                C
C                                                                      C
C     if MODE1<0 then the values returned are scaled by an exponential C
C                factor (dependent only on Z) to bring nearer unity    C
C                the functions for large |Z|, small |XN| < |Z|         C
C     So define SCALE = (  0        if MODE1 > 0                       C
C                       (  REAL(Z)  if MODE1 < 0                       C
C        then FI = EXP(-ABS(SCALE)) * I                                C
C             FIP= EXP(-ABS(SCALE)) * I'                               C
C         and FK = EXP(SCALE)       * K                                C
C             FKP= EXP(SCALE)       * K'                               C
C                                                                      C
C  Precision:  results to within 1-2 decimals of 'machine accuracy',   C
C              depending on the value of ACC in the calling sequence.  C
C   If ACC is too small or too large, a default ACCDEF is used.        C
C       The accuracy target actually used is returned in COMMON block  C
C         COMMON   /BSTEED/ ACCUR,NFP,NPQ(2),KASE                      C
C                                                                      C
C   BESSCC is coded for REAL*8 on IBM or equivalent  ACC > 2D-16       C
C                                                                      C
C   Use IMPLICIT COMPLEX*32 & REAL*16 on VS compiler   ACC > 1Q-24     C
C     (More GAM & CSC coefficients can be provided for ACC = 1Q-31)    C
C   For single precision CDC, CRAY etc reassign REAL*8=REAL etc.       C
C                                                                      C
C   IFAIL  in output:  -2 = argument out of range                      C
C                      -1 = one of the continued fractions failed,     C
C                           or arithmetic check before final recursion C
C                       0 = All Calculations satisfactory              C
C                      ge 0 : results available for orders up to & at  C
C                             position NL-IFAIL in the output arrays.  C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     Machine dependent parameters :                                   C
C                                                                      C
C     ACC      convergence criterion for continued fractions.          C
C              Except near zero-crossings of the functions,            C
C              the relative errors of the returned functions should be C
C              less than max(ACC, (100 + imag(Z))*unit-roundoff )      C
C                                                                      C
C     FPMAX    magnitude of largest floating point number * ACC        C
C     FPMIN    magnitude of smallest floating point number / ACC       C
C     FPLMIN   log(FPMIN)                                              C
C     FPHMIN   sqrt(FPMIN)                                             C
C     LIMIT    max. no. iterations for CF1, CF2 continued fractions    C
C              (If XMU > 0.35*|Z|,  then |Z| is limited to LIMIT).     C
C                                                                      C
C     GAM, CSC are the coefficients of the continued fraction form     C
C              of the diagonal Pade approximants for                   C
C              ln(Gamma(1+nu))/nu and  1/sin(pi.nu) - 1/(pi.nu) resp.  C
C              The number of terms required is a linear function of    C
C                  the number of digits accuracy, i.e. of log(ACC)  .  C
C      The given GAM & CSC parameters are sufficient for ACC > 1D-24   C
C      For fixed accuracy worse than this, expressions in the code     C
C      involving ACC may be pre-evaluated, and NGAM & NCSC reduced.    C
C                                                                      C
C     Associated routine  :   CF2E (appended)                          C
C                                                                      C
C     Intrinsic functions :   MIN, MAX, SQRT, REAL, IMAG, LOG, EXP,    C
C      (Generic names)        ABS, MOD, SIN, COS, INT                  C
C               Complex   :   LOG, EXP, SIN, SQRT, DCMPLX              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      PARAMETER(NGAM=33, NCSC=12)
      DIMENSION FI(NL),FK(NL),FIP(NL),FKP(NL)
      LOGICAL IFIP
      REAL*8 ACCUR,ACC,ACCDEF,ACCLOG,GAM(NGAM),CSC(NCSC),
     X       ZERO,ONE,TWO,HALF,PI,PI1,PI2,FPMAX,FPMIN,FPLMIN,FPHMIN,
     X       SCALE,ABSZ,XNU,XMU,XN0,XLL,YL0,C1,C2,G,G1,XM1,X2M,XX,K,EK,
     X       REAL,AIMAG,ABSC
      COMPLEX*16 TIMESI,CSINH,CTR,ZZ,Z,X,ZI,ALPHA,PQ1,PQ2,F2,
     X           FK0,FKP0,RAT,FIL,FIPL,SL,SL1,F1,D0,G2,SUM(2)
      PARAMETER(PI1 = 3217D+0/1024D+0)
C
      COMMON       /BSTEED/ ACCUR,NFP,NPQ(2),KASE
C***  common block is for information & storage only.
C     (it is not essential to working of the code)
C
      DATA ZERO,HALF,ONE,TWO,LIMIT /0.D+0, 0.5D+0, 1.D+0, 2.D+0, 20000/,
     X     FPMAX,FPMIN / 1D+60, 1D-60  /, FPHMIN,FPLMIN /1D-30,-140./,
     X     ACCDEF / 1D-6 /
      DATA GAM /
     X -5.772156649015328606065121D-01, 1.424886889659201737396080D+00,
     X -9.377115767249094049120554D-01,-9.773477272948880533901762D-02,
     X  6.130615858357687237828236D-01, 1.035548353597846452535594D-01,
     X  3.833453676357374951071144D-01, 1.676098533180169272436803D-01,
     X  3.423153070879739741315024D-01, 1.872096384664854899310796D-01,
     X  3.044176452057531884181363D-01, 2.038057743675276230967228D-01,
     X  3.036922239952957869745095D-01, 2.111473436553448271828676D-01,
     X  2.815492532986102124308247D-01, 2.186620399336151796598028D-01,
     X  2.887599926903399843284945D-01, 2.211365686356764978140098D-01,
     X  2.711931229606902430385447D-01, 2.278516401771275692711102D-01,
     X  2.798289016836091400090250D-01, 2.257719734937412432938853D-01,
     X  2.670653349467782427924255D-01, 2.338192475772592332514537D-01,
     X  2.720797129387430981577150D-01, 2.297304315464212589631377D-01,
     X  2.662328902179656243145202D-01, 2.355772901124245270936971D-01,
     X  2.665129963814811071666658D-01, 2.352423215988759927243035D-01,
     X  2.639101271705042366932099D-01, 2.349813957292258269500878D-01,
     X  2.659470304082729651745135D-01  /
      DATA CSC, PI2 /
     X  1.666666666666666666666667D-01,-1.166666666666666666666667D-01,
     X  1.122448979591836734693878D-02,-2.839620696763553906411049D-02,
     X  3.881831014317402702157693D-03,-1.256690168285914562406590D-02,
     X  1.959015840039126317490664D-03,-7.058674448331350123193777D-03,
     X  1.179797436463512584131754D-03,-4.514563690718488905367012D-03,
     X  7.880008568042515532566820D-04,-3.133992819253743972358894D-03,
     X -8.908910206761537356616720D-06 /
      CMPLX(C1,C2) = DCMPLX(C1,C2)
      AIMAG(X) = DIMAG(X)
      REAL(X) = DBLE(X)
      TIMESI(X) = CMPLX(-AIMAG(X),REAL(X))
      CTR(X,G) = CMPLX(REAL(X)*G,AIMAG(X)*G)
      ABSC(X) = ABS(REAL(X)) + ABS(AIMAG(X))
      CSINH(X) = - TIMESI(SIN(TIMESI(X)))
C
      MODE = ABS(MODE1)
      IFIP = MOD(MODE,2).EQ.1
      IFAIL = -2
      NFP   = 0
      KASE     = 0
      NPQ(1)   = 0
      NPQ(2)   = 0
      PI = PI1 + PI2
      ACCUR = ACC
      IF(ACC.GT.1D-4 .OR. ONE+ACC.EQ.ONE) ACCUR = ACCDEF
      ACCLOG = LOG(ACCUR)
C
      K = ZERO
      EK= ONE
      IF(REAL(ZZ).LT.ZERO.or.
     x   REAL(ZZ).eq.ZERO.and.AIMAG(ZZ).lt.ZERO) THEN
         EK = - ONE
         K  = SIGN(ONE,AIMAG(ZZ))
      ENDIF
      Z = CTR(ZZ,EK)
      X = TIMESI(Z)
      SCALE = ZERO
      IF(MODE1.LT.0) SCALE = REAL(Z)
C
      ABSZ = ABS(Z)
      IF(XNU.LE.-HALF .OR. ABSZ.GT.ONE/FPHMIN)  GO TO 180
      ZI  = ONE/Z
      IFAIL = -1
C
C ***       XNU  is initial nu value,
C ***       XMU  is final nu value,
C ***       XN0  is nu value near zero, and
C ***       XLL  is >> XMU, such that XLL-XN0 = 2m is even & > 0.
C          (The Y.. values are the X.. values + HALF, ie Coulomb Ls)
C
         XMU = XNU + NL - 1
         IST = - MAX(INT(XNU+0.4999),0)
           IF(ABSZ.LT.1D-3.AND. (XNU+IST)*LOG(ABSZ).GT.2.) IST = IST + 1
           IF(IST.GT.0) GO TO 180
         XN0 = XNU + IST
         YL0 = XN0 - HALF
         L0 = 1 + IST
C
C     Choose KASE of solution required :
C
      IF(ABSZ .LT. 3.0 .OR.
     X   ABS(REAL(Z)).LT.2.0 .AND. ABS(AIMAG(Z)).LT.20.) THEN
            KASE = 3
      ELSE IF(ABSZ .GT. 1.5*XMU  .AND. ABSZ.GT.25.) THEN
            KASE = 1
         ELSE
            KASE = 2
      ENDIF
C
C
C *** Normalisations for modified cylindrical Bessel functions
C
          ALPHA = CTR(ZI,HALF)
          FIL  = SQRT(CTR(ALPHA,PI))
          IF(REAL(FIL).LT.ZERO) FIL  = - FIL
C-----------------------------------------------------------------------
C
C  ***  Start Calculation with K(xn0)  (for xn0 near zero) :
C
      IF(KASE.LE.2) THEN
C
      PQ1 = CF2E(X,YL0,ONE,ACCUR,LIMIT,.TRUE.,F1,NPQ(1))
            C1 = SCALE - REAL(Z)
            G  = TWO * AINT( AIMAG(Z) / (PI+PI) )
            C2 = (AIMAG(Z) - G * PI1) - G *  PI2
         IF(C1 .LT. FPLMIN) GO TO 180
            G2 =  CTR(FIL * CMPLX(COS(C2),-SIN(C2)) , EXP(C1))
         FK0 = G2 * F1
         FKP0= (TIMESI(PQ1) - ALPHA) * FK0
C
         IF(KASE.EQ.1) THEN
C                           Find I(xn0) if OK to recur it upwards to XMU
C
           PQ2 = CF2E(X,YL0,-ONE,ACCUR,LIMIT,.FALSE.,F2,NPQ(2))
           F2 = CMPLX(ZERO,TWO) / (F1 * (PQ1-PQ2))
            SL = CMPLX(-TWO*REAL(Z),+YL0*PI-C2-C2)
              RAT = ZERO
            IF(REAL(SL).GT.ACCLOG-TWO) RAT = EXP(SL) * F1/F2
            D0  = ALPHA * F2 / G2
            SL1 = ONE - RAT
            FIL = D0 * SL1
            FIPL= (TIMESI(PQ2 - RAT*PQ1) - SL1*ALPHA) * D0
            GO TO 140
         ENDIF
      ENDIF
C
C     For KASE = 2 or 3, XMU too large to recur I up to it,
C                        so find XLL >> XMU, and recur downwards -stable
C                        NB. in these KASES until stmt 140, FIPL = I'/I
C
  10     I = MAX(NL+1,3+IST,INT(ABSZ))
         FIL = ZERO
         FIPL = ONE
            SL = CTR(ZI,TWO)
            G2  = CTR(SL,XNU+I-1)
            G = 10./ACCUR
            IF(KASE.EQ.2) G = SQRT(G)
         DO 15 NL1=I,LIMIT
            F1 = FIL - FIPL * G2
            IF(ABSC(F1) .GT. G) GO TO 16
               FIL = FIPL
               FIPL = F1
  15           G2 = G2 + SL
               GO TO 180
  16     NFP = NL1 - NL
         M = (NL1 - IST)/2
         IF(NL1-IST .EQ. M*2) NL1 = NL1 + 1
         XLL = XNU + NL1-1
      FIL = ONE
      FIPL = CTR(ZI,XLL) + CTR(Z,HALF/(XLL+ONE))
C
      SUM(1) = ONE
      SUM(2) = ONE
      C1 = ONE
      C2 = ONE
      YL0 = FPMAX*MIN(ONE,ABSZ**(XN0+ONE))
C
C *** downward recurrence  of I and I'/I to nu = XN0
C
      DO 25  L  = NL1-1,L0,-1
         SL = CTR(ZI,XNU+L)
         SL1 = SL - ZI
            IF(ABSC(FIL).GT.YL0) THEN
C                                       renormalise here & previously
               FIL = FIL * FPHMIN
               IF(KASE.EQ.3) THEN
                  SUM(1)= SUM(1)* FPHMIN
                  SUM(2)= SUM(2)* FPHMIN
                ENDIF
               DO 20 I=MAX(L+1,1),NL
               IF(ABSC(FI(I)).LT.FPHMIN) FI(I) = ZERO
  20            FI(I) = CTR(FI(I), FPHMIN)
            ENDIF
         D0 = SL + FIPL
         FIL = FIL * D0
         FIPL = SL1 + ONE/D0
           IF(L.GE.1.AND.L.LE.NL) THEN
            FI(L) =  FIL
            IF(IFIP) FIP(L)  = FIPL
           ENDIF
         IF(KASE.EQ.3 .AND. MOD(L-L0,2).EQ.0) THEN
                 XM1 = XN0 + (M-1)
                 X2M = XN0 + (M+M)
                 XX  = X2M
                 IF(M.GT.1) XX = XM1 * X2M / (M * (X2M-TWO))
                 C1 = - C1 / XX
             SUM(1) = SUM(1) + CTR(FIL,C1)
          IF(L.GE.L0+2 .AND. MODE.LE.2) THEN
                 C2 = C2 * (M-XN0) / (XX * (XN0+XM1))
             SUM(2) = SUM(2) + CTR(FIL,C2)
             ENDIF
          M = M - 1
          ENDIF
   25 CONTINUE
      IF(FIL .EQ. ZERO) FIL = + ACCUR
C
      IF(KASE.EQ.2) THEN
C                       Normalise I(nu=XN0) by Wronskian with K(nu)
C                       by determining RAT = I(XN0) / FIL
         RAT = ZI / ((FIPL*FK0 - FKP0) * FIL)
C
       ELSE
C                      Normalise I and K by I sums & Von Neumann series:
C          Begin by calculating G(xnu) = log(Gamma(1+xnu))/xnu  at XN0
      XM1 = -ACCLOG * .4343
      L   =  MIN(NGAM, INT(1.1+1.3*XM1)+1)
      G = GAM(L)
      DO 50 I=L-1,1,-1
  50  G = GAM(I) / (ONE + XN0*G)
C
C          Next calculate the first coefficient for the I sum :
      G2 = G + LOG(ZI+ZI)
      F2 = G2 * XN0
          IF(REAL(F2) .LT. FPLMIN) GO TO 180
      SL = EXP(F2)
      RAT = (C1 * EXP(-ABS(SCALE))) / (SL * SUM(1))
C
      IF(MODE.LE.2) THEN
C                             second coefficient for the K sum :
         SL1 = SL*SL * ((XN0+TWO)/(ONE-XN0))
C                             first coefficient for K sum
C                        G1 = - pi*half * (1/sin(pi*xn0) - 1/(pi*xn0))
            XX = PI * XN0
            L  = MIN(NCSC, INT(2+0.45*XM1)+1)
            G1 = CSC(L)
            DO 60 I=L-1,1,-1
  60        G1 = CSC(I) / (ONE + (XX*XX)*G1)
            G1 = - PI*HALF * XX * G1
C
            D0 = G2
            IF(XN0.NE.ZERO) D0 = G1 + SL * CSINH(F2)/XN0
         FK0 = (D0*(RAT*FIL) + (RAT*SL1)*SUM(2)*(ONE/C2))
     X             * EXP(SCALE+ABS(SCALE))
          IF(ABSC(FIPL).LT.30. .OR. ABSZ.LT.HALF) THEN
            FKP0 = FIPL * FK0 - ZI/(RAT*FIL)
          ELSE
C             Calculate K' independently of Wronskian I'K - K'I=1/z,
C                       as I(xn0) is near a zero.   Use CF2 to give K'/K
            PQ1 = CF2E(X,XN0-HALF,ONE,ACCUR,LIMIT,.FALSE.,F1,NPQ(1))
            FKP0 = (TIMESI(PQ1) - ALPHA) * FK0
            KASE = 4
          ENDIF
C
         ENDIF
      ENDIF
         YL0 = ABSC(RAT)
         RAT = RAT/YL0
C
C *** Upward recurrence of K from FK0,FKP0
C *** Upward recurrence of I from FIL,FIPL  if KASE=1
C ***   or   Renormalise FI,FIP at each nu  if KASE=2,3
C
  140 IFAIL = NL
         G1 = XN0 * K * PI
         F2= CMPLX(COS(G1),SIN(G1))
         G2 = CMPLX(ZERO,K*PI)
            G = EXP(MAX(-SCALE,FPLMIN))
            C2 = FPMIN / G
         C1 = ABSC(FIL) * 0.2
         J = 0
      SL = (XN0-ONE) * ZI
      DO 170 L = L0,NL
         SL1 = SL + ZI
         IF(KASE.GE.2) THEN
            IF(L.LT.1) GO TO 150
                     FIL = CTR( RAT* FI(L) ,YL0)
            IF(IFIP) FIPL= FIL*FIP(L)
C                    KASE=1 : I & I' at XN0 given, so recur upwards:
           ELSE IF(L.GT.L0) THEN
               F1 = FIPL - SL * FIL
               FIPL = FIL  - SL1*F1
               FIL  = F1
                  IF(ABSC(FIL) .LT. C1) THEN
                      J = J + 1
                  ELSE
                      J = 0
                  ENDIF
                  IF(J.GT.3) THEN
C                                LOSING ACCURACY in KASE = 1 !!!!!1
C                                Redo KASE=2 from large nu=XMU !
                      L0 = L-1
                      IST= L0-1
                      XN0 = XNU + IST
                      KASE = 2
                      GO TO 10
                   ENDIF
          ENDIF
                      IF(ABSC(FIL).LT.FPMIN) GO TO 180
          IF(L.GE.1) THEN
                     FI(L)  =     FIL  * F2
            IF(IFIP) FIP(L) = CTR(FIPL * F2,EK)
           ENDIF
  150  IF(MODE .GE. 3) GO TO 160
       IF(L.GT.L0) THEN
           F1 = SL * FK0 - FKP0
                         IF(ABSC(F1).GT.FPMAX*MIN(ABSZ,ONE)) GO TO 180
           FKP0 = - ( FK0 + SL1 * F1)
           FK0      = F1
        ENDIF
       IF(L.GE.1.AND.EK.GT.ZERO) THEN
                       FK(L)  = FK0
         IF(MODE.EQ.1) FKP(L) = FKP0
       ELSE IF(L.GE.1.AND.EK.LT.ZERO) THEN
         D0 = ZERO
         IF(G*ABSC(FK0).GT.C2) D0 = CMPLX(G*REAL(F2),-G*AIMAG(F2))
                       FK(L)  =    D0 * (G * FK0)  - G2 * FIL
         IF(MODE.EQ.1) FKP(L) = - (D0 * (G * FKP0) - G2 * FIPL)
       ENDIF
  160  IFAIL = MIN(IFAIL, NL - L)
      F2 = CTR(F2,EK)
  170  SL = SL1
C
  180 IF(IFAIL.GE.NL) IFAIL = -1
      RETURN
      END
      FUNCTION CF2E(X,XL,OMEGA,EPS,LIMIT,NORM,F20,N)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      LOGICAL NORM
      REAL*8 XL,AA,EPS,ERR,ASN,OMEGA,AIMAG,ABSC,ZERO,ONE,TWO
      DATA ZERO,ONE,TWO / 0D+0, 1D+0, 2D+0 /
      CMPLX(AA,ASN) = DCMPLX(AA,ASN)
      AIMAG(W) = IMAG(W)
      ABSC(W) = ABS(REAL(W)) + ABS(AIMAG(W))
C
C                                    (omega)        (omega)
C *** Evaluate  CF2  = p + PM.q  =  H   (0.0,X)' / H   (0.0,X)
C                                    XL             XL
C         where PM = omega.i
C
C     and evaluate H  (0.0,X) itself using the Steed version of Temme's
C                   XL                     algorithm.
C
C
      PQ = X
      AA = - XL * (XL + ONE)
      PM = CMPLX(ZERO,OMEGA)
      BB = TWO*(PQ + PM)
      DD = ONE/BB
      DL = AA * DD
           EM1 = ZERO
           EN  = PM
           N   = 1
           QQ = EN
           SN = ONE + QQ * DL
C          ASM = ABSC(SN)
   10 PQ    = PQ + DL
         N = N + 1
              EM2 = EM1
              EM1 = EN
            EN = CMPLX(ZERO,-AIMAG(PM)/N) * BB * EM1
     X                    - EM2 * (AA/((N-1)*N))
         AA = AA + (N+N-2)
         BB = BB + (PM+PM)
         DD = ONE/(AA * DD + BB)
         DL = DL * (BB * DD - ONE)
            ERR = ABSC(DL)/ABSC(PQ)
         IF(NORM) THEN
            QQ = QQ + EN
            DS = QQ * DL
            SN = SN + DS
               ASN = ABSC(SN)
               ERR = MAX(ERR,ABSC(DS)/ASN)
C              ASM = MAX(ASM,ASN)
          ENDIF
         IF(ERR.GE.EPS .AND. N.LE.LIMIT) GO TO 10
C
         CF2E  = (PQ + DL) * PM / X
       IF(.NOT.NORM) RETURN
         F20   = ONE / SN
C        ERR = ACCUR* ASM / ASN
      RETURN
      END
