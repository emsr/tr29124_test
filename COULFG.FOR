      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP, 
     *                  MODE1,KFN,IFAIL)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
C                                                                      C
C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
C  THIS VERSION WRITTEN UP       IN    CPC XX (1982) YYY-ZZZ           C
C                                                                      C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
C   THE DIRAC EQUATION ,ALSO SPHERICAL + CYLINDRICAL BESSEL EQUATIONS  C
C                                                                      C
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
C  STARTING ARRAY ELEMENT IS M1 = MAX0(  INT(XLMIN+ACCUR),0) + 1       C
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 C
C                                                                      C
C  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     C
C            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
C            = 3      F               CALL TO AT LEAST LENGTH (1)      C
C  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            C
C            = 1 SPHERICAL   BESSEL                                    C
C            = 2 CYLINDRICAL BESSEL                                    C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
C                                                                      C
C  PRECISION&  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
C   FOR MANTISSAS OF 56 + 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION    FC(1),GC(1),FCP(1),GCP(1)
      LOGICAL      ETANE0,XLTURN
      COMMON       /STEED/ PACCQ,NFP,NPQ,IEXP,M1
C***  COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
C***  COULFG HAS CALLS TO&  SQRT, ABS, AMOD ,  INT, SIGN, FLOAT, MIN1 
      DATA ZERO,ONE,TWO,TEN2,ABORT /0.0E0, 1.0E0, 2.0E0, 1.0E2, 2.0E4/
      DATA HALF,TM30 / 0.5E0, 1.0E-30 / 
      DATA RT2EPI /0.79788 45608 02865E0/
C *** THIS CONSTANT IS  DSQRT(TWO/PI)&  USE Q0 FOR IBM REAL*16& D0 FOR
C ***  REAL*8 + CDC DOUBLE P&  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.
C
                        ACCUR = R1MACH(4)
C ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
      MODE  = 1
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      ETA   = ETA1
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACC   = ACCUR 
      ACC4  = ACC*TEN2*TEN2
      ACCH  = SQRT(ACC)
C ***    TEST RANGE OF XX, EXIT IF.LE. SQRT(ACCUR) OR IF NEGATIVE
C
      IF(XX .LE. ACCH)                          GO TO 100
      X     = XX
      XLM   = XLMIN 
      IF(KFN .EQ. 2)  XLM = XLM - HALF
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 105
      E2MM1 = ETA*ETA + XLM*XLM + XLM
      XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM
      DELL  = XLMAX - XLMIN + ACC
      IF( ABS(MOD(DELL,ONE)) .GT. ACC) WRITE(6,2040)XLMAX,XLMIN,DELL
      LXTRA =   INT(DELL)
      XLL   = XLM +  FLOAT(LXTRA)
C ***       LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C ***       XLL  IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
C ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN
      M1  = MAX0(  INT(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
C
C ***    EVALUATE CF1  =  F   =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
C
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
    2 EK  = ETA / PK
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 =  PK + ONE
C ***   TEST ENSURES B1 .NE. ZERO FOR NEGATIVE ETA; FIXUP IS EXACT.
             IF( ABS(ETA*X + PK*PK1) .GT. ACC)  GO TO 3
             FCL  = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)
             PK   =  TWO + PK 
      GO TO 2
    3 D   =  ONE/((PK + PK1)*(XI + EK/PK1))
      DF  = -FCL*(ONE + EK*EK)*D
            IF(FCL .NE. ONE )  FCL = -ONE
            IF(D   .LT. ZERO)  FCL = -FCL
      F   =  F  + DF
C
C ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
C
      P     = ONE
    4 PK    = PK1
        PK1 = PK1 + ONE
        EK  = ETA / PK
        TK  = (PK + PK1)*(XI + EK/PK1)
        D   =  TK - D*(ONE + EK*EK)
              IF( ABS(D) .GT. ACCH)             GO TO 5
              WRITE (6,1000) D,DF,ACCH,PK,EK,ETA,X
              P = P  +   ONE
              IF( P .GT. TWO )                  GO TO 110
    5 D     = ONE/D 
              IF (D .LT. ZERO) FCL = -FCL
        DF  = DF*(D*TK - ONE) 
        F   = F  + DF
              IF(PK .GT. PX)                    GO TO 110
      IF( ABS(DF) .GE.  ABS(F)*ACC)             GO TO 4
                  NFP = PK - XLL - 1
      IF(LXTRA .EQ. 0)                          GO TO 7
C
C *** DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
C
      FCL = FCL*TM30
      FPL = FCL*F
      IF(MODE .EQ. 1) FCP(L1) = FPL
                      FC (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO 6  LP = 1,LXTRA
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL =  SQRT(ONE + EL*EL)
         SL    =  EL  + XL*XI 
         L     =  L1  - LP
         FCL1  = (FCL *SL + FPL)/RL
         FPL   =  FCL1*SL - FCL *RL
         FCL   =  FCL1
         FC(L) =  FCL
         IF(MODE .EQ. 1) FCP(L)  = FPL
         IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
    6 XL = XL - ONE 
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
C ***    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM 
C ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
C ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM 
C
    7 IF( XLTURN ) CALL JWKB(X,ETA,MAX(XLM,ZERO),FJWKB,GJWKB,IEXP)
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 9
          XLTURN = .FALSE.
      TA =  TWO*ABORT
      PK =  ZERO
      WI =  ETA + ETA
      P  =  ZERO
      Q  =  ONE - ETA*XI
      AR = -E2MM1
      AI =  ETA
      BR =  TWO*(X - ETA)
      BI =  TWO
      DR =  BR/(BR*BR + BI*BI)
      DI = -BI/(BR*BR + BI*BI)
      DP = -XI*(AR*DI + AI*DR)
      DQ =  XI*(AR*DR - AI*DI)
    8 P     = P  + DP
         Q  = Q  + DQ
         PK = PK + TWO
         AR = AR + PK
         AI = AI + WI
         BI = BI + TWO
         D  = AR*DR - AI*DI + BR
         DI = AI*DR + AR*DI + BI
         C  = ONE/(D*D + DI*DI)
         DR =  C*D
         DI = -C*DI 
         A  = BR*DR - BI*DI - ONE
         B  = BI*DR + BR*DI
         C  = DP*A  - DQ*B
         DQ = DP*B  + DQ*A
         DP = C
         IF(PK .GT. TA)                         GO TO 120
      IF( ABS(DP)+ ABS(DQ).GE.( ABS(P)+ ABS(Q))*ACC)   GO TO 8
                      NPQ   = PK/TWO
                      PACCQ = HALF*ACC/MIN( ABS(Q),ONE)
                      IF( ABS(P) .GT.  ABS(Q)) PACCQ = PACCQ* ABS(P)
C
C *** SOLVE FOR FCM = F AT LAMBDA = XLM,THEN FIND NORM FACTOR W=W/FCM 
C
      GAM = (F - P)/Q
            IF(Q .LE. ACC4* ABS(P))             GO TO 130
      W   = ONE/ SQRT((F - P)*GAM + Q)
            GO TO 10
C *** ARRIVE HERE IF G(XLM) .GT. 10**6 OR IEXP .GT. 250 + XLTURN = .TRUE.
    9 W   = FJWKB
      GAM = GJWKB*W 
      P   = F
      Q   = ONE
C
C *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C
   10                     ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .EQ. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .EQ. 2) BETA  =  SQRT(XI)*RT2EPI
      FCM  =  SIGN(W,FCL)*BETA
           FC(M1)  = FCM
                      IF(MODE .EQ. 3)           GO TO 11
           IF(.NOT. XLTURN)   GCL =  FCM*GAM
           IF(      XLTURN)   GCL =  GJWKB*BETA
           IF( KFN .NE. 0 )   GCL = -GCL
           GC(M1)  = GCL
           GPL =  GCL*(P - Q/GAM) - ALPHA*GCL
                      IF(MODE .EQ. 2)           GO TO 11
           GCP(M1) = GPL
           FCP(M1) = FCM*(F - ALPHA)
   11 IF(LXTRA .EQ. 0 ) RETURN
C *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
C *** RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
         W    = BETA*W/ ABS(FCL)
         MAXL = L1 - 1
      DO 12 L = M1,MAXL
                      IF(MODE .EQ. 3)           GO TO 12
                      XL = XL + ONE
         IF(ETANE0)   EL = ETA/XL
         IF(ETANE0)   RL = GC(L+1)
                      SL = EL + XL*XI
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1 
         GCL      = GCL1
         GC(L+1)  = GCL1
                      IF(MODE .EQ. 2)           GO TO 12
         GCP(L+1) = GPL
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
   12 FC(L+1)     = W* FC(L+1)
      RETURN
 1000 FORMAT(/' CF1 ACCURACY LOSS& D,DF,ACCH,K,ETA/K,ETA,X = ',1P7E9.2/)
C
C ***    ERROR MESSAGES
C
  100 IFAIL = -1
      WRITE(6,2000) XX,ACCH
 2000 FORMAT(' FOR XX = ',1PE12.3,' TRY SMALL-X  SOLUTIONS',
     *' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',E12.3/)
      RETURN
  105 IFAIL = -2
      WRITE (6,2005) XLMAX,XLMIN,XLM
 2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES&XLMAX,XLMIN,XLM = ',
     *1P3E15.6/)
      RETURN
  110 IFAIL =  1
      WRITE (6,2010) ABORT,F ,DF,PK,PX,ACC
 2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     *' F,DF,PK,PX,ACCUR =  ',1P5E12.3//)
      RETURN
  120 IFAIL =  2
      WRITE (6,2020) ABORT,P,Q,DP,DQ,ACC
 2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     *' P,Q,DP,DQ,ACCUR =  ',1P4E17.7,E12.3//)
      RETURN
  130 IFAIL =  3
      WRITE (6,2030) P,Q,ACC,DELL,LXTRA,M1
 2030 FORMAT(' FINAL Q.LE. ABS(P)*ACC*10**4 , P,Q,ACC = ',1P3E12.3,4X,
     *' DELL,LXTRA,M1 = ',E12.3,2I5 /)
      RETURN
 2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P3E20.10/)
      END 
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)
      REAL      XX,ETA1,XL,FJWKB,GJWKB, ZERO
C *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0
C *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C *** CALLS AMAX1,SQRT,ALOG,EXP,ATAN2,FLOAT,INT        BARNETT FEB 1981
      DATA   ZERO,HALF,ONE,SIX,TEN/ 0.0E0, 0.5E0, 1.0E0, 6.0E0, 10.0E0 /
      DATA  DZERO, RL35, ALOGE  /0.0E0, 35.0E0, 0.43429 45 E0 /
      X     = XX
      ETA   = ETA1
      GH2   = X*(ETA + ETA - X)
      XLL1  = MAX(XL*XL + XL,DZERO)
      IF(GH2 + XLL1 .LE. ZERO) RETURN
       HLL  = XLL1 + SIX/RL35 
       HL   = SQRT(HLL)
       SL   = ETA/HL + HL/X
       RL2  = ONE + ETA*ETA/HLL
       GH   = SQRT(GH2 + HLL)/X
       PHI  = X*GH - HALF*( HL*LOG((GH + SL)**2/RL2) - LOG(GH) )
          IF(ETA .NE. ZERO) PHI = PHI - ETA*ATAN2(X*GH,X - ETA)
      PHI10 = -PHI*ALOGE
      IEXP  =  INT(PHI10)
      IF(IEXP .GT. 250) GJWKB = TEN**(PHI10 - FLOAT(IEXP))
      IF(IEXP .LE. 250) GJWKB = EXP(-PHI)
      IF(IEXP .LE. 250) IEXP  = 0
      FJWKB = HALF/(GH*GJWKB) 
      RETURN
      END 
