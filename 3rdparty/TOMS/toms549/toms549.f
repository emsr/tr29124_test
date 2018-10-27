      COMPLEX Z, PLEM, PLEM1, CR, CI, P, P1                             MAN   10
      REAL W2, LC, ALPHA, FP, FP1, EPS, EPS1, A                         MAN   20
C                                                                       MAN   30
C     ******************************************************************MAN   40
C                                                                       MAN   50
C     THE PROGRAMS FOR WEIERSTRASS: P-FUNCTION AND ITS DERIVATIVE IN THEMAN   60
C     LEMNISCATIC CASE ARE TESTED FOR THE ARGUMENT VALUES               MAN   70
C     Z = 0.25, 0.5, 0.5 + 0.5*I, 0.25 + 0.25*I                         MAN   80
C     AND 10.5 - 7.5*I                                                  MAN   90
C     (SEE M. ABRAMOWITZ AND I. A. STEGUN& HANDBOOK OF MATHEMATICAL     MAN  100
C     FUNCTIONS, 18.14  FOR CORRESPONDING FUNCTION VALUES)              MAN  110
C     THE DRIVER PROGRAM PRINTS THE MAXIMAL RELATIVE ERRORS OF P AND P: MAN  120
C                                                                       MAN  130
C     ******************************************************************MAN  140
C                                                                       MAN  150
C                                                                       MAN  160
C     EVALUATION OF CONSTANTS                                           MAN  170
C     =======================                                           MAN  180
C                                                                       MAN  190
      LC = 2.6220575542921198104648396                                  MAN  200
C                                                                       MAN  210
C     LC IS THE LEMNISCATE CONSTANT                                     MAN  220
C                                                                       MAN  230
      W2 = 2.0*LC/SQRT(2.0)                                             MAN  240
      ALPHA = 1.0 + SQRT(2.0)                                           MAN  250
      CR = (1.0,0.0)                                                    MAN  260
      CI = (0.0,1.0)                                                    MAN  270
      FP = 1.0/(W2*W2)                                                  MAN  280
      FP1 = FP/W2                                                       MAN  290
C                                                                       MAN  300
C     * * * * * * * * * *                                               MAN  310
C                                                                       MAN  320
      Z = (0.25,0.0)                                                    MAN  330
      P = PLEM(Z)                                                       MAN  340
      P1 = PLEM1(Z)                                                     MAN  350
      EPS = CABS(P*FP-0.5*ALPHA*CR)/CABS(P)                             MAN  360
      EPS1 = CABS(P1*FP1+ALPHA*CR)/CABS(P1)                             MAN  370
C                                                                       MAN  380
C     * * * * * * * * * *                                               MAN  390
C                                                                       MAN  400
      Z = (0.5,0.0)                                                     MAN  410
      P = PLEM(Z)                                                       MAN  420
      P1 = PLEM1(Z)                                                     MAN  430
      A = CABS(P*FP-0.5*CR)/CABS(P)                                     MAN  440
      IF (A.GT.EPS) EPS = A                                             MAN  450
      A = CABS(P1*FP1)                                                  MAN  460
      IF (A.GT.EPS1) EPS1 = A                                           MAN  470
C                                                                       MAN  480
C     * * * * * * * * * *                                               MAN  490
C                                                                       MAN  500
      Z = (0.5,0.5)                                                     MAN  510
      P = PLEM(Z)                                                       MAN  520
      P1 = PLEM1(Z)                                                     MAN  530
      A = CABS(P*FP)                                                    MAN  540
      IF (A.GT.EPS) EPS = A                                             MAN  550
      A = CABS(P1*FP1)                                                  MAN  560
      IF (A.GT.EPS1) EPS1 = A                                           MAN  570
C                                                                       MAN  580
C     * * * * * * * * * *                                               MAN  590
C                                                                       MAN  600
      Z = (0.25,0.25)                                                   MAN  610
      P = PLEM(Z)                                                       MAN  620
      P1 = PLEM1(Z)                                                     MAN  630
      A = CABS(P*FP+0.5*CI)/CABS(P)                                     MAN  640
      IF (A.GT.EPS) EPS = A                                             MAN  650
      A = CABS(P1*FP1-(CR+CI)/SQRT(2.0))/CABS(P1)                       MAN  660
      IF (A.GT.EPS1) EPS1 = A                                           MAN  670
C                                                                       MAN  680
C     * * * * * * * * * *                                               MAN  690
C                                                                       MAN  700
      Z = (10.5,-7.5)                                                   MAN  710
C                                                                       MAN  720
C     THIS Z IS FOR TESTING THE REDUCTION TO FUNDAMENTAL PARALLELOGRAM. MAN  730
C     BECAUSE OF THIS REDUCTION, THE ERROR MIGHT BE AN ORDER OF         MAN  740
C     MAGNITUDE GREATER THAN THE BOUND GIVEN IN THE TEXT                MAN  750
C                                                                       MAN  760
      P = PLEM(Z)                                                       MAN  770
      P1 = PLEM1(Z)                                                     MAN  780
      A = CABS(P*FP)                                                    MAN  790
      IF (A.GT.EPS) EPS = A                                             MAN  800
      A = CABS(P1*FP1)                                                  MAN  810
      IF (A.GT.EPS1) EPS1 = A                                           MAN  820
C                                                                       MAN  830
C     * * * * * * * * * *                                               MAN  840
C     * * * * * * * * * *                                               MAN  850
C                                                                       MAN  860
      WRITE (6,99999)                                                   MAN  870
      WRITE (6,99998) EPS                                               MAN  880
      WRITE (6,99997) EPS1                                              MAN  890
      STOP                                                              MAN  900
99999 FORMAT (51H1   MAXIMAL RELATIVE ERROR FOR WEIERSTRASS P-FUNCTI,   MAN  910
     * 2HON/49H    IN THE LEMNISCATIC CASE AT FIVE SAMPLE POINTS/4X,    MAN  920
     * 50(1H=))                                                         MAN  930
99998 FORMAT (////25H    RELATIVE ERROR FOR P=, 5X, E10.3)              MAN  940
99997 FORMAT (//26H    RELATIVE ERROR FOR P1=, 4X, E10.3)               MAN  950
      END                                                               MAN  960
      COMPLEX FUNCTION PEQ(Z)                                           PEQ   10
C
C     WEIERSTRASS: P-FUNCTION IN THE EQUIANHARMONIC CASE
C     FOR COMPLEX ARGUMENT WITH UNIT PERIOD PARALLELOGRAM
C
      COMPLEX Z, Z2, Z4, Z6
      REAL ZR, ZI
      INTEGER M, N
C
C     REDUCTION TO FUNDAMENTAL PARALLELOGRAM
C
      ZI = 1.1547005383792515E0*AIMAG(Z) + 0.5E0
      M = INT(ZI)
      IF (ZI.LT.0E0) M = M - 1
      ZR = REAL(Z) - 0.5E0*FLOAT(M) + 0.5E0
      N = INT(ZR)
      IF (ZR.LT.0E0) N = N - 1
      Z2 = Z - FLOAT(N) - (0.5E0,0.86602540378443865E0)*FLOAT(M)
C
C     IF Z2=0 THEN Z COINCIDES WITH A LATTICE POINT.
C     SINCE P HAS POLES AT THE LATTICE POINTS,
C     A DIVISION ERROR WILL OCCUR
C
      Z2 = Z2*Z2
      Z4 = Z2*Z2
      Z6 = Z4*Z2
      PEQ = 1E0/Z2 + 6E0*Z4*(5E0+Z6)/(1E0-Z6)**2 + Z4*
     * (((((-2.6427662E-10*Z6+1.610954818E-8)*Z6+7.38610752879E-6)*
     * Z6+4.3991444671178E-4)*Z6+7.477288220490697E-2)*
     * Z6-6.8484153287299201E-1)/(((((6.2252191E-10*Z6+2.553314573E-7)*
     * Z6-2.619832920421E-5)*Z6-5.6444801847646E-4)*
     * Z6+4.565553484820106E-2)*Z6+1E0)
      RETURN
      END
      COMPLEX FUNCTION PEQ1(Z)                                          PEQ   10
C
C     FIRST DERIVATIVE OF WEIERSTRASS: P-FUNCTION IN THE
C     EQUIANHARMONIC CASE FOR COMPLEX ARGUMENT
C     WITH UNIT PERIOD PARALLELOGRAM
C
      COMPLEX Z, Z3, Z6
      REAL ZR, ZI
      INTEGER M, N
C
C     REDUCTION TO FUNDAMENTAL PARALLELOGRAM
C
      ZI = 1.1547005383792515E0*AIMAG(Z) + 0.5E0
      M = INT(ZI)
      IF (ZI.LT.0E0) M = M - 1
      ZR = REAL(Z) - 0.5E0*FLOAT(M) + 0.5E0
      N = INT(ZR)
      IF (ZR.LT.0E0) N = N - 1
      Z3 = Z - FLOAT(N) - (0.5E0,0.86602540378443865E0)*FLOAT(M)
C
C     IF Z3=0 THEN Z COINCIDES WITH A LATTICE POINT.
C     SINCE P: HAS POLES AT THE LATTICE POINTS,
C     A DIVISION ERROR WILL OCCUR
C
      Z3 = Z3*Z3*Z3
      Z6 = Z3*Z3
      PEQ1 = (((14E0*Z6+294E0)*Z6+126E0)*Z6-2E0)/(Z3*(1E0-Z6)**3) +
     * Z3*((((((-2.95539175E-9*Z6-2.6764693031E-7)*Z6+2.402192743346E-5)
     * *Z6+1.9656661451391E-4)*Z6+1.760135529461036E-2)*
     * Z6+8.1026243498822636E-1)*Z6-2.73936613149196804E0)/
     * ((((((4.6397763E-10*Z6+5.413482233E-8)*Z6-1.56293298374E-6)*
     * Z6-1.0393701076352E-4)*Z6+9.5553182532237E-4)*
     * Z6+9.131106969640212E-2)*Z6+1E0)
      RETURN
      END
      COMPLEX FUNCTION PLEM(Z)                                          PLE   10
C
C     WEIERSTRASS: P-FUNCTION IN THE LEMNISCATIC CASE
C     FOR COMPLEX ARGUMENT WITH UNIT PERIOD PARALLELOGRAM
C
      COMPLEX Z, Z2, Z4, Z6
      REAL ZR, ZI
      INTEGER M, N
C
C     REDUCTION TO FUNDAMENTAL PARALLELOGRAM
C
      ZR = REAL(Z) + 0.5E0
      ZI = AIMAG(Z) + 0.5E0
      M = INT(ZR)
      N = INT(ZI)
      IF (ZR.LT.0E0) M = M - 1
      IF (ZI.LT.0E0) N = N - 1
      Z2 = Z - FLOAT(M) - (0E0,1E0)*FLOAT(N)
C
C     IF Z2=0 THEN Z COINCIDES WITH A LATTICE POINT.
C     SINCE P HAS POLES AT THE LATTICE POINTS,
C     A DIVISION ERROR WILL OCCUR
C
      Z2 = Z2*Z2
      Z4 = Z2*Z2
      Z6 = Z4*Z2
      PLEM = 1E0/Z2 + 4E0*Z2*(3E0+Z4)/(1E0-Z4)**2 +
     * Z2*((((((((-7.233108E-11*Z4+1.714197273E-8)*Z4-2.5369036492E-7)*
     * Z4-7.98710206868E-6)*Z4+6.4850606909737E-4)*Z4+7.39624629362938E-
     * 3)*Z4+2.012382768497244E-2)*Z4+7.1177297543136598E-1)*
     * Z4-2.54636399353830738E0)/((((((((5.1161516E-10*Z4+6.61289408E-9)
     * *Z4+4.4618987048E-7)*Z4-8.42694918892E-6)*Z4+4.42886829095E-6)*
     * Z4-4.22629935217101E-3)*Z4+2.577496871700433E-2)*
     * Z4+4.2359940482277074E-1)*Z4+1E0)
      RETURN
      END
      COMPLEX FUNCTION PLEM1(Z)                                         PLE   10
C
C     FIRST DERIVATIVE OF WEIERSTRASS: P-FUNCTION IN THE
C     LEMNISCATIC CASE FOR COMPLEX ARGUMENT
C     WITH UNIT PERIOD PARALLELOGRAM
C
      COMPLEX Z, Z1, Z3, Z4
      REAL ZR, ZI
      INTEGER M, N
C
C     REDUCTION TO FUNDAMENTAL PARALLELOGRAM
C
      ZR = REAL(Z) + 0.5E0
      ZI = AIMAG(Z) + 0.5E0
      M = INT(ZR)
      N = INT(ZI)
      IF (ZR.LT.0E0) M = M - 1
      IF (ZI.LT.0E0) N = N - 1
      Z1 = Z - FLOAT(M) - (0E0,1E0)*FLOAT(N)
C
C     IF Z1=0 THEN Z COINCIDES WITH A LATTICE POINT.
C     SINCE P: HAS POLES AT THE LATTICE POINTS,
C     A DIVISION ERROR WILL OCCUR
C
      Z3 = Z1*Z1*Z1
      Z4 = Z3*Z1
      PLEM1 = (((1E1*Z4+9E1)*Z4+3E1)*Z4-2E0)/(Z1*(1E0-Z4))**3 +
     * Z1*((((((((((-3.9046302E-9*Z4-1.001487137E-8)*Z4+5.9573043092E-7)
     * *Z4-2.482518130524E-5)*Z4+1.4557266595395E-4)*
     * Z4+4.56633655643206E-3)*Z4+6.224782572111135E-2)*
     * Z4+1.038527937794269E-2)*Z4+1.19804620802637942E0)*
     * Z4+6.42791439683811718E0)*Z4-5.09272798707661477E0)/
     * ((((((((((4.726888E-11*Z4-3.0667983E-9)*Z4+1.0087596089E-7)*
     * Z4-8.060683451E-8)*Z4+1.184299251664E-5)*Z4-2.3096723361547E-4)*
     * Z4-2.90730903142055E-3)*Z4+1.338392411135511E-2)*
     * Z4+2.3098639320021426E-1)*Z4+8.4719880964554148E-1)*Z4+1E0)
      RETURN
      END
      COMPLEX Z, PEQ, PEQ1, CR, CI, P, P1                               MAN   10
      REAL W2, SQ3, FP, FP1, EPS, EPS1, A, C                            MAN   20
C                                                                       MAN   30
C     ******************************************************************MAN   40
C                                                                       MAN   50
C     THE PROGRAMS FOR WEIERSTRASS: P-FUNCTION AND ITS DERIVATIVE IN THEMAN   60
C     EQUIANHARMONIC CASE ARE TESTED FOR THE ARGUMENT VALUES            MAN   70
C     Z = 0.25, 0.5, 0.5 + 0.5/SQRT(3)*I, 0.25 + 0.25/SQRT(3)*I         MAN   80
C     AND 2.5 - 47/6*SQRT(3)*I                                          MAN   90
C     (SEE M. ABRAMOWITZ AND I. A. STEGUN& HANDBOOK OF MATHEMATICAL     MAN  100
C     FUNCTIONS, 18.13. FOR CORRESPONDING FUNCTION VALUES)              MAN  110
C     THE DRIVER PROGRAM PRINTS THE MAXIMAL RELATIVE ERRORS OF P AND P: MAN  120
C                                                                       MAN  130
C     ******************************************************************MAN  140
C                                                                       MAN  150
C                                                                       MAN  160
C     EVALUATION OF CONSTANTS                                           MAN  170
C                                                                       MAN  180
C     =======================                                           MAN  190
C                                                                       MAN  200
      SQ3 = SQRT(3.0)                                                   MAN  210
      CR = (1.0,0.0)                                                    MAN  220
      CI = (0.0,1.0)                                                    MAN  230
      C = 4.0**(-1.0/3.0)                                               MAN  240
      W2 = 1.52995403705719287491319417231                              MAN  250
      W2 = 2.0*W2                                                       MAN  260
      FP = 1.0/(W2*W2)                                                  MAN  270
      FP1 = FP/W2                                                       MAN  280
C                                                                       MAN  290
C     * * * * * * * * * *                                               MAN  300
C                                                                       MAN  310
      Z = (0.25,0.0)                                                    MAN  320
      P = PEQ(Z)                                                        MAN  330
      P1 = PEQ1(Z)                                                      MAN  340
      EPS = CABS(P*FP-C*(1.0+SQ3)*CR)/CABS(P)                           MAN  350
      EPS1 = CABS(P1*FP1+3.0**0.75*SQRT(2.0+SQ3)*CR)/CABS(P1)           MAN  360
C                                                                       MAN  370
C     * * * * * * * * * *                                               MAN  380
C                                                                       MAN  390
      Z = (0.5,0.0)                                                     MAN  400
      P = PEQ(Z)                                                        MAN  410
      P1 = PEQ1(Z)                                                      MAN  420
      A = CABS(P*FP-C*CR)/CABS(P)                                       MAN  430
      IF (A.GT.EPS) EPS = A                                             MAN  440
      A = CABS(P1*FP1)                                                  MAN  450
      IF (A.GT.EPS1) EPS1 = A                                           MAN  460
C                                                                       MAN  470
C     * * * * * * * * * *                                               MAN  480
C                                                                       MAN  490
      Z = 0.5*CR + 0.5/SQ3*CI                                           MAN  500
      P = PEQ(Z)                                                        MAN  510
      P1 = PEQ1(Z)                                                      MAN  520
      A = CABS(P*FP)                                                    MAN  530
      IF (A.GT.EPS) EPS = A                                             MAN  540
      A = CABS(P1*FP1-CI)/CABS(P1)                                      MAN  550
      IF (A.GT.EPS1) EPS1 = A                                           MAN  560
C                                                                       MAN  570
C     * * * * * * * * * *                                               MAN  580
C                                                                       MAN  590
      Z = Z*0.5                                                         MAN  600
      P = PEQ(Z)                                                        MAN  610
      P1 = PEQ1(Z)                                                      MAN  620
      A = CABS(P*FP+2.0**(1.0/3.0)*0.5*(-CR+SQ3*CI))/CABS(P)            MAN  630
      IF (A.GT.EPS) EPS = A                                             MAN  640
      A = CABS(P1*FP1-(0.0,3.0))/CABS(P1)                               MAN  650
      IF (A.GT.EPS1) EPS1 = A                                           MAN  660
C                                                                       MAN  670
C     * * * * * * * * * *                                               MAN  680
C                                                                       MAN  690
      Z = 2.5*CR - 47.0*SQ3/6.0*CI                                      MAN  700
C                                                                       MAN  710
C     THIS Z IS FOR TESTING THE REDUCTION TO FUNDAMENTAL PARALLELOGRAM. MAN  720
C     BECAUSE OF THIS REDUCTION, THE ERROR MIGHT BE AN ORDER OF         MAN  730
C     MAGNITUDE GREATER THAN THE BOUND GIVEN IN THE TEXT                MAN  740
C                                                                       MAN  750
      P = PEQ(Z)                                                        MAN  760
      P1 = PEQ1(Z)                                                      MAN  770
      A = CABS(P*FP)                                                    MAN  780
      IF (A.GT.EPS) EPS = A                                             MAN  790
      A = CABS(P1*FP1-CI)/CABS(P1)                                      MAN  800
      IF (A.GT.EPS1) EPS1 = A                                           MAN  810
C                                                                       MAN  820
C     * * * * * * * * * *                                               MAN  830
C     * * * * * * * * * *                                               MAN  840
C                                                                       MAN  850
      WRITE (6,99999)                                                   MAN  860
      WRITE (6,99998) EPS                                               MAN  870
      WRITE (6,99997) EPS1                                              MAN  880
      STOP                                                              MAN  890
99999 FORMAT (51H1   MAXIMAL RELATIVE ERROR FOR WEIERSTRASS P-FUNCTI,   MAN  900
     * 2HON/52H    IN THE EQUIANHARMONIC CASE AT FIVE SAMPLE POINTS/4X, MAN  910
     * 50(1H=))                                                         MAN  920
99998 FORMAT (////25H    RELATIVE ERROR FOR P=, 5X, E10.3)              MAN  930
99997 FORMAT (//26H    RELATIVE ERROR FOR P1=, 4X, E10.3)               MAN  940
      END                                                               MAN  950
