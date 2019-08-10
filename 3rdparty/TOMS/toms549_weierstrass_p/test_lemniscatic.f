      PROGRAM TEST_LEMNISCATIC

      COMPLEX Z, PLEM, PLEM1, CR, CI, P, P1
      REAL W2, LC, ALPHA, FP, FP1, EPS, EPS1, A
C
C     ******************************************************************
C
C     THE PROGRAMS FOR WEIERSTRASS: P-FUNCTION AND ITS DERIVATIVE IN THE
C     LEMNISCATIC CASE ARE TESTED FOR THE ARGUMENT VALUES
C     Z = 0.25, 0.5, 0.5 + 0.5*I, 0.25 + 0.25*I
C     AND 10.5 - 7.5*I
C     (SEE M. ABRAMOWITZ AND I. A. STEGUN& HANDBOOK OF MATHEMATICAL
C     FUNCTIONS, 18.14  FOR CORRESPONDING FUNCTION VALUES)
C     THE DRIVER PROGRAM PRINTS THE MAXIMAL RELATIVE ERRORS OF P AND P:
C
C     ******************************************************************
C
C
C     EVALUATION OF CONSTANTS
C     =======================
C
      LC = 2.6220575542921198104648396
C
C     LC IS THE LEMNISCATE CONSTANT
C
      W2 = 2.0*LC/SQRT(2.0)
      ALPHA = 1.0 + SQRT(2.0)
      CR = (1.0,0.0)
      CI = (0.0,1.0)
      FP = 1.0/(W2*W2)
      FP1 = FP/W2
C
C     * * * * * * * * * *
C
      Z = (0.25,0.0)
      P = PLEM(Z)
      P1 = PLEM1(Z)
      EPS = CABS(P*FP-0.5*ALPHA*CR)/CABS(P)
      EPS1 = CABS(P1*FP1+ALPHA*CR)/CABS(P1)
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C
      Z = (0.5,0.0)
      P = PLEM(Z)
      P1 = PLEM1(Z)
      A = CABS(P*FP-0.5*CR)/CABS(P)
      IF (A.GT.EPS) EPS = A
      A = CABS(P1*FP1)
      IF (A.GT.EPS1) EPS1 = A
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C
      Z = (0.5,0.5)
      P = PLEM(Z)
      P1 = PLEM1(Z)
      A = CABS(P*FP)
      IF (A.GT.EPS) EPS = A
      A = CABS(P1*FP1)
      IF (A.GT.EPS1) EPS1 = A
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C
      Z = (0.25,0.25)
      P = PLEM(Z)
      P1 = PLEM1(Z)
      A = CABS(P*FP+0.5*CI)/CABS(P)
      IF (A.GT.EPS) EPS = A
      A = CABS(P1*FP1-(CR+CI)/SQRT(2.0))/CABS(P1)
      IF (A.GT.EPS1) EPS1 = A
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C
      Z = (10.5,-7.5)
C
C     THIS Z IS FOR TESTING THE REDUCTION TO FUNDAMENTAL PARALLELOGRAM.
C     BECAUSE OF THIS REDUCTION, THE ERROR MIGHT BE AN ORDER OF
C     MAGNITUDE GREATER THAN THE BOUND GIVEN IN THE TEXT
C
      P = PLEM(Z)
      P1 = PLEM1(Z)
      A = CABS(P*FP)
      IF (A.GT.EPS) EPS = A
      A = CABS(P1*FP1)
      IF (A.GT.EPS1) EPS1 = A
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C     * * * * * * * * * *
C
      WRITE (6,99999)
      WRITE (6,99998) EPS
      WRITE (6,99997) EPS1
      STOP
99999 FORMAT (51H1   MAXIMAL RELATIVE ERROR FOR WEIERSTRASS P-FUNCTI,
     * 2HON/49H    IN THE LEMNISCATIC CASE AT FIVE SAMPLE POINTS/4X,
     * 50(1H=))
99998 FORMAT (////25H    RELATIVE ERROR FOR P=, 5X, E10.3)
99997 FORMAT (//26H    RELATIVE ERROR FOR P1=, 4X, E10.3)
      END PROGRAM
