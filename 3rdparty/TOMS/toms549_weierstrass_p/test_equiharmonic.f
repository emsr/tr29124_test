      PROGRAM TEST_EQUIANHARMONIC

      COMPLEX Z, PEQ, PEQ1, CR, CI, P, P1
      REAL W2, SQ3, FP, FP1, EPS, EPS1, A, C
C
C     ******************************************************************
C
C     THE PROGRAMS FOR WEIERSTRASS: P-FUNCTION AND ITS DERIVATIVE IN THE
C     EQUIANHARMONIC CASE ARE TESTED FOR THE ARGUMENT VALUES
C     Z = 0.25, 0.5, 0.5 + 0.5/SQRT(3)*I, 0.25 + 0.25/SQRT(3)*I
C     AND 2.5 - 47/6*SQRT(3)*I
C     (SEE M. ABRAMOWITZ AND I. A. STEGUN& HANDBOOK OF MATHEMATICAL
C     FUNCTIONS, 18.13. FOR CORRESPONDING FUNCTION VALUES)
C     THE DRIVER PROGRAM PRINTS THE MAXIMAL RELATIVE ERRORS OF P AND P:
C
C     ******************************************************************
C
C
C     EVALUATION OF CONSTANTS
C
C     =======================
C
      SQ3 = SQRT(3.0)
      CR = (1.0,0.0)
      CI = (0.0,1.0)
      C = 4.0**(-1.0/3.0)
      W2 = 1.52995403705719287491319417231
      W2 = 2.0*W2
      FP = 1.0/(W2*W2)
      FP1 = FP/W2
C
C     * * * * * * * * * *
C
      Z = (0.25,0.0)
      P = PEQ(Z)
      P1 = PEQ1(Z)
      EPS = CABS(P*FP-C*(1.0+SQ3)*CR)/CABS(P)
      EPS1 = CABS(P1*FP1+3.0**0.75*SQRT(2.0+SQ3)*CR)/CABS(P1)
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C
      Z = (0.5,0.0)
      P = PEQ(Z)
      P1 = PEQ1(Z)
      A = CABS(P*FP-C*CR)/CABS(P)
      IF (A.GT.EPS) EPS = A
      A = CABS(P1*FP1)
      IF (A.GT.EPS1) EPS1 = A
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C
      Z = 0.5*CR + 0.5/SQ3*CI
      P = PEQ(Z)
      P1 = PEQ1(Z)
      A = CABS(P*FP)
      IF (A.GT.EPS) EPS = A
      A = CABS(P1*FP1-CI)/CABS(P1)
      IF (A.GT.EPS1) EPS1 = A
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C
      Z = Z*0.5
      P = PEQ(Z)
      P1 = PEQ1(Z)
      A = CABS(P*FP+2.0**(1.0/3.0)*0.5*(-CR+SQ3*CI))/CABS(P)
      IF (A.GT.EPS) EPS = A
      A = CABS(P1*FP1-(0.0,3.0))/CABS(P1)
      IF (A.GT.EPS1) EPS1 = A
      WRITE(6,*)'Z=', Z, 'P=', P, 'P1=', P1
C
C     * * * * * * * * * *
C
      Z = 2.5*CR - 47.0*SQ3/6.0*CI
C
C     THIS Z IS FOR TESTING THE REDUCTION TO FUNDAMENTAL PARALLELOGRAM.
C     BECAUSE OF THIS REDUCTION, THE ERROR MIGHT BE AN ORDER OF
C     MAGNITUDE GREATER THAN THE BOUND GIVEN IN THE TEXT
C
      P = PEQ(Z)
      P1 = PEQ1(Z)
      A = CABS(P*FP)
      IF (A.GT.EPS) EPS = A
      A = CABS(P1*FP1-CI)/CABS(P1)
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
     * 2HON/52H    IN THE EQUIANHARMONIC CASE AT FIVE SAMPLE POINTS/4X,
     * 50(1H=))
99998 FORMAT (////25H    RELATIVE ERROR FOR P=, 5X, E10.3)
99997 FORMAT (//26H    RELATIVE ERROR FOR P1=, 4X, E10.3)
      END PROGRAM
