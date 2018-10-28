      SUBROUTINE SBESJ(X,LMAX,XJ,IFAIL)
C ***  REGULAR SPHERICAL BESSEL FUNCTIONS FROM L=0 TO L=LMAX
C ***  SPECIAL CASE OF COULOMB FUNCTIONS (ETA=0) SEE BARNETT CPC 1980
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XJ(1)
      DATA ZERO,ONE/ 0.0D0,1.0D0 /, ACCUR /1.0D-12/
      IFAIL= 0
      IF(DABS(X).LT.ACCUR) GO TO 5
      XI = ONE/X
      W  = XI + XI
      F  = ONE
      FP = DFLOAT(LMAX + 1)*XI
      B  = FP + FP + XI
      D  = ONE/B
      DEL= -D
C            DEL = D
      FP = FP + DEL
      END = B + 20000.0*W
    1 B  = B + W
      D  = ONE/(B - D )
C            D = ONE/(B + D)
      DEL= DEL*(B*D - ONE)
      FP = FP + DEL
      IF(D.LT.ZERO) F = -F
      IF(B.GT.END) GO TO 5
      IF(DABS(DEL).GE.DABS(FP)*ACCUR) GO TO 1
      FP = FP*F
      IF(LMAX.EQ.0) GO TO 3
      XJ(LMAX+1) = F
      XP2        = FP
C *** DOWNWARD RECURSION TO L=0 (N.B.  COULOMB FUNCTIONS)
      PL = DFLOAT(LMAX)*XI
      L  = LMAX
      DO 2 LP = 1,LMAX
      XJ(L) = PL*XJ(L+1) + XP2
      FP    = PL*XJ(L)   - XJ(L+1)
      XP2 = FP
      PL = PL - XI
    2 L  =  L - 1
      F  = XJ(1)
C *** SOLVE FOR THE L=0 COULOMB FUNCTIONS
    3 W  = XI/DSQRT(FP*FP + F*F)
      XJ(1) = W*F
      IF(LMAX.EQ.0) RETURN
      DO 4  L = 1,LMAX
    4 XJ(L+1) =  W*XJ(L+1)
      RETURN
    5 IFAIL   = 1
           WRITE(6,10) X
   10      FORMAT( 10H WITH X = ,1PD15.5,24H TRY ASYMPTOTIC SOLUTION)
      RETURN
      END
