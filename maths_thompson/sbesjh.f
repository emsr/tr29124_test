      SUBROUTINE SBESJH(X,LMAX,XJ,XJP,XH1,XH1P,IFAIL)
C ***                                                       I.J.Thompson
C ***                                                       31 May 1985.
C ***  COMPLEX SPHERICAL BESSEL FUNCTIONS from l=0 to l=LMAX
C ***    for X in the UPPER HALF PLANE ( Im(X) > -3)
C ***
C ***    XJ(l)   = j/l(x)          regular solution: XJ(0)=sin(x)/x
C ***    XJP(l)  = d/dx j/l(x)
C ***    XH1(l)  = h(1)/l(x)       irregular Hankel function:
C ***    XH1P(l) = d/dx h(1)/l(x)            XH1(0) = j0(x) + i. y0(x)
C ***                                               =(sin(x)-i.cos(x))/x
C ***                                               = -i.exp(i.x)/x
C ***  Using complex CF1, and trigonometric forms for l=0 solutions.
C ***
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      PARAMETER (LIMIT=20000)
      DIMENSION XJ(0:LMAX),XJP(0:LMAX),XH1(0:LMAX),XH1P(0:LMAX)
      REAL*8 ZERO,ONE,ACCUR,TM30,ABSC
      DATA ZERO,ONE/ 0.0D0,1.0D0 /, ACCUR /1.0D-12/, TM30 / 1D-30 /,
     #     CI / (0D0,1D0) /
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))
      IFAIL= -1
      IF(ABSC(X).LT.ACCUR .OR. IMAG(X).LT.-3.0) GO TO 5
      XI = ONE/X
      W  = XI + XI
      PL = LMAX*XI
      F = PL + XI
      B  = F + F + XI
      D  = ZERO
      C  = F
      DO 1 L=1,LIMIT
      D  = B - D
      C  = B - ONE/C
         IF(ABSC(D).LT.TM30) D = TM30
         IF(ABSC(C).LT.TM30) C = TM30
      D = ONE / D
      DEL= D * C
      F = F * DEL
      B = B + W
    1 IF(ABSC(DEL-ONE).LT.ACCUR) GO TO 2
        IFAIL = -2
        GO TO 5
C
    2 XJ(LMAX)   = TM30
      XJP(LMAX)  = F * XJ(LMAX)
C
C *** Downward recursion to l=0 (N.B.  Coulomb Functions)
C
      DO 3 L = LMAX-1,0,-1
      XJ(L) = PL*XJ(L+1) + XJP(L+1)
      XJP(L)= PL*XJ(L)   - XJ(L+1)
    3 PL = PL - XI
C *** Calculate the l=0 Bessel Functions
      XJ0  = XI * SIN(X)
      XH1(0) = EXP(CI*X) * XI * (-CI)
      XH1P(0)= XH1(0) * (CI - XI)
C
C *** Rescale XJ, XJP,  converting to spherical Bessels.
C *** Recur   XH1,XH1P             AS spherical Bessels.
C
        W = ONE/XJ(0)
         PL = XI
      DO 4  L = 0,LMAX
      XJ(L)  =  XJ0*(W*XJ(L))
      XJP(L) =  XJ0*(W*XJP(L)) - XI*XJ(L)
         IF(L.EQ.0) GO TO 4
      XH1(L) = (PL-XI) * XH1(L-1) - XH1P(L-1)
         PL = PL + XI
      XH1P(L)=- PL     * XH1(L)   + XH1(L-1)
   4  CONTINUE
      IFAIL = 0
      RETURN
    5      WRITE(6,10) IFAIL
   10      FORMAT( 'SBESJH : IFAIL = ',I4)
      RETURN
      END
