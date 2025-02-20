      SUBROUTINE DAIRY(DX,AI,AIP,BI,BIP)
C  FOR DOUBLE PRECISION ARGUMENTS, THIS ROUTINE CALCULATES THE AIRY
C     FUNCTION AI(X) AND ITS DERIVATIVE AIP(X).  IT ALSO FINDS
C     THE OTHER REAL LINEARLY INDEPENDENT SOLUTION BI(X) AND
C     ITS DERIVATIVE BIP(X).
C     THE DEFINITIONS AND NORMALIZATIONS ARE AS IN NBS HANDBOOK
C     OF MATHEMATICAL FUNCTIONS,P.446
C     THE METHODS USED ARE POWER SERIES EXPANSION FOR SMALL X
C     AND GAUSSIAN INTEGRATION FOR LARGE X
      DIMENSION X(16),W(16),XSQ(16)
      DOUBLE PRECISION DX,AI,AIP,BI,BIP
      DOUBLE PRECISION    XS ,XCUBE,AISUM,AIPSUM
      DOUBLE PRECISION DF,DFP,DG,DGP
      DOUBLE PRECISION FJM2,FJM1,FJ,FJP1,FJP2,FACTOR
      DOUBLE PRECISION C1,C2,ROOT3
      DOUBLE PRECISION DZETA,DARG,DROOTX
      DOUBLE PRECISION ROOT4X,S,CO,RATIO,EFAC,ZETASQ
      DOUBLE PRECISION SUMR,SUMI,SUMRP,SUMIP,TERMR,TERMI
      DOUBLE PRECISION DZERO,DA,DB,DEN,ONE
      DOUBLE PRECISION X,W,XSQ
      DOUBLE PRECISION RSQ,TEMP,   RTPI,RTPI2
      DOUBLE PRECISION TERMA,TERMB
      LOGICAL NEEDBI
      DATA DZERO,ONE /0.0D0,1.0D0/
      DATA ROOT3/1.732050807568877D0/
      DATA C1,C2 /.355028053887817D0, .258819403792807D0/
      DATA RTPI /.2820947917738781D0/
      DATA RTPI2/.5641895835477562D0/
      DATA RSQ,DEN,TERMI/3*0.0D0/
C *** INSERTED ARB 6.80
C  POSITIONS AND WEIGHTS FOR 10-TERM SUM FOR AIRY FUNCTIONS
       DATA W( 1) /  3.1542515762964787D-14/
       DATA W( 2) /  6.6394210819584921D-11/
       DATA W( 3) /  1.7583889061345669D-08/
       DATA W( 4) /  1.3712392370435815D-06/
       DATA W( 5) /  4.4350966639284350D-05/
       DATA W( 6) /  7.1555010917718255D-04/
       DATA W( 7) /  6.4889566103335381D-03/
       DATA W( 8) /  3.6440415875773282D-02/
       DATA W( 9) /  1.4399792418590999D-01/
       DATA W(10) /  8.1231141336261486D-01/
       DATA X( 1) /  1.4083081072180964D+01/
       DATA X( 2) /  1.0214885479197331D+01/
       DATA X( 3) /  7.4416018450450930D+00/
       DATA X( 4) /  5.3070943061781927D+00/
       DATA X( 5) /  3.6340135029132462D+00/
       DATA X( 6) /  2.3310652303052450D+00/
       DATA X( 7) /  1.3447970824609268D+00/
       DATA X( 8) /  6.4188858369567296D-01/
       DATA X( 9) /  2.0100345998121046D-01/
       DATA X(10) /  8.0594359172052833D-03/
      DATA XSQ( 1) /0.19833317248562170D 03/
      DATA XSQ( 2) /0.10434388535311650D 03/
      DATA XSQ( 3) /0.55377438020178170D 02/
      DATA XSQ( 4) /0.28165249974668990D 02/
      DATA XSQ( 5) /0.13206054139355800D 02/
      DATA XSQ( 6) /0.54338651079380440D 01/
      DATA XSQ( 7) /0.18084791929954200D 01/
      DATA XSQ( 8) /0.41202095387883690D 00/
      DATA XSQ( 9) /0.40402390924418070D-01/
      DATA XSQ(10) /0.64954507303538390D-04/
C  POSITIONS AND WEIGHTS FOR  4-TERM SUM FOR AIRY FUNCTIONS
       DATA W(11) /  4.7763903057577263D-05/
       DATA W(12) /  4.9914306432910959D-03/
       DATA W(13) /  8.6169846993840312D-02/
       DATA W(14) /  9.0879095845981102D-01/
       DATA X(11) /  3.9198329554455091D+00/
       DATA X(12) /  1.6915619004823504D+00/
       DATA X(13) /  5.0275532467263018D-01/
       DATA X(14) /  1.9247060562015692D-02/
      DATA XSQ(11) /0.15365090398596670D 02/
      DATA XSQ(12) /0.28613816631634610D 01/
      DATA XSQ(13) /0.25276291648668180D 00/
      DATA XSQ(14) /0.37044934027789980D-03/
C  POSITIONS AND WEIGHTS FOR  2-TERM SUM FOR AIRY FUNCTIONS
       DATA W(15) /  9.6807280595773604D-01/
       DATA W(16) /  3.1927194042263958D-02/
       DATA X(15) /  3.6800601866153044D-02/
       DATA X(16) /  1.0592469382112378D+00/
      DATA XSQ(15) /0.13542842977111070D-02/
      DATA XSQ(16) /0.11220040761098810D 01/
      IF(DX.LT.-5.0D0) GO TO 100
      NEEDBI=.FALSE.
      IF(DX.GT.3.7D0) GO TO 200
C     THIS ROUTE FOR SMALLX, USING POWER SERIES.
C     INITIALIZE
10    XS  = DX*DX
      XCUBE = XS *DX
      XS  = XS *0.5D0
      DF = C1
      DFP = C1*XS
      DG = C2*DX
      DGP = C2
      AISUM = DF - DG
      AIPSUM = DFP - DGP
      BI = DF + DG
      BIP = DFP + DGP
      FJM2=-2.0D0
20    FJM2=FJM2+3.0D0
      FJM1=FJM2+ONE
      FJ=FJM1+ONE
      FJP1=FJ+ONE
      FJP2=FJP1+ONE
      RATIO = XCUBE/FJ
      DF = DF*RATIO/FJM1
      DFP = DFP*RATIO/FJP2
      DG = DG*RATIO/FJP1
      DGP = DGP*RATIO/FJM2
      BI = BI + (DF+DG)
      BIP = BIP + (DFP+DGP)
      IF(NEEDBI) GO TO 80
      AISUM = AISUM + (DF-DG)
      AIPSUM = AIPSUM + (DFP-DGP)
C     CONVERGENCE TEST
80    IF(DABS(DF).GT.1.0D-16) GO TO 20
C     CONVERGENCE. COMPUTE FUNCTIONS
99    BI = ROOT3*BI
      BIP = ROOT3*BIP
C  THIS RETURNS IF X IS BETWEEN 3.7 AND 8.0, SINCE IN SUCH CASES MORE
C  ACCURATE VALUES OF AI AND AIP HAVE ALREADY BEEN FOUND BY GAUSSIAN
C  INTEGRATION
      IF(NEEDBI)RETURN
      AI = AISUM
      AIP = AIPSUM
      RETURN
C  GAUSSIAN INTEGRATION FOR LARGE NEGATIVE X
100   DROOTX = DSQRT(-DX)
      ROOT4X = DSQRT(DROOTX)
      DZETA = -.6666666666666667D0*DX*DROOTX
      DARG = DZETA - .7853981633974483D0
      SUMR = DZERO
      SUMI = DZERO
      SUMRP = DZERO
      SUMIP = DZERO
C  TEST TO SEE HOW MANY TERMS ARE NEEDED IN GAUSSIAN INTEGRATION
      IF(DX.LT.(-200.D0)) GO TO 140
      IF(DX.LT.(-15.D0)) GO TO 130
C  THIS CASE FOR DX BETWEEN -5.0 AND -15.0
      LIMLO=1
      LIMHI=10
      GO TO 149
C  THIS CASE FOR DX BETWEEN -15.0 AND -200.
130   LIMLO=11
      LIMHI=14
      GO TO 149
C  THIS CASE FOR DX.LT.-200.
140   LIMLO=15
      LIMHI=16
149   ZETASQ=DZETA**2
      DO 150 K=LIMLO,LIMHI
      TERMR=W(K)/((ZETASQ+XSQ(K))**2)
      SUMR = SUMR + TERMR
      TERMR=TERMR*X(K)
      SUMI=SUMI+TERMR
      TERMR=TERMR*X(K)
      SUMRP=SUMRP+TERMR
150   SUMIP=SUMIP+TERMR*X(K)
      SUMR=(SUMR*ZETASQ+SUMRP)*ZETASQ
      TEMP=SUMI*ZETASQ
      SUMI=(TEMP+SUMIP)*DZETA
      SUMRP=SUMRP*DZETA
      SUMIP=SUMIP-TEMP
C  FORM AIRY FUNCTIONS
196   S = DSIN(DARG)
      CO = DCOS(DARG)
      RATIO = RTPI2/ROOT4X
      AI = RATIO*(CO*SUMR + S*SUMI)
      BI = RATIO*(CO*SUMI - S*SUMR)
      SUMRP=SUMRP+SUMRP
      RATIO = -.25D0/DX
      FACTOR = -RTPI2*ROOT4X
      AIP = RATIO*AI - DROOTX*BI + FACTOR*(CO*SUMRP+S*SUMIP)
      BIP = RATIO*BI + DROOTX*AI + FACTOR*(CO*SUMIP-S*SUMRP)
      RETURN
C   GAUSSIAN INTEGRATION FOR LARGE POSITIVE X
200   DROOTX = DSQRT(DX)
      DZETA = .6666666666666667D0*DX*DROOTX
      EFAC = DEXP(-DZETA)
      ROOT4X = DSQRT(DROOTX)
      AI = DZERO
      BI = DZERO
      AIP = DZERO
      BIP = DZERO
      IF(DX.LT.8.0D0) NEEDBI=.TRUE.
C  TEST TO SEE HOW MANY TERMS ARE NEEDED IN GAUSSIAN INTEGRATION
      IF(DX.GT.15.0D0) GO TO 230
C  THIS CASE FOR DX BETWEEN 3.7 AND 15.
      LIMLO=1
      LIMHI=10
      GO TO 249
C  THIS CASE FOR DX GREATER THAN 15.
230   LIMLO=11
      LIMHI=14
249   DO 250 K=LIMLO,LIMHI
      DA=DZETA+X(K)
      TERMA = W(K)/DA
      AI = AI + TERMA
      AIP=AIP+TERMA*X(K)/DA
      IF(NEEDBI) GO TO 250
      DB=DZETA-X(K)
      TERMB = W(K)/DB
      BI = BI + TERMB
      BIP=BIP+TERMB*X(K)/DB
250   CONTINUE
C  FORM FUNCTIONS
      FACTOR=RTPI*DZETA/ROOT4X
      RATIO = 0.25D0/DX
      AI=AI*EFAC*FACTOR
      AIP=-(DROOTX+RATIO)*AI+RTPI*ROOT4X*EFAC*AIP
C  THIS IS SATISFIED ONLY FOR X BETWEEN 3.7 AND 8.0  IN THESE CASES
C  THE BI AND BIP ABOUT TO BE COMPUTED ARE NOT SUFFICIENTLY ACCURATE.
C  THUS RETURN TO POWER SERIES FOR BI AND BIP.
      IF(NEEDBI) GO TO 10
      FACTOR=FACTOR+FACTOR
      BI=BI*FACTOR/EFAC
      BIP=(DROOTX-RATIO)*BI-RTPI2*ROOT4X*BIP/EFAC
      RETURN
      END


