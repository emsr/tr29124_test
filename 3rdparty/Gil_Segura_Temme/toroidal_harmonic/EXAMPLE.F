
       PROGRAM EXAMPLE
       INTEGER NMAX,MODE,M,K,NEWN1,NEWN2,NEWN3,NSUP
       DOUBLE PRECISION EPS,PI,A,DL,DC,DCHA0,ALFA0,ALPHA,BETA,
     * DCHA,PHI,DFACTOR,QM10,PM10,PM1,DSUME,PRE,DDS,DFIELD
       DOUBLE PRECISION PL(0:2001),QL(0:2001),PL1(0:2001),
     * QL1(0:2001),PL2(0:2001),QL2(0:2001)       
       PARAMETER(EPS=1.D-8)     
       PI=3.14159265358979323D0
       NMAX=300
       MODE=0
       OPEN(UNIT=10,FILE='EXAMPLE.RES',STATUS='NEW')
       WRITE(*,*)'INTRODUCE GEOMETRIC PARAMETERS' 
       WRITE(*,*)'OF THE TORUS (A,L)'
1      READ(*,*)A,DL
       IF ((A.LT.0.D0).OR.(DL.LT.0.D0)) THEN
         WRITE(*,*)'IMPROPER ARGUMENTS FOR A OR L.
     *              MUST BE GREATER THAN ONE'
         GOTO 1
       END IF
       IF (DL.LT.A) THEN
         WRITE(*,*)'IMPROPER ARGUMENTS FOR A OR L.
     *               L MUST BE GREATER THAN A'
         GOTO 1
       END IF
       DC=DSQRT(DL*DL-A*A)
       DCHA0=DL/A
       ALFA0=DLOG(DCHA0+DSQRT(DCHA0*DCHA0-1.D0)) 
       WRITE(*,*)'..BEWARE: ALPHA MUST BE LOWER THAN', ALFA0
       WRITE(*,*)' '     
       WRITE(*,*)'INTRODUCE TOROIDAL COORDINATES (ALPHA, BETA)'
2      READ(*,*)ALPHA,BETA
       IF ((BETA.GT.PI).OR.(BETA.LE.-PI)) THEN
          WRITE(*,*)'IMPROPER ARGUMENT FOR BETA'
          GOTO 2
       END IF
       DCHA=DCOSH(ALPHA)
       IF (DCHA.GE.DCHA0) THEN
          WRITE(*,*)'I TOLD YOU THAT ALPHA MUST BE BELOW',ALFA0
          WRITE(*,*)' '
          WRITE(*,*)'..TRY AGAIN'
          GOTO 2
       END IF                  
       WRITE(*,*)'INTRODUCE M AND PHI'
       WRITE(*,*)'THE POTENTIAL IS GIVEN BY COS(M*PHI)' 
       READ(*,*)M,PHI
       DFACTOR=DCOS(M*PHI)/PI*DSQRT(2.D0*(DCHA-DCOS(BETA)))   
       CALL DTORH1(DCHA0,0,NMAX,PL,QL,NEWN1)
       QM10=QL(0)
       CALL DTORH1(DCHA0,M,NMAX,PL1,QL1,NEWN2)
       PM10=PL1(0)
       CALL DTORH1(DCHA,M,NMAX,PL2,QL2,NEWN3)
       PM1=PL2(0)
       NSUP=MIN(NEWN1,NEWN2,NEWN3) 
       DSUME=0.D0
       PRE=1.D0
       K=0
3      IF ((K.LE.NSUP).AND.(PRE.GT.EPS)) THEN
          K=K+1
          DSUME=DSUME+PL2(K)*QL(K)/PL1(K)*DCOS(K*BETA)
          PRE=ABS(DSUME-DDS)
          DDS=DSUME
          GOTO 3
       ENDIF
       DFIELD=DFACTOR*(PM1*QM10/PM10+2.D0*DSUME)
       WRITE(10,*)'NUMBER OF TERMS',K
       WRITE(10,*)'ELECTROSTATIC FIELD CREATED BY'
       WRITE(10,*)'A TOROIDAL CONDUCTOR OF GEOMETRIC'
       WRITE(10,*)'PARAMETERS (A,L)',A,DL
       WRITE(10,*)'AND POTENTIAL',DCOS(M*PHI)
       WRITE(10,*)'*******************************'
       WRITE(10,*) DFIELD
       END
   

