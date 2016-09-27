 
        SUBROUTINE DTORH3(Z,MDIM,NDIM,MMAX,NMAX,PL,QL,NEWM,NEWN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    VERSION 2.0 OF DTORH3  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  THIS ROUTINE EVALUATES TOROIDAL HARMONICS FOR TWO DIFFERENT CHOICES    C 
C  OF RELATIVE ACCURACY (PARAMETER IPRE, SEE BELOW) AND FOR TWO           C
C  DIFFERENT NORMALIZATIONS OF THE FUNCTIONS (PARAMETER MODE, SEE BELOW). C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C  INPUT :                                                                C  
C    Z        ARGUMENT OF THE FUNCTIONS                                   C
C    MDIM     M-DIMENSION OF THE ARRAYS. MDIM MUST BE GREATER THAN MMAX   C
C                                                                         C
C    NDIM     N-DIMENSION OF THE ARRAYS. NDIM MUST BE GREATER THAN NMAX   C
C                                                                         C
C    MMAX     MAXIMUM ORDER  OF THE FUNCTIONS :                           C           
C              FUNCTIONS OF ALL ORDERS BELOW MMAX ARE CALCULATED.         C            
C    NMAX     MAXIMUM DEGREE  OF THE FUNCTIONS :                          C           
C              FUNCTIONS OF ALL DEGREES BELOW MIN(NMAX,NEWN) ARE          C
C              CALCULATED.                                                C
C                                                                         C
C  PARAMETERS:                                                            C
C   IPRE:                                                                 C
C      *IF IPRE=1, PRECISION=10**(-12)                                    C
C      *IF IPRE=2, PRECISION=10**(-8)                                     C
C   MODE:                                                                 C
C      TWO DIFFERENT NORMALIZATIONS FOR THE OUTPUT FUNCTIONS CAN BE       C
C      CHOSEN, DEPENDING ON THE VALUE OF MODE (0 OR 1).                   C
C      FOR MODE=0, THE OUTPUTS ARE STANDARD TOROIDAL HARMONICS.           C
C      FOR MODE=1, A DIFFERENT NORMALIZATION IS CHOSEN. SEE BELOW.        C 
C   EPS:                                                                  C
C      REQUIRED ACCURACY FOR THE CONTINUED FRACTION (MODIFIED LENTZ) AND  C
C      SERIES. SHOULD BE CHOSEN SMALLER THAN THE PRECISION CORRESPONDING  C
C      TO THE VALUE OF IPRE.                                              C
C   TINY:                                                                 C
C      SMALL PARAMETER TO PREVENT OVERFLOWS (CLOSE TO THE                 C   
C      UNDERFLOW LIMIT)                                                   C 
C                                                                         C
C  OUTPUT :                                                               C
C   *IF MODE IS EQUAL TO 0:                                               C
C                                                                         C                                            
C    PL(M,N)                                                              C
C             THESE VALUES ARE KEPT IN AN ARRAY                           C
C    QL(M,N)                                                              C
C             THESE VALUES ARE KEPT IN AN ARRAY                           C
C                                                                         C
C    NEWM     MAXIMUM  ORDER OF FUNCTIONS CALCULATED WHEN                 C
C             QL (MMAX,0)   IS LARGER THAN 1/TINY                         C
C              (OVERFLOW LIMIT = 1/TINY).                                 C
C    NEWN     MAXIMUM  DEGREE OF FUNCTIONS CALCULATED WHEN                C          
C             PL (M,NMAX)   IS LARGER THAN 1/TINY  FOR SOME               C
C             M=0,...,NEWM                                                C
C              (OVERFLOW LIMIT = 1/TINY).                                 C
C                                                                         C
C   *IF MODE IS EQUAL TO 1:                                               C
C                                                                         C
C      THE SET OF FUNCTIONS EVALUATED IS:                                 C
C                 PL(M,N)/GAMMA(M+1/2),QL(M,N)/GAMMA(M+1/2),              C
C      WHICH ARE RESPECTIVELY STORED IN THE ARRAYS PL(M,N),QL(M,N)        C
C      NEWM AND NEWN REFER TO THIS NEW SET OF FUNCTIONS                   C
C                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    NOTE1: FOR A PRECISION OF 10**(-12):                                 C
C             IF Z>SQRT(2) THEN:                                          C
C             *PRIMAL ALGORITHM:                                          C
C               IF Z>5 AND (Z/MMAX)>0.22 THE CODE USES                    C
C               A SERIES EXPANSION FOR PL(MMAX,0).WHEN Z<20 AND           C
C               (Z/MMAX)<0.22, A CONTINUED FRACTION IS APPLIED.           C
C               IF Z>20 AND (Z/MMAX)<0.22 AN ASYMPTOTIC EXPANSION IS      C
C               USED.                                                     C
C             IF 1<Z<SQRT(2) THEN:                                        C
C               LET XL=Z/SQRT(Z*Z-1) THEN:                                C 
C             * DUAL ALGORITHM:                                           C
C                IF XL>5 AND (XL/NMAX)>0.22 THE CODE USES                 C
C                A SERIES EXPANSION FOR PL(NMAX,0)(XL).WHEN XL<20 AND     C   
C                (XL/NMAX)<0.22 A CONTINUED FRACTION IS APPLIED FOR THE   C
C                EVALUATION OF QL(0,NMAX)(Z).                             C 
C                IF XL>20 AND (XL/NMAX)<0.22 AN ASYMPTOTIC                C
C                EXPANSION IS USED FOR PL(NMAX,0)(XL).                    C                             
C    NOTE2: FOR A PRECISION OF 10**(-8):                                  C
C             IF Z>SQRT(2) THEN:                                          C
C             * PRIMAL ALGORITHM:                                         C
C                IF Z>5 AND (Z/MMAX)>0.12 THE CODE USES                   C
C                A SERIES EXPANSION FOR PL(MMAX,0).WHEN Z<20 AND          C
C                (Z/MMAX)<0.12, A CONTINUED FRACTION IS APPLIED.          C
C                IF Z>20 AND (Z/MMAX)<0.12 AN ASYMPTOTIC EXPANSION        C
C                IS USED.                                                 C
C             IF 1<Z<SQRT(2) THEN:                                        C
C               LET XL=Z/SQRT(Z*Z-1) THEN:                                C 
C             * DUAL ALGORITHM:                                           C
C                IF XL>5 AND (XL/NMAX)>0.12 THE CODE USES                 C
C                A SERIES EXPANSION FOR PL(NMAX,0)(XL).WHEN XL<20 AND     C   
C                (XL/NMAX)<0.12 A CONTINUED FRACTION IS APPLIED FOR THE   C
C                EVALUATION OF QL(0,NMAX)(Z).                             C 
C                IF XL>20 AND (XL/NMAX)<0.12 AN ASYMPTOTIC                C
C                EXPANSION IS USED FOR PL(NMAX,0)(XL).                    C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   DECLARATION OF VARIABLES   C             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IMPLICIT REAL*8 (A-H,O-Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   IN ORDER TO ENLARGE THE RANGES FOR C
C   THE ARGUMENTS AND IF THE           C
C   FORTRAN COMPILER SUPPORT QUADRUPLE C
C   ARITHMETICS, UNCOMMENT THE         C
C   FOLLOWING LINE:                    C
CCCCC       REAL*16 QZ,QDC1,QQ                 
       DIMENSION PL(0:MDIM,0:NDIM),QL(0:MDIM,0:NDIM)
       DIMENSION PR(2)
       PARAMETER(PI=3.14159265358979323D0,EPS=1.D-14,TINY=1.D-290,
     f           MODE=0,IPRE=1)
       IF (Z.LE.1.D0) THEN
          WRITE(*,*)'IMPROPER ARGUMENT.Z MUST BE GREATER THAN 1'
          STOP
       END IF
       IF ((IPRE.NE.1).AND.(IPRE.NE.2)) THEN
         WRITE(6,*)'IPRE MUST BE 1 OR 2'
         STOP
       END IF
       OVER=1.D0/TINY
       OVERQ=OVER*1.D-30
       UNDER=TINY
       TINYSQ=DSQRT(TINY)
       PR(1)=.22D0
       PR(2)=.12D0
       DSQ2=DSQRT(2.D0)
       PISQ=DSQRT(PI)
       DPPI=DSQRT(2.D0)/PISQ
       QZ=Z
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   IN ORDER TO ENLARGE THE RANGES FOR C
C   THE ARGUMENTS AND IF THE           C
C   FORTRAN COMPILER SUPPORT QUADRUPLE C
C   ARITHMETICS, UNCOMMENT THE         C
C   FOLLOWING LINE:                    C
CCCCC       QDC1=(QZ-1.Q0)*(QZ+1.Q0)
C   AND COMMENT THE FOLLOWING LINE:     
       QDC1=(QZ-1.D0)*(QZ+1.D0)
CCCCCCCCCCC END THE COMMENT  CCCCCCCCCCC
       QQ=SQRT(QDC1)
       QARGU=QZ/QQ
       XL=QARGU
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  WE USE PRIMAL OR DUAL ALGORITHM DEPENDING ON  C 
CC  THE VALUE OF Z:                    CCCCCCCCCCCC
CC    IF Z > SQRT(2) ----> PRIMAL      C 
CC    IF Z < SQRT(2) ----> DUAL        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
       IF (Z.LT.DSQ2) THEN       
CCCCCCCCCCCCCCCCCCCCCC
CC  DUAL ALGORITHM   C
CCCCCCCCCCCCCCCCCCCCCC    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   WE EVALUATE P^{0}_{-1/2},P^{0}_{+1/2}             C
CC   USING SLATEC ROUTINES FOR ELLIPTIC FUNCTIONS      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           DZZ=DSQRT(Z*Z-1.D0)
           ARGU1=DSQRT((Z-1.D0)/(Z+1.D0))
           ARGU2=DSQRT(2.D0*DZZ/(Z+DZZ))
           PL(0,0)=2.D0/PI*DSQRT(2.D0/(Z+1.D0))*ELLIP1(ARGU1)
           PL(0,1)=2.D0/PI*DSQRT(Z+DZZ)*ELLIP2(ARGU2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   WE APPLY FORWARD RECURRENCE IN N FOR P'S   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           NP=1
           IF (MODE.EQ.0) THEN
              DO WHILE ((NP.LE.NMAX).AND.(ABS(PL(0,NP)).LT.OVER))
                 PL(0,NP+1)=(2.D0*NP*Z*PL(0,NP)-(NP-0.5D0)*PL(0,NP-1))
     f                   /(NP+0.5D0)
                 NP=NP+1
              ENDDO       
              IF ((NP-1).NE.NMAX) NMAX=NP-1
              NEWN=NMAX
           ELSE
              PL(0,0)=PL(0,0)/PISQ
              PL(0,1)=PL(0,1)/PISQ
              DO WHILE ((NP.LE.NMAX).AND.(ABS(PL(0,NP)).LT.OVER))
                 PL(0,NP+1)=(2.D0*NP*Z*PL(0,NP)-(NP-0.5D0)*PL(0,NP-1))
     f                   /(NP+0.5D0)
                 NP=NP+1
              ENDDO
              IF ((NP-1).NE.NMAX) NMAX=NP-1
              NEWN=NMAX
           ENDIF 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  WE CHOOSE CF, SERIES OR ASYMPTOTIC EXPANSION    C
CC  DEPENDING ON THE VALUES OF XL AND NMAX          C                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           ICAL=1
           IF ((XL/NMAX).GT.PR(IPRE)) ICAL=2
           IF (XL.LT.5.D0) THEN
             ICAL=1
           ELSE IF (XL.GT.20.D0) THEN
             IF (ICAL.EQ.1) ICAL=0
           END IF
           IF (ICAL.EQ.1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  WE EVALUATE THE CONTINUED FRACTION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              CALL FRAC(Z,0,NMAX,EPS,TINYSQ,FC)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   WE EVALUATE QL(0,NMAX+1) AND QL(0,NMAX) USING :         C
C   THE WRONSKIAN : W{PL(0,NMAX),QL(0,NMAX)} =1./(1.-Z**2)  C
C   THE KNOWN VALUES OF PL(0,NMAX+1) AND PL(0,NMAX)         C
C   THE VALUE OF FC=QL(0,NMAX+1)/QL(0,NMAX)                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              IF (MODE.EQ.0) THEN
               QL(0,NMAX)=1.D0/(PL(0,NMAX+1)-FC*PL(0,NMAX))/(NMAX+0.5D0)
              ELSE
               QL(0,NMAX)=1.D0/(PL(0,NMAX+1)-FC*PL(0,NMAX))/
     f           ((NMAX+0.5D0)*PI) 
              ENDIF         
              QL(0,NMAX+1)=QL(0,NMAX)*FC
           ELSEIF (ICAL.EQ.2) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCC        
CC   WE CALCULATE SERIES  C
CCCCCCCCCCCCCCCCCCCCCCCCCCC  
 40           FL=NMAX/2.D0              
              CC=ABS(FLOAT(INT(FL))-FL)
              IF (CC.LT.0.4D0) THEN
                 AR=1.
              ELSE
                 AR=-1.
              END IF
              GAMMAPR=GAMMAH(NMAX,OVER)*AR
              GAMMA=GAMMAPR*PISQ
              IF (ABS(GAMMA).LT.TINY) THEN
                 NMAX=NMAX-5
                 GOTO 40
              ENDIF
              DFAC=PI**1.5/DSQ2*(XL*XL-1.D0)**0.25
              CALL EXPAN(XL,0,IPRE,OVER,Z,NMAX,PL0)
              QL(0,NMAX)=(PL0/GAMMA)*DFAC  
              IF (MODE.EQ.0) THEN   
                 QL(0,NMAX+1)=(PL(0,NMAX+1)*QL(0,NMAX)
     F                      -1.D0/(NMAX+0.5D0))/
     F                      PL(0,NMAX)
              ELSE
                 QL(0,NMAX)=QL(0,NMAX)/PISQ
                 QL(0,NMAX+1)=(PL(0,NMAX+1)*QL(0,NMAX)-
     F                      1.D0/(NMAX+0.5D0)/PI)/
     F                      PL(0,NMAX)
              ENDIF
           ELSEIF (ICAL.EQ.0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  WE USE THE ASYMPTOTIC EXPANSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 50           FL=NMAX/2.D0
              CC=ABS(FLOAT(INT(FL))-FL)
              IF (CC.LT.0.4D0) THEN
                 AR=1.
              ELSE
                 AR=-1.
              END IF
              GAMMAPR=GAMMAH(NMAX,OVER)*AR
              GAMMA=GAMMAPR*PISQ
              IF (ABS(GAMMA).LT.TINY) THEN
                 NMAX=NMAX-5
                 GOTO 50
              ENDIF
              DFAC=PI**1.5/DSQ2*(XL*XL-1.D0)**0.25
              CALL ASEXPAN(XL,NMAX,0,GAMMAPR,PL0)
              QL(0,NMAX)=(PL0/GAMMA)*DFAC  
              IF (MODE.EQ.0) THEN   
                 QL(0,NMAX+1)=(PL(0,NMAX+1)*QL(0,NMAX)
     F                      -1.D0/(NMAX+0.5D0))/
     F                      PL(0,NMAX)
              ELSE
                 QL(0,NMAX)=QL(0,NMAX)/PISQ
                 QL(0,NMAX+1)=(PL(0,NMAX+1)*QL(0,NMAX)-
     F                      1.D0/(NMAX+0.5D0)/PI)/
     F                      PL(0,NMAX)
              ENDIF
          ENDIF 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CFS FOR PS AND CALCULATION OF  C
C  PL(1,NMAX+1),PL(1,NMAX)        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          CALL FRACPS(XL,1,NMAX+1,EPS,TINYSQ,FCP1) 
          PL(1,NMAX+1)=FCP1*PL(0,NMAX+1)
          CALL FRACPS(XL,1,NMAX,EPS,TINYSQ,FCP2) 
          PL(1,NMAX)=FCP2*PL(0,NMAX)
          IF (MODE.EQ.1) THEN
             PL(1,NMAX+1)=PL(1,NMAX+1)*2.D0 
             PL(1,NMAX)=PL(1,NMAX)*2.D0
          ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  EVALUATION OF QL(1,NMAX+1),QL(1,NMAX)  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          DIZZ=-1.D0/DZZ
          IF (MODE.EQ.0) THEN 
             QL(1,NMAX+1)=(DIZZ+PL(1,NMAX+1)*QL(0,NMAX+1))/PL(0,NMAX+1)
             QL(1,NMAX)=(DIZZ+PL(1,NMAX)*QL(0,NMAX))/PL(0,NMAX) 
          ELSE
             DA=2.D0/PI
             QL(1,NMAX+1)=(DIZZ*DA+PL(1,NMAX+1)*QL(0,NMAX+1))
     f       /PL(0,NMAX+1)
             QL(1,NMAX)=(DIZZ*DA+PL(1,NMAX)*QL(0,NMAX))/PL(0,NMAX)  
          ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  EVALUATION OF THE SET {QL(M,N)}  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   WE APPLY FORWARD RECURRENCE IN M FOR Q'S   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IF (MODE.EQ.0) THEN
             MP=1
             DO WHILE ((MP.LE.MMAX).AND.(ABS(QL(MP,NMAX)).LT.OVERQ))
                QL(MP+1,NMAX)=-2.D0*MP*QARGU*QL(MP,NMAX) 
     f          -(MP-NMAX-0.5D0)*(MP+NMAX-0.5D0)*QL(MP-1,NMAX) 
                MP=MP+1
             ENDDO
             MMAX=MP-2
             MP=1
             DO WHILE ((MP.LE.MMAX).AND.(ABS(QL(MP,NMAX+1)).LT.OVERQ))
                QL(MP+1,NMAX+1)=-2.D0*MP*QARGU*QL(MP,NMAX+1) 
     f          -(MP-NMAX-1.5D0)*(MP+NMAX+0.5D0)*QL(MP-1,NMAX+1) 
                MP=MP+1
             ENDDO        
          ELSE
             MP=1
             DO WHILE ((MP.LE.MMAX).AND.(ABS(QL(MP,NMAX)).LT.OVER))
                D1=MP+0.5D0
                D2=MP-0.5D0
                QL(MP+1,NMAX)=-2.D0*MP*QARGU*QL(MP,NMAX)/D1 
     f          -(MP-NMAX-0.5D0)*(MP+NMAX-0.5D0)*QL(MP-1,NMAX)/(D1*D2) 
                MP=MP+1
             ENDDO
             MMAX=MP-2             
             MP=1
             DO WHILE ((MP.LE.MMAX).AND.(ABS(QL(MP,NMAX+1)).LT.OVER))
                D1=MP+0.5D0
                D2=MP-0.5D0
                QL(MP+1,NMAX+1)=-2.D0*MP*QARGU*QL(MP,NMAX+1)/D1 
     f          -(MP-NMAX-1.5D0)*(MP+NMAX+0.5D0)*QL(MP-1,NMAX+1)/(D1*D2) 
                MP=MP+1
             ENDDO        
          ENDIF    
          IF ((MP-1).LT.MMAX) MMAX=MP-1
          NEWM=MMAX+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FINALLY, FOR EACH M=0,...,MMAX APPLYING RECURRENCE   C
C  BACKWARDS, OBTAIN THE SET {QL(M,N)}                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          DO L=0,MMAX+1      
             DO I=1,NMAX
                NP=NMAX+1-I
                N=NP-1
                QL(L,N)=((NP+NP)*Z*QL(L,NP)-(NP-L+0.5D0)
     f          *QL(L,NP+1))/(NP+L-0.5D0)
             ENDDO
          ENDDO         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  EVALUATION OF THE SET {PL(M,N)}  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
CC CALCULATION OF THE CF FOR PS C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
           FL=MMAX/2.D0
           CC=ABS(FLOAT(INT(FL))-FL)
           IF (CC.LT.0.4D0) THEN
             AR=1.
           ELSE
             AR=-1.
           END IF
           IF (MODE.EQ.0) THEN
             GAMMA=GAMMAH(MMAX,OVER)*AR*PISQ               
           ELSE 
             GAMMA=AR
           END IF
           DFACQS=-(GAMMA/QL(MMAX+1,0))*GAMMA/PI
           CALL FRACPS(XL,MMAX+1,0,EPS,TINYSQ,FCP) 
           DD=MMAX+0.5D0
           IF (MODE.NE.0) THEN
            FCP=FCP/DD
            DFACQS=DFACQS/DD
           END IF
           PL(MMAX,0)=DFACQS/QQ
     f        /(1.D0-FCP*QL(MMAX,0)/QL(MMAX+1,0))
           PL(MMAX+1,0)=FCP*PL(MMAX,0)          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   EVALUATION OF PL(MMAX,1)     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           QM0=QL(MMAX,0)
           DFAC3=(GAMMA/QM0)*GAMMA/PI/(0.5D0-MMAX)
           CALL FRAC(Z,MMAX,0,EPS,TINYSQ,FC)
           PL(MMAX,1)=PL(MMAX,0)*FC+DFAC3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   EVALUATION OF PL(MMAX+1,1)      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          QM0=QL(MMAX+1,0)
          IF (MODE.EQ.0) THEN 
             DFACC=(0.5D0+MMAX)
          ELSE
             DFACC=1.D0/(0.5D0+MMAX)
          END IF      
          DFAC3=-(GAMMA/QM0)*GAMMA*DFACC/PI
          CALL FRAC(Z,MMAX+1,0,EPS,TINYSQ,FC2)
          PL(MMAX+1,1)=PL(MMAX+1,0)*FC2+DFAC3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  WE APPLY BACKWARD RECURRENCE OVER M TO GET THE  C
CC  SET PL(M,0),PL(M,1)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF (MODE.EQ.0) THEN
            DO I=1,MMAX
              MP=MMAX+1-I
              M=MP-1
              PL(M,0)=-(PL(MP+1,0)+2.D0*MP*QARGU*PL(MP,0))
     f              /((0.5D0-MP)*(0.5D0-MP))
            ENDDO
            DO I=1,MMAX
              MP=MMAX+1-I
              M=MP-1
              PL(M,1)=(PL(MP+1,1)+2.D0*MP*QARGU*PL(MP,1))
     f              /((1.5D0-MP)*(0.5D0+MP))
            ENDDO
         ELSE
            DO I=1,MMAX
              MP=MMAX+1-I
              M=MP-1
              PL(M,0)=((MP+0.5D0)*PL(MP+1,0)+2.D0*MP*QARGU*PL(MP,0))
     f              /(0.5D0-MP)
            ENDDO
            DO I=1,MMAX
              MP=MMAX+1-I
              M=MP-1
              PL(M,1)=((MP+0.5D0)*PL(MP+1,1)+2.D0*MP*QARGU*PL(MP,1))*
     f              (MP-0.5D0)/((1.5D0-MP)*(0.5D0+MP))
            ENDDO
         END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  NOW, WE PERFORM THE EVALUATION OVER N FOR EACH M  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC      WE START EVALUATING THE SET OF PL(M,N)        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO L=0,MMAX+1
            M=L
            FL=M/2.D0
            CC=ABS(FLOAT(INT(FL))-FL)
            IF (CC.LT.0.4D0) THEN
               AR=1.
            ELSE
               AR=-1.
            END IF
            IF (MODE.EQ.0) THEN
               GAMMA=GAMMAH(M,OVER)*AR*PISQ
            ELSE 
               GAMMA=AR
            END IF
            NP=1        
            DO WHILE ((NP.LE.NMAX).AND.
     f       (ABS((NP-M+0.5D0)*PL(L,NP)).LT.OVER))
                 PL(L,NP+1)=(2.D0*NP*Z*PL(L,NP)
     f           -(NP+M-0.5D0)*PL(L,NP-1))/(NP-M+0.5D0)
                 NP=NP+1
            ENDDO
            NMAX=NP-1
            NEWN=NMAX
         ENDDO
       ELSE
CCCCCCCCCCCCCCCCCCCCCCCC
CC  PRIMAL ALGORITHM   C
CCCCCCCCCCCCCCCCCCCCCCCC   
         ARGU1=DSQRT(2.D0/(Z+1.D0))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   WE EVALUATE Q^{0}_{-1/2},Q^{1}_{-1/2}           C
CC   USING SLATEC ROUTINES FOR ELLIPTIC FUNCTIONS    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          QL(0,0)=ARGU1*ELLIP1(ARGU1)
          QL(1,0)=-1.D0/DSQRT(2.D0*(Z-1.D0))
     f          *ELLIP2(ARGU1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    WE APPLY FORWARD RECURRENCE IN M FOR Q'S     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          MP=1
          IF (MODE.EQ.0) THEN
             DO WHILE ((MP.LE.MMAX).AND.(ABS(QL(MP,0)).LT.OVER))
               QL(MP+1,0)=-2.D0*MP*QARGU*QL(MP,0) 
     f          -(MP-0.5D0)*(MP-0.5D0)*QL(MP-1,0) 
               MP=MP+1
             ENDDO       
             IF ((MP-1).LT.MMAX) MMAX=MP-1
             NEWM=MMAX
          ELSE
             QL(0,0)=QL(0,0)/PISQ
             QL(1,0)=QL(1,0)*2.D0/PISQ
             DO WHILE ((MP.LE.MMAX).AND.(ABS(QL(MP,0)).LT.OVER))
               D1=MP+0.5D0
               QL(MP+1,0)=-2.D0*MP*QARGU*QL(MP,0)/D1
     f         -(MP-0.5D0)*QL(MP-1,0)/D1
               MP=MP+1
             ENDDO
             IF ((MP-1).LT.MMAX) MMAX=MP-1
             NEWM=MMAX
          END IF
          FL=MMAX/2.D0
          CC=ABS(FLOAT(INT(FL))-FL)
          IF (CC.LT.0.4D0) THEN
             AR=1.
          ELSE
             AR=-1.
          END IF
          IF (MODE.EQ.0) THEN
             GAMMAPR=GAMMAH(MMAX,OVER)*AR
             GAMMA=GAMMAPR*PISQ  
             IF (ABS(GAMMA).LT.TINY) THEN
                WRITE(6,*)'MMAX IS TOO LARGE FOR MODE=0'
                WRITE(6,*)'BETTER TRY MODE=1'
                STOP
             END IF              
          ELSE 
             GAMMAPR=AR
             GAMMA=AR
          END IF
          DFACQS=-(GAMMA/QL(MMAX+1,0))*GAMMA/PI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  EVALUATION OF PL(MMAX,0),PL(MMAX+1,0)   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  WE CHOOSE SERIES, CF OR ASYMPTOTIC EXPANSION  C
CC  FOR PL(MMAX,0)                                C
CC  DEPENDING ON THE VALUES OF Z AND MMAX         C                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ICAL=1
          IF ((Z/MMAX).GT.PR(IPRE)) ICAL=2
          IF (Z.LT.5.D0) THEN
             ICAL=1
          ELSE IF (Z.GT.20.D0) THEN
             IF ((MODE.NE.2).AND.(ICAL.EQ.1)) ICAL=0
          END IF
          IF (ICAL.EQ.1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCC        
CC  WE CALCULATE THE CF   C
CCCCCCCCCCCCCCCCCCCCCCCCCCC  
             CALL FRACPS(XL,MMAX+1,0,EPS,TINYSQ,FCP) 
             DD=MMAX+0.5D0
             IF (MODE.NE.0) THEN
               FCP=FCP/DD
               DFACQS=DFACQS/DD
             END IF
             PL(MMAX,0)=DFACQS/QQ
     f         /(1.D0-FCP*QL(MMAX,0)/QL(MMAX+1,0))
             PL(MMAX+1,0)=FCP*PL(MMAX,0)
          ELSEIF (ICAL.EQ.2) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCC        
CC   WE CALCULATE SERIES  C
CCCCCCCCCCCCCCCCCCCCCCCCCCC 
             CALL EXPAN(Z,MODE,IPRE,OVER,XL,MMAX,PL0)
             PL(MMAX,0)=PL0
             DD=MMAX+0.5D0
             IF (MODE.NE.0) THEN
                FCP=FCP/DD
                DFACQS=DFACQS/DD
             END IF
             PL(MMAX+1,0)=(QL(MMAX+1,0)/QL(MMAX,0))*
     f         (PL(MMAX,0)- DFACQS/QQ)
          ELSEIF (ICAL.EQ.0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
CC   WE CALCULATE ASYMPTOTIC   C
CC   EXPANSION                 C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             CALL ASEXPAN(Z,MMAX,MODE,GAMMAPR,PL0)
             PL(MMAX,0)=PL0
             DD=MMAX+0.5D0
             IF (MODE.NE.0) THEN
               FCP=FCP/DD
               DFACQS=DFACQS/DD
             END IF
             PL(MMAX+1,0)=(QL(MMAX+1,0)/QL(MMAX,0))*
     f         (PL(MMAX,0)-DFACQS/QQ)
          END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   EVALUATION OF PL(MMAX,1)   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          QM0=QL(MMAX,0)
          DFAC3=(GAMMA/QM0)*GAMMA/PI/(0.5D0-MMAX)
          CALL FRAC(Z,MMAX,0,EPS,TINYSQ,FC)
          PL(MMAX,1)=PL(MMAX,0)*FC+DFAC3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   EVALUATION OF PL(MMAX+1,1)    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          QM0=QL(MMAX+1,0)
          IF (MODE.EQ.0) THEN 
             DFACC=(0.5D0+MMAX)
          ELSE
             DFACC=1.D0/(0.5D0+MMAX)
          END IF      
          DFAC3=-(GAMMA/QM0)*GAMMA*DFACC/PI
          CALL FRAC(Z,MMAX+1,0,EPS,TINYSQ,FC2)
          PL(MMAX+1,1)=PL(MMAX+1,0)*FC2+DFAC3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  WE APPLY BACKWARD RECURRENCE OVER M TO GET THE  C
CC  SET PL(M,0),PL(M,1)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF (MODE.EQ.0) THEN
            DO I=1,MMAX
              MP=MMAX+1-I
              M=MP-1
              PL(M,0)=-(PL(MP+1,0)+2.D0*MP*QARGU*PL(MP,0))
     f              /((0.5D0-MP)*(0.5D0-MP))
            ENDDO
            DO I=1,MMAX
               MP=MMAX+1-I
               M=MP-1
               PL(M,1)=(PL(MP+1,1)+2.D0*MP*QARGU*PL(MP,1))
     f              /((1.5D0-MP)*(0.5D0+MP))
            ENDDO
         ELSE
            DO I=1,MMAX
               MP=MMAX+1-I
               M=MP-1
               PL(M,0)=((MP+0.5D0)*PL(MP+1,0)+2.D0*MP*QARGU*PL(MP,0))
     f              /(0.5D0-MP)
            ENDDO
            DO I=1,MMAX
              MP=MMAX+1-I
              M=MP-1
              PL(M,1)=((MP+0.5D0)*PL(MP+1,1)+2.D0*MP*QARGU*PL(MP,1))*
     f              (MP-0.5D0)/((1.5D0-MP)*(0.5D0+MP))
            ENDDO
         END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  NOW, WE PERFORM THE EVALUATION OVER N FOR EACH M  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC      WE START EVALUATING THE SET OF PL(M,N)        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO L=0,MMAX
            M=L
            FL=M/2.D0
            CC=ABS(FLOAT(INT(FL))-FL)
            IF (CC.LT.0.4D0) THEN
              AR=1.
            ELSE
              AR=-1.
            END IF
            IF (MODE.EQ.0) THEN
              GAMMA=GAMMAH(M,OVER)*AR*PISQ
            ELSE 
              GAMMA=AR
            END IF
            NP=1
            DO WHILE ((NP.LE.NMAX).AND.
     f        (ABS((NP-M+0.5D0)*PL(L,NP)).LT.OVER))
                 PL(L,NP+1)=(2.D0*NP*Z*PL(L,NP)
     f           -(NP+M-0.5D0)*PL(L,NP-1))/(NP-M+0.5D0)
                 NP=NP+1
            ENDDO
            NMAX=NP-1
            NEWN=NMAX
         ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   WE EVALUATE THE SET OF QL(M,N)        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
         DO L=0,MMAX
            M=L
            FL=M/2.D0
            CC=ABS(FLOAT(INT(FL))-FL)
            IF (CC.LT.0.4D0) THEN
              AR=1.
            ELSE
              AR=-1.
            END IF
            IF (MODE.EQ.0) THEN
               GAMMA=GAMMAH(M,OVER)*AR*PISQ
            ELSE 
               GAMMA=AR
            END IF     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   WE EVALUATE THE C.F. FOR Q'S USING LENTZ-THOMPSON      C
CC   IF M=0,1. IN OTHER CASE WE CALCULATE QL(M,NMAX) AND    C
CC   QL(M,NMAX+1) FROM THE RECURRENCE.                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            IF (M.EQ.2) THEN
            IF (MODE.EQ.0) THEN
              MP=1
              DO WHILE (MP.LE.MMAX)
                 QL(MP+1,NMAX)=-2.D0*MP*QARGU*QL(MP,NMAX) 
     f            -(MP-NMAX-0.5D0)*(MP+NMAX-0.5D0)*QL(MP-1,NMAX) 
                 MP=MP+1
              ENDDO
              MP=1
              DO WHILE (MP.LE.MMAX)
                 QL(MP+1,NMAX+1)=-2.D0*MP*QARGU*QL(MP,NMAX+1) 
     f           -(MP-NMAX-1.5D0)*(MP+NMAX+0.5D0)*QL(MP-1,NMAX+1) 
                 MP=MP+1
              ENDDO
            ELSE
              MP=1
              DO WHILE (MP.LE.MMAX)
                D1=MP+0.5D0
                D2=MP-0.5D0
                QL(MP+1,NMAX)=-2.D0*MP*QARGU*QL(MP,NMAX)/D1
     f           -(MP-NMAX-0.5D0)*(MP+NMAX-0.5D0)*QL(MP-1,NMAX)/(D1*D2)
                MP=MP+1
              ENDDO
              MP=1
              DO WHILE (MP.LE.MMAX)
                D1=MP+0.5D0
                D2=MP-0.5D0
                QL(MP+1,NMAX+1)=-2.D0*MP*QARGU*QL(MP,NMAX+1)/D1 
     f           -(MP-NMAX-1.5D0)*(MP+NMAX+0.5D0)
     f           *QL(MP-1,NMAX+1)/(D1*D2)
                MP=MP+1
              ENDDO
            END IF
           END IF
           IF ((M.EQ.0).OR.(M.EQ.1)) THEN
           CALL FRAC(Z,M,NMAX,EPS,TINYSQ,FC)
           DFAC4=(FACTCO(NMAX,PL(L,NMAX+1),M)*GAMMA)
     f        *GAMMA/PI/(NMAX+M+0.5D0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
CC     EVALUATION OF QL(L,NMAX+1) AND QL(L,NMAX) USING      C
CC     THE WRONSKIAN W{PL(L,NMAX),QL(L,NMAX)},              C
CC     THE KNOWN VALUES OF PL(L,NMAX+1) AND PL(L,NMAX)      C
CC     THE VALUE OF H = QL(L,NMAX+1)/QL(L,NMAX)             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           QL(L,NMAX)=DFAC4/(1.D0-FC*PL(L,NMAX)/PL(L,NMAX+1))
           QL(L,NMAX+1)=QL(L,NMAX)*FC
          END IF       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   WE USE THE BACKWARD RECURRENCE RELATION FOR Q'S          C                                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          DO I=1,NMAX
             NP=NMAX+1-I
             N=NP-1
             QL(L,N)=((NP+NP)*Z*QL(L,NP)-(NP-M+0.5D0)
     f       *QL(L,NP+1))/(NP+M-0.5D0)
          ENDDO
        END DO
        ENDIF
        RETURN
        END

 
