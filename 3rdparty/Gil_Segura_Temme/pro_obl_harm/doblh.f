         SUBROUTINE DOBLH(X,M,NMAX,MODE,RL,TL,NUEVO)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CALCULATION OF OBLATE SPHEROIDAL HARMONICS                             C
C                                                                         C
C  INPUT :                                                                C
C                                                                         C  
C    X        ARGUMENT OF THE FUNCTIONS                                   C
C                                                                         C
C    M        DEGREE OF THE SPHEROIDAL HARMONICS                          C
C                                                                         C
C    NMAX     MAXIMUM ORDER OF THE FUNCTIONS :                            C
C             WE SHALL GET  FUNCTIONS OF ALL THE ORDERS BELOW             C
C             MIN(NMAX,NUEVO). NUEVO IS DEFINED BELOW .                   C
C                                                                         C
C    MODE     MODE OF CALCULATION (SEE OUTPUT)                            C
C                                                                         C
C  OUTPUT :                                                               C
C   *IF MODE IS EQUAL TO 0:                                               C
C                                                                         C
C    RL(L+1)                                                              C
C             WE SHALL KEEP THESE VALUES IN AN ARRAY                      C
C    TL(L+1)                                                              C
C             WE SHALL KEEP THESE VALUES IN AN ARRAY                      C
C    NUEVO    MAXIMUM ORDER OF  FUNCTIONS CALCULATED WHEN                 C
C             RL (NMAX+1) IS LARGER THAN 1/TINY. TINY IS                  C
C             DEFINED BELOW.                                              C
C                                                                         C
C                                                                         C
C                                                                         C
C   *IF MODE IS EQUAL TO 1:                                               C
C                                                                         C
C    RL(L+1)/(2M-1)!!                                                     C
C             WE SHALL KEEP THESE VALUES IN AN ARRAY                      C
C    TL(L+1)/(2M)!!                                                       C
C             WE SHALL KEEP THESE VALUES IN AN ARRAY                      C
C    NUEVO    MAXIMUM ORDER OF  FUNCTIONS CALCULATED WHEN                 C
C             RL (NMAX+1,X)/(2M-1)!!  IS LARGER THAN 1/TINY.              C
C             TINY IS DEFINED BELOW.                                      C
C                                                                         C
C    *IF MODE IS EQUAL TO 2:                                              C
C                                                                         C
C    RL(L+1)*TINY/(2M-1)!!/(X*X+1)**(M/2)                                 C
C             WE SHALL KEEP THESE VALUES IN AN ARRAY                      C
C    TL(L+1)*(X*X+1)**(M/2)/(2M)!!/TINY                                   C
C             WE SHALL KEEP THESE VALUES IN AN ARRAY                      C
C    NUEVO    MAXIMUM ORDER OF  FUNCTIONS CALCULATED WHEN                 C
C             RL (NMAX+1,X)/(2M-1)!!  IS LARGER THAN 1/TINY.              C
C             TINY IS DEFINED BELOW.                                      C
C                                                                         C
C                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C   DECLARATION OF VARIABLES                                              C
C                                                                         C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION RL(0:NMAX+1),TL(0:NMAX+1)
       PARAMETER(PI=3.14159265358979323D0,EPS=1.D-16,TINY=1.D-290)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   EPS: REQUIRED ACCURACY FOR THE CONTINUED FRACTION (LENTZ-THOMPSON)   C
C   TINY: SMALL PARAMETER CLOSE TO THE UNDERFLOW LIMIT                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       OVER=1.D0/TINY
       TINYSQ=DSQRT(TINY)
       IF (X.LE.0.D0) THEN
         WRITE(6,*)'IMPROPER ARGUMENT. X MUST BE GREATER THAN 0'
         STOP
       END IF
       FL=M/2.
       IF (FLOAT(INT(FL)).EQ.FL) THEN
         AR=1.
       ELSE
         AR=-1.
       END IF      
       DZM=DSQRT(X*X+1.D0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C   WE USE THE CODE IF NMAX IS GREATER THAN OR EQUAL TO 2                 C
C                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       NMAXP=NMAX
       IF(NMAX.LT.2) NMAX=2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C   WE STORE THE VALUES OF RL(0) AND RL(1)                                C
C                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       RL0=1.D0
       IF (MODE.EQ.0) THEN
         IF (M.GT.0) THEN  
           J=1 
           DO I=1,M
             RL0=RL0*DZM*(2.*M-J)
             J=J+2
           END DO  
         END IF
         RL(0)=RL0
         RL(1)=X*(2.*M+1.D0)*RL(0)
       ELSE IF (MODE.EQ.1) THEN
         IF (M*DLOG(DZM).GT.(DLOG(OVER))) THEN
           WRITE(6,*)'BETTER TRY MODE=2'
           STOP
         END IF  
         IF (M.GT.0) THEN
           DO I=1,M
             RL0=RL0*DZM
           END DO 
         END IF
         RL(0)=RL0
         RL(1)=X*(2.*M+1.D0)*RL(0)
       ELSE
         RL(0)=RL0*TINY
         RL(1)=X*(2.*M+1.D0)*RL(0)                             
       END IF
       NP=1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C   WE USE THE RECURRENCE RELATIONS                                       C
C                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO WHILE ((NP.LE.NMAX).AND.(ABS((NP+1.D0)*RL(NP)).LT.OVER))
         RL(NP+1)=((2.D0*(NP+M)+1.D0)*X*RL(NP)+(NP+M+M)*
     &            RL(NP-1))/(NP+1.D0)
         NP=NP+1
       ENDDO
       NMAX=NP-1
       IF (MODE.EQ.0) THEN         
         FACTOR=FACTCO(NMAX,RL(NMAX+1),M,OVER)        
         DO WHILE ((FACTOR.EQ.0.D0).AND.(NMAX.GT.0)) 
           NMAX=NMAX-3
           FACTOR=FACTCO(NMAX,RL(NMAX+1),M,OVER)     
         END DO   
         IF (NMAX.LE.0.) THEN
           WRITE(6,*)'TRY ANOTHER M'
           STOP
         END IF
       ELSE
         FACTOR=FACTCONEW(NMAX,RL(NMAX+1),M)
       END IF         
       NUEVO=NMAX      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C   WE EVALUATE THE CONTINUED FRACTION USING                              C
C   LENTZ-THOMPSON ALGORITHM                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       N=NMAX
       MM=0
       B=(1.D0+(N+2.D0)/(M+M+N+1.D0) )*X
       A=1.D0
       FC=TINYSQ
       C0=FC
       D0=0.D0
10     D0=B+A*D0
       IF (DABS(D0).LT.TINYSQ) D0=TINYSQ
       C0=B+A/C0
       IF (DABS(C0).LT.TINYSQ) C0=TINYSQ
       D0=1.D0/D0
       DELTA=C0*D0
       FC=FC*DELTA
       MM=MM+1
       A=(N+MM+1.D0)/DFLOAT(N+M+M+MM) 
       B=(1.D0+(N+MM+2.D0)/(M+M+N+MM+1.D0))*X   
       IF (ABS(DELTA-1.D0).GT.EPS) GOTO 10
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C     WE EVALUATE TL(NMAX+1) AND TL(NMAX) USING :                         C
C     THE WRONSKIAN : W{RL(NMAX),TL(NMAX)} =1./(1.-X**2)                  C
C     THE KNOWN VALUES OF RL(NMAX+1) AND RL(NMAX)                         C
C     THE VALUE OF H = TL(NMAX+1)/TL(NMAX)                                C
C                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       TL(NMAX)=AR*FACTOR/(1.D0+FC*RL(NMAX)/RL(NMAX+1))
       TL(NMAX+1)=TL(NMAX)*FC
       DO I=1,NMAX
        NP=NMAX+1-I
        N=NP-1
        TL(N)=((NP+NP+M+M+1.D0)*X*TL(NP)+(NP+1.D0)*TL(NP+1))
     &        /DFLOAT(NP+M+M)
       ENDDO
       NMAX=NMAXP
       RETURN
       END


        FUNCTION FACTCO(I,PN,M,OVER)
        IMPLICIT REAL*8 (A-H,O-Z)
        FACTCO=1.D0/PN
        IF (M.GT.0) THEN
          J=M+M
          DO WHILE ((J.GT.1).AND.(FACTCO.LT.OVER))
            FACTCO=FACTCO*(I+J)
            J=J-1
          END DO  
          IF (J.GT.2) FACTCO=0.
        ELSE
          FACTCO=1.D0/(I+1.D0)/PN
        END IF
        RETURN          
        END


        FUNCTION FACTCONEW(N,PN,M)
        IMPLICIT REAL*8 (A-H,O-Z)
        FACTCONEW=1.D0/(N+1.D0)/PN
        IF (M.GT.0) THEN
         J=M+M
         DO L=1,N
           FACTCONEW=FACTCONEW*DFLOAT(J+L)/DFLOAT(L)
         END DO  
        END IF
        RETURN          
        END
       
 


      
       

       
      
       
       
