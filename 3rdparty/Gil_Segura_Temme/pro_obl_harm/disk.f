       PROGRAM DISK
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION PL(0:2001),QL(0:2001),POL(0:2001)	
       PARAMETER(EPS=1.D-8)     
       PI=3.14159265358979323D0
       NMAX=500
       WRITE(6,*)'INTRODUCE DISTANCE AND RADIUS' 
       WRITE(6,*)'L,A'
 1     READ(5,*)DL,A
       DLA=DL/A
       IF (DLA.LT.0.D0) THEN
         WRITE(6,*)'IMPROPER ARGUMENTS FOR THE RATIO.'
         WRITE(6,*)'MUST BE GREATER THAN ZERO'
         GOTO 1
       END IF
       DLA2=DLA*DLA
       A2=A*A
       WRITE(6,*)'INTRODUCE R AND CHARGE Q'
       WRITE(6,*)'R/A MUST BE LESS THAN ONE'
 2     READ(5,*)R,Q
       RA=R/A
       IF (RA.GE.1.D0) THEN
          WRITE(6,*)'IMPROPER ARGUMENTS FOR THE RATIO.'
          WRITE(6,*)'MUST BE LESS THAN ONE'
          GOTO 2
       END IF
       SINBET=RA
       COSBET=DSQRT(1.D0-SINBET*SINBET)
       DFACT=-Q/PI/A2  
       DFACT1=1.D0/DLA2*(1.D0+1.D0/DLA2*SINBET*SINBET)**(-1.5)
       DFACT2=2.D0/DSQRT(PI)/COSBET
       CALL DOBLH(DLA,0,NMAX,0,PL,QL,NUEVO)
       CALL POLIN(COSBET,NMAX,POL) 
       JJ=0
       NSUP=NUEVO
       OPEN(10,FILE='disk.dat',STATUS='NEW')       
       WRITE(10,*)'NSUP',NSUP        
       DSUME=0.D0
       DDS=1.D0
       PRE=1.D0
       K=0
       DO WHILE ((JJ.LE.NSUP).AND.(PRE.GT.EPS))
         JJ=2*K
         QLM=QL(JJ)
         DSUME=DSUME+FA(QLM,K)*(4.D0*K+1.D0)*2.D0**K
     &    *POL(JJ)/DSQRT(PI)
         PRE=ABS(1.D0-DSUME/DDS)
         DDS=DSUME
         K=K+1
       END DO
       TERM1=DFACT*DFACT1
       TERM2=DFACT*DFACT2*DSUME
       DENSITY=TERM1+TERM2
       WRITE(10,*)'NUMBER OF TERMS',JJ
       WRITE(10,*)'TERM1,TERM2',TERM1,TERM2
       WRITE(10,*)'SURFACE CHARGE DENSITY'
       WRITE(10,*) DENSITY
       CLOSE(10)
       END
   



       SUBROUTINE POLIN(X,NMAX,PL)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION PL(0:NMAX+1)
       PL(0)=1.D0
       PL(1)=X
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                         C
C   WE USE THE RECURRENCE RELATIONS                                       C
C                                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO NP=1,NMAX
         PL(NP+1)=((2.D0*NP+1.D0)*X*PL(NP)-NP*PL(NP-1))/(NP+1.D0)
       ENDDO
       RETURN
       END
 
      

       FUNCTION FA(QLM,N)
       IMPLICIT REAL*8 (A-H,O-Z)
       FA=QLM
        IF (N.GT.0) THEN
          J=1 
          DO I=1,N
            FA=FA*I/(2.D0*N-J)
            J=J+2
          END DO 
        END IF
       RETURN          
       END
