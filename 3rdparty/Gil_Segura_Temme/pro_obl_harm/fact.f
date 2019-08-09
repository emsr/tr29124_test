

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
      END FUNCTION

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
      END FUNCTION
