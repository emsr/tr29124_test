C Wynn epsilon algorithm.
C
      SUBROUTINE EPSAL(SOFN,N,E,LARRAY,ESTLIM)
        DIMENSION E(0:LARRAY)
        PARAMETER (HUGE = 1.E+60, TINY = 1.E-60, ZERO = 0.E0, ONE = 1.E0)
        E(N) = SOFN
        IF (N .EQ. 0) THEN
          ESTLIM = SOFN
        ELSE
          AUX2 = ZERO
          DO J = N,1,-1
            AUX1 = AUX2
            AUX2 = E(J-1)
            DIFF = E(J) - AUX2
            IF (ABS(DIFF) .LE. TINY) THEN
              E(J-1) = HUGE
            ELSE
              E(J-1) = AUX1 + ONE/DIFF
            END IF
          END DO
          IF ( MOD(N,2) .EQ. 0 ) THEN
            ESTLIM = E(0)
          ELSE
            ESTLIM = E(1)
          END IF
        END IF
        RETURN
      END SUBROUTINE
