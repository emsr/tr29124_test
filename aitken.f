C  Iterative Aitken Delta^2

      SUBROUTINE AITKEN(SOFN,N,A,LARRAY,ESTLIM)
        DIMENSION A(0:LARRAY)
        PARAMETER ( HUGE = 1.E+60 , TINY = 1.E-60 , TWO = 2.E0 )
        A(N) = SOFN
        IF (N .LT. 2) THEN
          ESTLIM = SOFN
        ELSE
          LOWMAX = N/2
          DO J = 1,LOWMAX
            M = N - 2*J
            DENOM = A(M+2) - TWO*A(M+1) + A(M)
            IF (ABS(DENOM) .LT. TINY) THEN
              A(M) = HUGE
            ELSE
              A(M) = A(M) - (A(M) - A(M+1))**2 / DENOM
            END IF
          ENd DO
          IF (MOD(N,2) .EQ. 0) THEN
            ESTLIM = A(0)
          ELSE
            ESTLIM = A(1)
          END IF
        END IF
        RETURN
      END SUBROUTINE
