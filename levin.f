      SUBROUTINE GLEVIN(SOFN,ROFN,BETA,N,ARUP,ARLO,LARRAY,ESTLIM)
        DIMENSION ARUP(0:LARRAY),ARLO(0:LARRAY)
        PARAMETER ( HUGE = 1.E+60 , TINY = 1.E-60 , ONE = 1.E0 )
        ARUP(N) = SOFN / ROFN
        ARLO(N) = ONE / ROFN
        IF ( N .GT. 0 ) THEN
          ARUP(N-1) = ARUP(N) - ARUP(N-1)
          ARLO(N-1) = ARLO(N) - ARLO(N-1)
          IF (N .GT. 1) THEN
            BN1 = BETA + FLOAT(N-1)
            BN2 = BETA + FLOAT(N)
            COEF = BN1 / BN2
            DO J = 2, N
              FACT = (BETA+FLOAT(N-J)) * COEF**(J-2) / BN2
              ARUP(N-J) = ARUP(N-J+1) - FACT * ARUP(N-J)
              ARLO(N-J) = ARLO(N-J+1) - FACT * ARLO(N-J)
            END DO
          END IF
        END IF
        IF (ABS(ARLO(0)) .LT. TINY) THEN
          ESTLIM = HUGE
        ELSE
          ESTLIM = ARUP(0) / ARLO(0)
        END IF
        RETURN
      END SUBROUTINE
