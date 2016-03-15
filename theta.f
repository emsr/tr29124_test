
      SUBROUTINE THETA(SOFN, N, A, B, LENGA, LENGB, ESTLIM)
        DIMENSION A(0:LENGA), B(0:LENGB)
        PARAMETER ( HUGE = 1.E+60 , TINY = 1.E-60 )
        PARAMETER ( ZERO = 0.E0 , ONE = 1.E0 , TWO = 2.E0 )
        JMAX = (2 * N + 1) / 3
        NMOD2 = MOD(N,2)
        IF (N .EQ. 0) THEN
          A(0) = SOFN
          ESTLIM = SOFN
          RETURN
        END IF
        IF ( NMOD2 .EQ. 0) THEN
          AUX2 = ZERO
          AUX1 = A(0)
          A(0) = SOFN
          DO J = 1, JMAX
            AUX3 = AUX2
            AUX2 = AUX1
            IF ( J .LT. JMAX ) THEN
              AUX1 = A(J)
            END IF
            IF ( MOD(J,2) .EQ. 0 ) THEN
              DENOM = A(J-1) - TWO * B(J-1) + AUX2
              IF ( ABS(DENOM) .LT. TINY ) THEN
                A(J) = HUGE
              ELSE
                A(J) = AUX3 + ( B(J-2) - AUX3 ) * ( A(J-1) - B(J-1) ) / DENOM
              END IF
            ELSE
              DIFF = A(J-1) - B(J-1)
              IF ( ABS(DIFF) .LT. TINY ) THEN
                A(J) = HUGE
              ELSE
                A(J) = AUX3 + ONE / DIFF
              END IF
            END IF
          END DO
          IF ( MOD(JMAX,2) .EQ. 0 ) THEN
            ESTLIM = A(JMAX)
          ELSE
            ESTLIM = A(JMAX-1)
          END IF
        ELSE
          AUX2 = ZERO
          AUX1 = B(0)
          B(0) = SOFN
          DO J = 1, JMAX
            AUX3 = AUX2
            AUX2 = AUX1
            IF ( J .LT. JMAX ) THEN
              AUX1 = B(J)
            END IF
            IF ( MOD(J,2) .EQ. 0 ) THEN
              DENOM = B(J-1) - TWO * A(J-1) + AUX2
              IF ( ABS(DENOM) .LT. TINY ) THEN
                B(J) = HUGE
              ELSE
                B(J) = AUX3 + ( A(J-2) - AUX3 ) * ( B(J-1) - A(J-1) ) / DENOM
              END IF
            ELSE
              DIFF = B(J-1) - A(J-1)
              IF ( ABS(DIFF) .LT. TINY ) THEN
                B(J) = HUGE
              ELSE
                B(J) = AUX3 + ONE / DIFF
              END IF
            END IF
          END DO
          IF ( MOD(JMAX,2) .EQ. 0 ) THEN
            ESTLIM = B(JMAX)
          ELSE
            ESTLIM = B(JMAX-1)
          END IF
        END IF
        RETURN
      END SUBROUTINE

      SUBROUTINE THETIT(SOFN, N, ARJ, LARRAY, ESTLIM)
        DIMENSION ARJ(0:LARRAY)
        PARAMETER ( HUGE = 1.E+60 , TINY = 1.E-60 )
        ARJ(N) = SOFN
        IF (N .LT. 3) THEN
          ESTLIM = SOFN
        ELSE
          LMAX = N/3
          M = N
          DO L = 1, LMAX
            M = M - 3
            DIFF0 = ARJ(M+1) - ARJ(M)
            DIFF1 = ARJ(M+2) - ARJ(M+1)
            DIFF2 = ARJ(M+3) - ARJ(M+2)
            DENOM = DIFF2 * (DIFF1 - DIFF0) - DIFF0 * (DIFF2 - DIFF1)
            IF (ABS(DENOM) .LT. TINY) THEN
              ARJ(M) = HUGE
            ELSE
              ARJ(M) = ARJ(M+1) - DIFF0 * DIFF1 * (DIFF2 - DIFF1) / DENOM
            END IF
          END DO
          ESTLIM = ARJ(MOD(N,3))
        END IF
        RETURN
      END SUBROUTINE
