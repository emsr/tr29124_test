
      PROGRAM SPHEROIDAL
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PL(0:2001),QL(0:2001)
      NMAX=2000
      M=2500
      MODE=2
      OPEN(10,FILE='data.dat',STATUS='NEW')
      WRITE(10,*)'PROLATE HARMONICS, M=2500, NMAX=2000, MODE 2'
      WRITE(10,*)
      WRITE(10,*)'PL2(N)=PL(N)/(X*X-1.)**M/2/(2M-1)!!*10**(-290)'
      WRITE(10,*)'QL2(N)=QL(N)*(X*X-1.)**M/2/(2M)!!*10**(290)'
      WRITE(10,*)
      WRITE(10,30)'X','N','PL2(N)',
     & 'QL2(N)','QL2(0)'
      WRITE(10,*)
      DO X=1.2D0,10.D0,1.6D0
         CALL DPROH(X,M,NMAX,MODE,PL,QL,NUEVO)
         WRITE(10,800) X,NUEVO,PL(NUEVO),QL(NUEVO),QL(0)
         NM=10
         WRITE(10,800) X,NM,PL(NM),QL(NM),QL(0)
      END DO
      M=120
      MODE=1
      WRITE(10,*)
      WRITE(10,*)'PROLATE HARMONICS, M=120, NMAX=2000, MODE 1'
      WRITE(10,*)
      WRITE(10,31)'X','N','PL(N)/(2M-1)!!',
     &             'QL(N)/(2M)!!','QL(0)/(2M)!!'
      WRITE(10,*)
      DO X=1.2D0,10.D0,1.6D0
         CALL DPROH(X,M,NMAX,MODE,PL,QL,NUEVO)
         WRITE(10,800) X,NUEVO,PL(NUEVO),QL(NUEVO),QL(0)
         NM=10
         WRITE(10,800) X,NM,PL(NM),QL(NM),QL(0)
      END DO
      WRITE(10,*)
      NMAX=2000
      M=50
      MODE=0
      WRITE(10,*)'OBLATE HARMONICS, M=50, NMAX=2000, MODE 0'
      WRITE(10,*)
      WRITE(10,32)'X','N', 'RL(N)', 'TL(N)', 'TL(0)'
      WRITE(10,*)
      DO X=0.2D0,10.D0,1.6D0
         CALL DOBLH(X,M,NMAX,MODE,PL,QL,NUEVO)
         WRITE(10,801) X,NUEVO,PL(NUEVO),QL(NUEVO),QL(0)
         NM=10
         WRITE(10,801) X,NM,PL(NM),QL(NM),QL(0)
      ENDDO
      WRITE(10,*)
   30 FORMAT (4X,A1,3X,A1,2X,3(8X,A7))
   31 FORMAT (4X,A1,3X,A1,5X,A14,2(2X,A14))
   32 FORMAT (4X,A1,3X,A1,3(9X,A6))
  800 FORMAT (F6.1,1X,I4,1X,3D16.9)
  801 FORMAT (F6.1,1X,I4,1X,3D16.9)
      CLOSE(10)
      END
