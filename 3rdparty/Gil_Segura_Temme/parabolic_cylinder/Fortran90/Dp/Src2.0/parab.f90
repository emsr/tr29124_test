  MODULE Parabolic
  IMPLICIT NONE
  INTEGER, PARAMETER  :: r8 = KIND(0.0d0)
  PRIVATE
  PUBLIC  :: parab
  CONTAINS       
    RECURSIVE SUBROUTINE parab(a,x,mode,uaxx,vaxx,ierr)
    ! ---------------------------------------------------------
    ! Calculation of the real parabolic cylinder functions
    !          U(a,x), V(a,x), x,a real and x>=0.
    ! These functions are independent solutions of the 
    ! differential equation
    !          d^2w/dx^2-(x^2/4+a)w=0
    ! The solutions U(a,x) and V(a,x), with initial conditions
    !  U(a,0)=pi**(1/2)/(2**(a/2+1/4)*Gamma(3/4+a/2)),  
    !  U'(a,0)=-pi**(1/2)/(2^{a/2-1/4}*Gamma (1/4+a/2)),
    !  V(a,0)=2**(a/2+1/4)*sin(pi*(3/4-a/2))/Gamma(3/4-a/2) 
    !  V'(a,0)=2**(a/2+1/4)*sin(pi*(1/4-a/2))/Gamma(1/4-a/2)
    !  constitute a numerically satisfactory pair for x>0.        
    ! ----------------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   uaxx , uaxx(1), U(a,x)
    !          uaxx(2), U'(a,x)
    !   vaxx , vaxx(1), V(a,x)
    !          vaxx(2), V'(a,x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems
    !          ierr=2, the argument x is out of range 
    ! -----------------------------------------------------------
    !           METHODS OF COMPUTATION
    ! -----------------------------------------------------------
    ! The present code uses different methods of computation
    ! depending on the values of a and x:
    !
    ! When (2.5< x<12.5) and (|a|<0.7) we compute U(a,x), V(a,x)
    !  and their derivatives using:
    !  - For U(a,x) (and its derivative), recurrence relations.   
    !  - For V(a,x) (and its derivative), series when x<=10.5 
    !        and recurrence relations when x>10.5
    ! In other cases:
    !
    ! 1) For a>0:
    !     When x<0.185781526149...
    !       For a>-0.23x*x+1.2x+18.72,
    !          For a<1.0/(x*x) and x<0.005 we compute: V(a,x) and its 
    !             derivative using series; U(a,x) and its
    !             derivative using asymptotic expansions in
    !             terms of elementary functions
    !          In other case, we use asymptotic expansions 
    !             in terms of elementary functions
    !       For a<-0.23x*x+1.2x+18.72, we use power series
    !     When 0.185781526149..<=x<3
    !        For a< 3.75/x-1.25, we use power series
    !        For 3.75/x-1.25 <= a<-0.23x*x+1.2x+18.72 
    !          For V(a,x) and its derivative, we use 
    !                     power series
    !          For U(a,x) and its derivative, we use 
    !                     recurrence relations         
    !        For a >= -0.23x*x+1.2x+18.72, we use asymptotic 
    !          expansions in terms of elementary functions 
    !     When 3<x<12
    !        For a< -0.23x*x+1.2x+18.72, 
    !          For V(a,x) and its derivative, we use power series
    !          For U(a,x) and its derivative, we use 
    !                     recurrence relations
    !        For a >= -0.23x*x+1.2x+18.72, we use asymptotic 
    !           expansions in terms of elementary functions   
    !     When x>=12
    !        For a < 2.5x-30 and a<150, we use asymptotic 
    !                                   expansions Poincare-type   
    !        For a >= 2.5x-30 or a>150, we use asymptotic expansions
    !            in terms of elementary functions
    ! 2) For a<0:
    !     When x<0.84483298487.. 
    !        For a <-0.21x*x-4.5x-40, we use asymptotic 
    !            expansions in terms of elementary functions 
    !        For a >=-0.21x*x-4.5x-40, we use power series
    !     When 0.84483298487..<= x<3
    !        For a < -0.21x*x-4.5x-40, we use asymptotic 
    !           expansions in terms of elementary functions 
    !        For -0.21x*x-4.5x-40<= a<-30/(x-0.3)+100/9, we use 
    !           quadrature 
    !        For a>=-30/(x-0.3)+100/9, we use power series
    !     When 3<=x<4 
    !        For a<-0.21x*x-4.5x-40, we use asymptotic 
    !           expansions in terms of elementary functions 
    !        For a>=-0.21x*x-4.5x-40, we use quadrature
    !     When 4<=x<43/35+3sqrt(155)/7 
    !        For a<-0.21x*x-4.5x-40, we use asymptotic 
    !            expansions in terms of elementary functions 
    !        For -0.21x*x-4.5x-40<=a< x-14, we use Airy-type 
    !            asymptotic expansions
    !        For a>=x-14, we use quadrature
    !     When 43/35+3sqrt(155)/7<= x<12
    !        For a< -0.21x*x-4.5x-40, we use asymptotic 
    !           expansions in terms of elementary functions
    !        For -0.21x*x-4.5x-40<=a<-7-0.14(x-4.8)**2, we use
    !           Airy-type asymptotic expansions 
    !        For a>= -7-0.14(x-4.8)*(x-4.8), we use quadrature
    !     When 12 <=x<-4/5+2/15*sqrt(18127)
    !        For a<-0.21x*x-4.5x-40, we use asymptotic 
    !           expansions in terms of elementary functions
    !        For -0.21x*x-4.5x-40<=a<-7-0.14(x-4.8)**2, we use 
    !           Airy-type asymptotic expansions 
    !        For -7-0.14(x-4.8)**2<=a<-0.23x*x+1.2x+18.72, we use
    !           quadrature
    !        For -0.23x*x+1.2x+18.72<=a<-2.5x+30, we use asymptotic 
    !          expansions in terms of elementary functions
    !        For a>=-2.5x+30, we use asymptotic expansions Poincare-type
    !     When -4/5+2/15*sqrt(18127)<=x<30
    !        For a<-0.21x*x-4.5x-40, we use asymptotic 
    !           expansions in terms of elementary functions
    !        For -0.21x*x-4.5x-40<=a<-0.23x*x+1.2x+18.72, we use 
    !           Airy-type asymptotic expansions
    !        For -0.23x*x+1.2x+18.72<=a<-2.5x+30, we use asymptotic 
    !           expansions in terms of elementary functions
    !        For a>=-2.5x+30, use asymptotic expansions Poincare-type
    !     When x >= 30
    !        For a< -0.295x*x+0.3x-107.5, we use asymptotic 
    !          expansions in terms of elementary functions
    !        For -0.295x*x+0.3x-107.5<=a<-0.1692*x*x, we use 
    !          Airy-type asymptotic expansions
    !        For -0.1692*x*x<=a<-2.5x+30, we use asymptotic 
    !          expansions in terms of elementary functions
    !        For a>=-2.5x+30, we use asymptotic expansions Poincare-type
    !          when |a|<150; in other cases, we use asymptotic 
    !          expansions in terms of elementary functions
    ! ----------------------------------------------------------------  
    !                  ACCURACY  
    !-----------------------------------------------------------------
    !  The aimed relative accuracy for scaled functions is better than 
    !  5.0e-14. Exceptions to this accuracy are the evaluation of the 
    !  function near their zeros and the error caused by the evaluation 
    !  of trigonometric functions of large arguments when |a|>>x.
    !  The routines always produce values for which
    !  the Wronskian relation is verified with a relative accuracy 
    !  better than 5.0e-14. The accuracy of the unscaled functions 
    !  is also better than 5.0e-14 for moderate values of x and a
    !  (except close to the zeros), while for large x and a
    !  the error is dominated by exponential and trigonometric 
    !  function evaluations.
    !  For IEEE standard double precision arithmetic, the
    !  accuracy is better than 5.0e-13 in the computable range of 
    !  unscaled PCFs (except close to the zeros). 
    ! ----------------------------------------------------------------
    ! Authors:
    !  Amparo Gil    (U. Cantabria, Santander, Spain)
    !                 e-mail: amparo.gil@unican.es
    !  Javier Segura (U. Cantabria, Santander, Spain)
    !                 e-mail: javier.segura@unican.es
    !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
    !                 e-mail: nico.temme@cwi.nl
    ! -------------------------------------------------------------
    !  References:
    !  1) Integral representations for computing real 
    !     parabolic cylinder functions.
    !     A. Gil, J. Segura, N.M. Temme         
    !     Numerische Mathematik 98(1) (2004) 105-134
    !  2) Numerical and asymptotic aspects of parabolic 
    !     cylinder functions.
    !     N.M. Temme  
    !     Journal of Computational and Applied Mathematics 
    !     (121), 221--246, 2000.
    ! -------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: mode
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT) :: uaxx(2)
    REAL(r8), INTENT(OUT) :: vaxx(2)   
    REAL(r8) :: uaxx1(2),vaxx1(2),x2,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10
    INTEGER,  INTENT(OUT) :: ierr
    ierr=0
    IF (x<0) THEN
      uaxx(1)=0.0_r8
      uaxx(2)=0.0_r8
      vaxx(1)=0.0_r8
      vaxx(2)=0.0_r8
      ierr=2
    ENDIF
    IF (ierr==0) THEN
      IF (((x>2.5_r8).AND.(x<12.5_r8)).AND.(abs(a)<0.7_r8)) THEN
        IF (x>10.5_r8) THEN
          ! Computation of V(a,x)
          CALL recurv(a,x,mode,vaxx,ierr)
        ELSE
          CALL series(a,x,mode,uaxx,vaxx,ierr)
        ENDIF
        CALL recuru(a,x,mode,uaxx,ierr)
      ELSE
        IF (a>=0) THEN
          IF (x<.1857815261497950467_r8) THEN
            x2=x*x
            f1=-0.23_r8*x2+1.2_r8*x+18.72_r8
            IF (a>f1) THEN
              IF ((a<1.0_r8/x2).AND.(x<0.005_r8)) THEN
                CALL series(a,x,mode,uaxx1,vaxx,ierr)
                CALL expaelem(a,x,mode,uaxx,vaxx1,ierr)
              ELSE
                CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
              ENDIF   
            ELSE
              CALL series(a,x,mode,uaxx,vaxx,ierr)
            ENDIF      
          ELSEIF (x<3.0_r8) THEN
            f1=-0.23_r8*x*x+1.2_r8*x+18.72_r8
            f2=3.75_r8/x-1.25_r8
            IF (a<f2) THEN
              CALL series(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF(a<f1) THEN
              ! For V(a,x)
              CALL series(a,x,mode,uaxx,vaxx,ierr)
              ! For U(a,x)
              CALL recuru(a,x,mode,uaxx,ierr)
            ELSE
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ENDIF
          ELSEIF (x<12.0_r8) THEN
            f1=-0.23_r8*x*x+1.2_r8*x+18.72_r8
            IF (a<f1) THEN
              ! For V(a,x)
              CALL series(a,x,mode,uaxx,vaxx,ierr)
              ! For U(a,x)
              CALL recuru(a,x,mode,uaxx,ierr)
            ELSE
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ENDIF        
          ELSE
            f7=2.5_r8*x-30.0_r8
            IF ((a<f7).AND.(a<150)) THEN
              CALL expax(a,x,mode,uaxx,vaxx,ierr)
            ELSE
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ENDIF
          ENDIF
        ENDIF
        IF (a<0) THEN
          IF (x<0.8448329848762434824327_r8) THEN
            f4=-0.21_r8*x*x-4.5_r8*x-40.0_r8
            IF (a<f4) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSE
              CALL series(a,x,mode,uaxx,vaxx,ierr)
            ENDIF      
          ELSEIF (x<3.0_r8) THEN
            f4=-0.21_r8*x*x-4.5_r8*x-40.0_r8
            f3=-30.0_r8/(x-0.3_r8)+100.0_r8/9.0_r8
            IF (a<f4) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f3) THEN
              CALL uvax(a,x,mode,uaxx,vaxx,ierr)
            ELSE
              CALL series(a,x,mode,uaxx,vaxx,ierr)  
            ENDIF
          ELSEIF (x<4.0_r8) THEN
            f4=-0.21_r8*x*x-4.5_r8*x-40.0_r8
            IF (a<f4) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSE
              CALL uvax(a,x,mode,uaxx,vaxx,ierr)  
            ENDIF
          ELSEIF (x<43.0_r8/35.0_r8+3.0_r8*sqrt(155.0_r8)/7.0_r8) THEN
            f4=-0.21_r8*x*x-4.5_r8*x-40.0_r8          
            f5=x-14.0_r8
            IF (a<f4) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f5) THEN
              CALL expair(a,x,mode,uaxx,vaxx,ierr)             
            ELSE
              CALL uvax(a,x,mode,uaxx,vaxx,ierr)
            ENDIF
          ELSEIF (x<12.0_r8) THEN
            f4=-0.21_r8*x*x-4.5_r8*x-40.0_r8   
            f6=-7.0_r8-0.14_r8*(x-4.8_r8)*(x-4.8_r8)
            IF (a<f4) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f6) THEN
              CALL expair(a,x,mode,uaxx,vaxx,ierr)
            ELSE
              CALL uvax(a,x,mode,uaxx,vaxx,ierr)
            ENDIF
          ELSEIF (x<-4.0_r8/5.0_r8+2.0_r8/15.0_r8*sqrt(18127.0_r8)) THEN
            x2=x*x 
            f4=-0.21_r8*x2-4.5_r8*x-40.0_r8   
            f6=-7.0_r8-0.14_r8*(x-4.8_r8)*(x-4.8_r8)
            f1=-0.23_r8*x2+1.2_r8*x+18.72_r8
            f8=-2.5_r8*x+30.0_r8
            IF (a<f4) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f6) THEN
              CALL expair(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f1) THEN
              CALL uvax(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f8) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSE
              CALL expax(a,x,mode,uaxx,vaxx,ierr)
            ENDIF
          ELSEIF (x<30.0_r8) THEN
            x2=x*x
            f4=-0.21_r8*x2-4.5_r8*x-40.0_r8
            f1=-0.23_r8*x2+1.2_r8*x+18.72_r8
            f8=-2.5_r8*x+30.0_r8
            IF (a<f4) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f1) THEN
              CALL expair(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f8) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSE
              CALL expax(a,x,mode,uaxx,vaxx,ierr)
            ENDIF
          ELSE
            x2=x*x
            f10=-0.295_r8*x2+0.3_r8*x-107.5_r8 
            f9=-0.1692_r8*x2
            f8=-2.5_r8*x+30.0_r8
            IF (a<f10) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f9) THEN
              CALL expair(a,x,mode,uaxx,vaxx,ierr)
            ELSEIF (a<f8) THEN
              CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
            ELSE
              IF (abs(a)<150) THEN
                CALL expax(a,x,mode,uaxx,vaxx,ierr)
              ELSE
                CALL expaelem(a,x,mode,uaxx,vaxx,ierr)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    END SUBROUTINE parab

    SUBROUTINE recuru(a,x,mode,uaxx,ierr)
    !--------------------------------------------------
    ! Application of recurrences for computing U(a,x)
    !--------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   uaxx , uaxx(1), U(a,x)
    !          uaxx(2), U'(a,x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    !----------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: mode
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT) :: uaxx(2)  
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: eps, uaxx1(2), vaxx1(2), uaxx2(2), vaxx2(2), &
                a1, a2, a3, xargu, xhalf, gammaa
    INTEGER :: ierr1,ierr2,mode2
    eps=epss
    ierr=0
    mode2=0
    xhalf=0.5_r8*x
    a1=-(int(a)-a)+21; a2=1+a1;
    CALL parab(a1,x,mode2,uaxx1,vaxx1,ierr1)
    CALL parab(a2,x,mode2,uaxx2,vaxx2,ierr2)
    a3=a1-1
    uaxx(1)=(x*uaxx1(1)+(a3+1.5_r8)*uaxx2(1))
    uaxx(2)=-xhalf*uaxx(1)-(a3+0.5_r8)*uaxx1(1)   
    DO WHILE(abs(a3-a)>eps*100)
      uaxx2(1)=uaxx1(1)  
      uaxx1(1)=uaxx(1)
      a2=a1;
      a1=a3;
      a3=a1-1; 
      uaxx(1)=(x*uaxx1(1)+(a3+1.5_r8)*uaxx2(1))
      uaxx(2)=-xhalf*uaxx(1)-(a3+0.5_r8)*uaxx1(1)   
    ENDDO
    IF (mode==1) THEN
      xargu=x*x*0.25_r8+a
      IF (xargu <0) THEN
        gammaa=exp(a*0.5_r8*(log(abs(a))-1.0_r8))
        uaxx(1)=uaxx(1)*gammaa
        uaxx(2)=uaxx(2)*gammaa
      ELSE
        gammaa=exp(a*log(x*0.5_r8+sqrt(xargu))+&
               x*0.5_r8*sqrt(xargu)-a*0.5_r8)
        uaxx(1)=uaxx(1)*gammaa
        uaxx(2)=uaxx(2)*gammaa
      ENDIF  
    ENDIF
    END SUBROUTINE recuru

    SUBROUTINE recurv(a,x,mode,vaxx,ierr)
    !-------------------------------------------------
    ! Application of recurrences for computing V(a,x)
    !-------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   vaxx , vaxx(1), V(a,x)
    !          vaxx(2), V'(a,x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    !----------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: mode
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT) :: vaxx(2)  
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: eps, uaxx1(2), vaxx1(2), uaxx2(2), vaxx2(2),&
                a1, a2, a3, xhalf, xargu, gammaa
    INTEGER :: ierr1,ierr2,mode2
    eps=epss
    ierr=0
    mode2=0
    xhalf=0.5_r8*x
    a1=-(int(a)-a)-20; a2=-1+a1;
    CALL parab(a1,x,mode2,uaxx1,vaxx1,ierr1)
    CALL parab(a2,x,mode2,uaxx2,vaxx2,ierr2)
    a3=a1+1
    vaxx(1)=x*vaxx1(1)+(a3-1.5_r8)*vaxx2(1)
    vaxx(2)=xhalf*vaxx(1)+(a3-0.5_r8)*vaxx1(1)
    DO WHILE(abs(a3-a)>eps*100)
      vaxx2(1)=vaxx1(1)  
      vaxx1(1)=vaxx(1)
      a2=a1;
      a1=a3;
      a3=a1+1;  
      vaxx(1)=x*vaxx1(1)+(a3-1.5_r8)*vaxx2(1)
      vaxx(2)=xhalf*vaxx(1)+(a3-0.5_r8)*vaxx1(1)
    ENDDO
    IF (mode==1) THEN
      xargu=x*x*0.25_r8+a
      IF (xargu <0) THEN
        gammaa=exp(a*0.5_r8*(log(abs(a))-1.0_r8))
        vaxx(1)=vaxx(1)/gammaa
        vaxx(2)=vaxx(2)/gammaa
      ELSE
        gammaa=exp(a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8)
        vaxx(1)=vaxx(1)/gammaa
        vaxx(2)=vaxx(2)/gammaa
      ENDIF
    ENDIF
    END SUBROUTINE recurv

    SUBROUTINE uvax(a,x,mode,uaxx,vaxx,ierr)
    ! --------------------------------------------------------
    ! Calculation of U(a,x), V(a,x) by using quadrature rules
    ! --------------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   uaxx , uaxx(1), U(a,x)
    !          uaxx(2), U'(a,x)
    !   vaxx , vaxx(1), V(a,x)
    !          vaxx(2), V'(a,x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! ---------------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: mode
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT) :: uaxx(2)
    REAL(r8), INTENT(OUT) :: vaxx(2)   
    REAL(r8) :: eps, t, aa, xargu, ffa, twoaeta, mu,& 
                a14, c, s, sq, w0, sqw0, sinpia, ijax(4), gh(6),&
                argum, a14sqrt2opi, dl, dnew, mub, mu2,&
                mu4, mus, smus, gmus, beta, gammaa2
    REAL(r8), DIMENSION(0:25) :: gk
    REAL(r8), DIMENSION(0:12) :: wk
    INTEGER k,ierr
    DATA gk/1.0_r8,.41666666666666666667e-1_r8,0,-.97463348765432098765e-2_r8,0,&
         .12318482904510826965e-1_r8,0,-.37822933917705539682e-1_r8,0,&
         .21514326767209896474_r8,0,-1.9630003732872175294_r8,0,&
         26.256830962378916652_r8,0,-484.19061617504532506_r8,0,&
         11773.948564802034554_r8,0,-365037.92569092371983_r8,0,&
         14054558.808383655048_r8,0,-657894020.31009296820_r8,0,&
         36795321248.737494263_r8/
    DATA wk/1.0_r8,-.17361111111111111111e-2_r8,.81219457304526748971e-3_r8,&
         -.11215312855682914598e-2_r8,.33920312789254653178e-2_r8,&
         -.18817620620360951108e-1_r8,.16870892343666095644_r8,&
         -2.2330644165264046334_r8,40.925670821571152909_r8,&
         -991.44221631790969480_r8,30664.092693023833555_r8,&
         -1178670.6504836921080_r8,55108658.068376820265_r8/
    eps=epss
    ierr=0
    IF (a > 0) THEN
      xargu=abs(x*x*0.25_r8+a)
      IF (mode==0) THEN
      ! ------------------------------------------------
      ! Check for possible overflow/underflow problems
      ! ------------------------------------------------
        IF (a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8 > log(giant)) THEN
          ierr=1
          uaxx(1)=dwarf
          uaxx(2)=dwarf
          vaxx(1)=giant 
          vaxx(2)=giant
        ENDIF
      ENDIF
      IF (ierr==0) THEN
        mu=sqrt(a); a14=sqrt(mu); 
        sinpia= sin(pi*a); 
        t=0.5*x/mu; 
        sq=sqrt(t*t+1.0_r8); w0=t+sq; sqw0=sqrt(w0)
        ! --------------------------------
        ! Calculation of Scaled functions
        ! --------------------------------
        CALL ij(a,x,eps,ijax)
        uaxx(1)=a14*ijax(1)/sqrttwopi         ! U(a,x)
        uaxx(2)=-mu*a14*ijax(3)/sqrttwopi     ! U'(a,x)
        IF (a>50) THEN
          ! -----------------------------
          ! Computation of G(mu) and S(mu)
          ! G(mu)=1/gmus
          !-----------------------------
          mub=sqrt(2*a)
          mu2=mub*mub
          mu4=mu2*mu2
          DO k=0,25
            IF (k==0) THEN
              mus=1.0_r8
              gmus=0.0_r8
            ELSE
              mus=mus*mu2
            ENDIF
            gmus=gmus+gk(k)/mus
          ENDDO
          DO k=0,12
            IF (k==0) THEN
              mus=1.0_r8
              smus=0.0_r8
            ELSE
              mus=mus*mu4
            ENDIF
            smus=smus+wk(k)/mus
          ENDDO
          dl=-2.0_r8*a*(t*sq+log(t+sq))
          IF (dl<log(dwarf)) THEN
            dnew=0.0_r8
          ELSE         
            dnew=exp(dl)*sqrttwopi*smus/(gmus*gmus)
          ENDIF
        ELSE
          ffa=exp(a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8)          
          dnew=(gamma(a+0.5_r8)/ffa)/ffa
        ENDIF
        ! V(a,x)
        vaxx(1)=a14*(sinpia*ijax(1)/sqrttwopi*dnew &
                       +sqw0*ijax(2))*onepi 
        ! V'(a,x)
        vaxx(2)=mu*a14*(-sinpia*ijax(3)/sqrttwopi*dnew&
                       +sqw0*ijax(4))*onepi 
        IF (mode==0) THEN  
          ffa=exp(a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8)
          uaxx(1)=uaxx(1)/ffa
          uaxx(2)=uaxx(2)/ffa
          vaxx(1)=ffa*vaxx(1)
          vaxx(2)=ffa*vaxx(2)
        ENDIF
      ENDIF
    ELSE 
      aa=-a 
      mu=sqrt(aa); a14=sqrt(mu) 
      t=0.5_r8*x/mu
      IF (t <= 1) THEN
        IF (mode==0) THEN
        ! ------------------------------------------------
        ! Check for possible overflow/underflow problems
        ! ------------------------------------------------
          IF (a*0.5_r8*(log(abs(a))-1.0_r8) <log(dwarf)) THEN
            ierr=1
            uaxx(1)=giant
            uaxx(2)=giant
            vaxx(1)=dwarf
            vaxx(2)=dwarf
          ENDIF
        ENDIF
        IF (ierr==0) THEN  
          ! ---------------------------------
          ! Calculation of Scaled Functions
          ! ---------------------------------
          sq=sqrt((1.0_r8-t)*(1.0_r8+t))
          a14sqrt2opi=sqrt2opi*a14   
          CALL gh12(a,x,eps,gh) 
          twoaeta=aa*(acos(t)-t*sq)
          argum=twoaeta-pikwart
          c=cos(argum) 
          s=sin(argum)
          uaxx(1)=a14sqrt2opi*(c*gh(1)-s*gh(2))        ! U(a,x)
          uaxx(2)=mu*a14sqrt2opi*(c*gh(3)-s*gh(4))     ! U'(a,x)
          IF (abs(a)>50) THEN
            mub=sqrt(-2.0_r8*a)
            mu2=mub*mub
            mu4=mu2*mu2
            DO k=0,25
              IF (k==0) THEN
                mus=1.0_r8
                gmus=0.0_r8
              ELSE
                mus=mus*mu2
              ENDIF
              gmus=gmus+gk(k)/mus
            ENDDO
            DO k=0,12
              IF (k==0) THEN
                mus=1.0_r8
                smus=0.0_r8
              ELSE
                mus=mus*mu4
              ENDIF
              smus=smus+wk(k)/mus
            ENDDO
            beta=gmus*gmus/(sqrttwopi*smus)
          ELSE
            beta=1.0_r8/(sqrttwopi*aux9(aa))
          ENDIF 
          vaxx(1)=a14sqrt2opi*beta*(-s*gh(1)-c*gh(2))     ! V(a,x)
          vaxx(2)=mu*a14sqrt2opi*beta*(-s*gh(3)-c*gh(4))  ! V'(a,x)
          IF (mode==0) THEN
           ! -----------------------------------
           ! Calculation of Unscaled Functions
           ! -----------------------------------
            gammaa2=exp(a*0.5_r8*(log(abs(a))-1.0_r8))
            uaxx(1)=uaxx(1)/gammaa2
            uaxx(2)=uaxx(2)/gammaa2
            vaxx(1)=vaxx(1)*gammaa2
            vaxx(2)=vaxx(2)*gammaa2
          ENDIF
        ENDIF       
      ELSE
        xargu=abs(x*x*0.25_r8+a)
        IF (mode==0) THEN
           ffa=a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8
        ! ------------------------------------------------
        ! Check for possible overflow/underflow problems
        ! ------------------------------------------------
          IF ((ffa < log(dwarf)).OR.(ffa>log(giant))) THEN
            ierr=1
            uaxx(1)=giant
            uaxx(2)=giant
            vaxx(1)=dwarf 
            vaxx(2)=dwarf
          ENDIF
        ENDIF
        IF (ierr==0) THEN 
        ! ---------------------------------
        ! Calculation of Scaled Functions
        ! ---------------------------------
          CALL gh123(a,x,eps,gh)         
          a14sqrt2opi=sqrt2opi*a14  
          uaxx(1)=a14sqrt2opi*gh(1)         ! U(a,x)
          uaxx(2)=mu*a14sqrt2opi*gh(3)      ! U'(a,x)
          IF (abs(a)>50) THEN
            mub=sqrt(-2*a)
            mu2=mub*mub
            mu4=mu2*mu2
            DO k=0,25
              IF (k==0) THEN
                mus=1.0_r8
                gmus=0.0_r8
              ELSE
                mus=mus*mu2
              ENDIF
              gmus=gmus+gk(k)/mus
            ENDDO
            DO k=0,12
              IF (k==0) THEN
                mus=1.0_r8
                smus=0.0_r8
              ELSE
                mus=mus*mu4
              ENDIF
              smus=smus+wk(k)/mus
            ENDDO
            beta=gmus*gmus/(sqrttwopi*smus)
          ELSE
            beta=1.0_r8/(sqrttwopi*aux9(aa))
          ENDIF
          dl=2.0_r8*a*(t*sqrt(t*t-1.0_r8)-log(t+sqrt(t*t-1.0_r8)))
          IF (dl<log(dwarf)) THEN
            dnew=0.0_r8
          ELSE         
            dnew=exp(dl)*beta
          ENDIF
          vaxx(1)=a14sqrt2opi*(dnew*gh(2)+beta*gh(5))     ! V(a,x)
          vaxx(2)=mu*a14sqrt2opi*(dnew*gh(4)+beta*gh(6))  ! V'(a,x)
          IF (mode==0) THEN
            ! -----------------------------------
            ! Calculation of Unscaled Functions
            ! -----------------------------------
            ffa=exp(a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8)
            uaxx(1)=uaxx(1)/ffa
            uaxx(2)=uaxx(2)/ffa
            vaxx(1)=ffa*vaxx(1)
            vaxx(2)=ffa*vaxx(2) 
          ENDIF
       ENDIF
      ENDIF
    ENDIF
    END SUBROUTINE uvax
                                          
    SUBROUTINE ij(a,x,eps,ijax)
    ! ----------------------------------------------------
    ! Auxiliary routine for uvax
    ! ----------------------------------------------------
    ! Computation of the functions I,J,Id,Jd of relations 
    ! of GST(2004) on PCF
    ! ----------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(IN) :: eps
    REAL(r8), INTENT(OUT):: ijax(4)
    REAL(r8) :: mu, t, delta, eps10, lneps10, a0, b0, u0, v0, q0, &
                q1, w0, w02, w02p1,up, h, h2, upoverw0, sigmap,&
                sigman, wr, sum(4), newsum(4)
    INTEGER nc
    LOGICAL stopl, stopr
    mu=sqrt(a); t=x*0.5_r8/mu; up=sqrt(t*t+1.0_r8); w0=t+up 
    upoverw0=up/w0; w02=w0*w0; w02p1=w0*w0+1.0_r8
    sigmap=sqrt(2.0_r8/(a*w02p1)); sigman=sqrt(2.0_r8/a)/w0
    eps10=eps*0.1_r8; lneps10=1.1_r8*log(eps10) 
    ! Upper limit for the integral
    b0=sqrt(-lneps10);
    a0=-b0
    v0=b0*sqrt(w02p1)/w0
    u0=b0
    q0=(1.0_r8+sigman*u0)*0.5_r8/u0; q1=(1.0_r8-sigman*u0)*0.5_r8/u0
    wr=-16.0_r8*b0/(twopi*pi)
    !here begins trap(sum)
    h=b0/16.0_r8 
    CALL intrho(0.0_r8,sum,stopr,stopl)
    newsum=0.0_r8 
    CALL summing(h,h,newsum)  
    sum=h*(sum*0.5_r8+newsum)
    delta=1
    nc=0 
    DO WHILE (((delta>eps).AND.(nc<=5)).OR.(nc<4))  
      nc=nc+1; h2=h*0.5_r8       
      newsum=0.0_r8
      CALL  summing(h,h2,newsum) 
      newsum=h*newsum
      sum=(sum+newsum)*0.5_r8
      h=h2
      delta=abs((sum(1)*sum(4)+sum(2)*sum(3))*wr-1.0_r8)
    ENDDO
    !here ends trap(sum)
    u0=4.0_r8*v0*sigmap*sqrt(w0/pi)
    v0=2.0_r8*sigman/sqrtpi 
    ijax(1)=u0*sum(1); ijax(3)=-u0*w0*sum(3)   
    ijax(2)=v0*sum(2); ijax(4)=-v0*w0*sum(4)
    CONTAINS
      SUBROUTINE summing(h,d,newsum)
      IMPLICIT NONE
      REAL(r8) :: r,h,d,s,newsum(4),fr(4)
      LOGICAL stopl, stopr 
      INTEGER j 
      r=d; s=1; stopl=.false.; stopr=.false.
      DO WHILE (((s>eps10).AND.(abs(r)<b0)).AND.((.NOT.stopr).OR.(.NOT.stopl)))
        CALL intrho(r,fr,stopr,stopl)
        s=0.0_r8
        DO j=1,4  
          newsum(j)=newsum(j)+fr(j); s=s+abs(fr(j)) 
        ENDDO
        r=r+h
      ENDDO
      END SUBROUTINE summing
      SUBROUTINE intrho(rho,intr,stopr,stopl)
      IMPLICIT NONE
      REAL(r8) :: rho, intr(4), dpdq, psir, errorf, rho2,&
                  v, intv, intvd, u, q, intp, intpd,&
                  intn, intnd
      LOGICAL stopr, stopl, posx
      posx=.true.   
      dpdq=0.0_r8
      IF (rho==0) THEN 
        intr(1)=1; intr(3)=-upoverw0; 
        intr(2)=2.0_r8/q0; intr(4)=-2.0_r8*upoverw0/q0  
      ELSE
        errorf=errorfunction(rho, .false., .false.) ! {erf(rho)}
        intr=0 
        rho2=rho*rho 
        IF (.NOT.stopr) THEN
          v=v0*errorf
          IF (v + eps10 < v0) THEN 
            CALL integrp(rho2, dpdq, v, psir, intv, intvd,.true.) 
            IF (psir>-lneps10) stopr=.true.
            intr(1)=intv; intr(3)=intvd
          ENDIF
        ENDIF
        intp=0.0_r8; intpd=0.0_r8; intn=0.0_r8; intnd=0.0_r8; q=0.0_r8
        IF ((.NOT.stopl).OR.(.NOT.stopr)) THEN
          q=errorf
          IF (q+eps10 < 1) THEN
            IF (.NOT.stopr) THEN
              u=1.0_r8/(q0+q1*q); dpdq=q0*u*u; u=q*u
              CALL integrp(rho2,dpdq,u, psir, intp,intpd,.false.) 
              IF (psir > -lneps10) stopr=.true.
            ENDIF               
            IF (.NOT.stopl) THEN
              q=-q; u=1.0_r8/(q0+q1*q); dpdq=q0*u*u; u=q*u 
              CALL  integrp(rho2,dpdq,u, psir, intn, intnd, .false.)
              IF (psir>-lneps10) stopl=.true.
            ENDIF
          ENDIF
          intr(2)=intp+intn; intr(4)=intpd+intnd
        ENDIF  
      ENDIF
      END SUBROUTINE intrho
      SUBROUTINE integrp(rho2,dpdq,p,psir,int,intd,posx)
      IMPLICIT NONE
      REAL(r8) :: rho2, dpdq, p,int,intd, psii, psir, r,&
                  s, c, chi, u, v, v2
      LOGICAL posx 
      IF (posx) THEN
        v=p*sigmap; v2=v*v 
        psir=rho2+p*p-a*v2*v2*aux7(v2)*0.25_r8
        s=-v2*v*aux8(v)/3.0_r8
        psii=a*s
        chi=-psii-(s+v)*0.5_r8
        c=cos(chi); s=sin(chi) 
        IF (psir > -lneps10) THEN
          r=0.0_r8
        ELSE
          r=exp(-psir)/sqrt(sqrt(1.0_r8+v2)) 
        ENDIF
        int=r*c; intd=-(upoverw0*c-s*v)*r 
      ELSE 
        u=p*sigman 
        psir=rho2+p*p+0.5*a*u*u*aux7(u) 
        IF (abs(u+1.0_r8) < 10*mactol) THEN
           u=u+mactol*10
        ENDIF 
        IF (psir > -lneps10) THEN
          r=0 
        ELSE
          r=dpdq*exp(-psir)/sqrt(1.0_r8+u)
        ENDIF 
        int=r; intd=-(upoverw0+u)*r 
      ENDIF
      END SUBROUTINE integrp  
    END SUBROUTINE ij
                                              
    RECURSIVE FUNCTION aux1(x) RESULT(res1)
    USE Someconstants
    ! x - sin x = (x^3/6) aux1(x)
    INTEGER k
    REAL(r8) :: x,  ax, y, e, t, res1
    ax=abs(x)
    IF (ax > 1.5) THEN
      y=6.0_r8*(x - sin(x))/(x*x*x) 
    ELSEIF (ax> 0.5) THEN
      e=aux1(x/3.0_r8);
      t=2.0_r8-x*x*e/27.0_r8;
      y=(e+t*(t*t))/9.0_r8
    ELSEIF (ax < 1.0e-10_r8) THEN
      y=1 
    ELSE
      ax=x*x; t=-ax/20.0_r8; y= 1.0_r8+t; k= 6
      DO WHILE (abs(t)>1.0e-20_r8)
        t=-t*ax/k; k=k+1; t=t/k; y=y+t; k=k+1 
      ENDDO
    ENDIF
    res1=y
    END FUNCTION aux1
                                            
    RECURSIVE FUNCTION aux2(x) RESULT(res2)
    USE Someconstants
    ! sinh(x) - x = (x^3/6) aux2(x)
    IMPLICIT NONE
    REAL(r8) :: x, ax, y, e, t, res2
    INTEGER k
    ax=abs(x)
    IF (ax > 1.5) THEN
      e=exp(x); 
      y=6.0_r8*((e-1.0_r8/e)*0.5_r8 - x)/(x*(x*x))
    ELSEIF (ax> 0.5) THEN
      e=aux2(x/3.0_r8); 
      t=2.0_r8+x*x*e/27.0_r8; 
      y=(e+t*(t*t))/9.0_r8
    ELSEIF (ax<1.0e-10_r8) THEN
      y=1 
    ELSE
      ax=x*x; t=ax/20.0_r8; y=1.0_r8+t; k=6
      DO WHILE (t > 1.0e-20_r8)
        t=t*ax/k; k=k+1; t=t/k; y=y+t; k=k+1 
      ENDDO
    ENDIF
    res2=y
    END FUNCTION aux2
                                    
    FUNCTION aux3(x)
    USE Someconstants
    ! 1 - x cot x = (x^2/3) aux3(x)
    IMPLICIT NONE
    REAL(r8) :: aux3, x, absx, xh, y
    absx=abs(x)
    IF (absx > 1) THEN
      y=3.0_r8*(1.0_r8-x*cos(x)/sin(x))/(x*x) 
    ELSEIF (absx < 1.0e-10_r8) THEN
      y=1.0_r8 
    ELSE
      xh=x*0.5_r8; y=sin(xh)/x
      y= 3.0_r8*(y-aux1(x)/(12.0_r8*y))/cos(xh)       
    ENDIF
    aux3=y
    END FUNCTION aux3
                                    
    FUNCTION aux4(x)
    USE Someconstants
    ! e^x - 1 = x aux4(x)
    IMPLICIT NONE
    REAL(r8) :: aux4, x, xh, y
    IF ((x < -0.7).OR.(x > 0.4)) THEN
      y=(exp(x)-1.0_r8)/x 
    ELSEIF (abs(x) < 1.0e-10_r8) THEN
      y=1.0_r8+x*0.5_r8
    ELSE
      xh=x*0.5_r8; y=exp(xh)*(1.0_r8+x*x*aux2(xh)/24.0_r8)
    ENDIF
    aux4=y
    END FUNCTION aux4
                                
    FUNCTION aux5(x)
    USE Someconstants
    ! ln(1+x) = x aux5(x)
    IMPLICIT NONE
    REAL(r8) :: aux5, x, r,  x1, y, y0
    x1=1.0_r8+x;
    IF (x1 <= 0) THEN
      y=log(dwarf) 
    ELSEIF ((x < -0.7).OR.(x > 1.36)) THEN
      y=log(x1)/x 
    ELSEIF (abs(x) < 1.0e-10_r8) THEN
      y=1.0_r8-x*0.5_r8 
    ELSE
      y0=log(x1); r=-x+x1*y0*aux4(-y0) 
      y=(y0-r)/x
    ENDIF
    aux5=y
    END FUNCTION aux5
                                      
    FUNCTION aux6(x)
    USE Someconstants
    ! ln(1+x)-x +1/2 x^2 = 1/3 x^3 aux6(x); expansion 4.1.29 of A&S, a=1
    IMPLICIT NONE
    REAL(r8) :: aux6, x, x1, y, r, r2, rk, tk
    INTEGER k
    x1=1+x
    IF (x1 <= 0) THEN
      y=log(dwarf) 
    ELSEIF ((x < -0.66).OR.(x > 2)) THEN 
      r=x*x; y=3.0_r8*(log(x1)-x+r*0.5_r8)/(x*r) 
    ELSEIF (abs(x) < 1.0e-10_r8) THEN
      y=1.0_r8-0.75_r8*x 
    ELSE
      r=x/(2.0_r8+x); r2=r*r; rk=1
      y=1.0_r8/5.0_r8; k=5; tk=1 
      DO WHILE (tk > mactol)
        k=k+2; rk=rk*r2; tk= rk/k; y= y+tk
      ENDDO    
      r=x+2.0_r8; r2=r*r; y=(2.0_r8+1.5_r8*r2+6*x*x*y/r2)/(r*r2)    
    ENDIF
    aux6=y
    END FUNCTION aux6
                                           
    FUNCTION aux7(x) 
    USE Someconstants
    ! ln(1+x)-x = -1/2 x^2  aux7(x)
    IMPLICIT NONE
    REAL(r8) :: aux7, x, x1, y 
    x1=1.0_r8+x
    IF (x1 <= 0) THEN
      y=log(dwarf) 
    ELSEIF ((x < -0.66).OR.(x > 2)) THEN
      y=-2.0_r8*(log(x1)-x)/(x*x)  
    ELSEIF (abs(x)<1.0e-10_r8) THEN
      y=1.0_r8-2.0_r8*x/3.0_r8 
    ELSE
      y=1.0_r8-2.0_r8*x*aux6(x)/3.0_r8
    ENDIF
    aux7=y
    END FUNCTION aux7
                                    
    FUNCTION aux8(x) 
    ! arctan(x)-x = -1/3 x^3  aux8(x) 
    IMPLICIT NONE
    REAL(r8) :: aux8, x, y, x2, u, w
    IF (abs(x)<1.0e-10_r8) THEN
      aux8=1.0_r8
    ELSE
      y=atan(x)
      IF (abs(x)>1) THEN 
        aux8=-3.0_r8*(y-x)/(x*x*x) 
      ELSE
        x2=x*x; w=sqrt(1.0_r8+x2); u=y/x
        aux8=u*(3.0_r8/(1.0_r8+w)-w*u*u*aux1(y)*0.5_r8)
      ENDIF
    ENDIF
    END FUNCTION aux8
                                 
    FUNCTION aux9(a) 
    USE Someconstants
    ! Gamma(a+1/2)*exp(a-a*ln(a))/(sqrt(2*pi) = aux9(a) 
    IMPLICIT NONE
    REAL(r8) :: aux9, a, z, z2, z2k, y, t,ck(0:10)
    INTEGER k
    IF (a == 0) THEN
      aux9=1.0_r8/sqrt2 
    ELSEIF (a == 1) THEN
      aux9=exp(1.0_r8)*0.5_r8/sqrt2
    ELSEIF (a < 8) THEN
      aux9=gamma(a+0.5_r8)*exp(a-a*log(a))/sqrttwopi 
    ELSE
      ck(0)=.0833333333333333333333333333333_r8
      ck(1)=-.00277777777777777777777777777778_r8
      ck(2)=.000793650793650793650793650793651_r8
      ck(3)=-.000595238095238095238095238095238_r8
      ck(4)=.000841750841750841750841750841751_r8
      ck(5)=-.00191752691752691752691752691753_r8
      ck(6)=.00641025641025641025641025641026_r8
      ck(7)= -.0295506535947712418300653594771_r8
      ck(8)=0.179644372368830573164938490016_r8
      ck(9)=-1.39243221690590111642743221691_r8
      ck(10)=13.4028640441683919944789510007_r8
      z=1.0_r8/(a+0.5_r8); z2= z*z; z2k= z2
      y=ck(0)+ck(1)*z2; k= 2; t=1.0_r8
      DO WHILE ((k < 11).AND.(abs(t) > 1.0e-20_r8))
        z2k=z2k*z2; t=ck(k)*z2k; y= y+t; k= k+1 
      ENDDO
      aux9=exp(y*z-aux7(0.5_r8/a)/(8.0_r8*a))
    ENDIF
    END FUNCTION aux9
                                            
    SUBROUTINE aux10(x, y, u, v, w) 
    USE Someconstants
    ! z=x+iy;  ln(1+z)-z =aux10(x, y) = u+ iv; 
    ! w = arctan(eta)-eta; eta = y/(1+x) 
    IMPLICIT NONE
    REAL(r8) :: x, y, u, v, w,x1, x2, y2, r2, xi, eta
    x2= x*x; y2= y*y; r2= x2+y2; x1= x+1.0_r8
    IF ((x1 == 0).AND.(y == 0)) THEN
      u=log(dwarf); v=0; w=0 
    ELSEIF ((abs(x) < mactol).AND.(abs(y) < mactol)) THEN 
      u=-(x2-y2)*0.5_r8; v=-x*y; w=-y2*y/3.0_r8  
    ELSE
      IF (r2 > 0.5) THEN
        IF (x1 == 0) THEN 
          u=-x + log(abs(y))
          IF (y > 0) THEN
            v=pihalf 
          ELSE
            v=-pihalf
          ENDIF
          w= v-y
        ELSE
          xi=1.0_r8+x+x+r2; eta= y/x1 
          IF (xi>0) THEN
            u= 0.5_r8*log(xi)-x 
          ELSE
            u= log(dwarf)
          ENDIF  
          IF (abs(eta) > 1) THEN
            v=phase(x1,y)-y; 
            w=atan(eta)-eta 
          ELSE
            w=-eta*eta*eta*aux8(eta)/3.0_r8; 
            v= w-eta*x
          ENDIF
        ENDIF
      ELSE
        xi=r2+x+x; eta=y/x1
        u=0.5*(r2-0.5_r8*xi*xi*aux7(xi))
        w=-eta*eta*eta*aux8(eta)/3.0_r8; v= w-eta*x
      ENDIF
    ENDIF
    END SUBROUTINE aux10
                                          
    FUNCTION hypergeom(a,c,x,eps)
    ! 1F1(a,c,x)
    IMPLICIT NONE
    REAL(r8) :: hypergeom, a, c, x, eps, del, s, t
    INTEGER k
    s= 1.0_r8; t= 1.0_r8; k= 0; del= 1.0_r8 
    DO WHILE ((del > eps).AND.(k<1000)) 
      t=t*x*(a+k)/(c+k); k= k+1; t= t/k
      s=s+t
      IF (s == 0) THEN
        del=1.0_r8 
      ELSE
        del=abs(t/s)
      ENDIF
    ENDDO
    hypergeom=s
    END FUNCTION hypergeom
                                        
    FUNCTION uaxhyp(a, z, eps)
    USE Someconstants
    USE AiryFunction, ONLY: xpowy
    ! U(a,x) by  1F1
    IMPLICIT NONE
    REAL(r8) :: uaxhyp, a, z, eps, s, t, y1, y2, u0, u0p
    s=0.25_r8+a*0.5_r8; t=s+0.5_r8 
    IF (t > 0) THEN 
      u0=sqrt(pi)*xpowy(2.0_r8,-a*0.5_r8-0.25_r8)/gamma(t) 
    ELSE
      u0=sin(pi*t)*gamma(1.0_r8-t)*xpowy(2.0_r8,-a*0.5_r8-0.25_r8)/sqrt(pi)
    ENDIF
    IF (s > 0) THEN 
      u0p=-sqrt(pi)*xpowy(2.0_r8,-a*0.5_r8+0.25_r8)/gamma(s) 
    ELSE
      u0p=-sin(pi*s)*gamma(1.0_r8-s)*xpowy(2.0_r8,-a*0.5_r8+0.25_r8)/sqrt(pi)
    ENDIF
    t= z*z*0.5_r8; s= exp(-t*0.5_r8)
    y1=s*hypergeom(0.5_r8*a+0.25_r8, 0.5_r8, t, eps)
    y2=z*s*hypergeom(0.5_r8*a+0.75_r8, 1.5_r8, t, eps)
    uaxhyp=u0*y1+u0p*y2
    END FUNCTION uaxhyp

    FUNCTION vaxhyp(a, z, eps)
    USE Someconstants
    USE AiryFunction, ONLY: xpowy
    ! V(a,x) by  1F1
    IMPLICIT NONE
    REAL(r8) :: vaxhyp, a, z, eps, s, t, y1, y2, v0, v0p
    s= 0.25_r8-a*0.5_r8; t= s+0.5_r8
    y1=sin(pi*t)
    y2=sin(pi*s)
    IF (y1 == 0) THEN
      v0=0.0_r8 
    ELSE
      IF (t > 0) THEN  
        v0=y1*xpowy(2.0_r8,a*0.5_r8+0.25_r8)/gamma(t)  
      ELSE
        v0=y1*y1*xpowy(2.0_r8,a*0.5_r8+0.25_r8)*gamma(1.0_r8-t)/pi
      ENDIF
    ENDIF
    IF (y2 == 0) THEN 
      v0p=0.0_r8 
    ELSE
      IF (s > 0) THEN 
        v0p= y2*xpowy(2.0_r8,a*0.5_r8+3.0_r8/4.0_r8)/gamma(s) 
      ELSE
        v0p= y2*y2*xpowy(2.0_r8,a*0.5_r8+3.0_r8/4.0_r8)*gamma(1.0_r8-s)/pi
      ENDIF
    ENDIF 
    t= z*z*0.5_r8; s= exp(-t*0.5_r8)
    y1= s*hypergeom(0.5_r8*a+0.25_r8, 0.5_r8, t, eps)
    y2= z*s*hypergeom(0.5_r8*a+0.75_r8, 1.5_r8, t, eps)
    vaxhyp= v0*y1+v0p*y2
    END FUNCTION vaxhyp
                                    
    FUNCTION uaxd(a,x,eps)
    ! U'(a,x) by  uaxhyp
    IMPLICIT NONE
    REAL(r8) :: uaxd, a, x, eps 
    IF (a > 0) THEN 
      uaxd= -x*0.5_r8*uaxhyp(a,x,eps)-(a+0.5_r8)*uaxhyp(a+1.0_r8,x,eps) 
    ELSE
      uaxd= x*0.5_r8*uaxhyp(a,x,eps)-uaxhyp(a-1.0_r8,x,eps)
    ENDIF
    END FUNCTION uaxd
                                     
    FUNCTION vaxd(a,x,eps)
    ! V'(a,x) by  vaxhyp
    IMPLICIT NONE
    REAL(r8) :: vaxd, a, x, eps 
    IF (a < 0) THEN
      vaxd=x*0.5_r8*vaxhyp(a,x,eps)+(a-0.5_r8)*vaxhyp(a-1.0_r8,x,eps) 
    ELSE
      vaxd=-x*0.5_r8*vaxhyp(a,x,eps)+vaxhyp(a+1.0_r8,x,eps)
    ENDIF
    END FUNCTION vaxd
                                                
    RECURSIVE FUNCTION phase(x, y) RESULT(phas)
    USE Someconstants
    ! The phase of the complex number z =  x + iy
    IMPLICIT NONE
    REAL(r8) :: x, y, phas
    IF ((x == 0).AND.(y == 0)) THEN
      phas= 0.0_r8 
    ELSEIF (y < 0) THEN
      phas=-phase(x,-y) 
    ELSEIF (x >= y) THEN
      phas= atan(y/x) 
    ELSEIF (x + y >= 0) THEN
      phas= pihalf - atan(x/y)  
    ELSE 
      phas= pi + atan(y/x)
    ENDIF  
    END FUNCTION phase
                                                 
    RECURSIVE FUNCTION errorfunction (x, erfc, expo) RESULT(errfu)
    ! coefficients are from Cody (1969), Math. Comp., 23, 631-637
    IMPLICIT NONE
    REAL(r8) :: x, y, z, r(0:8), s(0:8), errfu
    REAL(r8), PARAMETER :: oneoversqrtpi = 5.641895835477562869e-1_r8
    LOGICAL erfc, expo
    IF (erfc) THEN
      IF (x < -6.5) THEN
        y= 2.0_r8 
      ELSEIF (x < 0) THEN
        y= 2.0_r8 - errorfunction(-x, .true., .false.) 
      ELSEIF (x == 0) THEN
        y= 1.0_r8 
      ELSEIF (x < 0.5) THEN
        IF (expo) THEN
          y=exp(x*x)
        ELSE
          y=1.0_r8
        ENDIF
        y=y*(1.0_r8 - errorfunction(x, .false., .false.))
      ELSEIF (x < 4) THEN
        IF (expo) THEN
          y= 1.0_r8 
        ELSE
          y= exp(-x*x)
        ENDIF
        r(0)= 1.230339354797997253e3_r8
        r(1)= 2.051078377826071465e3_r8
        r(2)= 1.712047612634070583e3_r8
        r(3)= 8.819522212417690904e2_r8
        r(4)= 2.986351381974001311e2_r8
        r(5)= 6.611919063714162948e1_r8
        r(6)= 8.883149794388375941_r8
        r(7)= 5.641884969886700892e-1_r8
        r(8)= 2.153115354744038463e-8_r8
        s(0)= 1.230339354803749420e3_r8
        s(1)= 3.439367674143721637e3_r8
        s(2)= 4.362619090143247158e3_r8
        s(3)= 3.290799235733459627e3_r8
        s(4)= 1.621389574566690189e3_r8
        s(5)= 5.371811018620098575e2_r8
        s(6)= 1.176939508913124993e2_r8
        s(7)= 1.574492611070983473e1_r8
        y=y*fractio(x,8,r,s)
      ELSE
        z=x*x
        IF (expo) THEN
          y=1.0_r8 
        ELSE
          y= exp(-z)
        ENDIF
        z=1.0_r8/z
        r(0)=6.587491615298378032e-4_r8
        r(1)=1.608378514874227663e-2_r8
        r(2)=1.257817261112292462e-1_r8
        r(3)=3.603448999498044394e-1_r8
        r(4)=3.053266349612323440e-1_r8
        r(5)=1.631538713730209785e-2_r8
        s(0)=2.335204976268691854e-3_r8
        s(1)=6.051834131244131912e-2_r8
        s(2)=5.279051029514284122e-1_r8
        s(3)=1.872952849923460472_r8
        s(4)=2.568520192289822421_r8
        y=y*(oneoversqrtpi-z*fractio(z,5,r,s))/x
      ENDIF
      errfu=y
    ELSE
      IF (x == 0.0_r8) THEN 
        y=0
      ELSEIF (abs(x) > 6.5) THEN 
        y=x/abs(x)
      ELSEIF (x > 0.5) THEN
        y=1.0_r8 - errorfunction(x, .true., .false.) 
      ELSEIF (x < -0.5) THEN
        y=errorfunction(-x, .true., .false.)-1.0_r8
      ELSE
        r(0)=3.209377589138469473e3_r8
        r(1)=3.774852376853020208e2_r8
        r(2)=1.138641541510501556e2_r8
        r(3)=3.161123743870565597e0_r8
        r(4)=1.857777061846031527e-1_r8
        s(0)=2.844236833439170622e3_r8
        s(1)=1.282616526077372276e3_r8
        s(2)=2.440246379344441733e2_r8
        s(3)=2.360129095234412093e1_r8
        z=x*x
        y=x*fractio(z,4,r,s)
      ENDIF  
      errfu= y
    ENDIF        
    END FUNCTION errorfunction
                                   
    FUNCTION fractio(x,n,r,s)
    IMPLICIT NONE
    INTEGER n,k
    REAL(r8) :: x, fractio, r(0:8), s(0:8), a, b
    a=r(n); b=1
    DO k=n-1,0,-1 
      a=a*x+r(k); b=b*x+s(k) 
    ENDDO
    fractio=a/b
    END FUNCTION fractio
                                       
    FUNCTION gamma(x) 
    !---------------------------------------------------------------------
    ! This routine calculates the gamma function for real arguments x.
    ! This is a Fortran 90 version of the ACM TOMS routine for the
    ! gamma function by W.J. Cody which is included in:
    !   Algorithm 715: SPECFUN-- A Portable FORTRAN Package of Special 
    !                  Function Routines and Test Drivers. 
    !   ACM Trans. Math. Soft. 19(1) 22-30.  
    ! --------------------------------------------------------------------
    ! Explanation of machine-dependent constants.  Let
    !  beta   - radix for the floating-point representation
    !  maxexp - the smallest positive power of beta that overflows.
    ! Then the following machine-dependent constants must be declared
    ! in DATA statements:
    ! XBIG   - the largest argument for which GAMMA(X) is representable
    !          in the machine, i.e., the solution to the equation
    !                  GAMMA(XBIG) = beta**maxexp
    ! XINF   - the largest machine representable floating-point number;
    !          approximately beta**maxexp
    ! EPS    - the smallest positive floating-point number such that
    !          1.0+EPS .GT. 1.0
    ! XMININ - the smallest positive floating-point number such that
    !          1/XMININ is machine representable
    !
    !  We use the Fortran 90 intrinsic functions TINY, EPSILON, HUGE
    !  for the computation of XMININ, EPS and XINF, respectively.
    !  
    !   Approximate values for XBIG some important machines are:
    !                            beta       maxexp        XBIG
    ! CRAY-1         (S.P.)        2         8191        966.961
    ! Cyber 180/855
    !   under NOS    (S.P.)        2         1070        177.803
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (S.P.)        2          128        35.040
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (D.P.)        2         1024        171.624
    ! IBM 3033       (D.P.)       16           63        57.574
    ! VAX D-Format   (D.P.)        2          127        34.844
    ! VAX G-Format   (D.P.)        2         1023        171.489
    !--------------------------------------------------------------------
    ! Error returns
    !  The program returns the value XINF for singularities or
    !     when overflow would occur.  The computation is believed
    !     to be free of underflow and overflow.
    !  Intrinsic functions required are:
    !     INT, DBLE, EXP, LOG, REAL, SIN
    !
    ! References: "An Overview of Software Development for Special
    !              Functions", W. J. Cody, Lecture Notes in Mathematics,
    !              506, Numerical Analysis Dundee, 1975, G. A. Watson
    !              (ed.), Springer Verlag, Berlin, 1976.
    !              Computer Approximations, Hart, Et. Al., Wiley and
    !              sons, New York, 1968.
    !----------------------------------------------------------------------
    IMPLICIT NONE
    REAL (r8), INTENT(IN)  :: x
    REAL (r8)              :: gamma, fn_val
    ! Local variables
    INTEGER    :: i, n
    LOGICAL    :: parity
    REAL (r8)  :: fact, sum, xden, xnum, y, y1, ysq, z
    !-----------------------------------------------------
    !  Mathematical constants
    !-----------------------------------------------------
    REAL (r8), PARAMETER  :: one=1.0_r8, half=0.5_r8, twelve=12.0_r8,  &
                             two=2.0_r8, zero=0.0_r8,  &
                             sqrtpi=0.9189385332046727417803297_r8,  &
                             pi=3.1415926535897932384626434_r8
    !------------------------------------------------------
    !  Machine dependent parameters
    !------------------------------------------------------
    REAL (r8), PARAMETER  :: xbig = 171.624_r8, xminin=TINY(0.0_r8),   &
                             eps = EPSILON(0.0_r8), xinf =HUGE(0.0_r8)
    !----------------------------------------------------------------------
    !  Numerator and denominator coefficients for rational minimax
    !     approximation over (1,2).
    !----------------------------------------------------------------------
    REAL (r8), PARAMETER  :: P(8) =  &
           (/ -1.71618513886549492533811e+0_r8,  2.47656508055759199108314e+1_r8,  &
              -3.79804256470945635097577e+2_r8,  6.29331155312818442661052e+2_r8,  &
               8.66966202790413211295064e+2_r8, -3.14512729688483675254357e+4_r8,  &
              -3.61444134186911729807069e+4_r8,  6.64561438202405440627855e+4_r8 /)
    REAL (r8), PARAMETER  :: Q(8) =  &
           (/ -3.08402300119738975254353e+1_r8,  3.15350626979604161529144e+2_r8,  &
              -1.01515636749021914166146e+3_r8, -3.10777167157231109440444e+3_r8,  &
               2.25381184209801510330112e+4_r8,  4.75584627752788110767815e+3_r8,  &
              -1.34659959864969306392456e+5_r8, -1.15132259675553483497211e+5_r8 /)
    !------------------------------------------------------------------
    !  Coefficients for minimax approximation over (12, INF).
    !------------------------------------------------------------------
    REAL (r8), PARAMETER  :: c(7) =  &
              (/ -1.910444077728e-03_r8, 8.4171387781295e-04_r8,  &
              -5.952379913043012e-04_r8, 7.93650793500350248e-04_r8,  &
              -2.777777777777681622553e-03_r8, 8.333333333333333331554247e-02_r8,  &
               5.7083835261e-03_r8 /)
    INTEGER k
    !------------------------------------------------------------------
    parity= .false.
    fact=one
    n=0
    y=x
    k=0
    IF (y <= zero) THEN
   !---------------------------------------------------
   !  Argument is negative
   !---------------------------------------------------
      y=-x
      y1=AINT(y)
      fn_val=y-y1
      IF (fn_val /= zero) THEN
        IF (y1 /=aint(y1*half)*two) parity = .true.
        fact=-pi/sin(pi*fn_val)
        y=y+one
      ELSE
        fn_val=xinf
        k=1
      ENDIF
    ENDIF
    IF (k==0) THEN
    !------------------------------------------------
    !  Argument is positive
    !------------------------------------------------
       IF (y < eps) THEN
       !----------------------------------------------
       !  Argument < EPS
       !----------------------------------------------
         IF (y >= xminin) THEN
           fn_val=one/y
         ELSE
           fn_val=xinf
           k=1 
         ENDIF
       ELSEIF (y < twelve) THEN
         y1=y
         IF (y < one) THEN
         !---------------------------------------------
         !  0.0 < argument < 1.0
         !---------------------------------------------
           z=y
           y=y+one
         ELSE
         !-------------------------------------------------------
         !  1.0 < argument < 12.0, reduce argument if necessary
         !-------------------------------------------------------
           n=INT(y)-1
           y=y-n
           z=y-one
         ENDIF
         !-------------------------------------------------------
         !  Evaluate approximation for 1.0 < argument < 2.0
         !-------------------------------------------------------
         xnum=zero
         xden=one
         DO i=1,8
           xnum=(xnum+p(i))*z
           xden=xden*z+q(i)
         END DO
         fn_val=xnum/xden+one
         IF (y1 < y) THEN
         !-------------------------------------------------------
         !  Adjust result for case  0.0 < argument < 1.0
         !-------------------------------------------------------
           fn_val=fn_val/y1
         ELSEIF (y1 > y) THEN
         !---------------------------------------------------------
         !  Adjust result for case  2.0 < argument < 12.0
         !---------------------------------------------------------
           DO i=1,n
             fn_val=fn_val*y
             y=y+one
           ENDDO
         ENDIF
       ELSE
         !----------------------------------------------------
         !  Evaluate for argument .GE. 12.0,
         !----------------------------------------------------
         IF (y <= xbig) THEN
           ysq=y*y
           sum=c(7)
           DO i=1,6
             sum=sum/ysq+c(i)
           ENDDO
           sum=sum/y-y+sqrtpi
           sum=sum+(y-half)*log(y)
           fn_val=EXP(sum)
         ELSE
           fn_val=xinf
           k=1
         ENDIF
       ENDIF
    ENDIF
    !--------------------------------------------------------------
    !  Final adjustments and return
    !--------------------------------------------------------------
    IF (k==0) THEN
      IF (parity) fn_val=-fn_val
      IF (fact /= one) fn_val=fact/fn_val
    ENDIF
    gamma=fn_val
    END FUNCTION gamma
                                         
    SUBROUTINE gh12(aneg, x, eps, gh)
    USE Someconstants
    ! --------------------------------------------------------------
    ! For 0 <= x <= 2*sqrt(-a), a < 0; computes I_j, J_j; 
    ! see GST(2004) on PCF;
    ! the value of gh12(a,x,..) is the Wronskian relation 
    ! normalized to zero
    ! ---------------------------------------------------------------
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: aneg, x, eps
    REAL(r8), INTENT(OUT) :: gh(4)
    REAL(r8)  mu, a, t, teta, up, up2, vp, lambda, sigma, gammastar,&
              errorf, delta, eps10, h, lneps10, a0, b0, p0, p1, p2, p3,&
              p4, q1, q2, psir0, psir1, psir2, psii0, psii1, psii2, r0,&
              rho0, s0, s1, h2, sum(4), newsum(4)
    INTEGER nc
    LOGICAL tnearone
    a=-aneg; mu=sqrt(a); t=x/(2.0_r8*mu); teta=asin(t)
    up2=(1.0_r8-t)*(1.0_r8+t); up= sqrt(up2); vp= t;
    sigma=sqrt(2.0_r8/a); tnearone=(t >= 0.9)
    eps10=eps*0.1_r8; lneps10= 1.1_r8*log(eps10)
    r0=sqrt(-lneps10) 
    rho0=r0; a0=-rho0; b0=rho0
    gammastar=aux9(a)
    IF (tnearone) THEN
      q1=(1.0_r8+up)/vp; q2= 1.0_r8/(vp*vp) 
      psir0=up*(1.0_r8+up); psir1= q1; psir2= 0.5_r8*q2
      psii0=1.0_r8+up2; psii1=(up+up2+1.0_r8)/vp; psii2= up*q2
      r0=vp*sqrt(sqrt(up*up*(1.0_r8+up)*(1.0_r8+up)-2.0_r8*lneps10/a)-up*(1.0_r8+up))
      r0=r0/sigma
      p0=2.0_r8*vp*r0
      s1=sigma*r0; s0=(vp+s1)/p0; s1=(vp-s1)/p0 
    ELSE
      p0=up*vp; p1= up+up2; p2= p0*p1; p3= p0/p1 
      q1=2.0_r8*up*r0
      s1=sigma*r0; s0=(up+s1)/q1; s1=(up-s1)/q1 
    ENDIF
    lambda= 4.0_r8*a*sigma*sigma/(pi*pi*gammastar)  
    ! begin {main of trap}
    h = b0/32.0_r8; errorf=0.0_r8
    CALL intrho12(0.0_r8,errorf, sum) 
    sum=h*sum
    CALL summing12(h,h,newsum) 
    sum= sum + newsum
    nc=1
    delta=1.0_r8
    DO WHILE (((delta > eps).AND.(nc<10)).OR.(nc<5))
      nc= nc+1; h2= h*0.5_r8
      CALL summing12(h, h2,newsum)  
      sum=(sum + newsum)*0.5_r8
      delta=abs(lambda*(sum(3)*sum(2)-sum(1)*sum(4))-1.0_r8)
      h= h2
    ENDDO
    ! end (trap)
    gh= sum
    p0=2.0_r8*sigma/sqrtpi 
    gh=p0*gh
    ! G1=gh(1), G2=gh(2), H1=gh(3), H2=gh(4)
    CONTAINS
      SUBROUTINE summing12(h,d,newsum)
      IMPLICIT NONE
      REAL(r8) h, d, newsum(4), fr(4), errorf, r 
      newsum=0.0_r8
      r=d 
      DO WHILE ((r<b0).OR.(r<-a0)) 
        errorf= errorfunction(r, .false., .false.) ! {errorf = erf(r)}
        IF (r < b0) THEN 
          CALL intrho12(r, errorf, fr)
          newsum= newsum+fr
        ENDIF
        errorf= -errorf
        IF (r < -a0) THEN
          CALL intrho12(-r, errorf, fr)
          newsum=newsum+fr
        ENDIF
        r=r+h 
      ENDDO 
      newsum=h*newsum
      END SUBROUTINE summing12
      SUBROUTINE intrho12(rho, errorf, frho)
      IMPLICIT NONE
      REAL(r8) rho, errorf, frho(4), s, r, drds, p, q,&
               psir, psii, dpdq, dqdp, kappa, tau, alpha,&
               beta, gamma, chi, g1, g2, h1, h2 
      s= errorf 
      IF (((abs(s) > 1 -eps10).OR.(rho <a0)).OR.(rho > b0)) THEN 
        frho= 0.0_r8
      ELSE 
        r= 1.0_r8/(s0+s1*s); drds= s0*( r*r); r= s*r 
        IF (tnearone) THEN
          q= sigma*r; p= q*(q1+q2*q); dpdq= q1+2*q2*q; r= q*q
          kappa= p*up+q*vp; tau= q*up-p*vp
          CALL aux10(kappa, tau, alpha, beta, gamma) 
          psir= q2*r*(psir0+q*(psir1+psir2*q))-alpha
          psii= q2*q*r*(psii0+q*(psii1+psii2*q))/(1.0_r8+kappa)-gamma
        ELSE
          p= sigma*r 
          p4= 1.0_r8/(p1+p); dqdp= p2*p4*p4; p4= -p3*p4; q= p*(p3+p*p4) 
          kappa= p*up+q*vp; tau= q*up-p*vp
          q1= 1.0_r8+vp*p4*p 
          r= -q1*q1*(p*p4*p1-vp)/((1.0_r8+up)*(1.0_r8+kappa))
          CALL aux10(kappa, tau, alpha, beta, gamma) 
          psir=(p*p-q*q)*0.5_r8-alpha 
          psii=p*p*p*(r+p4*(2*up+p*p0*p4))-gamma
        ENDIF
        r=rho*rho+a*psir+(kappa+alpha)*0.5_r8
        IF (r > -lneps10) THEN 
          frho=0.0_r8
          IF (rho>0) THEN
            b0=rho 
          ELSE
            a0=rho
          ENDIF
        ELSE
          r= exp(-r)*drds; chi= a*psii+(teta+tau+beta)*0.5_r8
          alpha=cos(chi); beta=sin(chi)
          IF (tnearone) THEN
            g1=alpha*dpdq+beta; g2=beta*dpdq-alpha 
          ELSE
            g1=alpha+beta*dqdp; g2=beta-alpha*dqdp 
          ENDIF
          h1=-g1*q+g2*(p+up); h2=-q*g2-g1*(p+up)
          frho(1)=g1; frho(2)=g2; frho(3)=h1; frho(4)=h2 
          frho= r*frho
        ENDIF
      ENDIF
      END SUBROUTINE intrho12
      END SUBROUTINE gh12  
                                              
      SUBROUTINE gh123(aneg, x, eps, ghuv)
      USE Someconstants
      ! ----------------------------------------------------
      ! For x >= 0, aneg < 0; computes integrals I, J 
      ! and derivatives for t > 1. See GST(2004) on PCFs.
      ! ----------------------------------------------------
      IMPLICIT NONE
      REAL(r8), INTENT(IN) :: aneg, x, eps
      REAL(r8), INTENT(OUT):: ghuv(6)
      REAL(r8) :: a, t, mu, a2xi, ep, em, em2, eps10, lneps10, sigmau, &
                  sigmav, p0, q1, q2, a0, b0, cu, cv, rho0, vp, vm, sq,&
                  vpsq, vpsq2, vmsq, errorf, delta, gammastar, wronskian,&
                  sum(6), newsum(6), v, p, q, h, h2
      INTEGER nc
      ! begin {main of gh123}
       a= -aneg; mu=sqrt(a); t= 0.5_r8*x/mu; sq=sqrt((t-1.0_r8)*(t+1.0_r8))
       vp= t+sq; vm= 1.0_r8/vp
       vpsq= vp*sq; vpsq2= 2.0_r8*vpsq; vmsq= vm*sq
       sigmau= sqrt(2.0_r8/a); sigmav= 1.0_r8/x   ! sigmau= 1
       eps10=eps*0.1_r8; lneps10=1.1_r8*log(eps10) 
       a0=vm*vm; b0= -2.0_r8*a0*lneps10/a
       b0=b0-(b0-a0*log(1.0_r8+b0)-b0)*(1.0_r8+b0)/(1.0_r8+b0-a0)
       p0=sqrt(b0)/(sigmau*vm)
       rho0=sqrt(-lneps10)
       a0=-rho0; b0=rho0
       q1=2.0_r8*sq/sigmav; q2=vm/sigmav
       cu=1.0_r8/sqrtpi*p0*sigmau/sqrt(vp)
       cv=1.0_r8/sqrtpi*sigmav
       gammastar=aux9(a); wronskian=pi*gammastar/a
       a2xi= a*(t*sq-log(vp))
       IF (a2xi>log(giant)) THEN
         ep=giant
         em=0.0_r8; em2=0.0_r8   
       ELSE
         ep=exp(a2xi); 
         em=1.0_r8/ep; em2=em*em
       ENDIF
       !  begin trap(ghuv);
       h=b0/16.0_r8;
       errorf=0.5_r8; 
       CALL intrho123(0.0_r8,errorf,sum) 
       sum= h*sum
       q=10.0_r8; nc=1
       CALL summing123(h,h,newsum) 
       sum=sum+newsum
       delta=1.0_r8
       DO WHILE (((delta > eps).AND.(nc <= 5)).OR.(nc<4))
         nc=nc+1; h2=h*0.5_r8
         CALL summing123(h, h2,newsum) 
         sum=(sum+newsum)*0.5_r8
         p=sum(1)
         delta=abs(em2*(sum(1)*sum(4)-sum(2)*sum(3))&
               +sum(1)*sum(6)-sum(3)*sum(5)-wronskian)
         h=h2;v=abs(p-q);q=p
       ENDDO 
       ghuv=sum
     CONTAINS
      SUBROUTINE summing123(h,d,newsum)
      IMPLICIT NONE
      REAL(r8) h, d, newsum(6), fr(6), rho, errorf
      newsum=0.0_r8
      rho=d 
      DO WHILE ((rho<b0).OR.(rho<-a0)) 
        errorf= 0.5_r8*errorfunction(rho, .true., .false.) !  {0.5*erfc(rho)}
        IF (rho < -a0) THEN 
          CALL intrho123(-rho, errorf, fr)
          newsum= newsum+fr
        ENDIF
        errorf=1.0_r8-errorf
        IF (rho < b0) THEN 
          CALL intrho123(rho, errorf, fr)
          newsum=newsum+fr
        ENDIF
        rho=rho+h 
      ENDDO 
      newsum=h*newsum
      END SUBROUTINE summing123
                                           
      SUBROUTINE intrho123(rho,errorf,fuv)
      IMPLICIT NONE
      REAL(r8) rho, errorf, fuv(6), rho2, r,&
               s, fu(4), fv(6)
      INTEGER j
      rho2=rho*rho
      r=errorf
      CALL intu(r,rho2,fu)
      CALL intv(r,rho2,fv)
      DO j= 1,4 
        fuv(j)=fu(j)
      ENDDO
      DO j= 5,6 
        fuv(j)=fv(j)
      ENDDO
      s=0.0_r8
      DO j= 1,6 
        s=s+abs(fuv(j))
      ENDDO
      IF (s < eps10) THEN
        IF (rho > 0) THEN
          b0=rho 
        ELSE
          a0=rho
        ENDIF
      ENDIF
      END SUBROUTINE intrho123
                                       
      SUBROUTINE intu(r,rho2,fu)
      IMPLICIT NONE
      REAL(r8) r, rho2, fu(4), u, p, q, qq, c, s, s2, chi,&
               psir, psii, gh(4)  
      p=r*p0 
      u=sigmau*p 
      q=u*vm; qq=q*q 
      s=p*vm; s2=s*s
      psir=rho2+s2*(vpsq2+s2*aux7(qq)/a) !psir= rho*rho+a*(u*u/2-0.5*ln(1+qq));
      IF (psir > -lneps10) THEN 
        fu= 0.0_r8
      ELSE
        s=-q*qq*aux8(q)/3.0_r8 !s = arctan(q)-q
        psii=a*s
        chi=-psii+(s+q)*0.5_r8
        c=cos(chi); s=sin(chi)
        gh(1)=c; gh(2)=s
        gh(3)=-sq*c-u*s; gh(4)=u*c -sq*s 
        q=cu*exp(-psir)/sqrt(sqrt(1.0_r8+qq))
        fu=q*gh
      ENDIF
      END SUBROUTINE intu
                                                  
      SUBROUTINE intv(r,rho2, fv)
      IMPLICIT NONE
      REAL(r8) r, rho2, fv(6), p, q, s, s2,&
               psi1, psi2, gh(4)
      INTEGER j
      q=r*q1; p=sigmav*q; s=p*vp; s2=s*s 
      psi1=rho2+a*s2*(vmsq-s*aux6(s)/3.0_r8)
      IF (psi1 > -lneps10) THEN
        DO j= 1,2 
          fv(j)= 0.0_r8
        ENDDO
      ELSE
        gh(1)=1.0_r8; gh(2)=sq-p 
        s=cv*q1*exp(-psi1)/sqrt(p+vm) 
        fv(1)= s; fv(2)= s*gh(2)
      ENDIF
      q=r*q2; p=-sigmav*q; s=p*vp; s2=s*s 
      IF (s+1>eps10) THEN 
        psi2=rho2+a*s2*(vmsq-s*aux6(s)/3.0_r8) 
      ELSE
        psi2=-2.0_r8*lneps10
      ENDIF 
      IF ((psi2 > -lneps10).OR.(p+vm <= 0)) THEN
        DO j= 3,4 
          fv(j)=0.0_r8
        ENDDO
      ELSE
        s= cv*q2*exp(-psi2)/sqrt(p+vm)
        gh(3)= 1.0_r8; gh(4)= sq-p
        fv(3)=s; fv(4)= s*gh(4)   
      ENDIF
      fv(5)= fv(1)+fv(3)
      fv(6)= fv(2)+fv(4)
      fv(3)= 0.0_r8; fv(4)= 0.0_r8  
      END SUBROUTINE intv
    END SUBROUTINE gh123
                      
     SUBROUTINE series(a,x,mode,uax,vax,ierr)
    ! -------------------------------------------------
    ! Calculation of U(a,x), V(a,x) by using McLaurin 
    ! series
    ! -------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   uaxx , uaxx(1), U(a,x)
    !          uaxx(2), U'(a,x)
    !   vaxx , vaxx(1), V(a,x)
    !          vaxx(2), V'(a,x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! -----------------------------------------------------
    USE Someconstants
    USE AiryFunction, ONLY: xpowy
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT):: uax(2)
    REAL(r8), INTENT(OUT):: vax(2)
    REAL(r8) eps, ah, ahp14, ahp34, zerk, f1, f2, err1, err2, x2, am1,&
             a0, bm1, b0, y1, y2, x2k, x2km, x2kd, facto1, facto2,&
             a2, b2, acof, bcof, twok, argu, y1p, y2p,&
             facto1d, facto2d, acofd, bcofd, err1p, err2p,&
             t, t2, ta, xargu, gammaa
    REAL(r8) acofn, bcofn, acofdn, bcofdn
    REAL(r8) mubA,mu2A,mu4A,mubB,mu2B,mu4B,musA,gmusA,&
             musB,gmusB,smusA,smusB
    REAL(r8) betaA,betaB,sa,coefs(0:8),yk,y,ha,dd
    REAL(r8), DIMENSION(0:25) :: gk
    REAL(r8), DIMENSION(0:12) :: wk
    INTEGER, INTENT(IN):: mode
    INTEGER, INTENT(OUT):: ierr
    INTEGER k,i,li
    DATA gk/1.0_r8,.41666666666666666667e-1_r8,0,&
            -.97463348765432098765e-2_r8,0,&
            .12318482904510826965e-1_r8,0,-.37822933917705539682e-1_r8,0,&
            .21514326767209896474_r8,0,-1.9630003732872175294_r8,&
            0,26.256830962378916652_r8,0,-484.19061617504532506_r8,0,&
            11773.948564802034554_r8,0,-365037.92569092371983_r8,&
            0,14054558.808383655048_r8,0,-657894020.31009296820_r8,&
            0,36795321248.737494263_r8/
    DATA wk/1.0_r8,-.17361111111111111111e-2_r8,.81219457304526748971e-3_r8,&
            -.11215312855682914598e-2_r8,.33920312789254653178e-2_r8,&
            -.18817620620360951108e-1_r8,.16870892343666095644_r8,&
            -2.2330644165264046334_r8,40.925670821571152909_r8,&
            -991.44221631790969480_r8,30664.092693023833555_r8,&
            -1178670.6504836921080_r8,55108658.068376820265_r8/
    eps=epss
    ierr=0
    ah=a*0.5_r8
    ahp14=0.25_r8+ah
    ahp34=0.5_r8+ahp14
    err1=1.0_r8; err2=1.0_r8
    err1p=1.0_r8; err2p=1.0_r8
    x2=x*x;
    am1=1.0_r8; a0=a;
    bm1=1.0_r8; b0=a;
    y1=am1+a0*x2*0.5_r8;
    y1p=a0*x;
    y2=(bm1+b0*x2/6.0_r8)*x;
    y2p=y1;
    k=2
    x2k=x2; x2km=x2*x;
    x2kd=x; 
    facto1=0.5_r8; facto2=1.0_r8/6.0_r8
    DO WHILE ((err1 > eps).OR.(err2 > eps).OR.(err1p>eps).OR.(err2p>eps))
      acofn=0; bcofn=0;      
      acofdn=0; bcofdn=0
      DO li=1,2
        twok=2.0_r8*k
        zerk=0.5_r8*(k-1.0_r8)
        a2=a*a0+zerk*(twok-3.0_r8)*am1;
        b2=a*b0+zerk*(twok-1.0_r8)*bm1;
        x2k=x2k*x2; x2km=x2km*x2;
        x2kd=x2kd*x2;
        facto1d=facto1/(twok-1.0_r8);
        facto2d=facto2/twok;
        facto1=facto1d/twok;
        facto2=facto2d/(twok+1.0_r8);
        acof=a2*(x2k*facto1); bcof=b2*(x2km*facto2);
        acofd=a2*(x2kd*facto1d);
        bcofd=b2*(x2k*facto2d);
        y1=y1+acof; y2=y2+bcof;
        y1p=y1p+acofd; y2p=y2p+bcofd;
        k=k+1;  
        am1=a0; a0=a2;
        bm1=b0; b0=b2;
        acofn=acofn+acof; bcofn=bcofn+bcof
        acofdn=acofdn+acofd; bcofdn=bcofdn+bcofd
      ENDDO
    err1=abs(acofn/y1); err2=abs(bcofn/y2)
    err1p=abs(acofdn/y1p); err2p=abs(bcofdn/y2p)
    ENDDO
    IF ((mode==1).AND.(a>0)) THEN
      t=x/(2.0_r8*sqrt(a))
      IF (a>50) THEN
        mubA=sqrt(2.0_r8*(a*0.5_r8+0.25_r8))
        mu2A=mubA*mubA
        mu4A=mu2A*mu2A
        mubB=sqrt(2.0_r8*a)
        mu2B=mubB*mubB
        mu4B=mu2B*mu2B
        DO k=0,25
          IF (k==0) THEN
            musA=1.0_r8
            gmusA=0.0_r8
            musB=1.0_r8
            gmusB=0.0_r8
          ELSE
            musA=musA*mu2A
            musB=musB*mu2B
          ENDIF
          gmusA=gmusA+gk(k)/musA
          gmusB=gmusB+gk(k)/musB
        ENDDO
        DO k=0,12
          IF (k==0) THEN
            musA=1.0_r8
            smusA=0.0_r8
            musB=1.0_r8
            smusB=0.0_r8
          ELSE
            musA=musA*mu4A
            musB=musB*mu4B
          ENDIF
          smusA=smusA+wk(k)/musA
          smusB=smusB+wk(k)/musB
        ENDDO
        betaA=gmusA*gmusA/(sqrttwopi*smusA)
        betaB=gmusB*gmusB/(sqrttwopi*smusB)
      ELSE
        betaA=1.0_r8/(sqrttwopi*aux9(a*0.5_r8+0.25_r8))
        betaB=1.0_r8/(sqrttwopi*aux9(a))
      ENDIF 
      IF (a<30.0_r8) THEN
        y=0.5_r8/a
        sa=exp(1.0_r8/y*log(1.0_r8+y))/exp(1.0_r8)
      ELSE
        coefs(0)=1.0_r8; coefs(1)=-0.5_r8; coefs(2)=11.0_r8/24.0_r8
        coefs(3)=-7.0_r8/16.0_r8; coefs(4)=2447.0_r8/5760.0_r8
        coefs(5)=-959.0_r8/2304.0_r8; coefs(6)=238043.0_r8/580608.0_r8
        coefs(7)=-67223.0_r8/165888.0_r8; coefs(8)=559440199.0_r8/1393459200.0_r8
        sa=0.0_r8
        yk=1.0_r8
        y=0.5_r8/a
        i=0
        DO WHILE ((i<=8).AND.(yk*y>dwarf))
          sa=sa+coefs(i)*yk 
          yk=yk*y
          i=i+1
        ENDDO
      ENDIF
      ha=betaA/(xpowy(sa,0.25_r8)*xpowy(a+0.5_r8,0.25_r8))
      t2=t*t
      IF (a>200) THEN
        coefs(1)=2.0_r8
        coefs(2)=1.0_r8/3.0_r8
        coefs(3)=-1.0_r8/20.0_r8
        coefs(4)=1.0_r8/56.0_r8
        coefs(5)=-5.0_r8/576.0_r8
        coefs(6)=7.0_r8/1408.0_r8
        ta=coefs(1)*t
        yk=t
        y=t2
        i=2
        DO WHILE ((i<=6).AND.(yk*y>dwarf))
          yk=yk*y
          ta=ta+coefs(i)*yk 
          i=i+1
        ENDDO
        dd=exp(a*ta)
      ELSE
        dd=exp(a*(log(t+sqrt(t2+1.0_r8))+t*sqrt(t2+1.0_r8)))
      ENDIF
      f1=sqrtpi*ha*dd
      f2=-1.0_r8/(sqrt2*ha)*betaB*dd
      uax(1)=f1*y1+f2*y2
      uax(2)=f1*y1p+f2*y2p
      dd=1.0_r8/dd
      argu=pi*(0.25_r8-a*0.5_r8)
      f1=2.0_r8/sqrtpi*ha/betaB*dd*(1.0_r8+sin(pi*a))*0.5_r8
      f2=sqrt2/(pi*ha)*dd*(1.0_r8-sin(pi*a))*0.5_r8
      vax(1)=(f1*y1+f2*y2)     
      vax(2)=(f1*y1p+f2*y2p)
    ELSE
      IF (a>0) THEN
        f1=sqrtpi/(xpowy(2.0_r8,ahp14)*gamma(ahp34))
        f2=-sqrtpi/(xpowy(2.0_r8,ah-0.25_r8)*gamma(ahp14))
      ELSE
        argu=pi*ahp14
        f1=gamma(0.25_r8-ah)*cos(argu)/(sqrtpi*xpowy(2.0_r8,ahp14))
        f2=-gamma(0.75_r8-ah)*sin(argu)/(sqrtpi*xpowy(2.0_r8,ah-0.25_r8))
      ENDIF   
      uax(1)=f1*y1+f2*y2
      uax(2)=f1*y1p+f2*y2p
      IF (a >0) THEN  
        f1=xpowy(2.0_r8,ahp14)/pi*gamma(ahp14)*&
           sin(pi*(3.0_r8/4.0_r8-ah))*sin(pi*(3.0_r8/4.0_r8-ah))
        f2=xpowy(2.0_r8,ahp14+0.5_r8)/pi*gamma(ah+3.0_r8/4.0_r8)*&
           sin(pi*(0.25_r8-ah))*sin(pi*(0.25_r8-ah))
        vax(1)=(f1*y1+f2*y2)     
        vax(2)=(f1*y1p+f2*y2p)
      ELSE
        f1=xpowy(2.0_r8,ahp14)*sin(pi*(3.0_r8/4.0_r8-ah))/&
           gamma(3.0_r8/4.0_r8-ah)       
        f2=xpowy(2.0_r8,ahp14+0.5_r8)*sin(pi*(0.25_r8-ah))/&
           gamma(0.25_r8-ah)
        vax(1)=(f1*y1+f2*y2)
        vax(2)=(f1*y1p+f2*y2p)
        IF (mode==1) THEN
          xargu=x*x*0.25_r8+a
          IF (xargu <=0) THEN
            IF ((x<dwarf).AND.(abs(a)<dwarf)) THEN
              gammaa=1.0_r8
            ELSE
              xargu=x*x*0.25_r8+a
              gammaa=exp(a*0.5_r8*(log(abs(a))-1.0_r8))
            ENDIF
            uax(1)=uax(1)*gammaa
            uax(2)=uax(2)*gammaa
            vax(1)=vax(1)/gammaa
            vax(2)=vax(2)/gammaa
          ELSE
            gammaa=exp(a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8)
            uax(1)=uax(1)*gammaa
            uax(2)=uax(2)*gammaa
            vax(1)=vax(1)/gammaa
            vax(2)=vax(2)/gammaa
          ENDIF  
        ENDIF
      ENDIF
    ENDIF
    END SUBROUTINE series
                         
    SUBROUTINE expair(a,x,mode,uaxx,vaxx,ierr)
    ! ------------------------------------------------------------
    !            Airy-type asymptotic expansions
    ! ------------------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   uaxx , uaxx(1), U(a,x)
    !          uaxx(2), U'(a,x)
    !   vaxx , vaxx(1), V(a,x)
    !          vaxx(2), V'(a,x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! --------------------------------------------------------------
    !  The Airy-type asymptotic expansions are used for a<<0 and 
    !  around the turning point x*x/4+a=0
    ! --------------------------------------------------------------
    USE Someconstants
    USE AiryFunction
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT) :: uaxx(2)
    REAL(r8), INTENT(OUT) :: vaxx(2)
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: eps, mu, mu2, mu4, mu8, mu13, mu16, mu23, mu43, mu83, t, t2, zeta,&
                y, zmasf, psi, aa, sq, etta, eta, argu, phis,&
                chi, sas, sbs, scs, sds, as, bs, asp, bso, bsp, bspo, etal,&
                twom13, hmu, mus, gmus, smus, gmu, mu4k, f2, air,&
                dair, bir, dbir, dfacu, dfacup, dfacv, dfacvp, ffa, xargu
    REAL(r8), DIMENSION(0:10,0:40) :: ast, bst
    REAL(r8), DIMENSION(0:30) :: fik,chik
    REAL(r8), DIMENSION(0:25) :: gk
    REAL(r8), DIMENSION(0:12) :: wk
    REAL(r8), DIMENSION(0:18) :: cozmas
    INTEGER, DIMENSION(0:6) :: inda, indb
    INTEGER mode, k, j, l
    DATA inda/40,30,30,30,20,20,20/
    DATA indb/37,30,30,30,20,20,20/
    DATA cozmas/1.885618083164126731735585_r8,.2828427124746190097603378_r8,&
         -.2525381361380526872860159e-1_r8,.4910463758239913363894753e-2_r8,&
         -.1255516301822705121450363e-2_r8,.3718259816936472859679921e-3_r8,&
         -.1208434440504353679395974e-3_r8,.4188900896706268006309575e-4_r8,&
         -.1522610358835666495714500e-4_r8,.5739999368626520519558632e-5_r8,&
         -.2227369320217030245089599e-5_r8,.8848730844862201973674137e-6_r8,&
         -.3584555319099271632854106e-6_r8,.1476133191751092628648807e-6_r8,&
         -.6164726751264643754437699e-7_r8,.2605998126670963041648663e-7_r8,&
         -.1113366163939335549490076e-7_r8,.4801280953394988359287492e-8_r8,&
         -.2087736482939914810074796e-8_r8/
    DATA ast(1,0:39)/-.86458333333333333333e-2_r8,.11118506493506493506e-1_r8,&
       -.11338461538461538462e-1_r8,.10485072951739618406e-1_r8,-.91750414291590762179e-2_r8,&
        .77467458168452212003e-2_r8,-.63785507498774352182e-2_r8,.51549769720023356846e-2_r8,&
       -.41065415400340848415e-2_r8,.32340443813071477202e-2_r8,-.25232263480572751191e-2_r8,&
        .19534081514864215152e-2_r8,-.15023865064299121639e-2_r8,.11490357422627664756e-2_r8,&
       -.87453393997191988595e-3_r8,.66279160411710427603e-3_r8,-.50044124329147178642e-3_r8,&
        .37660561306010746994e-3_r8,-.28257320676903403086e-3_r8,.21145378792045091763e-3_r8,&
       -.15785247292738707697e-3_r8,.11758031706072545009e-3_r8,-.87407457034219086216e-4_r8,&
        .64858335619867878134e-4_r8,-.48045267802574547378e-4_r8,.35535265384237208329e-4_r8,&
       -.26244779856397880734e-4_r8,.19357337318355832102e-4_r8,-.14259616234056385015e-4_r8,&
        .10492181059381593327e-4_r8,-.77117350891421712751e-5_r8,&
        .56623460684682441023388000908606482749399795844647e-5_r8,&
        -.41536087769992473560528795985505977591447090572397e-5_r8,&
        .30441367115808179827127378360809570747317003238776e-5_r8,&
         -.22291254179152811330836670701699386921930200692429e-5_r8,&
        .16310120592387743503983919363881550173527463656721e-5_r8,&
       -.11924797593173726521905994130534586329421615084496e-5_r8,&
         .87123116433670338207099105124965843641863370154514e-6_r8,&
        -.63609261679946347627956411620578214436890625764449e-6_r8,&
        .46411616429358644143012047921855806058625426820438e-6_r8/
    DATA ast(2,0:30)/.56266224308317252266e-2_r8,-.11647291855662719633e-1_r8,&
        .17438935550993679015e-1_r8,-.22248420786275842558e-1_r8,.25681645587606477421e-1_r8,&
       -.27651752813946224558e-1_r8,.28276480334473051867e-1_r8,-.27784472881795235575e-1_r8,&
        .26445276180024289855e-1_r8,-.24523529386345315127e-1_r8,.22252947171162342216e-1_r8,&
       -.19824623260148005537e-1_r8,.17384726703684670525e-1_r8,-.15037729762938088493e-1_r8,&
        .12852408813094941275e-1_r8,-.10868797174027507674e-1_r8,.91049867228358693020e-2_r8,&
       -.75631853088958125607e-2_r8,.62347763146567497230e-2_r8,-.51043369269941506053e-2_r8,&
        .41526903371370890122e-2_r8,-.33591243802274860175e-2_r8,.27029276999493624052e-2_r8,&
       -.21643904392277104271e-2_r8,.17254005528927170693e-2_r8,-.13697460566867932515e-2_r8,&
        .10832120482574887245e-2_r8,-.85354146066446264128e-3_r8,.67031131907963448661e-3_r8,&
       -.52476209811926895365e-3_r8,.40960652212630493071e-3_r8/
    DATA ast(3,0:30)/-.11580906180050927561e-1_r8,.31661887057447902497e-1_r8,&
       -.60537519726707374400e-1_r8,.96022104837560466155e-1_r8,&
       -.13486858218773139158_r8,.17360451334378319335_r8,&
       -.20913158018379167172_r8,.23907654949776246859_r8,&
       -.26192548161477102970_r8,.27699864052741948107_r8,&
       -.28432695853487095482_r8,.28448186552194680157_r8,&
       -.27839658337617157169_r8,.26720327422032545919_r8,&
       -.25209898969048886827_r8,.23424493493810997748_r8,&
       -.21469801828258118088_r8,.19437049669898026028_r8,&
       -.17401212913960327669_r8,.15420903436758014570_r8,&
       -.13539394342452291497_r8,.11786338661988688362_r8,&
       -.10179832039934106068_r8,.87285631081354269540e-1_r8,&
       -.74338768749844687828e-1_r8,.62916431192973831951e-1_r8,&
       -.52938730452196001459e-1_r8,.44300646423185523639e-1_r8,&
       -.36882824571960708148e-1_r8,.30559932289668292664e-1_r8,&
       -.25206873839963280882e-1_r8/
    DATA ast(4,0:30)/.49393729669138401882e-1_r8,-.16460159056974091096_r8,&
        .37671106920022966088_r8,-.70446154818843445310_r8,&
        1.1517696715633945107_r8,-1.7071311993240408422_r8,&
        2.3458575390237919962_r8,-3.0341470817455636226_r8,3.7338648123687769017_r8,&
       -4.4071025283752592817_r8,5.0199043219968235593_r8,-5.5448551987388746055_r8,&
        5.9624855820903817327_r8,-6.2616193802379769230_r8,6.4388904169360639081_r8,&
       -6.4976858367265413542_r8,6.4467642666502533981_r8,-6.2987589550594794269_r8,&
        6.0687264338988964841_r8,-5.7728499953831028737_r8,5.4273610811194061563_r8,&
       -5.0477040039578825218_r8,4.6479413723592056878_r8,-4.2403787420074181144_r8,&
        3.8353760898208164707_r8,-3.4413090294005663100_r8,3.0646425650883640354_r8,&
       -2.7100830787497833682_r8,2.3807788873442699130_r8,-2.0785451219757146248_r8,&
        1.8040941616937906026_r8/
    DATA ast(5,0:28)/-.36006972101648543486_r8,1.3987393559898456251_r8,&
       -3.6892362137755083799_r8,7.8737303343225028183_r8,&
       -14.569386493592020890_r8,24.261162048097651407_r8,&
       -37.213581231065030993_r8,53.417291724674799405_r8,&
       -72.576260862138309707_r8,94.133040133549441048_r8,&
       -117.32368949349208833_r8,141.25103106722442453_r8,&
       -164.96454465710299113_r8,187.53671898263579271_r8,-208.12824122199390069_r8,&
        226.03734163390257598_r8,-240.73138542763382112_r8,251.86109378733459125_r8,&
       -259.25942968790957949_r8,262.92818670198517035_r8,-263.01574432667643324_r8,&
        259.78942292465782198_r8,-253.60552167078665437_r8,244.87958311821722765_r8,&
       -234.05880699965059900_r8,221.59791652287623744_r8,-207.93921917262608293_r8,&
        193.49713382040369616_r8,-178.64709066614237649_r8/
    DATA ast(6,0:25)/4.0084065348046763094_r8,-17.642621241600336767_r8,&
        52.311391233887190108_r8,-124.66354447894727130_r8,256.05306717453359351_r8,&
       -470.80863230664034499_r8,793.63941447186478508_r8,-1246.6101080073623060_r8,&
        1846.1783018847890331_r8,-2600.7449029420965852_r8,3509.0494384975085126_r8,&
       -4559.5844595960940674_r8,5731.0427527316521420_r8,-6993.6748735645957022_r8,&
        8311.3386629781516967_r8,-9643.9722024686523983_r8,10950.214002409855203_r8,&
       -12189.920564766787596_r8,13326.380783686538323_r8,-14328.087748340125799_r8,&
        15169.991649192512856_r8,-15834.215245207554837_r8,16310.260903736659595_r8,&
       -16594.773183346838010_r8,16690.942899479998399_r8,-16607.648660184592203_r8/
    DATA ast(7,0:22)/-63.290222545472642845_r8,309.49063765889768949_r8,&
       -1013.7206553394096370_r8,2655.5227786862968090_r8,-5969.3059766708673407_r8,&
        11964.743360104188382_r8,-21906.959390128289116_r8,37252.392923816787966_r8,&
       -59544.188680047082730_r8,90278.109242533664710_r8,-130754.87888773039296_r8,&
        181936.88302880890521_r8,-244326.24121616523617_r8,317877.96339339844357_r8,&
       -401957.01966033293854_r8,495342.63338040766157_r8,-596277.79843283893916_r8,&
        702557.58142432196808_r8,-811646.60776554538511_r8,920814.40627216124150_r8,&
       -1027276.9472424641737_r8,1128333.5453824016662_r8,-1221490.0112516386311_r8/
    DATA bst(0,0:36)/-.32142857142857142857e-1_r8,.15555555555555555556e-1_r8,&
       -.10085343228200371058e-1_r8,.68923076923076923077e-2_r8,&
       -.47973770749280953363e-2_r8,.33672793089263677499e-2_r8,&
       -.23740671688092001398e-2_r8,.16782309959031527659e-2_r8,&
       -.11883320683635189568e-2_r8,.84238779139280588983e-3_r8,&
       -.59762438660929409956e-3_r8,.42422138636439454428e-3_r8,&
       -.30126023626926855833e-3_r8,.21400917293042855689e-3_r8,&
       -.15206634775554541912e-3_r8,.10807399278785359819e-3_r8,&
       -.76820938020122548526e-4_r8,.54612918802955620593e-4_r8,&
       -.38829207461698758135e-4_r8,.27609666542412784727e-4_r8,&
       -.19633469331004862369e-4_r8,.13962435063909612290e-4_r8,&
       -.99300043010631579079e-5_r8,.70625008660394052299e-5_r8,&
       -.50232596984969135714e-5_r8,.35729626884250873769e-5_r8,&
       -.25414708010873361644e-5_r8,.18078147381341295781e-5_r8,&
       -.12859778274667428677e-5_r8,.91479251956680900289e-6_r8,&
       -.65075912613628042224e-6_r8,.46294093040524422199e-6_r8,&
       -.32933491162045547170e-6_r8,.23429130240268590734e-6_r8,&
       -.16667872738882887922e-6_r8,.11857940849768171016e-6_r8,&
       -.84361252507627603914e-7_r8/
    DATA bst(1,0:30)/.11616363324175824176e-1_r8,-.10738690547767928720e-1_r8,&
        .11273908748814210999e-1_r8,-.11330200650226966016e-1_r8,&
        .10882662128717442695e-1_r8,-.10073431983586506837e-1_r8,&
        .90528041420056262576e-2_r8,-.79437161494535131346e-2_r8,&
        .68354205437865891183e-2_r8,-.57866877852097172840e-2_r8,&
        .48319188837933835377e-2_r8,-.39875048709666930997e-2_r8,&
        .32573797517277923555e-2_r8,-.26374472093943733107e-2_r8,&
        .21188963477292799694e-2_r8,-.16905597364664082581e-2_r8,&
        .13405067192519113713e-2_r8,-.10570579364055743550e-2_r8,&
        .82938060092290405418e-3_r8,-.64779233018712147103e-3_r8,&
        .50387108159099364466e-3_r8,-.39044280607171355606e-3_r8,&
        .30149757745508396450e-3_r8,-.23206892067503710671e-3_r8,&
        .17809916581662199089e-3_r8,-.13630510145467900797e-3_r8,&
        .10405223392871875512e-3_r8,-.79241928092168523010e-4_r8,&
        .60213081558437480541e-4_r8,-.45658357387345341779e-4_r8,&
        .34554059520366810843e-4_r8/
    DATA bst(2,0:30)/-.17619791271984698754e-1_r8,.22460738436828023930e-1_r8,&
       -.31095536623705537180e-1_r8,.39843610073496917522e-1_r8,&
       -.47521689852425072253e-1_r8,.53475980779992135034e-1_r8,&
       -.57414884418372632722e-1_r8,.59320052235418060306e-1_r8,&
       -.59362779310300406734e-1_r8,.57828321677111474388e-1_r8,&
       -.55053943676791592811e-1_r8,.51382942898551505119e-1_r8,&
       -.47133919046907858690e-1_r8,.42582911713303221181e-1_r8,&
       -.37955441784732551387e-1_r8,.33425554856704672102e-1_r8,&
       -.29119366448832191622e-1_r8,.25121136014826860211e-1_r8,&
       -.21480423732279927589e-1_r8,.18219346128000877782e-1_r8,&
       -.15339318528807313005e-1_r8,.12826952106907783895e-1_r8,&
       -.10658970903423837305e-1_r8,.88061445510104236352e-2_r8,&
       -.72363110956985624537e-2_r8,.59166055051532985471e-2_r8,&
       -.48150249286128209578e-2_r8,.39014607485083487893e-2_r8,&
       -.31483167600335778722e-2_r8,.25308172470515622861e-2_r8,&
       -.20270915098905906085e-2_r8/
    DATA bst(3,0:30)/.60909763139637582786e-1_r8,-.96541486771214866841e-1_r8,&
        .16264377675453827437_r8,-.24915885452321001414_r8,&
        .35009667669419877168_r8,-.45836730339191387087_r8,&
        .56648575100308131913_r8,-.66748754950937166485_r8,&
        .75561704712993985977_r8,-.82671466079887381067_r8,&
        .87832844401420206021_r8,-.90961364830930523024_r8,&
        .92109362061208216666_r8,-.91434938542624464396_r8,&
        .89169172213458344119_r8,-.85585378144126278784_r8,&
        .80972748910911579700_r8,-.75615478142621795241_r8,&
        .69777564969133675222_r8,-.63692892519825209993_r8,&
        .57559825199333248396_r8,-.51539418379669499028_r8,&
        .45756322033392356537_r8,-.40301535937866445856_r8,&
        .35236298156274064514_r8,-.30596531074994938002_r8,&
        .26397410487303626816_r8,-.22637751012529200117_r8,&
        .19304009378463816988_r8,-.16373793768410157069_r8,&
        .13818833224511412946_r8/
    DATA bst(4,0:29)/-.37829872479673768094_r8,.70699348356478890943_r8,&
        -1.3865797541213790568_r8,2.4460312321724494685_r8,&
        -3.9207588293947925720_r8,5.8080313668527459768_r8,&
        -8.0630182513775208787_r8,10.603669134030102425_r8,&
       -13.320620699455056397_r8,16.089517395703050953_r8,&
       -18.783542996117974847_r8,21.284502840495644149_r8,&
       -23.491416999286173127_r8,25.326169364855214008_r8,&
       -26.736231662339052310_r8,27.694813083767226802_r8,&
       -28.198977054536741489_r8,28.266337819667546672_r8,&
       -27.930931277990815536_r8,27.238778053930476484_r8,&
       -26.243549745521884151_r8,25.002633225307177664_r8,&
       -23.573777915301600688_r8,22.012416302016647166_r8,&
       -20.369672976431342515_r8,18.691023002223376311_r8,&
       -17.015524914036864532_r8,15.375534273694109775_r8,&
       -13.796797055459840802_r8,12.298824764252335256_r8/
    DATA bst(5,0:26)/3.7008098833796096974_r8,-7.8943179119715393505_r8,&
        17.523666638205992170_r8,-34.732110691221732984_r8,&
        62.145391921905206236_r8,-102.17014739446708717_r8,&
        156.60003201213783498_r8,-226.31321551629452367_r8,311.08916186164296289_r8,&
       -409.56133862117968200_r8,519.30074205116124335_r8,-637.00727007547042847_r8,&
        758.77505311063660762_r8,-880.39404139944726093_r8,997.65225418279181001_r8,&
       -1106.6093103403454609_r8,1203.8202784542850274_r8,-1286.4978451886305591_r8,&
        1352.6090567441321809_r8,-1400.9096495972083876_r8,1430.9238827471423665_r8,&
       -1442.8807756116574511_r8,1437.6189354564596975_r8,-1416.4720513211302631_r8,&
        1381.1460162370242338_r8,-1333.5968879195104838_r8,1275.9168367360627308_r8/
    DATA bst(6,0:22)/-52.440232872505911846_r8,124.91492555422010956_r8,&
       -307.80505874513857089_r8,673.65965906811191089_r8,-1324.7001416829727542_r8,&
        2383.2377235388771970_r8,-3981.7258172782825530_r8,6249.8655385949204720_r8,&
       -9300.3111310734301499_r8,13214.918162292874992_r8,-18033.350332471047341_r8,&
        23745.434126825075803_r8,-30288.059814249641659_r8,37546.796886266045900_r8,&
       -45361.823914328331689_r8,53537.336176016916081_r8,-61853.323266034541002_r8,&
        70078.506118453102599_r8,-77983.267923314792836_r8,85351.571541509186577_r8,&
       -91991.086631042572063_r8,97741.013229975714811_r8,-102477.35110602911554_r8/
    DATA fik/1.0_r8,-.1_r8,.42142857142857142857e-1_r8,-.21912698412698412698e-1_r8,&
         .12448755411255411255e-1_r8,-.74288549743906886764e-2_r8,&
         .45746350641201831678e-2_r8,-.28791700676332804184e-2_r8,&
         .18414150026089362744e-2_r8,-.11923004331977092900e-2_r8,&
         .77957025714134722714e-3_r8,-.51376250832318276629e-3_r8,&
         .34081277524518474479e-3_r8,-.22733442456475481892e-3_r8,&
         .15235591836762461067e-3_r8,-.10252250859554482318e-3_r8,&
         .69234005415816302727e-4_r8,-.46900261834713094097e-4_r8,&
         .31859086498603514106e-4_r8,-.21695264658176838875e-4_r8,&
         .14806797927678100534e-4_r8,-.10125785547514577093e-4_r8,&
         .69372508648302511184e-5_r8,-.47606638519588906843e-5_r8,&
         .32719619838504728166e-5_r8,-.22519366067543023254e-5_r8,&
         .15519011644445928413e-5_r8,-.10707548626041065912e-5_r8,&
         .73959982335703643779e-6_r8,-.51138836489547516290e-6_r8,&
         .35393386504561175831e-6_r8/
    DATA chik/-.1_r8,.74285714285714285714e-1_r8,-.54095238095238095238e-1_r8,&
         .39063615749330035044e-1_r8,-.28085509411223696938e-1_r8,&
         .20139984242514854760e-1_r8,-.14417859514984444956e-1_r8,&
         .10309442176430897606e-1_r8,-.73655122348323755815e-2_r8,&
         .52589148288898677708e-2_r8,-.37529918053932661301e-2_r8,&
         .26772689420573039019e-2_r8,-.19092899409669881562e-2_r8,&
         .13612619157861153785e-2_r8,-.97033149280715834007e-3_r8,&
         .69154709536062995295e-3_r8,-.49278576070368617280e-3_r8,&
         .35110626274417083999e-3_r8,-.25013277029410498759e-3_r8,&
         .17818059964612554284e-3_r8,-.12691507492553136865e-3_r8,&
         .90392689282224654705e-4_r8,-.64376056356811587852e-4_r8,&
         .45844739898716106234e-4_r8,-.32646108214340078356e-4_r8,&
         .23246222222742878232e-4_r8,-.16552150417025056676e-4_r8,&
         .11785263078514715153e-4_r8,-.83908989293805415224e-5_r8,&
         .59739749862631200831e-5_r8,-.42530962002148905661e-5_r8/
    DATA gk/1.0_r8,.41666666666666666667e-1_r8,0,-.97463348765432098765e-2_r8,0,&
         .12318482904510826965e-1_r8,0,-.37822933917705539682e-1_r8,&
         0,.21514326767209896474_r8,0,-1.9630003732872175294_r8,&
         0,26.256830962378916652_r8,0,-484.19061617504532506_r8,0,&
         11773.948564802034554_r8,0,-365037.92569092371983_r8,0,&
         14054558.808383655048_r8,0,-657894020.31009296820_r8,0,&
         36795321248.737494263_r8/
    DATA wk/1.0_r8,-.17361111111111111111e-2_r8,.81219457304526748971e-3_r8,&
         -.11215312855682914598e-2_r8,.33920312789254653178e-2_r8,&
         -.18817620620360951108e-1_r8,.16870892343666095644_r8,&
         -2.2330644165264046334_r8,40.925670821571152909_r8,&
         -991.44221631790969480_r8,30664.092693023833555_r8,&
         -1178670.6504836921080_r8,55108658.068376820265_r8/
    eps=epss
    ierr=0
    mu=sqrt(-2.0_r8*a); t=abs(x)/(mu*sqrt2)
    mu2=mu*mu
    IF (mode==0) THEN   
      IF (t <= 1) THEN
        ! ------------------------------------------------
        ! Check for possible overflow/underflow problems
        ! ------------------------------------------------
        IF (a*0.5_r8*(log(abs(a))-1.0_r8) <log(dwarf)) THEN
          ierr=1
          uaxx(1)=giant
          uaxx(2)=giant
          vaxx(1)=dwarf
          vaxx(2)=dwarf
        ENDIF
      ELSE
        xargu=abs(x*x*0.25_r8+a)
        ffa=a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8
        ! ------------------------------------------------
        ! Check for possible overflow/underflow problems
        ! ------------------------------------------------
        IF ((ffa < log(dwarf)).OR.(ffa>log(giant))) THEN
          ierr=1
          uaxx(1)=giant
          uaxx(2)=giant
          vaxx(1)=dwarf 
          vaxx(2)=dwarf
        ENDIF
      ENDIF
    ENDIF
    IF (ierr==0) THEN
      mu4=mu2*mu2
      mu8=mu4*mu4
      mu13=xpowy(mu,onethird)
      mu16=xpowy(mu,onesix)
      mu23=mu13*mu13
      mu43=xpowy(mu4,onethird)
      mu83=xpowy(mu8,onethird)
      t2=t*t
      IF (abs(1.0_r8-t)<eps*100) THEN
        psi=0.0_r8
        zeta=0.0_r8
      ELSEIF (t > 1.1_r8) THEN
        psi=0.5_r8*(t*sqrt(t2-1.0_r8)-log(t+sqrt(t2-1.0_r8)))
        zeta=xpowy(3*psi*0.5_r8,twothird)
      ELSEIF (t < 0.9_r8) THEN
        sq=sqrt((1.0_r8-t)*(1.0_r8+t))
        etta=0.5_r8*(acos(t)-t*sq)
        zeta=-xpowy(3*etta*0.5_r8,twothird)
      ELSE
        ! Taylor series
        y=t-1.0_r8
        zmasf=cozmas(12)
        DO k=0,11
          j=11-k
          zmasf=cozmas(j)+y*zmasf
        ENDDO
        zmasf=0.5_r8*xpowy(abs(y),1.5_r8)*zmasf
        IF (t < 1) THEN
          zeta=-xpowy(1.5_r8*zmasf,twothird)
        ELSE
          zeta=xpowy(1.5_r8*zmasf,twothird)
        ENDIF
      ENDIF
      argu=zeta*mu43
      ! Computation of g(mu)
      DO k=0,25
        IF (k==0) THEN
          mus=1.0_r8
          gmus=0.0_r8
        ELSE
          mus=mus*mu2
        ENDIF
        gmus=gmus+gk(k)/mus
      ENDDO
      DO k=0,12
        IF (k==0) THEN
          mus=1.0_r8
          smus=0.0_r8
        ELSE
          mus=mus*mu4
        ENDIF
        smus=smus+wk(k)/mus
      ENDDO
      IF ((t > 0.9_r8).AND.(t < 1.1_r8)) THEN
      ! Taylor series for the functions phi, chi
        phis=0.0_r8
        chi=0.0_r8
        twom13=xpowy(2.0_r8,-onethird)
        eta =zeta*twom13
        etal=eta
        DO k=0,20
          IF (k==0) THEN
            etal=1.0_r8
            phis=0.0_r8
            chi=0.0_r8
          ELSE
            etal=etal*eta
          ENDIF
          phis=phis+fik(k)*etal
          chi=chi+chik(k)*etal
        ENDDO
        phis=xpowy(2.0_r8,-1.0_r8/6.0_r8)*phis
        chi=twom13*chi
      ELSE
      ! Exact values
        phis=xpowy(zeta/((t-1.0_r8)*(t+1.0_r8)),0.25_r8)
        chi=0.25_r8*(1.0_r8-2.0_r8*t*xpowy(phis,6.0_r8))/zeta
      ENDIF
      sas=1.0_r8
      sbs=0.0_r8
      scs=0.0_r8
      sds=1.0_r8
      bs=0.0_r8
      bsp=0.0_r8
      twom13=xpowy(2.0_r8,-onethird)
      eta =zeta*twom13
      DO k=0,6
        bso=bs
        bspo=bsp
        as=0.0_r8
        bs=0.0_r8
        asp=0.0_r8
        bsp=0.0_r8
        IF (k > 0) THEN
          mu4k=mu4k*mu4
        ELSE
          mu4k=1.0_r8
        ENDIF
        f2=1.0_r8/mu4k
        IF (k <= 5) THEN
          DO l=0,40
            IF (l.EQ.0) THEN
              etal=1.0_r8
            ELSE
              etal=etal*eta
            ENDIF
            IF (k >0 ) THEN
              IF (l < inda(k)) THEN
                as=as+ast(k,l)*etal
              ENDIF
              IF ((l+1) < inda(k)) THEN
                asp=asp+(l+1)*ast(k,l+1)*etal
              ENDIF
            ENDIF
            IF (l < indb(k)) THEN
              bs=bs+bst(k,l)*etal
            ENDIF
            IF ((l+1) < indb(k)) THEN
              bsp=bsp+(l+1)*bst(k,l+1)*etal
            ENDIF
          ENDDO
          bs=bs/twom13
          asp=asp*twom13
          sas=sas+as*f2
          sbs=sbs+bs*f2    
        ENDIF   
        IF (k>=1) THEN
          sds=sds+(as+chi*bso+bspo)*f2
          scs=scs+(chi*as+asp+zeta*bs)*f2
        ELSE
          scs=scs+(chi+zeta*bs)*f2
        ENDIF
      ENDDO
      IF (t <=1) THEN
        CALL aibi(argu,air,bir)
        CALL aibip(argu,dair,dbir)
        aa=abs(a)
        hmu=xpowy(aa,-0.25_r8)/sqrt2
        gmu=hmu/gmus
        dfacu=2*sqrtpi*mu13*gmu*phis
        dfacup=sqrt2*sqrtpi*mu23*gmu/phis
        dfacv=phis/(mu23*gmu*smus)
        dfacvp=1.0_r8/(sqrt2*mu13*gmu*phis*smus)  
      ELSE 
        CALL airsca(argu,air,bir,dair,dbir)
        dfacu=twoexp14*sqrt2*sqrtpi/mu16/gmus*phis
        dfacup=twoexp14*sqrtpi*mu16/(gmus*phis)
        dfacv=twoexp14/(mu16)*phis*gmus/smus
        dfacvp=1.0_r8/twoexp14*mu16*gmus/(smus*phis)
      ENDIF
      IF (mode==0) THEN
        IF (t>1) THEN
         ! -----------------------------------
         ! Calculation of Unscaled Functions
         ! -----------------------------------
          xargu=abs(x*x*0.25_r8+a)
          ffa=exp(a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8)
          dfacu=dfacu/ffa
          dfacup=dfacup/ffa
          dfacv=ffa*dfacv
          dfacvp=ffa*dfacvp
        ELSE
          ! -----------------------------------
          ! Calculation of Unscaled Functions
          ! ----------------------------------- 
          ffa=exp(a*0.5_r8*(log(abs(a))-1.0_r8))
          dfacu=dfacu/ffa
          dfacup=dfacup/ffa
          dfacv=dfacv*ffa
          dfacvp=dfacvp*ffa
        ENDIF
      ENDIF
      uaxx(1)=dfacu*(air*sas+dair*sbs/mu83)  
      uaxx(2)=dfacup*(air*scs/mu43+dair*sds)
      vaxx(1)=dfacv*(bir*sas+dbir*sbs/mu83)  
      vaxx(2)=dfacvp*(bir*scs/mu43+dbir*sds)
    ENDIF    
    END SUBROUTINE expair

    SUBROUTINE expaelem(a,x,mode,uaxx,vaxx,ierr)
    ! --------------------------------------------------
    ! Calculation of U(a,x), V(a,x) by using asymptotic
    ! expansions in terms of elementary functions
    ! --------------------------------------------------
    ! Inputs:
    !   a ,    order of the functions
    !   x ,    argument of the functions
    !   mode , mode=0, unscaled functions
    !          mode=1, scaled functions
    ! Outputs:
    !   uaxx , uaxx(1), U(a,x)
    !          uaxx(2), U'(a,x)
    !   vaxx , vaxx(1), V(a,x)
    !          vaxx(2), V'(a,x)
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow or/and underflow problems 
    ! ---------------------------------------------------------
    USE Someconstants
    USE AiryFunction
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: mode
    REAL(r8), INTENT(OUT) :: uaxx(2)
    REAL(r8), INTENT(OUT) :: vaxx(2)
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) eps, mu, mu2, t, t2p1, facu, t2, t3, t4, t5, t6, t7,&
             t8, t9, t10, t11, t12, t13, t14, t15, t16, t17,&
             t18, t19, t20, t21, t22, t23, t24, t25, t26, t27,&
             t28, t29, t30, t31, t32, t33, t34, t35, t36, t37,&
             t38, t39, t40, t41, tau, tau2, tau3, tau4, tau5,&
             tau6, tau7, tau8, tau9, tau10, tau11, tau12, tau13,&
             tau14, tau15, tau16, tau17, tau18, tau19, tau20,&
             tau21, tau22, tau23, tau24, tau25, tau26, tau27,&
             tau28, tau29, tau30, tau31, tau32, tau33, tau34,&
             tau35, tau36,& 
             ffa, x2, xargu, sinpia, smus, musp, musn, sphip, schip,&
             sphin, schin, gammaa2, e1, e2, e3, e4, delta1,&
             delta2, delta3, delta4, maxi, mu4, tts, ttss, &
             eta, mus, gmus, s1, s2, s3, s4, mus1, mus2
    REAL(r8) ::argu3, dl,dnew, beta
    REAL(r8), DIMENSION(0:15) :: phi,chi,u,v
    REAL(r8), DIMENSION(0:25) :: gk
    REAL(r8), DIMENSION(0:12) :: wk
    INTEGER k,itwok
    DATA gk/1.0_r8,.41666666666666666667e-1_r8,0,&
            -.97463348765432098765e-2_r8,0,.12318482904510826965e-1_r8,&
            0,-.37822933917705539682e-1_r8,0,.21514326767209896474_r8,&
            0,-1.9630003732872175294_r8,0,26.256830962378916652_r8,0,&
            -484.19061617504532506_r8,0,11773.948564802034554_r8,0,&
            -365037.92569092371983_r8,0,14054558.808383655048_r8,0,&
            -657894020.31009296820_r8,0,36795321248.737494263_r8/
    DATA wk/1.0_r8,-.17361111111111111111e-2_r8,.81219457304526748971e-3_r8,&
            -.11215312855682914598e-2_r8,.33920312789254653178e-2_r8,&
            -.18817620620360951108e-1_r8,.16870892343666095644_r8,&
            -2.2330644165264046334_r8,40.925670821571152909_r8,&
            -991.44221631790969480_r8,30664.092693023833555_r8,&
            -1178670.6504836921080_r8,55108658.068376820265_r8/
    eps=epss
    ierr=0
    mu=sqrt(2.0_r8*abs(a))
    t=abs(x)/(mu*sqrt2)
    IF ((a<0).AND.(t<1)) THEN
      IF (mode==0) THEN
      ! ------------------------------------------------
      ! Check for possible overflow/underflow problems
      ! ------------------------------------------------
        IF (a*0.5_r8*(log(abs(a))-1.0_r8) <log(dwarf)) THEN
          ierr=1
          uaxx(1)=giant
          uaxx(2)=giant
          vaxx(1)=dwarf
          vaxx(2)=dwarf
        ENDIF
      ENDIF
      IF (ierr==0) THEN
        mu2=mu*mu
        t2=t*t; t3=t2*t; t4=t3*t; t5=t4*t; t6=t5*t; t7=t6*t; t8=t7*t; t9=t8*t;
        t10=t9*t; t11=t10*t; t12=t11*t; t13=t12*t; t14=t13*t; t15=t14*t;
        t16=t15*t; t17=t16*t; t18=t17*t; t19=t18*t; t20=t19*t; t21=t20*t;
        t22=t21*t; t23=t22*t; t24=t23*t; t25=t24*t; t26=t25*t; t27=t26*t;
        t28=t27*t; t29=t28*t; t30=t29*t; t31=t30*t; t32=t31*t; t33=t32*t;
        t34=t33*t; t35=t34*t; t36=t35*t; t37=t36*t; t38=t37*t; t39=t38*t;
        t40=t39*t; t41=t40*t  
        u(0)=1.0_r8
        u(1)=.4166666666666666666666667e-1_r8*t*(t2-6.0_r8)
        u(2)=.1258680555555555555555556_r8+.2161458333333333333333333_r8*t2-&
             .7812500000000000000000000e-2_r8*t4
        u(3)=-.6252170138888888888888889_r8*t-.3665002893518518518518519_r8*t3&
             -.6820746527777777777777778e-1_r8*t5+.4385850694444444444444444e-1_r8*t7&
             -.9746334876543209876543210e-2_r8*t9
        u(4)=-.3892736866640946502057613e-2_r8*t6+.3208267274707433127572016_r8&
             -.8071183569637345679012346e-2_r8*t8+.1827437789351851851851852e-2_r8*t10&
             +3.079461293161651234567901_r8*t2+1.279432885440779320987654_r8*t4
        u(5)=.1495080253037428963354889e-9_r8*t*(-34009066266.0_r8-119582875013.0_r8*t2&
             +1994971575.0_r8*t10-3630137104.0_r8*t8+4433574213.0_r8*t6&
             -37370295816.0_r8*t4+82393456.0_r8*t14-617950920.0_r8*t12)
        u(6)=25.84962798338368105790252_r8*t6+.1713039028908536874755046e-1_r8*t14&
            +2.581943865065753536501370_r8-.2309715544595780055849500e-2_r8*t16&
            -.5076538089737495643136799e-1_r8*t12-.3186668406760671571869489_r8*t8&
            +.8601473529456107027269498e-1_r8*t10+62.68920892021549979121503_r8*t2&
            +121.7179460820865799300044*t4;
        u(7)=-.7865531634245733182633045e-15_r8*t*(110065274285017326.0_r8&
            -125776738623632286.0_r8*t10+839442893976461885.0_r8*t2-51589093643118259.0_r8*t14&
            +96739097673864045.0_r8*t12+1014872932457342145.0_r8*t4+18208957214179542.0_r8*t16&
            +80153092162171710.0_r8*t6+114221755795701455.0_r8*t8-3833677186738392.0_r8*t18&
            +365112113022704.0_r8*t20)
        u(8)=5596.879044604074283137768_r8*t6+13.69521272837774834409345_r8*t14&
            +2.638371331228334872355069_r8*t18-7.390056136076261345369326_r8*t16&
            -17.49598518388745077640624_r8*t12+679.5130553150720497530461_r8*t8&
            +.5384626640674465938317402e-1_r8*t22+11.60187790297422517485544_r8*t10&
            -.5608986084035902019080627_r8*t20+1769.339898073369191761939_r8*t2&
            +7072.617883848053965523396_r8*t4+46.00893448266345421442011_r8
        u(9)=.1716369641152583388897846e-22_r8*t*(-.1359234799748243655770368e27_r8&
            +196970983448449506582784.0_r8*t26+.3935158819364260070148776e27_r8*t10&
            -.4824119188126012644493997e27_r8*t12+.4501074357324417329611653e27_r8*t14&
            -.5072207498713186429103949e27_r8*t8-.4586684500765682622356092e28_r8*t4&
            -.2551586343346702497275826e28_r8*t6+.1671990097243978109434037e27_r8*t18&
            -.3176839480721397029728588e27_r8*t16+.1661834830959361591291224e26_r8*t22&
            -.6369835812416149521724992e26_r8*t20-.1825949193391480686660700e28_r8*t2&
            -2659108276554068338867584.0_r8*t24)
        u(10)=914112.2670988274244565704_r8*t6+1457.987476485281604308530_r8*t14&
             +981.4809913805934941962403_r8*t18-52.76053458609560685621252_r8*t24&
             -.6338906553354126033547987_r8*t28+8.504699625750119095010216_r8*t26&
             +1213.983451915325067883233_r8-1375.290270384664889535434_r8*t16&
             -1323.383956656255885190354_r8*t12+415376.7168317029950197364_r8*t8&
             +200.5566577594646000950658_r8*t22+35768.04278735927027177930_r8*t10&
             -521.6682753837636432341761_r8*t20+75167.67289016707645386544_r8*t2&
             +532260.3452982493777355059_r8*t4
        u(11)=-.3359006667469412593393634e-31_r8*t*(.2788947144378624557826141e37_r8+&
             .1212394967494477547839171e34_r8*t32-.2000451696365887953934633e35_r8*t30&
             +.1550285194341124782611328e36_r8*t28-.7492595928004916139721089e36_r8*t26&
             +.2528614856038361343385493e37_r8*t24-.6321380226416724568027542e37_r8*t22&
             +.1211628603883129655678535e38_r8*t20-.1817616151670370488438873e38_r8*t18&
             -.2039447456964012213382125e38_r8*t14+.1526729786684424672364467e38_r8*t12&
             -.4167726542712734865896683e35_r8*t10+.1286775691925214910694996e39_r8*t8&
             +.3375397423010140569245605e39_r8*t6+.2668957183771871703360651e39_r8*t4&
             +.6033537049761229856043034e38_r8*t2+.2158825430879550563836032e38_r8*t16)
        u(12)=153053249.2910557053308385_r8*t6+78424.56411652519477233020_r8*t14+&
             128182.4650008160699499318_r8*t18-38558.08568947489765851297_r8*t24&
             -4637.584888683471835001935_r8*t28+15541.01589815635919400013_r8*t26&
             -119592.7457128667753060960_r8*t16+2857220.354347191785055925_r8*t12&
             +151197224.7859743079793903_r8*t8+965.8070562838161807622995_r8*t30&
             +73303.13746459791327933258_r8*t22+46150931.05589019507584725_r8*t10&
             +7.635830211413084690217302_r8*t34-125.3548793040314736644007_r8*t32&
             -108987.9643972328459888073_r8*t20+4457999.843158161718417669_r8*t2&
             +50231562.58692156905275991_r8*t4+48082.01508844808488819284_r8
        u(13)=.7703374808710407147877427e-40_r8*t*(-.6926180733786124297783135e47_r8&
             -.3490732449110820673555558e50_r8*t6+.8082605068869249980049792e48_r8*t22&
             +.2682009600079171280174385e48_r8*t26-.5172558239822205556645540e48_r8*t24&
             -.1109805333478410099533049e48_r8*t28+.1487945383913154715160031e46_r8*t34&
             +.3580113360465377112428930e47_r8*t30-.8679377325149432331738492e46_r8*t32&
             +.8249366361784484160838564e43_r8*t38-.1608626440547974411363520e45_r8*t36&
             -.9370275546601561256767938e48_r8*t16+.6651854621331235420882608e48_r8*t14&
             -.7836351143702021345383874e48_r8*t12-.7007020293314608957292234e49_r8*t10&
             -.2782407158660428165686684e50_r8*t8-.1553433992040236087729803e50_r8*t4&
             -.2221765150562660149162442e49_r8*t2-.1032888153750967506466895e49_r8*t20&
             +.1084709216891677773429753e49_r8*t18)
        u(14)=29777856511.31303518396996_r8*t6+405725357.0120253215335119_r8*t14&
             +12581472.10822014648877845_r8*t18-11128041.92193217407836428_r8*t24&
             -3746199.514655045078233992_r8*t28+7174530.658451246530535730_r8*t26&
             -9965848.257226286974264551_r8*t16+2313.542955855017656723246_r8*t38&
             +7954048017.175341662787203_r8*t12+52707258002.79912630559216_r8*t8&
             +1560495.910682062257811549_r8*t30+14108384.19513404699749891_r8*t22&
             +35147200819.52248373734473_r8*t10+123538.5886859896169363174_r8*t34&
             -506576.1722262477285863913_r8*t32-21294.73752084727438922035_r8*t36&
             -14687754.72489477423506592_r8*t20+358625041.4599066211836996_r8*t2&
             +6078341745.370762847643443_r8*t4-119.1524269109880338226565_r8*t40&
             +2674542.439978280740741820_r8
        u(15)=-414790432.6994505167392778_r8*t+6945.991168730899603944917_r8+&
             510652.2760188471673646379_r8*t2+3506139.062413368413555602_r8*t18&
             -4424351.480098610189178738_r8*t20+4507540.556879283080170755_r8*t22&
             +2493813.214756674570007723_r8*t26-3725064.570331952034347715_r8*t24&
             +5707621.898105787921162792_r8*t4-19329184319.93965811896088_r8*t3&
             -204822062695.9835996571353_r8*t5-746830449214.8311385338279_r8*t7&
             -1071759390145.197565881790_r8*t9+21556032.59225403242059432_r8*t6&
             +35040737.47008682938237320_r8*t8+24773160.94128538674006004_r8*t10&
             +1461033.332711555204742717_r8*t14+5971318.466177208703908447_r8*t12&
             -2216567.774820235631058887_r8*t16-1343045.511157374057238458_r8*t28&
             -190443.7698884165190432579_r8*t32+574046.2785186452537735790_r8*t30&
             +47314.98337952692225675031_r8*t34+913.5019396509082593070329_r8*t38&
             -8288.707145438840911983684_r8*t36-47.66097076439521352906259_r8*t40&
             +2815957.942354381680464703_r8*t29-4949452.062439341149009949_r8*t27&
             +6883478.737523361808313990_r8*t25-122150797756.7143021925426_r8*t13&
             +16322778.36996041226643312_r8*t17-5600428212.110544205675050_r8*t15&
             -610079501215.4117897742389_r8*t11-4100118.049694015514682209_r8*t19&
             +6373075.827367999613237065_r8*t21-7547338.888154381034617815_r8*t23&
             -110849.0793687075287393032_r8*t35+432724.8534694151550145610_r8*t33&
             -1257756.708472409889487035_r8*t31+119.1524269109880338226565_r8*t41&
             -2244.037373490274636993364_r8*t39+19926.89640129648768991938_r8*t37
        v(0)=u(0)
        v(1)=u(1)+0.5_r8*t*u(0)
        v(2)=-.1241319444444444444444444_r8-.2838541666666666666666667_r8*t2&
            +.1302083333333333333333333e-1_r8*t4
        v(3)=.6252170138888888888888889_r8*t+.5749059606481481481481482_r8*t3&
            -.8773871527777777777777778e-1_r8*t5+.4385850694444444444444444e-1_r8*t7&
            -.9746334876543209876543210e-2_r8*t9
        v(4)=-.3816782005529835390946500e-2_r8*t6-.3043902864181455761316872_r8&
            +.1385806990258487654320987e-1_r8*t8-.3045729648919753086419753e-2_r8*t10&
            -3.334384192949459876543209_r8*t2-1.443856321735146604938272_r8*t4
        v(5)=5.084628339853796939300410_r8*t+19.57347561662252711334020_r8*t3&
            +.3028328551887274489271016_r8*t11-.5607805781707374951009209_r8*t9&
            +.5729826674329610798794825_r8*t7+5.264663972579893261316874_r8*t5&
            +.1231848290451082696453067e-1_r8*t15-.9238862178383120223397999e-1_r8*t13
        v(6)=-28.17555842815603248376414_r8*t6-.2906392060283023236943954e-1_r8*t14&
            -2.502684474788043402799041_r8+.3849525907659633426415834e-2_r8*t16&
            +.9037170913188460136551619e-1_r8*t12+.4309883571133621996252204_r8*t8&
            -.1608534918423846448703581_r8*t10-64.67370051767834022936110_r8*t2&
            -129.7003433719719481799922_r8*t4
        v(7)=86.57218967207372000385268_r8*t+860.9772677404372680912355_r8*t3&
            +1121.756470592261700060500_r8*t5+257.5263468684699442083095_r8*t7&
            -92.32733782717116555772998_r8*t9-76.18059281400228878625290_r8*t13&
            -14.32808718833711404666971_r8*t17+40.61769611078856415222867_r8*t15&
            +99.05203232887895800013400_r8*t11+3.015390918777700925457745_r8*t19&
            -.2871800875026381833769281_r8*t21
        v(8)=-8015.583638506729061155720_r8*t6-24.41008329244471177080736_r8*t14&
            -4.519609259635658228333184_r8*t18+12.88479579364600159304798_r8*t16&
            +32.09421420855616969862385_r8*t12-1043.207983739144800727558_r8*t8&
            -.8974377734457443230529005e-1_r8*t22-23.18565586367109152015149_r8*t10&
            +.9467968509852602608208101_r8*t20-1986.189381518536462868531_r8*t2&
            -9133.569273415523777237354_r8*t4-40.56325518941026578943259_r8
        v(9)=2332.949345485996505888585_r8*t+36760.58162380117473168501_r8*t3&
            +99845.48057991644426620303_r8*t5+65183.71340269782868490893_r8*t7&
            +2622.201227310636710039871_r8*t9+3.380750161788867217892261_r8*t27&
            -45.64012718414970744154552_r8*t25-8306.938556402440897450785_r8*t13&
            -5467.955632719681768619915_r8*t17+7751.066214952751247768557_r8*t15&
            +6712.462538057135629060264_r8*t11+2875.684484702232084974429_r8*t19&
            -1094.638706632625876620833_r8*t21+285.3669009128752146814545_r8*t23
        v(10)=-80675.59761993737478733866_r8*t2-1745.835002504339507610403_r8*t18&
             +911.9476960925337440830110_r8*t20-345.5461064752209510674120_r8*t22&
             -14.31536396632473462576254_r8*t26+89.78156942102429575363385_r8*t24&
             -629192.8778816414355475387_r8*t4-1163750.663354504353033275_r8*t6&
             -560765.7238305186542428340_r8*t8-50992.88361699267893018953_r8*t10&
             -2718.991768401117794132920_r8*t14+2481.071341994842760629068_r8*t12&
             +2498.827701212957733409790_r8*t16+1.056484425559021005591331_r8*t28&
             -1118.965893570671438005353_r8
        v(11)=93680.92053187578373908477_r8*t+2244650.940705589178033889_r8*t3+&
             10624527.77751646968478233_r8*t5+14380781.82932476042057132_r8*t7&
             +5381368.427872546425567057_r8*t9-5209.003030909290742154185_r8*t29&
             +25188.20188878916819094802_r8*t27-85058.16289349076452263097_r8*t25&
             -516900.3318170604402267402_r8*t13-727425.9526984207442109454_r8*t17&
             +687103.9096254180006068751_r8*t15+718764.2117855061019497234_r8*t11&
             +612373.6241206691026703661_r8*t19-408052.1263148301297946975_r8*t21&
             +212774.0620423199828384262_r8*t23-40.72442779420311834782560_r8*t33&
             +671.9530586043514527391227_r8*t31
        v(12)=-4666637.427699311288632484_r8*t2-236576.3442879792221912621_r8*t18&
             +197095.9817518425408469117_r8*t20-130314.1457376034450561984_r8*t22&
             -26890.53644971165967586539_r8*t26+67556.99854178264515084666_r8*t24&
             -56407070.38512263316820071_r8*t4-181816409.6548781888410872_r8*t6&
             -188159903.2775835541339353_r8*t8-59729728.73424203726696958_r8*t10&
             -158804.8694791222203286819_r8*t14-3778065.140617472679475600_r8*t12&
             +226299.4976083004679687682_r8*t16+7937.202199806361407713609_r8*t28&
             +210.6216499981442527051607_r8*t32-1637.030496200573732230987_r8*t30&
             -12.72638368568847448369550_r8*t34-45598.90544342769885089189_r8
       v(13)=5335496.618522339267627189_r8*t+183576348.0194920448725676_r8*t3&
             +1354176851.693565444970589_r8*t5+3188659860.422165550420750_r8*t7&
             +2627237695.599145398984883_r8*t9-8559691.708355332039465864_r8*t29&
             +20693853.48001686619594154_r8*t27-39924064.49900019391213026_r8*t25&
             +10733492.58653441250421101_r8*t13-72327111.51346468735941796_r8*t17&
             +51092798.08238289191185946_r8*t15+717509281.1168328033139334_r8*t11&
             +83750832.87636387843170424_r8*t19-79752878.22653923708993263_r8*t21&
             +62400837.71297169795531112_r8*t23+114641.0994472638480012588_r8*t35&
             -668909.4451481646322592866_r8*t33+2760161.897061770530576782_r8*t31&
             +635.4796101919361803875012_r8*t39-12391.85239874275551755627_r8*t37
       v(14)=-2660954.178544058526885368_r8-360244269.1046558724915919_r8*t2&
             -23905166.63593070422718515_r8*t18+27364022.18954367225325487_r8*t20&
             -25828929.45718207123186239_r8*t22-12760534.41026642324470170_r8*t26&
             +20063135.60588227514388882_r8*t24-6152008132.743840646273932_r8*t4&
             -30330498547.79309275046998_r8*t6-53977141029.83300025224457_r8*t8&
             -36163818479.23891864661839_r8*t10-425356977.9744537092972964_r8*t14&
             -8214360062.046573670722601_r8*t12+18349861.97314849584742475_r8*t16&
             +6582218.353750901954628165_r8*t28+871571.7482998278264261389_r8*t32&
             -2711831.071120029682914462_r8*t30-210623.3475234238287272353_r8*t34&
             -3882.383243516360102054890_r8*t38+36005.76814847969026409722_r8*t36&
             +198.5873781849800563710940_r8*t40
       v(15)=416127703.9194396571096487_r8*t+6945.991168730899603944917_r8&
             +510652.2760188471673646379_r8*t2+3506139.062413368413555602_r8*t18&
             -4424351.480098610189178738_r8*t20+4507540.556879283080170755_r8*t22&
             +2493813.214756674570007723_r8*t26-3725064.570331952034347715_r8*t24&
             +5707621.898105787921162792_r8*t4+19508496840.66961142955273_r8*t3&
             +207861233568.6689810809570_r8*t5+761719377470.4876561258129_r8*t7&
             +1098113019146.597129034586_r8*t9+21556032.59225403242059432_r8*t6&
             +35040737.47008682938237320_r8*t8+24773160.94128538674006004_r8*t10&
             +1461033.332711555204742717_r8*t14+5971318.466177208703908447_r8*t12&
             -2216567.774820235631058887_r8*t16-1343045.511157374057238458_r8*t28&
             -190443.7698884165190432579_r8*t32+574046.2785186452537735790_r8*t30&
             +47314.98337952692225675031_r8*t34+913.5019396509082593070329_r8*t38&
             -8288.707145438840911983684_r8*t36-47.66097076439521352906259_r8*t40&
             -4689057.699681904219581699_r8*t29+8536717.391664964414277814_r8*t27&
             -12447499.69848944884749613_r8*t25+126127821765.3019730239362_r8*t13&
             -21305702.49857355575356540_r8*t17+5803290890.616556866441806_r8*t15&
             +627653101625.1730316429113_r8*t11+10390854.10380408875907143_r8*t19&
             -13716953.18981538673077002_r8*t21+14601530.98572140453336727_r8*t23&
             +172618.3737117023372074619_r8*t35-686012.9395825390193077566_r8*t33&
             +2038004.663813441018392810_r8*t31-178.7286403664820507339848_r8*t41&
             +3400.808851417783465354987_r8*t39-30574.26516172012488452956_r8*t37
        t2p1=1.0_r8-t2
        mu4=mu2*mu2
        tts=t2p1*t2p1*t2p1
        ttss=sqrt(tts)
        eta=0.5_r8*(acos(t)-t*sqrt(1.0_r8-t2))
        ! ------------------------------
        ! Computation of G(mu) and S(mu)
        ! G(mu)=1/gmus
        ! -----------------------------
        DO k=0,25
          IF (k==0) THEN
            mus=1.0_r8
            gmus=0.0_r8
          ELSE
            mus=mus*mu2
          ENDIF
          gmus=gmus+gk(k)/mus
        ENDDO
        DO k=0,12
          IF (k==0) THEN
            mus=1.0_r8
            smus=0.0_r8
          ELSE
            mus=mus*mu4
          ENDIF
          smus=smus+wk(k)/mus
        ENDDO
        maxi=1.0_r8
        k=0
        DO WHILE ((2*k+1<=15).AND.(maxi>eps))
          IF (k==0) THEN
            mus1=1.0_r8
            mus2=mu2*ttss
            s1=0.0_r8
            s2=0.0_r8
            s3=0.0_r8
            s4=0.0_r8
          ELSE
            mus1=-mus1*mu4*tts
            mus2=-mus2*mu4*tts
          ENDIF
          itwok=2*k
          delta1=u(itwok)/mus1
          delta2=u(itwok+1)/mus2
          delta3=v(itwok)/mus1
          delta4=v(itwok+1)/mus2
          s1=s1+delta1
          s2=s2+delta2
          s3=s3+delta3
          s4=s4+delta4
          e1=abs(delta1/s1)
          e2=abs(delta2/s2)
          e3=abs(delta3/s3)
          e4=abs(delta4/s4)
          k=k+1
          maxi=max(e1,e2,e3,e4)
        ENDDO        
        xargu=abs(x*x*0.25_r8+a)
        facu=xpowy(xargu,0.25_r8)
        argu3=2.0_r8*a*eta+pi*0.25_r8
        uaxx(1)=sqrt2/(gmus*facu)*(cos(argu3)*s1+sin(argu3)*s2)
        uaxx(2)=sqrt2*facu/gmus*(-sin(argu3)*s3+cos(argu3)*s4)
        IF (abs(a)>100) THEN
          vaxx(1)=gmus/(sqrtpi*smus*facu)*&
                  (sin(argu3)*s1-cos(argu3)*s2)
          vaxx(2)=facu*gmus/(sqrtpi*smus)*&
                 (cos(argu3)*s3+sin(argu3)*s4)
        ELSE
          beta=1.0_r8/(sqrttwopi*aux9(abs(a)))
          vaxx(1)=sqrt2/facu*beta/gmus*&
                   (sin(argu3)*s1-cos(argu3)*s2)
          vaxx(2)=sqrt2*facu*beta/gmus*&
                   (cos(argu3)*s3+sin(argu3)*s4)   
        ENDIF
        IF (mode==0) THEN
           gammaa2=exp(a*0.5_r8*(log(abs(a))-1.0_r8))
           uaxx(1)=uaxx(1)/gammaa2
           uaxx(2)=uaxx(2)/gammaa2
           vaxx(1)=vaxx(1)*gammaa2
           vaxx(2)=vaxx(2)*gammaa2
        ENDIF
      ENDIF
    ELSE
      xargu=abs(x*x*0.25_r8+a)
      IF (mode==0) THEN
      ! ------------------------------------------------
      ! Check for possible overflow/underflow problems
      ! ------------------------------------------------
        IF (a>0) THEN
        ! For a>0
          IF (a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8 > log(giant)) THEN
            ierr=1
            uaxx(1)=dwarf
            uaxx(2)=dwarf
            vaxx(1)=giant 
            vaxx(2)=giant
          ENDIF
        ELSE
          ffa=a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8
          IF ((ffa> log(giant)).OR.(ffa<log(dwarf))) THEN
            ierr=1
            uaxx(1)=giant
            uaxx(2)=giant
            vaxx(1)=dwarf 
            vaxx(2)=dwarf
          ENDIF
        ENDIF
      ENDIF
      IF (ierr==0) THEN
        x2=x*x
        tau=-a*0.5_r8/(x2*0.25_r8+a+x*0.5_r8*sqrt(x2*0.25_r8+a))
        tau2=tau*tau; tau3=tau2*tau; tau4=tau3*tau; tau5=tau4*tau;
        tau6=tau5*tau; tau7=tau6*tau; tau8=tau7*tau; tau9=tau8*tau;
        tau10=tau9*tau; tau11=tau10*tau; tau12=tau11*tau; tau13=tau12*tau;
        tau14=tau13*tau; tau15=tau14*tau; tau16=tau15*tau;
        tau17=tau16*tau; tau18=tau17*tau; tau19=tau18*tau; 
        tau20=tau19*tau; tau21=tau20*tau; tau22=tau21*tau;
        tau23=tau22*tau; tau24=tau23*tau; tau25=tau24*tau;
        tau26=tau25*tau; tau27=tau26*tau; tau28=tau27*tau;
        tau29=tau28*tau; tau30=tau29*tau; tau31=tau30*tau;
        tau32=tau31*tau; tau33=tau32*tau; tau34=tau33*tau;
        tau35=tau34*tau; tau36=tau35*tau; 
        phi(0)=1.0_r8
        phi(1)=-1.666666666666666666666667_r8*tau3-2.5_r8*tau2-.75_r8*tau
        phi(2)=21.38888888888888888888889_r8*tau6+64.16666666666666666666667_r8*tau5&
              +67.375_r8*tau4+27.875_r8*tau3+3.28125_r8*tau2
        phi(3)=-525.2160493827160493827160_r8*tau9-2363.472222222222222222222_r8*tau8&
              -4254.25_r8*tau7-3861.229166666666666666667_r8*tau6-1814.5125_r8*tau5&
              -396.328125_r8*tau4-27.0703125_r8*tau3
        phi(4)=19126.61779835390946502058_r8*tau12+114759.7067901234567901235_r8*tau11&
              +292637.2523148148148148148_r8*tau10+411244.1666666666666666667_r8*tau9&
              +344694.4244791666666666667_r8*tau8+173584.125_r8*tau7&
              +49877.23984375_r8*tau6+7077.45703125_r8*tau5+329.91943359375_r8*tau4
        phi(5)=-5328.1988525390625_r8*tau5-1496065.721609933035714286_r8*tau7&
              -153266.177490234375_r8*tau6-7447882.350390625_r8*tau8&
              -22100856.07764756944444444_r8*tau9-42015571.25390625_r8*tau10&
              -52708852.66637731481481481_r8*tau11-22880216.54128086419753086_r8*tau13&
              -43565653.69020061728395062_r8*tau12-6933398.951903292181069959_r8*tau14&
              -924453.1935871056241426612_r8*tau15
        phi(6)=107230.0019073486328125_r8*tau6+3911743.44635009765625_r8*tau7&
              +49592558.34827532087053571_r8*tau8+327327256.2444545200892857_r8*tau9&
              +3526512665.739431423611111_r8*tau11+1320314667.667214471726190_r8*tau10&
              +6506075517.774877025462963_r8*tau12+7726873009.8828125_r8*tau14&
              +8449099135.763165509259259_r8*tau13+4881459532.087512860082305_r8*tau15&
              +501515857.5210048010973937_r8*tau17+2031139222.960069444444444_r8*tau16&
              +55723984.16900053345526597_r8*tau18
        phi(7)=-2585008.974552154541015625_r8*tau7-115109601.2886428833007812_r8*tau8&
              -1814958496.40869140625_r8*tau9-78228459666.49902913411458_r8*tau11&
              -15118535470.70343279157366_r8*tau10-273482436440.6177463107639_r8*tau12&
              -1223138659577.850206163194_r8*tau14-677676841014.5960557725694_r8*tau13&
              -1629074225394.340964988426_r8*tau15-1152270884857.328880529835_r8*tau17&
              -1602934313802.891348379630_r8*tau16-588917533590.4779128086420_r8*tau18&
              -42266641992.18690462581923_r8*tau20-202879881562.4971422039323_r8*tau19&
              -4025394475.446371869125641_r8*tau21
        phi(8)=338971759786.5465644792884_r8*tau24+72622595.87882459163665771_r8*tau8&
              +3839296326.682519912719727_r8*tau9+741010050728.2100511823382_r8*tau11&
              +72972869796.83574485778809_r8*tau10+4728818470068.397703588577_r8*tau12&
              +65097477007191.15089355469_r8*tau14+20667683341938.66336140951_r8*tau13&
              +152479523376491.2169460720_r8*tau15+367482071900019.2868381076_r8*tau17&
              +270645546827389.1750793457_r8*tau16+382220257297374.9787138846_r8*tau18&
              +179020321106247.5114691394_r8*tau20+302446999586890.6026957948_r8*tau19&
              +76811168492401.26177650130_r8*tau21+4067661117438.558773751461_r8*tau23&
              +22575519201784.00119432061_r8*tau22
        phi(9)=-10820487030026246.15802560_r8*tau24-2329974951.112288981676102_r8*tau9&
              -3204421618368.187572338364_r8*tau11-143186772419.2560920119286_r8*tau10&
              -38658804680937.31829602378_r8*tau12-1564351256647306.276325171_r8*tau14&
              -295696972623483.3523212160_r8*tau13-6033601440981970.440630913_r8*tau15&
              -39340249994654496.14279334_r8*tau17-17548357451443068.50438477_r8*tau16&
              -68913028780415733.17775943_r8*tau18-103213237253410851.2086602_r8*tau20&
              -95001751263219690.44791428_r8*tau19-87982837747558695.44900832_r8*tau21&
              -29268899384432143.34313245_r8*tau23-58193404452593856.87994868_r8*tau22&
              -440154830082830.7139763560_r8*tau26-2772975429521833.498051043_r8*tau25&
              -32604061487617.08992417452_r8*tau27
        phi(10)=20356704334716668369.03719_r8*tau24+3526672650910581.893464877_r8*tau30&
              +5906050907256.692442968488_r8*tau11+84053846361.37582501396537_r8*tau10&
              +152794280429199.7814958678_r8*tau12+19259274438496755.70484603_r8*tau14&
              +2146857055341595.890335684_r8*tau13+120401447638376916.9311805_r8*tau15&
              +1937739604392837793.628372_r8*tau17+553463759919030815.3303338_r8*tau16&
              +5291532561100032069.755284_r8*tau18+19849372125387206085.07073_r8*tau20&
              +11452797885661395365.50547_r8*tau19+27701776760098776464.83803_r8*tau21&
              +28202657288836392386.51597_r8*tau23+31172857703295082254.68334_r8*tau22&
              +5055359283148227235.199949_r8*tau26+11563080908824599399.17167_r8*tau25&
              +1641667477501437855.453314_r8*tau27+52900089763658728.40197315_r8*tau29&
              +372945632833794035.2339107_r8*tau28
        phi(11)=-10421163556087996323925.44_r8*tau24-267313411091476194892.6407_r8*tau30&
              -3367884798525.126806809567_r8*tau11-266994698245891.9833631720_r8*tau12&
              -126736372928041779.7376528_r8*tau14-7868339035378478.878049733_r8*tau13&
              -1310856141939127651.253943_r8*tau15-51013867921762764148.51193_r8*tau17&
              -9505095396241888043.965687_r8*tau16-210101665156711621956.3567_r8*tau18&
              -1766396419169469155581.958_r8*tau20-680809067349466450984.1967_r8*tau19&
              -3714461044542513155277.256_r8*tau21-8996413456126446920199.301_r8*tau23&
              -6381319935521667021236.195_r8*tau22-7669111682905600049152.052_r8*tau26&
              -9899229839297131395937.954_r8*tau25-4798668661333286324214.086_r8*tau27&
              -923291756397106785011.7944_r8*tau29-2388134461189160429354.980_r8*tau28&
              -6991628530430228603.794118_r8*tau32-54534702537355783109.59412_r8*tau31&
              -423735062450316885.0784314_r8*tau33
        phi(12)=2265181499981701606442225.0_r8*tau24+1161521061286761581713892.0_r8*tau30&
              +148397423935013.3999250465_r8*tau12+435451384030404205.1201922_r8*tau14&
              +13128818258663568.25439738_r8*tau13+7937408485723114496.434515_r8*tau15&
              +773568243048000449287.6993_r8*tau17+93358506268269984746.98223_r8*tau16&
              +4769376606728078799223.182_r8*tau18+85555881143101015387051.24_r8*tau20&
              +22697751563763377591312.77_r8*tau19+260197218851336240835649.8_r8*tau21&
              +1328296900615223882220191.0_r8*tau23+647062456346154782161169.5_r8*tau22&
              +3822476128794570206783031.0_r8*tau26+3220716113222632840960304.0_r8*tau25&
              +3781841350960097251345735.0_r8*tau27+2103782484922223264132877.0_r8*tau29&
              +3106620169009037645045708.0_r8*tau28+178206675733161601307968.3_r8*tau32&
              +514292924731074141830449.8_r8*tau31+46546067778486203906898.96_r8*tau33&
              +1007853846038078711159.049_r8*tau35+8617150383625572980409.870_r8*tau34&
              +55991880335448817286.61384_r8*tau36
        chi(0)=1.0_r8
        chi(1)=2.333333333333333333333333_r8*tau3+3.5_r8*tau2+1.25_r8*tau
        chi(2)=-25.27777777777777777777778_r8*tau6-75.83333333333333333333333_r8*tau5&
              -79.95833333333333333333333_r8*tau4-33.625_r8*tau3-4.21875_r8*tau2
        chi(3)=587.0061728395061728395062_r8*tau9+2641.527777777777777777778_r8*tau8&
              +4759.027777777777777777778_r8*tau7+4330.520833333333333333333_r8*tau6&
              +2047.1125_r8*tau5+453.109375_r8*tau4+31.9921875*tau3
        chi(4)=-7859.49609375_r8*tau5-373.90869140625_r8*tau4-54833.22265625_r8*tau6&
              -189793.2208333333333333333_r8*tau7-375736.0005208333333333333_r8*tau8&
              -447530.4166666666666666667_r8*tau9-124738.8117283950617283951_r8*tau11&
              -318189.0131172839506172840_r8*tau10-20789.80195473251028806584_r8*tau12
        chi(5)=5889.0618896484375_r8*tau5+167081.378173828125_r8*tau6+&
              1618334.674874441964285714_r8*tau7&
              +8018883.473046875_r8*tau8+23725519.15985243055555556_r8*tau9&
              +56413247.88802083333333333_r8*tau11+45020597.45963541666666667_r8*tau10&
              +46595060.47183641975308642_r8*tau12+7411564.396862139917695473_r8*tau14&
              +24461987.83320473251028807_r8*tau13+988208.5862482853223593964_r8*tau15
        chi(6)=-116554.3498992919921875_r8*tau6-4209790.52947998046875_r8*tau7&
              -53046844.35861642020089286_r8*tau8-348738575.260791015625_r8*tau9&
              -3740348588.554578993055556_r8*tau11-1402928860.327645399305556_r8*tau10&
              -6891971112.248734085648148_r8*tau12-8173643819.459394290123457_r8*tau14&
              -8942654096.485291280864198_r8*tau13-5161879657.283629115226337_r8*tau15&
              -530173906.5222050754458162_r8*tau17-2147389212.053647976680384_r8*tau16&
              -58908211.83580056393842402_r8*tau18
        chi(7)=2776491.120815277099609375_r8*tau7+122708978.7942123413085938_r8*tau8&
              +1925312041.345973423549107_r8*tau9+82508725486.56345738002232_r8*tau11&
              +15983788916.83610055106027_r8*tau10+287952809783.2672276475694_r8*tau12&
              +1285107966580.498679832176_r8*tau14+712650296975.9644274450231_r8*tau13&
              +1710504064563.104781539352_r8*tau15+1208933528167.414560828189_r8*tau17&
              +1682287248047.253970550412_r8*tau16+617751113141.8105388374486_r8*tau18&
              +44328429406.43992436366408_r8*tau20+212787605947.7454370522786_r8*tau19&
              +4221755181.565707082253721_r8*tau21
        chi(8)=-353396089990.2293970103219_r8*tau24-77307924.64520037174224854_r8*tau8&
              -4062968417.365064620971680_r8*tau9-778781085218.4028647286551_r8*tau11&
              -76913846048.55869483947754_r8*tau10-4959356633780.345308757964_r8*tau12&
              -68082808844916.25546820747_r8*tau14-21641132068083.78743394640_r8*tau13&
              -159329371436569.0371229384_r8*tau15-383535732852395.0897533275_r8*tau17&
              -282612499203202.3440902145_r8*tau16-398767118978192.2197332605_r8*tau18&
              -186682859410003.5980701839_r8*tau20-315453595729638.2593476723_r8*tau19&
              -80088585766097.86195177945_r8*tau21-4240753079882.752764123863_r8*tau23&
              -23536984672244.36711526127_r8*tau22
       chi(9)=11229654792749015.22701717_r8*tau24+2463116376.890134066343307_r8*tau9&
              +3360186490494.481565139510_r8*tau11+150652583603.0128768086433_r8*tau10&
              +40437382791704.19147014618_r8*tau12+1630902159338730.110003549_r8*tau14&
              +308720787738688.0545137042_r8*tau13+6283277020019534.490934813_r8*tau15&
              +40904941176357109.48990767_r8*tau17+18258673913072200.31158583_r8*tau16&
              +71616091309247803.17800093_r8*tau18+107183846394401245.5057093_r8*tau20&
              +98687634471617789.29578043_r8*tau19+91345622223600985.09434355_r8*tau21&
              +30378295319311232.27795892_r8*tau23+60406850807305980.38153996_r8*tau22&
              +456764446312371.4956358411_r8*tau26+2877683806119897.731818695_r8*tau25&
              +33834403430546.03671376601_r8*tau27
       chi(10)=-21053597837396406673.0365_r8*tau24-3646220876365177.889853517_r8*tau30&
              -6184760673789.150449261069_r8*tau11-88364300020.93355963006616_r8*tau10&
              -159549747541818.1452511725_r8*tau12-20033854775307941.92047071_r8*tau14&
              -2236938226816488.952569230_r8*tau13-125079628228935167.7962316_r8*tau15&
              -2009277564547413368.390572_r8*tau17-574371677952986224.0321184_r8*tau16&
              -5483258169122136680.737239_r8*tau18-20548963029842481258.16041_r8*tau20&
              -11861476998383736648.72274_r8*tau19-28668438826053413980.12674_r8*tau21&
              -29172870122949767579.89601_r8*tau23-32252013255464105175.54592_r8*tau22&
              -5227324779598662081.545248_r8*tau26-11957474746245008204.2076_r8*tau25&
              -1697408519086965241.663463_r8*tau27-54693313145477668.34780275_r8*tau29&
              -385594378487915085.2699942_r8*tau28
       chi(11)=10749119360152201659095.30_r8*tau24+275549372270078767749.3455_r8*tau30&
              +3524530603107.690844335593_r8*tau11+278502821903198.7234563101_r8*tau12&
              +131627099146616557.0788647_r8*tau14+8187462234898173.302466749_r8*tau13&
              +1359407981935523997.779927_r8*tau15+52789252821617587071.26402_r8*tau17&
              +9845374898471241989.909180_r8*tau16+217243684855444804880.3725_r8*tau18&
              +1824328572800792450321.910_r8*tau20+703503512901659205452.3662_r8*tau19&
              +3834651403368051540134.680_r8*tau21+9281659895577361070060.539_r8*tau23&
              +6585513896734531838538.180_r8*tau22+7907844921517711462070.132_r8*tau26&
              +10208879014255827826055.24_r8*tau25+4947481686715095751572.891_r8*tau27&
              +951780479796384621502.0100_r8*tau29+2461979171794923193819.241_r8*tau28&
              +7206755562135774099.295476_r8*tau32+56213398719189220090.88340_r8*tau31&
              +436773064371865096.9269985_r8*tau33
       chi(12)=-2331551239435886577025091.0_r8*tau24-1194449360583363232233128.0_r8*tau30&
              -154712207932248.0126878145_r8*tau12-451767349373865480.3252141_r8*tau14&
              -13649617202765419.19403282_r8*tau13-8221062051469232071.896864_r8*tau15&
              -799266188640123245085.8777_r8*tau17-96565089437574073764.48836_r8*tau16&
              -4923418321891009149545.962_r8*tau18-88200546293836821900475.57_r8*tau20&
              -23413619423670496869405.09_r8*tau19-268104858375645682223269.2_r8*tau21&
              -1367612737884965288180114.0_r8*tau23-666447368381836165309217.2_r8*tau22&
              -3932747815369540922524828.0_r8*tau26-3314285837643629657858286.0_r8*tau25&
              -3890307474197961385216392.0_r8*tau27-2163604453746493613327353.0_r8*tau29&
              -3195297616564396799146471.0_r8*tau28-183238172201351976450037.1_r8*tau32&
              -528838703410831864182048.1_r8*tau31-47858665362437589747700.5_r8*tau33&
              -1036244095222249942459.304_r8*tau35-8859971761162727071404.065_r8*tau34&
              -57569116401236107914.40578_r8*tau36
       mu2=2.0_r8*a
       maxi=1.0_r8
       k=0
       DO WHILE ((k <= 12).AND.(maxi > eps))
         IF (k==0) THEN
           musp=1.0_r8
           musn=1.0_r8
           sphip=0.0_r8
           schip=0.0_r8
           sphin=0.0_r8
           schin=0.0_r8
         ELSE
           musn=musn*mu2
           musp=-musp*mu2
         ENDIF
         delta1=phi(k)/musp
         delta2=chi(k)/musp
         delta3=phi(k)/musn
         delta4=chi(k)/musn
         sphip=sphip+delta1
         schip=schip+delta2
         sphin=sphin+delta3
         schin=schin+delta4
         e1=abs(delta1/sphip)
         e2=abs(delta2/schip)
         e3=abs(delta3/sphin)
         e4=abs(delta4/schin)
         k=k+1
         maxi=max(e1,e2,e3,e4)
       ENDDO
       xargu=x*x*0.25_r8+a
       facu=xpowy(xargu,0.25_r8)
       uaxx(1)=1.0_r8/(facu*sqrt2)*sphip
       uaxx(2)=-facu/sqrt2*schip
       IF (a<0) THEN
         vaxx(1)=1.0_r8/(sqrt(pi)*facu)*sphin
         vaxx(2)=facu/sqrt(pi)*schin
         IF (mode==0) THEN
           ffa=exp(a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8)
           uaxx(1)=uaxx(1)/ffa
           uaxx(2)=uaxx(2)/ffa
           vaxx(1)=ffa*vaxx(1)
           vaxx(2)=ffa*vaxx(2)
         ENDIF
       ELSE
         sinpia=sin(pi*a)
        ! -----------------------------
        ! Computation of G(mu) and S(mu)
        ! G(mu)=1/gmus
        !-----------------------------
         mu2=mu*mu
         mu4=mu2*mu2
         DO k=0,25
           IF (k==0) THEN
             mus=1.0_r8
             gmus=0.0_r8
           ELSE
             mus=mus*mu2
           ENDIF
           gmus=gmus+gk(k)/mus
         ENDDO
         DO k=0,12
           IF (k==0) THEN
             mus=1.0_r8
             smus=0.0_r8
           ELSE
             mus=mus*mu4
           ENDIF
           smus=smus+wk(k)/mus
         ENDDO
         dl=-2.0_r8*a*(t*sqrt(t*t+1.0_r8)+log(t+sqrt(t*t+1.0_r8)))
         IF (dl<log(dwarf)) THEN
           dnew=0.0_r8
         ELSE         
           dnew=exp(dl)*sqrttwopi*smus/(gmus*gmus)
         ENDIF
         vaxx(1)=1.0_r8/facu*(sinpia*dnew*sphip/(sqrt2*pi)+sphin/sqrt(pi))
         vaxx(2)=-facu*(sinpia*dnew*schip/(sqrt2*pi)-schin/sqrt(pi))
         IF (mode==0) THEN
           ffa=exp(a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8)
           uaxx(1)=uaxx(1)/ffa
           uaxx(2)=uaxx(2)/ffa
           vaxx(1)=vaxx(1)*ffa
           vaxx(2)=vaxx(2)*ffa
         ENDIF
       ENDIF
     ENDIF
   ENDIF 
   END SUBROUTINE expaelem

   SUBROUTINE expax(a,x,mode,uax,vax,ierr)
   ! -------------------------------------------------
   ! Calculation of U(a,x), V(a,x) by using 
   ! asymptotic expansions Poincare-type
   ! -------------------------------------------------
   ! Inputs:
   !   a ,    order of the functions
   !   x ,    argument of the functions
   !   mode , mode=0, unscaled functions
   !          mode=1, scaled functions
   ! Outputs:
   !   uaxx , uaxx(1), U(a,x)
   !          uaxx(2), U'(a,x)
   !   vaxx , vaxx(1), V(a,x)
   !          vaxx(2), V'(a,x)
   !   ierr , error flag
   !          ierr=0, computation succesful
   !          ierr=1, overflow or/and underflow problems 
   ! ----------------------------------------------------
   USE Someconstants
   USE AiryFunction, ONLY: xpowy
   IMPLICIT NONE
   REAL(r8), INTENT(IN) :: a
   REAL(r8), INTENT(IN) :: x
   REAL(r8), INTENT(OUT):: uax(2)
   REAL(r8), INTENT(OUT):: vax(2)
   REAL(r8) eps, aph, amh, aphl, amhl, err1, err2, err1p, err2p, &
            x2, a2, x2k, b2, y1, y2, acof, bcof, facto1, facto2, &
            y1p, y2p, a2p, b2p, acofd, bcofd, lambda, facto1p, &
            facto2p, sqrtx, facu, afacu, phiax, over
   REAL(r8) xargu,ffa
   INTEGER, INTENT(IN) :: mode
   INTEGER, INTENT(OUT):: ierr
   INTEGER k,l
   INTEGER m
   eps=epss
   over=giant*1.e-5_r8
   m=0
   ierr=0
   aph=a+0.5_r8
   amh=a-0.5_r8
   x2=x*x
   sqrtx=sqrt(x)
   xargu=x*x*0.25_r8+a 
   ffa=a*log(x*0.5_r8+sqrt(xargu))+x*0.5_r8*sqrt(xargu)-a*0.5_r8
   IF (mode==0) THEN 
     ! -----------------------------------------------
     ! Checking possible underflow/overflow problems
     ! -----------------------------------------------    
      IF ((ffa> log(giant)).OR.(ffa<log(dwarf))) THEN
        uax(1)=dwarf
        uax(2)=dwarf
        vax(1)=giant
        vax(2)=giant
        ierr=1
      ENDIF
    ENDIF
    IF (ierr==0) THEN
     IF (mode==1)THEN
       lambda=1.0_r8+4.0_r8*a/(x*x)
       facu=1.0_r8+sqrt(lambda)
       afacu=a/(x*facu)
       phiax=xpowy(2.0_r8/facu,a)*exp(2.0_r8*afacu*afacu)
     ELSE
       phiax=exp(0.25_r8*x*x+a*log(x))
     ENDIF
     facto1=1.0_r8/(sqrtx*phiax)
     facto1p=-sqrtx*0.5_r8/phiax
     facto2=sqrt2/(sqrtpi*sqrtx)*phiax
     facto2p=sqrtx/(sqrtpi*sqrt2)*phiax
     err1=1.0_r8; err2=1.0_r8
     err1p=1.0_r8;
     err2p=1.0_r8
     y1=1.0_r8; 
     y2=1.0_r8;
     y1p=1.0_r8;
     y2p=1.0_r8;
     k=1
     a2=1.0_r8
     b2=1.0_r8
     x2k=1.0_r8  
     DO WHILE (((err1 > eps).OR.(err2 > eps).OR.(err1p >eps).OR.&
              (err2p > eps)).AND.(m<500))
       l=2*k
       aphl=aph+l
       amhl=amh-l
       a2p=-amhl*(amh+l-1.0_r8)*a2
       a2=-(aphl-2.0_r8)*(aphl-1.0_r8)*a2
       b2=(amhl+2.0_r8)*(amhl+1.0_r8)*b2
       IF (abs(aph-l)< dwarf) THEN
         b2p=0.0_r8 
       ELSE
         b2p=aphl/(aph-l)*b2
       ENDIF
       x2k=l*x2*x2k
       acof=a2/x2k; bcof=b2/x2k;
       acofd=a2p/x2k;
       bcofd=b2p/x2k;
       y1=y1+acof; y2=y2+bcof;
       y1p=y1p+acofd
       y2p=y2p+bcofd    
       k=k+1;  
       err1=abs(acof/y1); err2=abs(bcof/y2)
       err1p=abs(acofd/y1p); 
       err2p=abs(bcofd/y2p)
       IF (abs(a2)>over) m=500
       IF (abs(b2)>over) m=501
     ENDDO
     IF (m==500) THEN
       ierr=1
       uax(1)=giant
       uax(2)=giant
     ELSE
       uax(1)=facto1*y1
       uax(2)=facto1p*y1p;
     ENDIF
     IF (m==501) THEN
       ierr=1
       vax(1)=giant        
       vax(2)=giant
     ELSE
       vax(1)=facto2*y2   
       vax(2)=facto2p*y2p;
     ENDIF
   ENDIF     
   END SUBROUTINE expax
 END MODULE Parabolic 
