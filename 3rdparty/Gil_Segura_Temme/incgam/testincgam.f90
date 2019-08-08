  PROGRAM testincgam
! -------------------------------------------------------
! This is a test program for the Fortran 90 
! module incgamFI. This module includes the
! public routines:
!  a) Routine incgam: computation of the
!     incomplete gamma function rations P(a,x), Q(a,x)  
!  b) Routine invincgam: computation of x in the equations
!     P(a,x)=p, Q(a,x)=q, where p, q and a are provided.
! ------------------------------------------------------- 
  USE IncgamFI
  USE Someconstants
  IMPLICIT NONE
  REAL(r8) :: eps, a, x, p, q, xr, delta, d0, erxr
  INTEGER :: ierr,ierri,i,j
! ----------------------------------------------------
! 1) Test of the incgam routine: the function 
!    checkincgam checks the relative accuracy 
!    in the recursions 
!       Q(a+1,x)=Q(a,x)+x^a*exp(-x)/Gamma(a+1)
!       P(a+1,x)=P(a,x)-x^a*exp(-x)/Gamma(a+1)
! ----------------------------------------------------
  eps=0.5e-17_r8;
  d0=-1; 
  DO i=1,100
    a=i*0.1_r8;
    DO j=1,100 
      x=j*0.1_r8;
      delta=abs(checkincgam(a,x,eps)); 
      IF (delta>d0) THEN
        d0=delta;
      ENDIF
    ENDDO
  ENDDO 
  print*,'Max. error (direct computation): ',d0
! ----------------------------------------------------
! 2) Test of the invincgam routine: the composition of
!   the functions with their inverse is the identity.
! ----------------------------------------------------
  d0=-1; 
  DO i=1,100
    a=i*0.1_r8;
    DO j=1,100 
      x=j*0.1_r8;
      CALL incgam(a,x,p,q,ierr)
      CALL invincgam(a,p,q,xr,ierri)
      erxr=abs(1.0_r8-x/xr)
      IF (erxr>d0) THEN
        d0=erxr;
      ENDIF
    ENDDO
  ENDDO
  print*,'Max. error (inversion): ',d0
  END PROGRAM testincgam

     


