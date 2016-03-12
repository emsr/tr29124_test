ADGMBESSIMIN, BESSIMSE.  A CODE TO EVALUATE MODIFIED BESSEL FUNCTIONS   ADGM0000
1   BASED ON THE CONTINUED FRACTION METHOD.  J. SEGURA,                 ADGM0000
2   P. FERNANDEZ DE CORDOBA, YU. L. RATIS.                              ADGM0000
REF. IN COMP. PHYS. COMMUN. 105 (1997) 263                              ADGM0000
CPC FREE FORMAT DATA 735 CARDS                                                  
cc       bessimse.f                                                             
      SUBROUTINE bessimse (x,Nmax,bi,bk,nuevo)                                  
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
cc Modified Bessel functions of the first and second kinds (half-integral)  cc  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c                                                                            c  
cc INPUT:  x-->argument / Nmax-->maximum order of the Bessel functions      cc  
c                                                                            c  
cc OUTPUT: bi(l+1)-->mbf of the first kind of order l (array)               cc  
c          bk(l+1)-->mbf of the second kind of order l  (array)              c  
c        nuevo--> max. order of bf's calculated when bk(Nmax+1,x)>10**expma  c  
c                                                                            c  
c  PARAMETERS: eps->control of accuracy in the continued fraction            c  
c              mode->when mode=1 (0) exp(+-x) is (not) factored out          c  
c              expma-> overflow/underflow numbers = 10**(+-expma)            c  
c                                                                            c  
c CALLS: slimit(x,expma,mode,Nmax,nuevo)--> evaluates the max. order (nuevo) c  
c                                                                            c  
cc NOTE: this code is for real positive x and Nmax less than 2000           cc  
c        Double precission accuracy                                          c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
                                                                                
      IMPLICIT REAL*8 (a-h,o-z)                                                 
      INTEGER MODE                                                              
      DIMENSION bi(2001),bk(2001)                                               
                                                                                
       PARAMETER(pi=3.14159265358979323D0,mode=0,eps=1.D-16,                    
     + expma=300.)                                                              
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccc                                              
c   we use the code if Nmax>=2   c                                              
cccccccccccccccccccccccccccccccccc                                              
                                                                                
      nuevo = Nmax                                                              
      nmaxp=Nmax                                                                
      IF(Nmax.LT.2) Nmax=2                                                      
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccccccc                                          
c slimit evaluates the maximum order c                                          
c that can be calculated, avoiding   c                                          
c overflows.                         c                                          
c comment the following line         c                                          
c when not needed                    c                                          
cccccccccccccccccccccccccccccccccccccc                                          
                                                                                
      CALL slimit(x,expma,mode,Nmax,nuevo)                                      
                                                                                
                                                                                
                                                                                
      Nmax=nuevo                                                                
                                                                                
ccccccccccccccccccccccccccccccccc                                               
c   values of bk(1) and bk(2)   c                                               
ccccccccccccccccccccccccccccccccc                                               
                                                                                
      xinv = 1./x                                                               
      xinv2=xinv*xinv                                                           
      IF(MODE.eq.0)then                                                         
          xk0 = xinv*DEXP(-x)                                                   
          xk1 = (xinv + xinv2 )*DEXP(-x)                                        
      ELSE                                                                      
          xk0=xinv                                                              
          xk1= (xinv+xinv2)                                                     
      ENDIF                                                                     
      bk(1) = xk0*0.5D0*pi                                                      
      bk(2) = xk1*0.5D0*pi                                                      
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                
c   upward recurrence relations for bk(l+1) l= 2,3,4,...,Nmax  c                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                
                                                                                
      DO 7 np=2,Nmax                                                            
          n=np-1                                                                
          bk(np+1)=DFLOAT(n+n+1)*bk(np)*xinv+bk(n)                              
   7  CONTINUE                                                                  
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c   c.f.  ( Steed's algorithm ) to calculate h = bi(Nmax+1)/bi(Nmax)     c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
                                                                                
      n=Nmax                                                                    
      d=x/DFLOAT(n+n+1)                                                         
      del=d                                                                     
      h=del                                                                     
      b=DFLOAT(n+n+3)*xinv                                                      
   8  d=1.D0/(b+d)                                                              
      del=(b*d-1.D0)*del                                                        
      h=h+del                                                                   
      b=b+2.d0*xinv                                                             
      IF(ABS(del/h).GT.eps) GO TO 8                                             
                                                                                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     we evaluate bi(Nmax+1) and bi(Nmax) using :                         c     
c     the Wronskian : w{bi(Nmax),bi(Nmax)} =-0.5*pi/x**2                  c     
c     the known values of bi(Nmax+1) and bk(Nmax)                         c     
c     the value of h = bi(Nmax+1)/bi(Nmax)                                c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
                                                                                
      bi(Nmax)=0.5D0*pi/(x*x*(h*bk(Nmax)+bk(Nmax+1)))                           
      bi(Nmax+1)=h*bi(Nmax)                                                     
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c downward recurrence relations for bi(l+1) l= Nmax-2,Nmax-3,..,2,1,0  c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
                                                                                
      DO 9 i=2,Nmax                                                             
          np=Nmax+2-i                                                           
          n=np-1                                                                
          bi(n)=DFLOAT(n+n+1)*bi(np)*xinv+bi(np+1)                              
   9  CONTINUE                                                                  
      Nmax=nmaxp                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
                                                                                
       SUBROUTINE slimit(x,expma,mode,Nmax,nuevo)                               
       IMPLICIT REAL*8 (a-h,o-z)                                                
       INTEGER nuevo                                                            
cccccccccccccccccccccccccccccccccccccccccc                                      
c   limits to the value of Nmax = nuevo  c                                      
c   bk(nuevo+1,x)<10**expma              c                                      
cccccccccccccccccccccccccccccccccccccccccc                                      
      fcmax=expma*LOG(10.)                                                      
      fc=fcmax+1.                                                               
      fn=DFLOAT(Nmax)+0.5                                                       
      DO WHILE(fc.GT.fcmax)                                                     
          alfn=LOG(2.*fn)                                                       
          IF(mode.EQ.0)THEN                                                     
              fc=(fn-0.5)*alfn-fn-(fn+0.5)*LOG(x)                               
          ELSE                                                                  
              fnx=SQRT(fn*fn+x*x)                                               
              fc=fn*LOG(fn+fnx)-0.5*alfn-fnx-(fn+0.5)*LOG(x)+x                  
          ENDIF                                                                 
                                                                                
          fn=fn-3.                                                              
      ENDDO                                                                     
      fn=fn+3.                                                                  
      nuevo = INT(fn)                                                           
      RETURN                                                                    
      END                                                                       
c----------                                                                     
cc       bessimin.f                                                             
      SUBROUTINE bessimin(x,Nmax,bi,bk,nuevo)                                   
                                                                                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
cc   Modified Bessel Functions of the first and second kind (integer)      cc   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
cc INPUT:  x-->argument / Nmax--> maximum order of the  mbf's              cc   
c             eps-->control of the accuracy in the c.f.                     c   
c                                                                           c   
cc OUTPUT: bi(l+1)-->mbf of the first kind of order l (array)              cc   
c          bk(l+1)-->mbf of the second kind of order l (array               c   
c         nuevo-->maxi. order mbf's calculated when bk(Nmax+1,x)>10**expma  c   
c                                                                           c   
c  PARAMETERS: eps-->control of the accuracy in the c.f.                    c   
c              mode-->when mode=1 (0) exp(+-x) is (not) factored out        c   
c              expma-->overflow/undeflow numbers=10**(+-expma)              c   
c                                                                           c   
c  CALLS: limit(x,expma,mode,Nmax,nuevo)-> evaluate the max. order (nuevo)  c   
c          bessim01(x,mode,xk0,xk1)-> evaluates K(0) and K(1)               c   
c                                                                           c   
cc NOTE :this version of the code is for real x>0 and Nmax less than 2000  cc   
c        the precission is limited by the accuracy in the calculation of    c   
c        K(0) (bk(1)) and K(1) (bk(2)); prec= 2.d-7 for bessim01            c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
                                                                                
                                                                                
      IMPLICIT REAL*8 (a-h,o-z)                                                 
      INTEGER MODE                                                              
      DIMENSION bi(2001),bk(2001)                                               
      PARAMETER(pi=3.14159265358979323D0,mode=0,eps=1.D-8,                      
     + expma=300.D0)                                                            
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccc                                              
c   we use the code if Nmax>=2   c                                              
cccccccccccccccccccccccccccccccccc                                              
                                                                                
      nuevo = Nmax                                                              
      nmaxp=Nmax                                                                
      IF(Nmax.LT.2) Nmax=2                                                      
                                                                                
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccccccc                                          
c limit evaluates the maximum order  c                                          
c that can be calculated, avoiding   c                                          
c overflows.                         c                                          
c comment the following line         c                                          
c when not needed                    c                                          
cccccccccccccccccccccccccccccccccccccc                                          
                                                                                
                                                                                
                                                                                
      CALL limit(x,expma,mode,Nmax,nuevo)                                       
                                                                                
      Nmax=nuevo                                                                
                                                                                
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccc                                              
c   values of bk(1) and bk(2)    c                                              
cccccccccccccccccccccccccccccccccc                                              
                                                                                
      CALL  bessim01(x,mode,xk0,xk1)                                            
                                                                                
      bk(1) = xk0                                                               
      bk(2) = xk1                                                               
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc              
c   upward recurrence relations for  bk(l+1) l= 2,3,4,...,Nmax   c              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc              
                                                                                
      xinv=1.0D0/x                                                              
      DO 7 np=2,Nmax                                                            
          n=np-1                                                                
          bk(np+1)=DFLOAT(n+n)*bk(np)*xinv+bk(n)                                
   7  CONTINUE                                                                  
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                      
c   c.f. (Steed's) to calculate h = bi(Nmax+1)/bi(Nmax)  c                      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                      
                                                                                
      n=Nmax                                                                    
      d=x/DFLOAT(n+n)                                                           
      del=d                                                                     
      h=del                                                                     
      b=2.d0*DFLOAT(n+1)*xinv                                                   
  8   d=1.d0/(b+d)                                                              
      del=(b*d-1.D0)*del                                                        
      h=h+del                                                                   
      b=b+2.D0*xinv                                                             
      IF(ABS(del/h).GT.eps) GO TO 8                                             
                                                                                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               
c     we evaluate bi(Nmax+1) and bi(Nmax) using :               c               
c     the Wronskian : W{bi(Nmax),bk(Nmax)} = 2./pi/x            c               
c     the known values of bk(Nmax+1) and bk(Nmax)               c               
c     the value of h = bi(Nmax+1)/bi(Nmax)                      c               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               
                                                                                
      bi(Nmax+1)=xinv/(bk(Nmax)+bk(Nmax+1)/h)                                   
      bi(Nmax)=bi(Nmax+1)/h                                                     
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c   downward recurrence relations for bi(l+1) l= Nmax-2,Nmax-3,..,2,1,0    c    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
                                                                                
      DO 9 i=2,Nmax                                                             
          np=Nmax+2-i                                                           
          n=np-1                                                                
          bi(n)=DFLOAT(n+n)*bi(np)*xinv+bi(np+1)                                
   9  CONTINUE                                                                  
      Nmax=nmaxp                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
                                                                                
       SUBROUTINE limit(x,expma,mode,Nmax,nuevo)                                
       IMPLICIT REAL*8 (a-h,o-z)                                                
       INTEGER nuevo                                                            
cccccccccccccccccccccccccccccccccccccccccc                                      
c   limits to the value of Nmax = nuevo  c                                      
c   bk(nuevo+1,x)<10**expma              c                                      
cccccccccccccccccccccccccccccccccccccccccc                                      
      fcmax=expma*LOG(10.)                                                      
      fc=fcmax+1.                                                               
      fn=DFLOAT(Nmax)                                                           
      DO WHILE(fc.GT.fcmax)                                                     
          alfn=LOG(2.*fn)                                                       
          IF(mode.EQ.0)THEN                                                     
              fc=(fn-0.5)*alfn-fn-fn*LOG(x)                                     
          ELSE                                                                  
              fnx=SQRT(fn*fn+x*x)                                               
              fc=fn*LOG(fn+fnx)-0.5*alfn-fnx-fn*LOG(x)+x                        
          ENDIF                                                                 
                                                                                
          fn=fn-3.                                                              
      ENDDO                                                                     
      fn=fn+3.                                                                  
      nuevo = INT(fn)                                                           
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
                                                                                
      SUBROUTINE bessim01(x,mode,k0,k1)                                         
      IMPLICIT REAL*8 (a-h,o-z)                                                 
      REAL*8 k0,k1,i0,i1                                                        
      SAVE a1,a2,a3,a4,a5,a6,a7                                                 
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7                            
      SAVE c1,c2,c3,c4,c5,c6,c7                                                 
      SAVE r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7                            
      DATA a1,a2,a3,a4,a5,a6,a7/1.0D0,3.5156229D0,3.0899424D0,                  
     + 1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/                          
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566D0,0.42278420D0,                     
     + 0.23069756D0,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/                
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414D0,-0.7832358D-1,                     
     + 0.2189568D-1,-0.1062446D-1,0.587872D-2,-0.251540D-2,                     
     + 0.53208D-3/                                                              
      DATA c1,c2,c3,c4,c5,c6,c7/0.5D0,0.87890594D0,0.51498869D0,                
     + 0.15084934D0,0.2658733D-1,0.301532D-1,0.32411D-3/                        
      DATA r1,r2,r3,r4,r5,r6,r7/1.0D0,0.15443144D0,-0.67278579D0,               
     + -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/                     
      DATA s1,s2,s3,s4,s5,s6,s7/1.25331414D0,0.23498619D0,                      
     + -0.3655620D-1,0.1504268D-1,-0.780353D-2,0.325614D-2,                     
     + -0.68245D-3/                                                             
                                                                                
                                                                                
       IF(x.LE.2.0)THEN                                                         
           y=(x/3.75)**2                                                        
           i0=a1+y*(a2+y*(a3+y*(a4+y*(a5+y*(a6+y*a7)))))                        
           i1=x*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))                    
           y=x*x/4.0                                                            
           fac2=LOG(x/2.0)                                                      
           IF(mode.EQ.1)THEN                                                    
              fac3=EXP(x)                                                       
           ELSE                                                                 
              fac3=1.                                                           
           ENDIF                                                                
           k0=fac3*((-fac2*i0)+(p1+y*(p2+y*(p3+y*(p4+y*(                        
     +     p5+y*(p6+y*p7)))))))                                                 
           k1=fac3*((fac2*i1)+(1.0/x)*(r1+y*(r2+y*(r3+y*(                       
     +     r4+y*(r5+y*(r6+y*r7)))))))                                           
       ELSE                                                                     
           y=2.0/x                                                              
           IF(mode.EQ.0)THEN                                                    
               fac3=EXP(-x)/SQRT(x)                                             
           ELSE                                                                 
               fac3=1/SQRT(x)                                                   
           ENDIF                                                                
           k0=fac3*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(                              
     +     q6+y*q7))))))                                                        
           k1=fac3*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*(                              
     +     s6+y*s7))))))                                                        
       ENDIF                                                                    
                                                                                
                                                                                
       RETURN                                                                   
       END                                                                      
c----------                                                                     
cc      testpro.f                                                               
       PROGRAM test                                                             
       IMPLICIT REAL*8 (a-h,o-z)                                                
       DIMENSION bi(2001),bk(2001)                                              
                                                                                
                                                                                
       WRITE(4,29)'SEMI-INTEGRAL ORDER,NMAX=1000,MODE 0'                        
       Nmax=1000                                                                
       WRITE(4,30)                                                              
       WRITE(4,31)'x','n','i(n)','k(n)','i(0)'                                  
       DO x=1,5                                                                 
          CALL bessimse(x,Nmax,bi,bk,nuevo)                                     
          l=nuevo+1                                                             
          WRITE(4,32) x, l-1, bi(l),bk(l),bi(1)                                 
          l=21                                                                  
          WRITE(4,32) x, l-1, bi(l),bk(l),bi(1)                                 
       ENDDO                                                                    
                                                                                
                                                                                
       Nmax=2000                                                                
                                                                                
       WRITE(4,30)                                                              
       WRITE(4,30),'INT. ORDER, NMAX=2000, MODE 1'                              
       WRITE(4,30)                                                              
       WRITE(4,31)'x','n','I(n)','K(n)','I(0)'                                  
       DO ix=0,3,1                                                              
           x=10.D0**ix                                                          
           CALL bessimin(x,Nmax,bi,bk,nuevo)                                    
           l=nuevo+1                                                            
           WRITE(4,33) x, l-1,bi(l),bk(l),bi(1)                                 
           l=101                                                                
           WRITE(4,33) x, l-1,bi(l),bk(l),bi(1)                                 
       ENDDO                                                                    
                                                                                
  29   FORMAT (A37)                                                             
  30   FORMAT (A30)                                                             
  31   FORMAT (4X,A1,3X,A1,(10X,A4),2(12X,A4))                                  
  32   FORMAT (F6.0,I4,1x,3D17.11)                                              
  33   FORMAT (F6.0,I4,1x,3D16.9)                                               
                                                                                
       END                                                                      
                                                                                
                                                                                
      SUBROUTINE bessimse (x,Nmax,bi,bk,nuevo)                                  
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
cc Modified Bessel functions of the first and second kinds (half-integral)  cc  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c                                                                            c  
cc INPUT:  x-->argument / Nmax-->maximum order of the Bessel functions      cc  
c                                                                            c  
cc OUTPUT: bi(l+1)-->mbf of the first kind of order l (array)               cc  
c          bk(l+1)-->mbf of the second kind of order l  (array)              c  
c        nuevo--> max. order of bf's calculated when bk(Nmax+1,x)>10**expma  c  
c                                                                            c  
c  PARAMETERS: eps->control of accuracy in the continued fraction            c  
c              mode->when mode=1 (0) exp(+-x) is (not) factored out          c  
c              expma-> overflow/underflow numbers = 10**(+-expma)            c  
c                                                                            c  
c CALLS: slimit(x,expma,mode,Nmax,nuevo)--> evaluates the max. order (nuevo) c  
c                                                                            c  
cc NOTE: this code is for real positive x and Nmax less than 2000           cc  
c        Double precission accuracy                                          c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
                                                                                
      IMPLICIT REAL*8 (a-h,o-z)                                                 
      INTEGER MODE                                                              
      DIMENSION bi(2001),bk(2001)                                               
                                                                                
       PARAMETER(pi=3.14159265358979323D0,mode=0,eps=1.D-16,                    
     + expma=300.)                                                              
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccc                                              
c   we use the code if Nmax>=2   c                                              
cccccccccccccccccccccccccccccccccc                                              
                                                                                
      nuevo = Nmax                                                              
      nmaxp=Nmax                                                                
      IF(Nmax.LT.2) Nmax=2                                                      
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccccccc                                          
c slimit evaluates the maximum order c                                          
c that can be calculated, avoiding   c                                          
c overflows.                         c                                          
c comment the following line         c                                          
c when not needed                    c                                          
cccccccccccccccccccccccccccccccccccccc                                          
                                                                                
      CALL slimit(x,expma,mode,Nmax,nuevo)                                      
                                                                                
                                                                                
                                                                                
      Nmax=nuevo                                                                
                                                                                
ccccccccccccccccccccccccccccccccc                                               
c   values of bk(1) and bk(2)   c                                               
ccccccccccccccccccccccccccccccccc                                               
                                                                                
      xinv = 1./x                                                               
      xinv2=xinv*xinv                                                           
      IF(MODE.eq.0)then                                                         
          xk0 = xinv*DEXP(-x)                                                   
          xk1 = (xinv + xinv2 )*DEXP(-x)                                        
      ELSE                                                                      
          xk0=xinv                                                              
          xk1= (xinv+xinv2)                                                     
      ENDIF                                                                     
      bk(1) = xk0*0.5D0*pi                                                      
      bk(2) = xk1*0.5D0*pi                                                      
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                
c   upward recurrence relations for bk(l+1) l= 2,3,4,...,Nmax  c                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                
                                                                                
      DO 7 np=2,Nmax                                                            
          n=np-1                                                                
          bk(np+1)=DFLOAT(n+n+1)*bk(np)*xinv+bk(n)                              
   7  CONTINUE                                                                  
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c   c.f.  ( Steed's algorithm ) to calculate h = bi(Nmax+1)/bi(Nmax)     c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
                                                                                
      n=Nmax                                                                    
      d=x/DFLOAT(n+n+1)                                                         
      del=d                                                                     
      h=del                                                                     
      b=DFLOAT(n+n+3)*xinv                                                      
   8  d=1.D0/(b+d)                                                              
      del=(b*d-1.D0)*del                                                        
      h=h+del                                                                   
      b=b+2.d0*xinv                                                             
      IF(ABS(del/h).GT.eps) GO TO 8                                             
                                                                                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     we evaluate bi(Nmax+1) and bi(Nmax) using :                         c     
c     the Wronskian : w{bi(Nmax),bi(Nmax)} =-0.5*pi/x**2                  c     
c     the known values of bi(Nmax+1) and bk(Nmax)                         c     
c     the value of h = bi(Nmax+1)/bi(Nmax)                                c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
                                                                                
      bi(Nmax)=0.5D0*pi/(x*x*(h*bk(Nmax)+bk(Nmax+1)))                           
      bi(Nmax+1)=h*bi(Nmax)                                                     
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c downward recurrence relations for bi(l+1) l= Nmax-2,Nmax-3,..,2,1,0  c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
                                                                                
      DO 9 i=2,Nmax                                                             
          np=Nmax+2-i                                                           
          n=np-1                                                                
          bi(n)=DFLOAT(n+n+1)*bi(np)*xinv+bi(np+1)                              
   9  CONTINUE                                                                  
      Nmax=nmaxp                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
                                                                                
       SUBROUTINE slimit(x,expma,mode,Nmax,nuevo)                               
       IMPLICIT REAL*8 (a-h,o-z)                                                
       INTEGER nuevo                                                            
cccccccccccccccccccccccccccccccccccccccccc                                      
c   limits to the value of Nmax = nuevo  c                                      
c   bk(nuevo+1,x)<10**expma              c                                      
cccccccccccccccccccccccccccccccccccccccccc                                      
      fcmax=expma*LOG(10.)                                                      
      fc=fcmax+1.                                                               
      fn=DFLOAT(Nmax)+0.5                                                       
      DO WHILE(fc.GT.fcmax)                                                     
          alfn=LOG(2.*fn)                                                       
          IF(mode.EQ.0)THEN                                                     
              fc=(fn-0.5)*alfn-fn-(fn+0.5)*LOG(x)                               
          ELSE                                                                  
              fnx=SQRT(fn*fn+x*x)                                               
              fc=fn*LOG(fn+fnx)-0.5*alfn-fnx-(fn+0.5)*LOG(x)+x                  
          ENDIF                                                                 
                                                                                
          fn=fn-3.                                                              
      ENDDO                                                                     
      fn=fn+3.                                                                  
      nuevo = INT(fn)                                                           
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      SUBROUTINE bessimin(x,Nmax,bi,bk,nuevo)                                   
                                                                                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
cc   Modified Bessel Functions of the first and second kind (integer)      cc   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
cc INPUT:  x-->argument / Nmax--> maximum order of the  mbf's              cc   
c             eps-->control of the accuracy in the c.f.                     c   
c                                                                           c   
cc OUTPUT: bi(l+1)-->mbf of the first kind of order l (array)              cc   
c          bk(l+1)-->mbf of the second kind of order l (array               c   
c         nuevo-->maxi. order mbf's calculated when bk(Nmax+1,x)>10**expma  c   
c                                                                           c   
c  PARAMETERS: eps-->control of the accuracy in the c.f.                    c   
c              mode-->when mode=1 (0) exp(+-x) is (not) factored out        c   
c              expma-->overflow/undeflow numbers=10**(+-expma)              c   
c                                                                           c   
c  CALLS: limit(x,expma,mode,Nmax,nuevo)-> evaluate the max. order (nuevo)  c   
c          bessim01(x,mode,xk0,xk1)-> evaluates K(0) and K(1)               c   
c                                                                           c   
cc NOTE :this version of the code is for real x>0 and Nmax less than 2000  cc   
c        the precission is limited by the accuracy in the calculation of    c   
c        K(0) (bk(1)) and K(1) (bk(2)); prec= 2.d-7 for bessim01            c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
                                                                                
                                                                                
      IMPLICIT REAL*8 (a-h,o-z)                                                 
      INTEGER MODE                                                              
      DIMENSION bi(2001),bk(2001)                                               
      PARAMETER(pi=3.14159265358979323D0,mode=1,eps=1.D-8,                      
     + expma=300.D0)                                                            
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccc                                              
c   we use the code if Nmax>=2   c                                              
cccccccccccccccccccccccccccccccccc                                              
                                                                                
      nuevo = Nmax                                                              
      nmaxp=Nmax                                                                
      IF(Nmax.LT.2) Nmax=2                                                      
                                                                                
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccccccc                                          
c limit evaluates the maximum order  c                                          
c that can be calculated, avoiding   c                                          
c overflows.                         c                                          
c comment the following line         c                                          
c when not needed                    c                                          
cccccccccccccccccccccccccccccccccccccc                                          
                                                                                
                                                                                
                                                                                
      CALL limit(x,expma,mode,Nmax,nuevo)                                       
                                                                                
      Nmax=nuevo                                                                
                                                                                
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccc                                              
c   values of bk(1) and bk(2)    c                                              
cccccccccccccccccccccccccccccccccc                                              
                                                                                
      CALL  bessim01(x,mode,xk0,xk1)                                            
                                                                                
      bk(1) = xk0                                                               
      bk(2) = xk1                                                               
                                                                                
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc              
c   upward recurrence relations for  bk(l+1) l= 2,3,4,...,Nmax   c              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc              
                                                                                
      xinv=1.0D0/x                                                              
      DO 7 np=2,Nmax                                                            
          n=np-1                                                                
          bk(np+1)=DFLOAT(n+n)*bk(np)*xinv+bk(n)                                
   7  CONTINUE                                                                  
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                      
c   c.f. (Steed's) to calculate h = bi(Nmax+1)/bi(Nmax)  c                      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                      
                                                                                
      n=Nmax                                                                    
      d=x/DFLOAT(n+n)                                                           
      del=d                                                                     
      h=del                                                                     
      b=2.d0*DFLOAT(n+1)*xinv                                                   
  8   d=1.d0/(b+d)                                                              
      del=(b*d-1.D0)*del                                                        
      h=h+del                                                                   
      b=b+2.D0*xinv                                                             
      IF(ABS(del/h).GT.eps) GO TO 8                                             
                                                                                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               
c     we evaluate bi(Nmax+1) and bi(Nmax) using :               c               
c     the Wronskian : W{bi(Nmax),bk(Nmax)} = 2./pi/x            c               
c     the known values of bk(Nmax+1) and bk(Nmax)               c               
c     the value of h = bi(Nmax+1)/bi(Nmax)                      c               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc               
                                                                                
      bi(Nmax+1)=xinv/(bk(Nmax)+bk(Nmax+1)/h)                                   
      bi(Nmax)=bi(Nmax+1)/h                                                     
                                                                                
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c   downward recurrence relations for bi(l+1) l= Nmax-2,Nmax-3,..,2,1,0    c    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
                                                                                
      DO 9 i=2,Nmax                                                             
          np=Nmax+2-i                                                           
          n=np-1                                                                
          bi(n)=DFLOAT(n+n)*bi(np)*xinv+bi(np+1)                                
   9  CONTINUE                                                                  
      Nmax=nmaxp                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
                                                                                
       SUBROUTINE limit(x,expma,mode,Nmax,nuevo)                                
       IMPLICIT REAL*8 (a-h,o-z)                                                
       INTEGER nuevo                                                            
cccccccccccccccccccccccccccccccccccccccccc                                      
c   limits to the value of Nmax = nuevo  c                                      
c   bk(nuevo+1,x)<10**expma              c                                      
cccccccccccccccccccccccccccccccccccccccccc                                      
      fcmax=expma*LOG(10.)                                                      
      fc=fcmax+1.                                                               
      fn=DFLOAT(Nmax)                                                           
      DO WHILE(fc.GT.fcmax)                                                     
          alfn=LOG(2.*fn)                                                       
          IF(mode.EQ.0)THEN                                                     
              fc=(fn-0.5)*alfn-fn-fn*LOG(x)                                     
          ELSE                                                                  
              fnx=SQRT(fn*fn+x*x)                                               
              fc=fn*LOG(fn+fnx)-0.5*alfn-fnx-fn*LOG(x)+x                        
          ENDIF                                                                 
                                                                                
          fn=fn-3.                                                              
      ENDDO                                                                     
      fn=fn+3.                                                                  
      nuevo = INT(fn)                                                           
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
                                                                                
      SUBROUTINE bessim01(x,mode,k0,k1)                                         
      IMPLICIT REAL*8 (a-h,o-z)                                                 
      REAL*8 k0,k1,i0,i1                                                        
      SAVE a1,a2,a3,a4,a5,a6,a7                                                 
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7                            
      SAVE c1,c2,c3,c4,c5,c6,c7                                                 
      SAVE r1,r2,r3,r4,r5,r6,r7,s1,s2,s3,s4,s5,s6,s7                            
      DATA a1,a2,a3,a4,a5,a6,a7/1.0D0,3.5156229D0,3.0899424D0,                  
     + 1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/                          
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566D0,0.42278420D0,                     
     + 0.23069756D0,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/                
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414D0,-0.7832358D-1,                     
     + 0.2189568D-1,-0.1062446D-1,0.587872D-2,-0.251540D-2,                     
     + 0.53208D-3/                                                              
      DATA c1,c2,c3,c4,c5,c6,c7/0.5D0,0.87890594D0,0.51498869D0,                
     + 0.15084934D0,0.2658733D-1,0.301532D-1,0.32411D-3/                        
      DATA r1,r2,r3,r4,r5,r6,r7/1.0D0,0.15443144D0,-0.67278579D0,               
     + -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/                     
      DATA s1,s2,s3,s4,s5,s6,s7/1.25331414D0,0.23498619D0,                      
     + -0.3655620D-1,0.1504268D-1,-0.780353D-2,0.325614D-2,                     
     + -0.68245D-3/                                                             
                                                                                
                                                                                
       IF(x.LE.2.0)THEN                                                         
           y=(x/3.75)**2                                                        
           i0=a1+y*(a2+y*(a3+y*(a4+y*(a5+y*(a6+y*a7)))))                        
           i1=x*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))                    
           y=x*x/4.0                                                            
           fac2=LOG(x/2.0)                                                      
           IF(mode.EQ.1)THEN                                                    
              fac3=EXP(x)                                                       
           ELSE                                                                 
              fac3=1.                                                           
           ENDIF                                                                
           k0=fac3*((-fac2*i0)+(p1+y*(p2+y*(p3+y*(p4+y*(                        
     +     p5+y*(p6+y*p7)))))))                                                 
           k1=fac3*((fac2*i1)+(1.0/x)*(r1+y*(r2+y*(r3+y*(                       
     +     r4+y*(r5+y*(r6+y*r7)))))))                                           
       ELSE                                                                     
           y=2.0/x                                                              
           IF(mode.EQ.0)THEN                                                    
               fac3=EXP(-x)/SQRT(x)                                             
           ELSE                                                                 
               fac3=1/SQRT(x)                                                   
           ENDIF                                                                
           k0=fac3*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(                              
     +     q6+y*q7))))))                                                        
           k1=fac3*(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*(                              
     +     s6+y*s7))))))                                                        
       ENDIF                                                                    
                                                                                
                                                                                
       RETURN                                                                   
       END                                                                      
                                                                        ADGM****
