#include   <iostream>
#include   <cstdlib>
#include   <cmath>

/*
  The Author would be happy to help in any way he can
  in the use of these routins.

  Contact Address:
  CERI, KACST, P.O. Box 6086,
  Riyadh 11442,
  Saudi Arabia,
  Fax:966+1+4813764
  Email: alhargan@kacst.edu.sa

  Legal Matters:
  The program is here as is, no guarantees are given, stated or implied,
  you may use it at your own risk.

  last modified Oct 1998
*/
/*
  If you do not understand the routines do not worry, all you need is
  the use of the functions
  Coefficients(...) for first kind functions
  CoefficientSec(...) for second kind circumferential functions
  These take care of evaluating the MCNs of the other routines.

  For theoretical derivation see:
  The companion paper
  or
  F.A. Alhargan, "A complete method for the computations of Mathieu
  characteristic numbers of integer orders," SIAM Review, vol-38, No. 2,
  pp. 239-255, June 1996.

*/

const double pi=3.141592653589793;
inline double abs(double z)
     {return  fabs(z);}
inline double acosh(double z)
       {return log(z+sqrt(z*z-1));}
inline double asinh(double z)
       {return log(z+sqrt(z*z+1));}
inline double atanh(double z)
       {return 0.5*log((1+z)/(1-z));}

/*--------------------------------------------------------------*/
/*                     MCNR header file                         */
/*--------------------------------------------------------------*/


     double Estimatmcn(char typ,int n, double h);
     double Estimatmcnl(char typ,int n, double t);
     double nChain(char typ, int n, double t);


     double MCNRoot(double a, char typ, int n, double h, double accr, double &acc);
     double NewtonRaphsons(double a, char typ, int n, double h, double &acp);

     double Asympn(char typ, int n, double h);
     double Mcapn(int n, double h);
     double Approxn(char typ, int n, double  h);

   double Coefficients(char typ,int n, double h, double Bm[], int CDIM, double &mcnr, double &Norm);
   int CoefficientSec(char typ, int n, double h, double mcn, double Bm[], int CDIM, double Gm[], double &Norm);

/*========================================================================*/

/*--------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
double MCNRoot(double a, char typ, int n, double h,double accr, double &acc)
{
     /*
       inputs:
       a    = the estimated MCN
       typ  = type of MCN e or o
       n    = order of the MCN
       h    = the parameter
       accr = the required accuracy


       ouput:
       The function returns the improved MCN with accuracy acc
       acc  = the obtained accuracy


     */

    int i,aprxt;
    double df,tst,c,est;
    double MCNimpro,t1,act;

    act=accr;
    aprxt=0;
    MCNimpro=1.0*n*n;

   /* if h is sufficiently small then approx is accurate enough */
    if ( fabs(MCNimpro-a)<1e-16 && fabs(h)<0.25*n ) aprxt=1;


    if(aprxt)
      {
	 MCNimpro=a;
	 NewtonRaphsons(MCNimpro,typ,n,h,acc);
      }


     if(!aprxt)
     {
       c=1;MCNimpro=a;i=1;t1=0;
       t1=NewtonRaphsons(MCNimpro,typ,n,h,c);MCNimpro=t1;est=c;
       while (c>act && i<60 && acc!=0.0)  /* Start Newton iteration */
	{
	 t1=NewtonRaphsons(MCNimpro,typ,n,h,c);        //   c=abs(f(x)/f'(x))
	 acc=fabs(MCNimpro-t1);
	 MCNimpro=t1;
	 i=i+1;
	}
       if (c>1e-6 && c>act) std::cout <<"MCNRoot: Newton Raphsons has not converged Acc=" <<c;
       df=fabs(a-MCNimpro);
       tst=df/(fabs(a)+1);
       if (tst>0.15)  /* Assuming that the guess value is close to root */
	{
	 if (fabs(h)>fabs(n/5.0))
	 {
	   std::cout <<"\n";
	   std::cout << "MCNRoot:";
	   std::cout << "Warning  mcn not accurate enough, Type= " <<typ <<", n=" << n <<", h=" << h <<"\n";
	   std::cout << "Estimated Accrcy=" << est << ", Improved Accrcy=" << c <<"\n";
	   std::cout << "Estimated MCN ="<<a <<", Improved root="<< MCNimpro <<"\n";
	   std::cout << "Test=" <<tst <<"\n";
	   std::cout << "Estimated was taken to be the correct MCN as the Improved is in error";
	   std::cout <<"\n";
	 }
	 MCNimpro=a;
	}
       acc=c;
     }


     return MCNimpro;
}
/*-----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* This algorithem is from  J. Wimp p83 */

double NewtonRaphsons(double a, char typ, int n,double h, double &acp)
{
    /*
      This algorithem is from  J. Wimp p83
      also see equations (47)-(52).
      This computes the a-f(a)/f'(a) for the MCNs
      inputs:
       a=estimated MCN, typ='e' or 'o'
       n=order, h=parameter
      outputs:
       improved MCN
       accuracy of the computation acp
    */
    int p,ex=0;
    double  r,M=1.2*n+15;
    double  d,ih2,h2,Y0,Y1,Yr,Yr1,Yr2,Yd0,Yd1,Ydr,Ydr1,Ydr2,f,fd;
    double  nr,dm,mr;

    M=int(M);
    p=n%2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"NewtonRaphsons: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"\n" <<"NewtonRaphsons: n can not be zero for odd functions :STOP" <<"\n";
       exit(0);
     }

    if(h==0) {acp=0; return 1.0*n*n;}  // when h is zero, MCN=n*n

    d=a;
    if ( p==0 && typ=='o' ) ex=1;
    h2=0.5*h;h2=h2*h2;
    ih2=1.0/h2;
    Ydr2=0;Ydr1=0;Yr2=0;Yr1=1;
    for ( r=M ; r>=ex; r--)
     {
       mr=(2.0*r+2.0+p);
       mr=mr*mr;
       dm=(d-mr);
       dm=dm/h2;
       Yr=dm*Yr1-Yr2;
       Ydr=Yr1*ih2+dm*Ydr1-Ydr2;
       Yr2=Yr1;Yr1=Yr;
       Ydr2=Ydr1;Ydr1=Ydr;
       nr=fabs(Yr);
       if(nr>1e4)  /* Normalize to avoid overflow */
	   {
	     nr=1/nr;
	     Yr1=Yr1*nr;
	     Yr2=Yr2*nr;
	     Ydr1=Ydr1*nr;
	     Ydr2=Ydr2*nr;
	   }

     }

    Y0=Yr1;Y1=Yr2;
    Yd0=Ydr1;Yd1=Ydr2;

     if ( typ=='e'  )
      {
       if ( p==0 ) {f= 0.5*d*Y0-h2*Y1; fd=0.5*d*Yd0-h2*Yd1+0.5*Y0;}
       if ( p==1 ) {f= (d-1-h2)*Y0-h2*Y1; fd=(d-1-h2)*Yd0-h2*Yd1+Y0;}
      }
     else
      {
       if ( p==1 ) {f=(d-1+h2)*Y0-h2*Y1; fd=(d-1+h2)*Yd0-h2*Yd1+Y0;}
       if ( p==0 ) {f=(d-4)*Y0-h2*Y1; fd=(d-4)*Yd0-h2*Yd1+Y0;}
      }

   fd=fd+1e-90;  /* in case fd=zero so it does not crash */
   acp=fabs(f);
   r=fabs(f/fd);
   if (acp>r) acp=r;     //This value is returned f(x)/f'(x)

   return d-f/fd;
}
/*=======================================================================*/
/*-----------------------------------------------------------------------*/
double Mcapn(int n, double h)
{
   double R,n2,hp2,hp4,hp8,hp12,m3,m5;
    /*
       McLachlan (p17) Approximation for |h|~<1.4n for n>3
       Approximation valid for both an and bn.
       Note for higher values of n (n>100) the limit is |h|~<1.1n
       and may be less as n becomes large.

       Computes McLachlan estimate of MCN for n>3
       inputs:
       n=order, h=parameter
       output:
       R=estimated MCN


    */
   if(h==0) return 0;
   hp2=0.25*h*h;
   hp4=hp2*hp2;
   hp8=hp4*hp4;
   hp12=hp8*hp4;
   n2=n*n;
   m3=(n2-1)*(n2-1)*(n2-1);
   m5=m3*(n2-1)*(n2-1);
   R=n2+hp4/(2.0*(n2-1))+(5*n+7)*hp8/(32.0*m3*(n2-4));
   if(h>0.2*n) R=R+(9.0*n2*n2+58.0*n2+29)*hp12/(64.0*m5*(n2-4.0)*(n2-9));

   return R;

}

/*------------------------------------------------------------------*/
double Asympn(char typ, int n, double h)
{
/*
  These functions give the asymptotic values for a(n) and b(n)
  in the region  2n<=~ h <infinity  for a(n)   n>0, n=0  3<=h<infinity
  in the region  2n<=~ h <infinity for b(n)  n>1, n=1  3<=h<infinity
  See McLachlan p. 232
  inputs:
    typ='e' or 'o', n=order, h=parameter
  output:
    R=estimated MCN

*/

    double m,m2,m3,m4,m6;
    double tw10,tw14,tw16,tw20;
    double ih,ih2,ih3,ih4,ih5;
    double R;

    if ( typ=='e' ) m=2*n+1; else m=2*n-1;

    m2=m*m;m3=m2*m;m4=m*m3;m6=m4*m2;
    ih=1.0/h;
    ih2=ih*ih;
    ih3=ih2*ih;ih4=ih3*ih;ih5=ih4*ih;
    tw10=1024.0;tw10=1.0/tw10; // pow(2,-10)
    tw14=tw10*0.0625;
    tw16=tw14*0.25;
    tw20=tw16*0.0625;

    R= m*h - 0.5*h*h-0.125*(m2+1)-0.015625*(m3+3*m)*ih;
    R=R-(tw10)*(5*m4+34*m2+9)*(ih2);
    R=R-(tw14)*m*(33*m4+410*m2+405)*(ih3)-(tw16)*(63*m6+1260*m4+2943*m2+486)*(ih4);
    R=R-(tw20)*m*(527*m6+15617*m4+69001*m2+41607)*(ih5);

     return R;

}
/*--------------------------------------------------------------------*/

double Approxn(char typ, int n, double  h)
{
/*
   This function gives the approximate value for a(n) in the region
   of 0<h<0.8n for n>3, for n<3 the region is 0<h<4.
   Also it gives the approximate value for b(n) in the region
   of 0<h<0.8n for n>3, for n<3 the region is 0<h<4
   These are new approximations see
     F.A. Alhargan, "A complete method for the computations of Mathieu
     characteristic numbers of integer orders," SIAM Review, vol-38, No. 2,
     pp. 239-255, June 1996.
   They give good results for n<=3. However, for n>3 McLachlan's Appr.
   give better results and greater region for h.

  inputs:
    typ='e' or 'o', n=order, h=parameter
  output:
    R=approximate MCN

*/
   double  Q,R,T,S,c2,c1,c0,t1,M,th,w;

   if (h==0)  return n*n;
   if ( (typ=='e') )
     {
     if ( n==0 )
	   {
	     T=h*h*0.25;
	     R= 2-sqrt(4.0+2.0*T*T);  // equ (20b)
	     return R;
	   }

     if ( n==1)
	   {
	     T=h*h*0.25;
	     S=5.0*T*T-16.0*T+64.0;
	     Q=sqrt(S);
	     R=5.0+0.5*(T-Q);  // equ (21a)
	     R=R;
	     return  R;
	   }
     if ( n==2 )   // equ (20d)
	{
	 T=h*h*0.25;
	 c0=20.0*T*T;
	 c1=-48.0-3.0*T*T;
	 c2=-8;
       }

     if ( n==3 )   // equ (21b)
	{
	 T=h*h*0.25;
	 c0=T*T*(T+8.0);
	 c1=16.*T-128.0-2.*T*T;
	 c2=-T-8;
	}

     }

   if ( typ=='o' )
    {
     if ( n==1)
	  {
	    T=h*h*0.25;
	    R=5-0.5*(T+sqrt(64.0+16.0*T+5.0*T*T));  // equ (21c)
	    return R;
	  }

     if ( n==2 )
	  {
	    T=h*h*0.25;
	    R=10-sqrt(36+T*T);  // equ (21d)
	    return R;
	  }

     if ( n==3 )  // equ (21e)
      {
	T=h*h*0.25;
	c0=T*T*(8-T);
	c1=-16.0*T-128-2.0*T*T;
	c2=T-8;
      }

     }

     /* This approximation developed for general n>3 in similar manner
	as the above approximation, but was not mentioned in the paper
	as it is not accurate enough, development shown in my PhD thesis */
     if ( n>3 )
       {
	T=h*h*0.25;
	c0=8.0*T*T;
	c1=(16-2*T*T-16.0*n*n);
	c2=-8;
       }

     /*  Solution of the cubic x^3+c2 x^2+c1 x+c0=0  */
     Q=(3*c1-c2*c2)/9 ;R=(9*c2*c1-27*c0-2*c2*c2*c2)/54;
     w=Q*Q*Q+R*R;
     if ( w>=0 )
      {
       t1=R+sqrt(w);
       S=abs(t1)/t1*pow(abs(t1),1/3.);
       t1=R-sqrt(w);
       T=abs(t1)/t1*pow(abs(t1),1/3.);
      }
     else
      {
       M=pow(R*R+fabs(w),1/6.);
       th=atan(sqrt(abs(w))/R)/3.0+4.0*pi/3.0;
       S=M*cos(th);T=M*cos(th);
      }
     R=(3*S+3*T-c2)/3.0;

     if ( R<0 && h>0.2*n )
       {
	 R=Asympn(typ,n,h);  /* in this case return asymptotic value */
       }
     else
       R=n*n+fabs(R);

   return R;
}

/*====================================================================*/
/*---------------------------------------------------------------------*/
double Estimatmcn(char typ, int n, double h)
{
  /*
    This is the main routine for estimating MCNs
    it combines all other routines (asymptotic,approximate,Chaining) to
    obtain an estimate of an MCN for  the given inputs.

  inputs:
    typ='e' or 'o', n=order, h=parameter
  output:
    estimated MCN

  */

  double ps[4];   /* Contains positions for regions where
		 approximation is valid see Table 3 for n<4*/
  double MCNE,Yn,Xn,Av,t;


  if(h==0) return n*n;

  if (n<4)
   {
   if ( typ=='e' )
     { ps[0]=4;ps[1]=4;ps[2]=3.5;ps[3]=5;}
    else
     {ps[1]=4;ps[2]=4.5;ps[3]=5.25;}

    if ( fabs(h)<=ps[n]   )   return Approxn(typ,n,h);
    if ( fabs(h)>ps[n]    )    return Asympn(typ,n,h);
   }


   t=fabs(h)/(1.0*n);
   Xn=Mcapn(n,h);
   Yn=Asympn(typ,n,h);
   Av=0.5*(Xn+Yn);

  if(Yn>Xn && t>=2.0) return Yn;
  if(Xn>Yn && t<=1.0)  return Xn;

  if (n<19)
   {
    if (t<=1.4)  return Xn;
    if (t>1.4)   return Yn;
   }


  if(n>18 && n<70)        /* note even can be  to n<79 & odd  to n<122 */
   {
     if(Xn>Yn && t<1.4) return Xn;
     if(Xn>Yn && t>1.7) return Yn;
     return Av;
   }


   MCNE=nChain(typ,n,t);       /* if the above does not cover the rang,
				  then  chain in n-direction */

   return MCNE;

}
/*----------------------------------------------*/
double nChain(char typ, int n, double t)
{
   double rt,ac,co,cn,h;
   int r,ro,m,p;
   p=n%2;
   m=50+p;
   ro=m+10;
   rt=1.0*ro/m;
   co=rt*rt*Estimatmcnl(typ,m,t);
   cn=co;

  for(r=m+50;r<n;r+=50)  /* it works for steps up to 50 */
   {
     rt=1.0*r/ro;
     h=t*ro;
     cn=rt*rt*MCNRoot(co,typ,ro,h,0.001,ac);
     co=cn;
     ro=r;
   }
     rt=1.0*n/ro;
     h=t*ro;
     cn=rt*rt*MCNRoot(co,typ,ro,h,0.001,ac);


     return cn;

}
/*--------------------------------------------------------------------*/
double Estimatmcnl(char typ, int n, double t)
{
  double MCNE,Yn,Xn,Av,h;

  if(n<20) {std::cout <<"Estimatmcnl: n can not be less than 20, n="<<n; exit(0);}
  if(n>60) {std::cout <<"Estimatmcnl: n can not be greater than 60, n="<<n; exit(0);}

   h=t*n;

   Xn=Mcapn(n,h);
   Yn=Asympn(typ,n,h);
   Av=0.5*(Xn+Yn);

   if(Xn>Yn && t<1.6) MCNE=Xn;
   if(Xn>Yn && t>1.6) MCNE=Yn;

   if(t>1.7) MCNE=Yn;
   if(t<1.5) MCNE=Xn;

   if(t>=1.5 && t<=1.7) MCNE=Av;

   MCNE=MCNRoot(MCNE,typ,n,h,0.001,Av);

   return MCNE;

}
/*---------------------------------------------------------------------*/

/*=====================================================================*/
/*============== Coeffecients For Mathieu Functions  ==================*/
/*
  This calculates the coefficients acoording to the standard porcedure
  outlined by McLachlan, with some modifications.
  The MCN is estimated then improved using MCNRoot.
  Since coefficients are sensitive to the accuracy of the MCN
  The main causes of non-convergence are
  1) Inaccurrate estimation of MCN
  2) The MCN is very close to a singularity
  3) problem 2) is overcome using forward and backward computation
  The  printed warnings often not very serious,
       however results must be checked
  The coefficient are normalized according to Morse

  Input: typ,n,h
  Output: Cofficients and Normlization constant M and
  accuracy of computation is returned


*/

/*-------- Coeffecients For Even Mathieu functions Of odd Order  B(2m+1) ----*/

double Coefficients(char typ,int n, double h, double Bm[], int CDIM, double &mcnr, double &Norm)
{
/*
       inputs:
       typ  = type of the coefficiens 'e' or 'o'
       n    = order of  the function
       h    = the parameter
       CDIM = number of the coefficients array Bm[]


       ouput:
       An array of the computed coefficients  Bm[]
       The function returns the accuracy of the computation
       mcnr = Obtained Mathieu characteristic Number
       Norm = Normalization constant
    


*/
  int MD,m,p,ni,ex=0;
  double acc,aprtest;
  double *V,Nrm,v0,Ti,Vni,mcn,h2,ih2;


    p=n%2;
    h2=0.25*h*h;

    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"Coefficients: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     {
      // std::cout <<"\n" <<"Coefficients: n can not be zero for odd functions :STOP" <<"\n"; exit(0);
	  for(m=0;m<=MD;m++)  Bm[2*m+p]=0;
	  mcnr=0;
	  Norm=1e250;   // infinity
	 return 0;
     }

     V= new double [CDIM];
    if (typ=='o' && p==0)  {ex=1;Bm[0]=0;}    else ex=0;

    MD=CDIM-1;

    mcn=Estimatmcn(typ,n,h);      //Estimate MCN={an,bn}
    acc=1e-8;
    mcn=MCNRoot(mcn,typ,n,h,1e-16,acc); // Improve the Accuracy of MCN
    mcnr=mcn;

      if (abs(h)==0)
      {
	for(m=0;m<=MD;m++)
	 {Bm[2*m+p]=0;}
	if ( typ=='e' )
	  {
	     // make sure that n<2*CDIM+p;
	    if(n<2*MD+p) Bm[n]=1;
	    if ( n==0 ) Nrm=2*pi; else Nrm=pi;
	   }
	else
	  { Bm[n]=1./n;Nrm=pi*(Bm[n]*Bm[n]); }

       Norm=Nrm;
       delete[] V;
       return 1e-8;
     }

    aprtest=fabs(mcn-n*n);

    if(aprtest<1e-7 && abs(h)<0.4*(n+1))   // Use approximation
      {
	for(m=0;m<=MD;m++) {Bm[2*m+p]=0;}

	if ( typ=='e' )
	  {
	   Bm[n]=1;
	   if ( n==0 ) Nrm=2*pi; else Nrm=pi;
	  }
	else
	 {
	   Bm[n]=1.0/n;
	   Nrm=pi*Bm[n]*Bm[n];
	 }

       Norm=Nrm;
       delete[] V;
       return 1e-6;
     }

  ih2=1.0/h2;
  if (n==0)  /* the n=0 is a special case where forward is not needed */
     {
      V[MD]=-h2/((2*MD+2+p)*(2*MD+2+p)-mcn);
      for(m=MD;m>=2-p;m--)
	 {
	  Ti=1.0/((2.0*m+p)*(2.0*m+p));
	  V[m-1]=-Ti*h2/(1-mcn*Ti+Ti*h2*V[m]);
	 }
      V[0]=-.5*h2/(1-mcn/4+h2*V[1]/4);v0=mcn*ih2;
      acc=abs(V[0]/v0-1);
    }
  else
    {   /* Use Forward & Backward  */

       ni=(n-p)/2;    /* Work out the point where f&b computation should meet */
       if ( ni>=MD ) ni=MD;

		     /* initialized for forward computation */
       if ( typ=='e'  )
       {
	if ( p==0 ) {V[0]= mcn*ih2;V[1]=(mcn-4)*ih2-2/V[0];}
	if ( p==1 ) {V[0]=((mcn-1)*ih2-1);V[1]=(mcn-9)*ih2-1/V[0];}
       }

       if ( typ=='o'  )
       {
	if ( p==0 ) {V[0]=0;V[1]=(mcn-4)*ih2;}
	if ( p==1 ) {V[0]=(1+(mcn-1)*ih2);V[1]=(mcn-9)*ih2-1/V[0];}
       }

       for(m=2;m<=ni;m++)  /* Forward Computation */
	 {V[m]=(mcn-(2.0*m+p)*(2.0*m+p))*ih2-1.0/V[m-1];}

      Vni=V[ni];  /* Store last value of forward computation */

		 /* start backward computation */
      V[MD]=-h2/((2*MD+2+p)*(2*MD+2+p)-mcn);  /* The previous V assumed to be zero */
      for(m=MD;m>ni;m--)
	{
	 Ti=1.0/((2.0*m+p)*(2.0*m+p));
	 V[m-1]=-Ti*h2/(1-mcn*Ti+Ti*h2*V[m]);
	}

      acc=abs(V[ni]/Vni-1);     /*Check that forward & backward agree*/
    }

      Bm[p+2*ex]=1;
      for(m=0+ex;m<MD;m++)
       {Bm[2*m+2+p]=V[m]*Bm[2*m+p];}

  /*- Normalization is different for the even and the odd Coeff. ----*/

    if ( typ=='e'  )
     {
       Nrm=1;
       for(m=1;m<=MD;m++)
	 {Nrm=Nrm+(Bm[2*m+p]);}

       Bm[p]=(1.0/Nrm);
       if ( p==0 ) Nrm=2*Bm[0]*Bm[0]; else Nrm=Bm[1]*Bm[1];
     }

   if ( typ=='o'  )
    {
       if ( p==0 ) Nrm=2; else Nrm=1;
       for(m=2-p;m<=MD;m++)
	  {Nrm=Nrm+(2.0*m+p)*Bm[2*m+p];}

       Bm[2-p]=(1.0/Nrm);
       Nrm=(Bm[2-p])*(Bm[2-p]);
    }

    if (acc>1e-5){ std::cout <<"\n";
		 std::cout << "Coefficients: Warning Coeffecient  Have not converged" <<"\n";
		 std::cout << "type is " << typ;
		 std::cout << " n=" << n << " h=" <<h;
		 std::cout << " acc=" << acc <<"\n";
	       }
   /*-------------------------------*/

    for(m=1+ex;m<=MD;m++)
      {
	Bm[2*m+p]=Bm[2*m+p]*Bm[p+2*ex];
	Nrm=Nrm+Bm[2*m+p]*Bm[2*m+p];
      }

    Norm=pi*Nrm;   //Normalization Constant

    delete[] V;
    return acc;

}

/*--------------------------------------------------------------------*/

int CoefficientSec(char typ, int n, double h, double mcn, double Bm[], int CDIM, double Gm[], double &Norm)
{
/*
   Calculation of the Second Kind Mathieu Coefficients 
   Inputs:
   Bm(),mcn,n,h,typ
   Bm= Coefficients of the first kind Bem or Bom
   mcn= Mathieu Characteristic Number an or bn
   n=order h=paramtere , typ= 'e' for even typ='o' for odd

   Outputs:
   Gm()=Second Kind Mathieu Coefficients,  Dem() or Dom()
   Norm=Second Kind Normalization constant
   returns 1 if computation has converged or 0 if not
*/

  int i,s,p,r,sign=1,conv=1;
  double  *w,h2,ih2,th,Nrm,acc;

    p=n%2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"CoefficientSec: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     {
      // std::cout <<"\n" <<"CoefficientSec: n can not be zero for odd functions :STOP" <<"\n";exit(0);
	   for(r=0;r<=CDIM-1;r++)
	      Gm[2*r+p]=0;
       Norm=1e250;  // infinity
     }

    w= new double [2*CDIM+3];

   s=CDIM-1;h2=0.25*h*h;ih2=1.0/h2;

   if(typ=='o') sign=-1;
   w[2*s+p]=0;
   w[2*s+p-2]=-2.0*(2.0*s+p)*sign*Bm[2*s+p]*ih2;

   for(r=s-1;r>=1;r--)
     w[2*r+p-2]=((mcn-(2.0*r+p)*(2.0*r+p))*w[2*r+p]-2.0*(2.0*r+p)*sign*Bm[2*r+p])*ih2-w[2*r+2+p];

   if(typ=='e')
    {
     if ( p==0 )  th=((mcn-4.0)*w[2]*ih2-w[4])/(2.0*Bm[0])-2.0*mcn/(h2*h2);
     if ( p==1 )  th=((mcn-1+h2)*w[1]*ih2-w[3])/(2.0*Bm[1])-ih2;
    }

   if ( typ=='o' )
    {
     if ( p==0 )  {w[0]=(4.0*Bm[2]+(mcn-4.0)*w[2])/(2*h2)-0.5*w[4];th=(w[2]-mcn*w[0]*ih2)/Bm[2];}
     if ( p==1 )  th=(w[3]-(mcn-1-h2)*w[1]*ih2)/(2.0*Bm[1])-ih2;
    }

   for(r=0;r<=s;r++)
      Gm[2*r+p]=w[2*r+p]-th*Bm[2*r+p];

  if(typ=='e')
    {
     Gm[0]=0;
     Nrm=1;
     for(i=0;i<=s;i++)
      {Nrm=Nrm+(2.0*i+p)*Gm[2*i+p];}

     Nrm=1.0/Nrm;
    }

  if ( typ=='o' )
    {
     Nrm=Gm[p];
     for (i=1; i<=s; i++)
      { Nrm=Nrm+Gm[2*i+p];}

     Nrm=1/Nrm;
    }
    r=4;
    acc=(mcn-(2.0*r+p)*(2.0*r+p))*Gm[2*r+p]-h2*h2*(Gm[2*r+2+p]+Gm[2*r-2+p]);
    if(typ=='e') acc=acc-2.0*(2.0*r+p)*Bm[2*r+p];
	else acc=fabs(acc+2.0*(2.0*r+p)*Bm[2*r+p]);

    Norm=Nrm;
    if (acc>1e-5){ std::cout <<"\n";
		 std::cout << "CoefficientSec: Warning Coeffecient  Have not converged" <<"\n";
		 std::cout << "type is " << typ;
		 std::cout << " n=" << n << " h=" <<h;
		 std::cout << " acc=" << acc <<"\n";
		 conv=0;
	       }

    delete[] w;
    return conv;

}
/*==================================================================*/

/*====================================================================*/

