#include <iostream>
#include <cstdlib>
#include <cmath>

/*========================================================================*/
/*  		         Modified Mathieu Functions                       */
/*========================================================================*/

/*
  The Author would be happy to help in any way he can
  in the use of these routins.

  CERI, KACST, P.O. Box 6086, Riydh, Saudi Arabia,
  Fax:966+1+4813764
  Email: alhargan@kacst.edu.sa


*/



/*
  Modified Mathieu functions of real argument
  Qn=(Qe,Qo,Ee,Eo) and MZ=(Ie,Io,Ke,Ko)

  last modified  Oct. 1998

*/
#define CDIMs 60  /* Important: CIDMs must be greater than the order n, i.e. CDIMs=n+15 */

const double pi=3.141592653589793;

inline double abs(double z)
     {return  fabs(z);}
inline double acosh(double z)
       {return log(z+sqrt(z*z-1));}
inline double asinh(double z)
       {return log(z+sqrt(z*z+1));}
inline double atanh(double z)
       {return 0.5*log((1+z)/(1-z));}
inline int signf(int m)
{
  if(m%2==0) return 1;
   else return -1;
}

#define h_type double
	/*
	  For complex h replace double by complex, routins in this file
	  do not need modifing. However, you will not be able to run the
	  complex  part unless you have the complex routins for computing
	  Mathieu  Charactristic Numbers (MCNs) mcnc.cpp.
	  The one I have just developed  has limited range for h complex,
	  I am working in its extension.
	  Though it will not be  available  for sometime.
	*/

     /* routines required from file mcnr.cpp  */
   h_type Estimatmcn(char typ,int n, h_type h);
   h_type MCNRoot(h_type a, char typ, int n, h_type h, double accr, double &acc);
   double Coefficients(char typ,int n, h_type h, h_type Bm[], int CDIM, h_type &mcnr, h_type &Norm);
   int CoefficientSec(char typ, int n, h_type h, h_type mcn, h_type Bm[], int CDIM, h_type Dm[], h_type &Norm);

   /* routines required from file bsslr.cpp */
       double BesselJ(double z, int M, double bJn[]);
       double BesselJ(double z, int M, double bJn[], double JDn[]);

       double BesselY(double z, int M, double Yn[]);
       double BesselY(double z, int M, double Yn[], double YDn[]);
   /* routines required from file mbsslr.cpp */
	  void BesselI(double z, int M, double bIn[]);
	  void BesselI(double z, int M, double bIn[], double IDn[]);
	  void BesselK(double z, int M, double bKn[]);
	  void BesselK(double z, int M, double bKn[], double KDn[]);

   /* routines required from file mathur.cpp */
   h_type  MathuSn(char typ, int n, double v, h_type Bm[], int CDIM);
   h_type  MathuSDn(char typ, int n, double v, h_type Bm[], int CDIM);
   h_type  MathuFn(char typ, int n, double v, h_type Dm[], h_type Norm, h_type Bm[], int CDIM);
   h_type  MathuFDn(char typ, int n, double v, h_type Dm[], h_type Norm, h_type Bm[], int CDIM);

   /*   Routines available in this file  */
   h_type MathuMZn(char typ,int n, h_type h, double u, short kind);
   h_type MathuMZDn(char typ,int n, h_type h, double u, short kind);

   h_type MathuQn(char typ,int n, h_type h, double v, short kind);
   h_type MathuQDn(char typ,int n, h_type h, double v, short kind);

   h_type  MathuIn(char typ, int n, h_type Am[], h_type Ims[], h_type Zm[], int CDIM);
   h_type  MathuKn(char typ, int n, h_type Am[], h_type Kms[], h_type Zm[], int CDIM);
   h_type  MathuIDn(char typ, int n, h_type h, double u, h_type Am[], h_type Ims[], h_type IDms[], h_type Zm[], h_type ZDm[], int CDIM);
   h_type  MathuKDn(char typ, int n, h_type h, double u, h_type Am[], h_type Ims[], h_type IDms[], h_type Zm[], h_type ZDm[], int CDIM);

   h_type  MathuQn(char typ, int n, double v, h_type Am[], int CDIM);
   h_type  MathuQDn(char typ, int n, double v, h_type Am[], int CDIM);
   h_type  IQn(char typ,int ni, double vi, double ti, double ui, h_type Am[], int CDIM);
   h_type  MathuEn(char typ, int n, double v, h_type Cm[], h_type Norm, h_type Am[], int CDIM);
   h_type  MathuEDn(char typ, int n, double v, h_type Cm[], h_type Norm, h_type Am[], int CDIM);


   double MCoefficients(char typ,int n, h_type h, h_type Am[], int CDIM, h_type &mcnr, h_type &Norm);
   int MCoefficientSec(char typ, int n, h_type h, h_type mcn, h_type Am[], int CDIM, h_type Cm[], h_type &Norm);

/* ====================================================================*/
/* ================== Mathieu Functions ===============================*/
/*
  Important Convensions
1. CDIM is the number of terms summed in the series, which also determins
   the size of the arrays.
2. The value p is used to identify if n is even (p=0) or odd (p=1).
3. D in any function means differential e.g JDn(x)= d(Jn(x))/dx=J'n(x)
4. The series used are that which are convergent in all the z-plane.

5. The Coefficients of the Even functions is assumed be stored in
   one array Be() for both even and odd orders, similarly for odd functions Bo().
6. The Coefficients Be() and Bo() have to be calculated before calling the routine
   for Mathieu functions.
7. In the case of Jen and Jon the Bessel functions has to be calculated first
   Also for Yen and Yon.
8. Warning are not serious, but care needs to be taken as to validity of the
   results, in most cases the results are correct.
9. In the case of a serious error the programs haults with a messag

10.  Bm stands for Bem or Bom     1st kind coefficients
     Dm stands for Dem or Dom     2nd kind coefficients

11. All arrays are of size 2*CDIM, i.e Bm[0]...Bm[2*CDIM-1]
    some arrays are of size CDIM, e.g. Zm[0]...Zm[CDIM-1]

*/
/*------------------------------------------------------------------*/
/*---------------- Easy Functions -----------------------------------*/
/*
  When computing Mathieu functions, one usually would be requiring several
  different functions of the same order, in this case it is best to use
  the functions which take the coefficients and arrays of Bessel functions
  as inputs.
  On the other hand, one may want to compute a single function of one order
  in this case, 'The Easy Functions' coming shortly may be better and easily
  used.
*/
double MathuMZn(char typ,int n, double h, double u, short kind)
{
/*
   inputs:
   typ ='e' or 'o'
   n  = order
   h  = parameter
   u  = argument 
   kind = 1 or 2

   output:
   The function returns the value of the required modified radial function


*/

  int CDIM=CDIMs-1;
  double Am[2*CDIMs],Ims[2*CDIMs],Zm[2*CDIMs];
  double mcnr,Norm,Vr;
  double xp,xs;
  

  if(kind<1 || kind>2)
     {std::cout <<"MathuMZn: The kind has not been specified correctly:STOP"; exit(0);}


  xp=0.5*h*exp(u);
  xs=0.5*h*exp(-u);

  MCoefficients(typ,n,h,Am,CDIM,mcnr,Norm);

  BesselI(xs,2*CDIM,Ims);
  if(kind==1)
    { BesselI(xp,2*CDIM,Zm);
      Vr=MathuIn(typ,n,Am,Ims,Zm,CDIM);
    }
  if(kind==2)
     {  BesselK(xp,2*CDIM,Zm);
	Vr=MathuKn(typ,n,Am,Ims,Zm,CDIM);
     }

      return Vr;


}
double MathuMZDn(char typ,int n, double h, double u, short kind)
{
  int CDIM=CDIMs-1;
  double Am[2*CDIMs],Ims[2*CDIMs],Zm[2*CDIMs];
  double IDms[2*CDIMs],ZDm[2*CDIMs];
  double mcnr,Norm,Vr;
  double xp,xs;
  

  if(kind<1 || kind>2)
     {std::cout <<"MathuMZn: The kind has not been specified correctly:STOP"; exit(0);}


  xp=0.5*h*exp(u);
  xs=0.5*h*exp(-u);


  MCoefficients(typ,n,h,Am,CDIM,mcnr,Norm);

  BesselI(xs,2*CDIM,Ims,IDms);
  if(kind==1)
    { BesselI(xp,2*CDIM,Zm,ZDm);
       Vr=MathuIDn(typ,n,h,u,Am,Ims,IDms,Zm,ZDm,CDIM);
    }
  if(kind==2)
     {  BesselK(xp,2*CDIM,Zm,ZDm);
	Vr=MathuKDn(typ,n,h,u,Am,Ims,IDms,Zm,ZDm,CDIM);
     }

    return Vr;
}
/*------------------------------------------------------*/
double MathuQn(char typ,int n, double h, double v, short kind)
{
/*
   inputs:
   typ ='e' or 'o'
   n  = order
   h  = parameter
   v  = argument in radians
   kind = 1 or 2

   output:
   The function returns the value of the required modified circumferential function


*/

  int CDIM=CDIMs-1;
  h_type Am[2*CDIMs],Cm[2*CDIMs];
  h_type Vr,mcnr,Norm1,Norm2;


  MCoefficients(typ,n,h,Am,CDIM,mcnr,Norm1);

  if(kind==1) Vr=MathuQn(typ,n,v,Am,CDIM);
  if(kind==2)
  {
   MCoefficientSec(typ,n,h,mcnr,Am,CDIM,Cm,Norm2);
   Vr=MathuEn(typ,n,v,Cm,Norm2,Am,CDIM);
  }

   return Vr;

}
double MathuQDn(char typ,int n, double h, double v, short kind)
{
  int CDIM=CDIMs;
  h_type Am[2*CDIMs],Cm[2*CDIMs];
  h_type Vr,mcnr,Norm1,Norm2;


  MCoefficients(typ,n,h,Am,CDIM,mcnr,Norm1);

  if(kind==1) Vr=MathuQDn(typ,n,v,Am,CDIM);
  if(kind==2)
  {
   MCoefficientSec(typ,n,h,mcnr,Am,CDIM,Cm,Norm2);
   Vr=MathuEDn(typ,n,v,Cm,Norm2,Am,CDIM);
  }

   return Vr;

}
/*=======================================================================*/
/*====================================================================*/
/*-------- Integrals of Qen && Qon For Ports away from the center -----*/
/*
     This uses Romberg integration Routine
     The routine is general only the function call is chaged
*/
h_type IQn(char typ,int ni, double vi, double ti, double ui, h_type Am[], int CDIM)
{
  int K,M,N,N1,N2,TMP;
  double a,b,ACC,DX,R,X;
  double sui,sna,snb;
  h_type T[2][11],INTG,Z;

  a=vi-ti;b=vi+ti;    // Limits
  sui=sinh(ui);
  sui=sui*sui;
  sna=sin(a); sna=sna*sna;
  snb=sin(b);snb=snb*snb;
  INTG =MathuQn(typ,ni,a,Am,CDIM)*sqrt(sui+sna);
  INTG+=MathuQn(typ,ni,b,Am,CDIM)*sqrt(sui+snb);      //Function CALL
  T[0][0]=0.5*(b-a)*INTG;

  N1=0;N2=1;K=0;ACC=1;
  N=1;
  while((ACC>1e-6 && K<10) || K<4)    // Max NUMBER OF ITERATIONS
   {
    K++;
    N=2*N; DX=(b-a)/N;Z=0.0;

   for (R=1;R<=N;R+=2)
	{
	 X=a+(R*DX);
	 Z +=MathuQn(typ,ni,X,Am,CDIM)*pow(pow(sinh(ui),2)+pow(sin(X),2),.5);      // Function CALL
	}

   T[N2][0]=.5*T[N1][0]+DX*Z;   // VALUE OF T(K,0)

   for ( M=1 ; M<=K; M++)
    { T[N2][M]=(pow(4.0,M)*T[N2][M-1]-T[N1][M-1])/(pow(4.0,M)-1);}

   ACC=abs(T[N2][K]-T[N1][K-1]);
   TMP=N1;N1=N2;N2=TMP; //SWITCH ARRAY
  }

  INTG=T[N1][K];
  if (ACC>0.05) std::cout <<"INTERGRATION DID NOT CONVERGE";


   return INTG;

}
/*--------------------------------------------------------------------*/
/*====================================================================*/
/*---Modified Circumferential First kind Even && Odd Mathieu  Function ---*/
double  MathuQn(char typ, int n, double v, double Am[],  int CDIM)
{

 double MathuQ;
		/*
		  Make sure that the coefficients correspond correctly to
		  the type and the order.
		  i.e the coefficients passed are assumed to be
		  Modified Mathieu Coefficients;
		*/

      MathuQ=MathuSn(typ,n,v,Am,CDIM);

       return MathuQ;

}
/*--------------------------------------------------------------------*/
/*-Dif. of Modified Circumferential First kind Even && Odd Mathieu  Function-*/
double  MathuQDn(char typ, int n, double v, double Am[], int CDIM)
{

 double MathuQ;

	MathuQ=MathuSDn(typ,n,v,Am,CDIM);

       return MathuQ;

}
/*--- Circumferential Second Kind Mathieu Function (Non-periodic) ---*/
/*---------- Fen() && Fon()  ------*/
/*
   Bm is Bem or Bom 1st kind Mathieu coefficients
   Dm is Dem or Dom 2nd kind Mathieu coefficients
   Norm=second kind normalization constant
    not finished yet
*/
h_type  MathuEn(char typ, int n, double v, h_type Cm[], h_type Norm, h_type Am[], int CDIM)
{

 h_type MathuE;

	  /*   assuming that Cm,Am and Norm are the modified elements  */


     MathuE=MathuFn(typ,n,v,Cm,Norm,Am,CDIM);

    return MathuE;
}
/*----------------------------------------------------------------------*/
/*- Dif. of Circumferential Second Kind Mathieu Function (Non-periodic)--*/
/*------- Fe'n() and F'n() -------*/
h_type  MathuEDn(char typ, int n, double v, h_type Cm[], h_type Norm, h_type Am[], int CDIM)
{
 
 h_type MathuED;

	  /*   assuming that Cm,Am and Norm are the modified Coefficients  */

     MathuED=MathuFDn(typ,n,v,Cm,Norm,Am,CDIM);

    return MathuED;
}
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*------ Modified Radial Mahtieu function  -----------------*/
h_type  MathuIn(char typ, int n, h_type Am[], h_type Ims[], h_type Zm[], int CDIM)
{
 int  m,ni,p,s;
 h_type Sum,MathZ,Tmp;

	 /* assuming that  Am's are the Modified Mathieu Coefficients */

    p=n%2;
    ni=(n-p)/2;

    s=ni;
    if(s>=CDIM) s=CDIM-1;
    Sum=0;
    if (typ=='e')
    {
	    Tmp=1/(Am[2*s+p]+1e-50);
	    if (s==0 && n!=1) Tmp=0.5/(Am[2*s+p]+1e-50);

		for (m=s;m<CDIM-s-p;m++)
		    {Sum +=Am[2*m+p]*(Ims[m-s]*Zm[m+s+p]+Ims[m+s+p]*Zm[m-s]);}
		for (m=0;m<s;m++)
		    {Sum +=Am[2*m+p]*(Ims[s-m]*Zm[m+s+p]+Ims[m+s+p]*Zm[s-m]);}
     }

     if ( typ=='o' )
	{
	   Tmp=1/(Am[2*s+p]+1e-50);
		for (m=s;m<CDIM-s-p;m++)
		    {Sum +=Am[2*m+p]*(Ims[m-s]*Zm[m+s+p]-Ims[m+s+p]*Zm[m-s]);}
		for (m=0;m<s;m++)
		    {Sum +=Am[2*m+p]*(Ims[s-m]*Zm[m+s+p]-Ims[m+s+p]*Zm[s-m]);}
     }
	  MathZ=Sum*Tmp;
	   return MathZ;
}
/*------------------------------------------------------------------------------*/
h_type  MathuKn(char typ, int n, h_type Am[], h_type Ims[], h_type Zm[], int CDIM)
{
 int  m,ni,p,s;
 h_type Sum,MathZ,Tmp;
 double sign=1.0;

    p=n%2;
    ni=(n-p)/2;
    if(p==1) sign=-1.0;
    s=ni;
    if(s>=CDIM) s=CDIM-1;
    Sum=0;
    if ( typ=='e' )
    {
	    Tmp=1.0/(Am[2*s+p]+1e-50);
	    if (s==0 && n!=1) Tmp=0.5/(Am[2*s+p]+1e-50);

		for (m=s;m<CDIM-s-p;m++)
		    {Sum +=signf(m)*Am[2*m+p]*(Ims[m-s]*Zm[m+s+p]
						+sign*Ims[m+s+p]*Zm[m-s]);}
		for (m=0;m<s;m++)
		    {Sum +=signf(m)*Am[2*m+p]*(Ims[s-m]*Zm[m+s+p]
						+sign*Ims[m+s+p]*Zm[s-m]);}
     }

     if ( typ=='o' )
	{
	   Tmp=1/(Am[2*s+p]+1e-50);
		for (m=s;m<CDIM-s-p;m++)
		    {Sum +=signf(m)*Am[2*m+p]*(Ims[m-s]*Zm[m+s+p]
						-sign*Ims[m+s+p]*Zm[m-s]);}
		for (m=0;m<s;m++)
		    {Sum +=signf(m)*Am[2*m+p]*(Ims[s-m]*Zm[m+s+p]
						-sign*Ims[m+s+p]*Zm[s-m]);}
     }
	  MathZ=Sum*Tmp*signf(s);

   return MathZ;
}
/*--------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*- Differential of the Radial Second Kind Even and Odd Mathieu Function -*/
/*
  This Function takes the coefficients Am[], the Bessel functions
  Ims[]=Im(0.5h exp(-u)) and IDms[]. Zm[]=Zm(0.5h exp(u)) and ZDm[].
   where Zm=Im
  sums the series giving  Ze'n(h,cosh u) if typ='e'
  or Zo'n(h, cosh u) if typ='o'
  Note that differentiation is wrt u.
*/
h_type  MathuIDn(char typ, int n, h_type h, double u, h_type Am[], h_type Ims[], h_type IDms[], h_type Zm[], h_type ZDm[], int CDIM)
{

 int  m,ni,p,s;
 h_type Sum,MathZD,Tmp;
 double un,up;

	/* assuming that Am's are the Modified Mathieu Coefficients */
    un=exp(-u);
    up=exp(u);

    p=n%2;
    ni=(n-p)/2;

    s=ni;
    if(s>=CDIM) s=CDIM-1;
    Sum=0;
    if ( typ=='e' )
     {
	Tmp=1/(Am[2*s+p]+1e-50);
	if (s==0 && n!=1) Tmp=0.5/(Am[2*s+p]+1e-50);

	 for(m=s;m<CDIM-s-p;m++)
	  {Sum +=Am[2*m+p]*(up*ZDm[m+s+p]*Ims[m-s]-un*Zm[m+s+p]*IDms[m-s]
			   +up*ZDm[m-s]*Ims[m+s+p]-un*Zm[m-s]*IDms[m+s+p]);}

	 for(m=0;m<s;m++)
	  {Sum +=Am[2*m+p]*(up*ZDm[m+s+p]*Ims[s-m]-un*Zm[m+s+p]*IDms[s-m]
			   +up*ZDm[s-m]*Ims[m+s+p]-un*Zm[s-m]*IDms[m+s+p]);}
      }

    if ( typ=='o' )
     {
	Tmp=1/(Am[2*s+p]+1e-50);

	 for(m=s;m<CDIM-s-p;m++)
	  {Sum +=Am[2*m+p]*(up*ZDm[m+s+p]*Ims[m-s]-un*Zm[m+s+p]*IDms[m-s]
			   -up*ZDm[m-s]*Ims[m+s+p]+un*Zm[m-s]*IDms[m+s+p]);}

	 for(m=0;m<s;m++)
	  {Sum +=Am[2*m+p]*(up*ZDm[m+s+p]*Ims[s-m]-un*Zm[m+s+p]*IDms[s-m]
			    -up*ZDm[s-m]*Ims[m+s+p]+un*Zm[s-m]*IDms[m+s+p]);}
      }

	MathZD=Sum*.5*h*Tmp;

     return MathZD;
}
/*------------------------------------------------------------------*/
h_type  MathuKDn(char typ, int n, h_type h, double u, h_type Am[], h_type Ims[], h_type IDms[], h_type Zm[], h_type ZDm[], int CDIM)
{

 int  m,ni,p,s;
 h_type Sum,MathZD,Tmp;
 double un,up,sign=1.0;

    un=exp(-u);
    up=exp(u);

    p=n%2;
    if(p==1) sign=-1.0;
    ni=(n-p)/2;
    s=ni;
    if(s>=CDIM) s=CDIM-1;
    Sum=0;
    if ( typ=='e' )
     {
	Tmp=1.0/(Am[2*s+p]+1e-50);
	if (s==0 && n!=1) Tmp=0.5/(Am[2*s+p]+1e-50);

	 for(m=s;m<CDIM-s-p;m++)
	  {Sum +=signf(m)*Am[2*m+p]*(up*ZDm[m+s+p]*Ims[m-s]-un*Zm[m+s+p]*IDms[m-s]
			       +sign*(up*ZDm[m-s]*Ims[m+s+p]-un*Zm[m-s]*IDms[m+s+p]));}

	 for(m=0;m<s;m++)
	  {Sum +=signf(m)*Am[2*m+p]*(up*ZDm[m+s+p]*Ims[s-m]-un*Zm[m+s+p]*IDms[s-m]
			       +sign*(up*ZDm[s-m]*Ims[m+s+p]-un*Zm[s-m]*IDms[m+s+p]));}

      }

    if ( typ=='o' )
     {
	Tmp=1.0/(Am[2*s+p]+1e-50);

	 for(m=s;m<CDIM-s-p;m++)
	  {Sum +=signf(m)*Am[2*m+p]*(up*ZDm[m+s+p]*Ims[m-s]-un*Zm[m+s+p]*IDms[m-s]
			       -sign*(up*ZDm[m-s]*Ims[m+s+p]-un*Zm[m-s]*IDms[m+s+p]));}

	 for(m=0;m<s;m++)
	  {Sum +=signf(m)*Am[2*m+p]*(up*ZDm[m+s+p]*Ims[s-m]-un*Zm[m+s+p]*IDms[s-m]
			       -sign*(up*ZDm[s-m]*Ims[m+s+p]-un*Zm[s-m]*IDms[m+s+p]));}

      }

      MathZD=Sum*.5*h*Tmp*signf(s);

     return MathZD;
}
/*=======================================================================*/
/*----------  Modified Mathieu Coefficients of the first kind ----------*/
double MCoefficients(char typ,int n, h_type h, h_type Am[], int CDIM, h_type &mcnr, h_type &Norm)
{
/*
   Calculation of the Modified First Kind  Mathieu Coefficients
   Inputs:
   n,h,typ
   n=order, h=paramtere , typ= 'e' for even typ='o' for odd

   Outputs:
   Am= First Kind Modified Mathieu  Coefficients  Aem or Aom
   mcn= Modified Mathieu Characteristic Number 
   Norm=First Kind Modified Normalization constant
   returns the accuracy of the computation acc


*/
  int p,m;
  char ntyp;
  double acc,Tmp;

  p=n%2;
  ntyp=typ;
  if(p==1) {if (typ=='e') ntyp='o'; else ntyp='e';}
  acc=Coefficients(ntyp,n,h,Am,CDIM,mcnr,Norm);

    if (typ=='e') Tmp=1.0/MathuSn(ntyp,n,0.5*pi,Am,CDIM);
    if (typ=='o') Tmp=1.0/MathuSDn(ntyp,n,0.5*pi,Am,CDIM);
    if (typ=='o' && p==1) Tmp=-Tmp;

    for(m=0;m<CDIM;m++)
      {
	   Am[2*m+p]=signf(m)*Am[2*m+p]*Tmp;
      }

    Norm=Norm*Tmp;

  return acc;
}
/*---------  Modified Second Kind Coeffecients --------------------*/
int MCoefficientSec(char typ, int n, h_type h, h_type mcn, h_type Am[], int CDIM, h_type Cm[], h_type &Norm)
{

	  /* The same routin for computing the standard 2nd Kind Coef
	     is used for computing the Modified 2nd Kind Coef
	     with small modifications
	     1. Bm replace by modified coeficients Am
	     2. h2 is replaced by -h2
	     3. for (p==1) mcn is interchanged between the even and the odd
	     Input:typ,n,h,mcn,Am[],CDIM
	     Output: Cm[],Norm
	   */

  int i,s,p,r,sign=1,conv=1;
  h_type  *w,h2,ih2,th,Nrm;

    w= new h_type [2*CDIM+3];
    p=n%2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"MCoefficientSec: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"\n" <<"MCoefficientSec: n can not be zero for odd functions :STOP" <<"\n";
     exit(0);
     }

     s=CDIM-1;
     h2=-0.25*h*h;  // it is negative for the modified
     ih2=1.0/h2;

   if(typ=='o') sign=-1;
   w[2*s+p]=0;
   w[2*s+p-2]=-2.0*(2.0*s+p)*sign*Am[2*s+p]*ih2;

   for(r=s-1;r>=1;r--)
     w[2*r+p-2]=((mcn-(2.0*r+p)*(2.0*r+p))*w[2*r+p]-2.0*(2.0*r+p)*sign*Am[2*r+p])*ih2-w[2*r+2+p];

   if(typ=='e')
    {
     if ( p==0 )  th=((mcn-4.0)*w[2]*ih2-w[4])/(2.0*Am[0])-2.0*mcn/(h2*h2);
     if ( p==1 )  th=((mcn-1.0+h2)*w[1]*ih2-w[3])/(2*Am[1])-ih2;
    }

   if ( typ=='o' )
    {
     if ( p==0 )  {w[0]=(4*Am[2]+(mcn-4)*w[2])/(2*h2)-w[4]/2;th=(w[2]-mcn*w[0]*ih2)/Am[2];}
     if ( p==1 )  th=(w[3]-(mcn-1-h2)*w[1]*ih2)/(2*Am[1])-ih2;
    }

   for(r=0;r<=s;r++)
      Cm[2*r+p]=w[2*r+p]-th*Am[2*r+p];

  if(typ=='e')
    {
     Cm[0]=0;
     Nrm=1;
     for(i=0;i<=s;i++)
      {Nrm=Nrm+(2.0*i+p)*Cm[2*i+p];}

     Nrm=1/Nrm;
    }

  if ( typ=='o' )
    {
     Nrm=Cm[p];
     for (i=1; i<=s; i++)
      { Nrm=Nrm+Cm[2*i+p];}

     Nrm=1/Nrm;
    }

    Norm=Nrm;
    delete[] w;
    return conv;

}
/*==================================================================*/
#undef h_type
