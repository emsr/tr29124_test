/* standard include files*/
#include <iostream>
#include <cstdlib>
#include <cmath>

#define CDIMS 60  /* Important: CIDMs must be greater than the order n, i.e. CDIMs=n+15 */
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

*/
/*
  For theoretical derivation see:
  The companion paper
  or
  F.A. Alhargan, "A complete method for the computations of Mathieu
  characteristic numbers of integer orders," SIAM Review, vol-38, No. 2,
  pp. 239-255, June 1996.

*/
/*
  Mathieu functions of real argument
  last modified 20-12-1998

*/
/*
  To be able to run any of these routines you need the following programs:
  bsslr.cpp
  mcnr.cpp

*/

#define CDIMs 60    /* this should be higher than the order n, i.e CDIMs=n+15 */


const  double pi=3.141592653589793;

inline double abs(double z)
     {return  fabs(z);}
inline double acosh(double z)
       {return log(z+sqrt(z*z-1));}
inline double asinh(double z)
       {return log(z+sqrt(z*z+1));}
inline double atanh(double z)
       {return 0.5*log((1+z)/(1-z));}
inline int signf(int m)     /* sign function */
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

   /*   Routines available in this file  */
   h_type MathuZn(char typ,int n, h_type h, double u, short kind);
   h_type MathuZDn(char typ,int n, h_type h, double u, short kind);
   h_type MathuSn(char typ,int n, h_type h, double v, short kind);
   h_type MathuSDn(char typ,int n, h_type h, double v, short kind);

   h_type  MathuZn(char typ, int n, h_type Bm[], h_type Jms[], h_type Zm[], int CDIM);
   h_type  MathuZDn(char typ, int n, h_type h, double u, h_type Bm[], h_type Jms[], h_type JDms[], h_type Zm[], h_type ZDm[], int CDIM);

   h_type  MathuSn(char typ, int n, double v, h_type Bm[], int CDIM);
   h_type  MathuSDn(char typ, int n, double v, h_type Bm[], int CDIM);
   h_type  ISn(char typ,int ni, double vi, double ti, double ui, h_type Bm[], int CDIM);
   h_type  MathuFn(char typ, int n, double v, h_type Dm[], h_type Norm, h_type Bm[], int CDIM);
   h_type  MathuFDn(char typ, int n, double v, h_type Dm[], h_type Norm, h_type Bm[], int CDIM);

     /* The following do not converge well for all u, better to use MathuZn */
   h_type  MathuJn(char typ, int n, double u, h_type Bm[], h_type Jm[], int CDIM);
   h_type  MathuJDn(char typ, int n, h_type h, double u, h_type Bm[], h_type Jm[], h_type JDm[], int CDIM);
   h_type  MathuYn(char typ, int n, h_type Bm[], h_type Jms[], h_type Ym[], int CDIM);
   h_type  MathuYDn(char typ, int n, h_type h, double u, h_type Bm[], h_type Jms[], h_type JDms[], h_type Ym[], h_type YDm[], int CDIM);


/* ====================================================================*/
/* ================== Mathieu Functions ===============================*/
/*
  Important Convensions
1. CDIM is the number of terms summed in the series, which also determins
   the size of the arrays.
2. The value p is used to identify if n is even (p=0) or odd (p=1).
3. D in any function means differential e.g JDn(x)= d(Jn(x))/dx=J'n(x)
4. The series used are that which are convergent in all the z-plane.
   This is MahtuZm.
5. The Coefficients of the Even functions is assumed to be stored in  one array
   Be() for both even and odd orders, similarly for odd functions in Bo().
6. The Coefficients Be() and Bo() have to be calculated before calling the
   routines  for Mathieu functions.
7. In the case of Jen and Jon the Bessel functions have to be calculated first,
   also for Yen and Yon.
8. Warning are not serious, but care needs to be taken as to the validity of
   the results, in most cases the results are correct.
9. In the case of a serious error the programs haults with a message

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
h_type MathuZn(char typ,int n, h_type h, double u, short kind)
{
/*
   inputs:
   typ ='e' or 'o'
   n  = order
   h  = parameter
   u  = argument 
   kind = 1 or 2

   output:
   The function returns the value of the required  radial function


*/

  int CDIM=CDIMs-1;
  h_type Bm[2*CDIMs],Jms[2*CDIMs],Zm[2*CDIMs];
  h_type mcnr,Norm;
  h_type xp,xs;

  if(kind<1 || kind>2)
     {std::cout <<"MathuZn: The kind has not been specified correctly:STOP"; exit(0);}

  xp=0.5*h*exp(u);
  xs=0.5*h*exp(-u);
  Coefficients(typ,n,h,Bm,CDIM,mcnr,Norm);

  BesselJ(xs,2*CDIM,Jms);
  if(kind==1) BesselJ(xp,2*CDIM,Zm);
  if(kind==2) BesselY(xp,2*CDIM,Zm);

   return MathuZn(typ,n,Bm,Jms,Zm,CDIM);

}
h_type MathuZDn(char typ,int n, h_type h, double u, short kind)
{
  int CDIM=CDIMs-1;
  h_type Bm[2*CDIMs];
  h_type Jms[2*CDIMs],Zm[2*CDIMs];
  h_type JDms[2*CDIMs],ZDm[2*CDIMs];
  h_type mcnr,Norm;
  h_type xp,xs;

  if(kind<1 || kind>2)
     {std::cout <<"MathuZn: The kind has not been specified correctly:STOP"; exit(0);}

  xp=0.5*h*exp(u);
  xs=0.5*h*exp(-u);
  Coefficients(typ,n,h,Bm,CDIM,mcnr,Norm);

  BesselJ(xs,2*CDIM,Jms,JDms);
  if(kind==1) BesselJ(xp,2*CDIM,Zm,ZDm);
  if(kind==2) BesselY(xp,2*CDIM,Zm,ZDm);

   return MathuZDn(typ,n,h,u,Bm,Jms,JDms,Zm,ZDm,CDIM);

}
h_type MathuSn(char typ,int n, h_type h, double v, short kind)
{
/*
   inputs:
   typ ='e' or 'o'
   n  = order  
   h  = parameter
   v  = argument in radians
   kind = 1 or 2
  
   output:
   The function returns the value of the required circumferential function


*/
  int CDIM=CDIMs-1;
  h_type Bm[2*CDIMs],Dm[2*CDIMs];
  h_type Vr,mcnr,Norm1,Norm2;

  if(kind<1 || kind>2)
     {std::cout <<"MathuSn: The kind has not been specified correctly:STOP"; exit(0);}

  Coefficients(typ,n,h,Bm,CDIM,mcnr,Norm1);

  if(kind==1)  Vr=MathuSn(typ,n,v,Bm,CDIM);
  if(kind==2)
    {
     CoefficientSec(typ,n,h,mcnr,Bm,CDIM,Dm,Norm2);
     Vr=MathuFn(typ,n,v,Dm,Norm2,Bm,CDIM);
    }

   return Vr;

}
h_type MathuSDn(char typ,int n, h_type h, double v, short kind)
{
  int CDIM=CDIMs-1;
  h_type Bm[2*CDIMs],Dm[2*CDIMs];
  h_type Vr,mcnr,Norm1,Norm2;

  if(kind<1 || kind>2)
     {std::cout <<"MathuSn: The kind has not been specified correctly:STOP"; exit(0);}

  Coefficients(typ,n,h,Bm,CDIM,mcnr,Norm1);

  if(kind==1)  Vr=MathuSDn(typ,n,v,Bm,CDIM);
  if(kind==2)
    {
     CoefficientSec(typ,n,h,mcnr,Bm,CDIM,Dm,Norm2);
     Vr=MathuFDn(typ,n,v,Dm,Norm2,Bm,CDIM);
    }

   return Vr;

}
/*==============================================================*/
/*----------------- General way of computing Zen and Zon ------------*/
/*
   This Function takes the coefficients Bm[], the Bessel functions
   Jms[]=Jm(0.5h exp(-u))  Zm[]=Zm(0.5h exp(u)) all of size CDIM
   Z={Jn,Yn,Hn}
   This is a general Mathieu function which uses the best convergent series
   for Mathieu functions. This function is suitable for computing;
   first kind, second kind, third kind and fourth kind radial Mathieu functions
   of course Zm will have to changed as appropriate
   for Jen or Jon, Zm[]=Jm[]
       Yen or Yon  Zm[]=Ym[]
       He(1)n      Zm[]=H(1)m[]  etc
  Note for Hen some variables need to be changed to complex.

  Make sure when typ='e' to pass the Coefficients Bem[]
  and when  typ='o' pass the Coefficients Bom[]
*/
h_type  MathuZn(char typ, int n, h_type Bm[], h_type Jms[], h_type Zm[], int CDIM)
{
 int  m,ni,p,s;
 h_type Sum,MathZ,Tmp;

    p=n%2;
    ni=(n-p)/2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"MathuZn: Type has not been specified correctly:STOP";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"MathuZn: n can not be zero for odd functions :STOP";
	exit(0);
     }
    s=ni;
    if(s>=CDIM) s=CDIM-1;
    Sum=0;
    if ( typ=='e' )
    {
	    Tmp=1/(Bm[2*s+p]+1e-50);
	    if (s==0 && n!=1) Tmp=0.5/(Bm[2*s+p]+1e-50);

		for (m=s;m<CDIM-s-p;m++)
		    {Sum +=signf(ni-m)*Bm[2*m+p]*(Jms[m-s]*Zm[m+s+p]+Jms[m+s+p]*Zm[m-s]);}
		for (m=0;m<s;m++)
		    {Sum +=signf(ni-s)*Bm[2*m+p]*(Jms[s-m]*Zm[m+s+p]+Jms[m+s+p]*Zm[s-m]);}

     }

     if ( typ=='o' )
	{
	   Tmp=1/(Bm[2*s+p]+1e-50);
		for (m=s;m<CDIM-s-p;m++)
		    {Sum +=signf(ni-m)*Bm[2*m+p]*(Jms[m-s]*Zm[m+s+p]-Jms[m+s+p]*Zm[m-s]);}
		for (m=0;m<s;m++)
		    {Sum +=signf(ni-s)*Bm[2*m+p]*(Jms[s-m]*Zm[m+s+p]-Jms[m+s+p]*Zm[s-m]);}

     }
	  MathZ=pow(0.5*pi,0.5)*Sum*Tmp;
	   return MathZ;
}
/*-----------------------------------------------------------------------*/
/*- Differential of the Radial Second Kind Even and Odd Mathieu Function -*/
/*
  This Function takes the coefficients Bm[], the Bessel functions
  Jms[]=Jm(0.5h exp(-u)) and JDms[]. Zm[]=Zm(0.5h exp(u)) and ZDm[].
   where Zm={Jm,Ym,Hm}
  sums the series giving  Ze'n(h,cosh u) if typ='e'
  or Zo'n(h, cosh u) if typ='o'
  Note that differentiation is wrt u.
*/
h_type  MathuZDn(char typ, int n, h_type h, double u, h_type Bm[], h_type Jms[], h_type JDms[], h_type Zm[], h_type ZDm[], int CDIM)
{

 int  m,ni,p,s;
 h_type Sum,MathZD,Tmp;
 double un,up;

    un=exp(-u);
    up=exp(u);

    p=n%2;
    ni=(n-p)/2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"MathuZDn: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"\n"<<"MathuZDn: n can not be zero for odd functions :STOP"<<"\n";
     exit(0);
     }

    s=ni;
    if(s>=CDIM) s=CDIM-1;
    Sum=0;
    if ( typ=='e' )
     {
	Tmp=1/(Bm[2*s+p]+1e-50);
	if (s==0 && n!=1) Tmp=0.5/(Bm[2*s+p]+1e-50);

	 for(m=s;m<CDIM-s-p;m++)
	  {Sum +=signf(ni-m)*Bm[2*m+p]*(up*ZDm[m+s+p]*Jms[m-s]-un*Zm[m+s+p]*JDms[m-s]+up*ZDm[m-s]*Jms[m+s+p]-un*Zm[m-s]*JDms[m+s+p]);}

	 for(m=0;m<s;m++)
	  {Sum +=signf(ni-s)*Bm[2*m+p]*(up*ZDm[m+s+p]*Jms[s-m]-un*Zm[m+s+p]*JDms[s-m]+up*ZDm[s-m]*Jms[m+s+p]-un*Zm[s-m]*JDms[m+s+p]);}


      }

    if ( typ=='o' )
     {
	Tmp=1/(Bm[2*s+p]+1e-50);

	 for(m=s;m<CDIM-s-p;m++)
	  {Sum +=signf(ni-m)*Bm[2*m+p]*(up*ZDm[m+s+p]*Jms[m-s]-un*Zm[m+s+p]*JDms[m-s]-up*ZDm[m-s]*Jms[m+s+p]+un*Zm[m-s]*JDms[m+s+p]);}

	 for(m=0;m<s;m++)
	  {Sum +=signf(ni-s)*Bm[2*m+p]*(up*ZDm[m+s+p]*Jms[s-m]-un*Zm[m+s+p]*JDms[s-m]-up*ZDm[s-m]*Jms[m+s+p]+un*Zm[s-m]*JDms[m+s+p]);}

      }
	MathZD=pow(0.5*pi,0.5)*Sum*0.5*h*Tmp;
     return MathZD;
}
/*=======================================================================*/
/*------- Circumferential First kind Even and Odd Mathieu  Function ------*/
h_type  MathuSn(char typ, int n, double v, h_type Bm[], int CDIM)
{
 int  m,p;
 h_type Sum,MathuS;

    p=n%2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"MathuSn: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"\n"<<"MathuSn: n can not be zero for odd functions :STOP"<<"\n";
     exit(0);
     }

    if ( typ=='e' )
       {
	  Sum=Bm[p]*cos(p*v);
	  for(m=1;m<CDIM;m++)
	     { Sum=Sum+Bm[2*m+p]*cos((2*m+p)*v);}

	  MathuS=Sum;
       }
    if ( typ=='o'  )
       {
	  Sum=Bm[p]*sin(p*v);
	  for(m=1;m<CDIM;m++)
	      {Sum=Sum+Bm[2*m+p]*sin((2*m+p)*v);}

	   MathuS=Sum;
       }

       return MathuS;

}
/*--------------------------------------------------------------------*/

/*--Dif. of Circumferential First kind Even and Odd Mathieu  Function --*/
h_type  MathuSDn(char typ, int n, double v, h_type Bm[], int CDIM)
{
 int  m,p;
 h_type Sum,MathuSD;

    p=n%2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"MathuSDn: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"\n"<<"MathuSDn: n can not be zero for odd functions :STOP"<<"\n";
     exit(0);
     }

    if (typ=='e')
       {
	  Sum=-Bm[p]*p*sin(p*v);
	  for(m=1;m<CDIM;m++)
	      {Sum=Sum-(2*m+p)*Bm[2*m+p]*sin((2*m+p)*v);}

	   MathuSD=Sum;
	}

     if ( typ=='o' )
	{
	  Sum=Bm[p]*p*cos(p*v);
	  for(m=1;m<CDIM;m++)
	      { Sum=Sum+(2*m+p)*Bm[2*m+p]*cos((2*m+p)*v);}

	  MathuSD=Sum;
	}

	return MathuSD;

}
/*--------------------------------------------------------------------*/
/*====================================================================*/
/*--- Circumferential Second Kind Mathieu Function (Non-periodic) ---*/
/*---------- Fen() and Fon()  ------*/
/*
   Bm is Bem or Bom 1st kind Mathieu coefficients
   Dm is Dem or Dom 2nd kind Mathieu coefficients
   Norm second kind normalization constant

*/
h_type  MathuFn(char typ, int n, double v, h_type Dm[], h_type Norm, h_type Bm[], int CDIM)
{
 int  m,p;
 h_type Sum,MathuF;

    p=n%2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"MathuFn: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"\n" <<"MathuFn: n can not be zero for odd functions :STOP"<<"\n";
     exit(0);
     }

    if ( typ=='e' )
      {
       Sum=Dm[p]*sin(p*v);
       for(m=1;m<CDIM;m++)
	   {Sum=Sum+Dm[2*m+p]*sin((2*m+p)*v);}
      }

    if ( typ=='o'  )
      {
       Sum=Dm[p]*cos(p*v);
       for (m=1;m<CDIM;m++)
	   {Sum=Sum+Dm[2*m+p]*cos((2*m+p)*v);}
       }

     MathuF=Norm*(v*MathuSn(typ,n,v,Bm,CDIM)+Sum);

    return MathuF;
}
/*----------------------------------------------------------------------*/
/*- Dif. of Circumferential Second Kind Mathieu Function (Non-periodic)--*/
/*------- Fe'n() and F'n() -------*/
h_type  MathuFDn(char typ, int n, double v, h_type Dm[], h_type Norm, h_type Bm[], int CDIM)
{
 int  m,p;
 h_type Sum,MathuFD;

    p=n%2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"MathuFDn: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"\n"<<"MathuFDn: n can not be zero for odd functions :STOP"<<"\n";
     exit(0);
     }


    if ( typ=='e' )
     {
      Sum=p*Dm[p]*cos(p*v);
      for(m=1;m<CDIM;m++)
       {Sum=Sum+Dm[2*m+p]*(2*m+p)*cos((2*m+p)*v);}
     }

    if ( typ=='o'  )
     { Sum=-p*Dm[p]*sin(p*v);
      for (m=1;m<CDIM;m++)
	{Sum=Sum-Dm[2*m+p]*(2*m+p)*sin((2*m+p)*v);}
    }

    MathuFD=Norm*(MathuSn(typ,n,v,Bm,CDIM)+v*MathuSDn(typ,n,v,Bm,CDIM)+Sum);

    return MathuFD;
}
/*--------------------------------------------------------------------*/
/*------------ Radial First Kind Even and Odd Mathieu Function ------*/
/*
  THESE SERIES ARE NOT CONVERGENT FOR ALL Z-PLANE
  BUT ARE LEFT HERE FOR TESTING PURPOSES.

 This function takes the coefficients Bm[] and the Bessel functions
 Jm(h cosh u)  of size 2*CDIM
 then sums the series giving Jen(h,cosh u) if typ='e'
 or  Jon(h, cosh u) if typ='o'
 Make sure when typ='e' to pass the Coefficients Bem[]
 and when  typ='o' pass the Coefficients Bom[]
*/

h_type  MathuJn(char typ, int n, double u, h_type Bm[], h_type Jm[], int CDIM)
{
 int  m,ni,p;
 h_type Sum,MathuJ;


    p=n%2;
    ni=(n-p)/2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"MathuJn: Type has not been specified correctly:STOP";
      exit(1);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"MathuJn: n can not be zero for odd functions :STOP";
     exit(2);
     }

    if ( typ=='e' )
     {
	    Sum=signf(ni)*Bm[p]*Jm[p];
	    for (m=1; m<CDIM; m++)
		{ Sum=Sum+signf(ni-m)*Bm[2*m+p]*Jm[2*m+p]; }

	    MathuJ=pow(pi/2.,.5)*Sum;
     }
    if ( typ=='o' )
      {
	    Sum=signf(ni)*p*Bm[p]*Jm[p];
	    for ( m=1; m<CDIM; m++)
		{ Sum=Sum+signf(ni-m)*Bm[2*m+p]*Jm[2*m+p]*(2*m+p); }

	    MathuJ=pow(pi/2.,.5)*Sum*tanh(u);
      }

       return MathuJ;
}
/*-----------------------------------------------------------------------*/
/*- Differential of the Radial First Kind Even and  Odd Mathieu Function-*/
/*
  This Function takes the coefficients Hm(), the Bessel functions Jm
  and the Dif  of Bessel functions Jm'(h cosh u) and n,h,u.
  sums the series giving  Je'n(h,cosh u) if typ='e'
  or Jo'n(h, cosh u) if typ='o'
  Note that differentiation is wrt u.
*/
h_type  MathuJDn(char typ, int n, h_type h, double u, h_type Bm[], h_type Jm[], h_type JDm[], int CDIM)
{
 int  m,ni,p;
 h_type Sum,MathuJD;

    p=n%2;
    ni=(n-p)/2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"MathuJDn: Type has not been specified correctly:STOP";
      exit(1);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"MathuJDn: n can not be zero for odd functions :STOP";
     exit(2);
     }

    if ( typ=='e' )
     {
	    Sum=signf(ni)*Bm[p]*JDm[p];
	    for ( m=1; m<CDIM; m++)
		{ Sum=Sum+signf(ni-m)*Bm[2*m+p]*JDm[2*m+p];}

	    MathuJD=pow(pi/2.,.5)*Sum*h*sinh(u);
     }
    if ( typ=='o' )
      {
	    Sum=signf(ni)*p*Bm[p]*(tanh(u)*sinh(u)*h*JDm[p]+Jm[p]/pow(cosh(u),2));
	    for ( m=1 ; m<CDIM; m++)
	     { Sum=Sum+signf(ni-m)*Bm[2*m+p]*(tanh(u)*sinh(u)*h*JDm[2*m+p]+Jm[2*m+p]/pow(cosh(u),2))*(2.*m+p);}

	     MathuJD=pow(pi/2.,.5)*Sum;
       }

       return MathuJD;

}
/*-----------------------------------------------------------------------*/
/*=======================================================================*/
/*------------ Radial Second Kind Even and Odd Mathieu Function ----------*/
/*
  This function takes the coefficients Hm() && the Bessel functions
  Jms(0.5h exp(-u)) and Ym(0.5h exp(u)).
  Then sums the series giving Yen(h,cosh u) if typ='e'
  or Yon(h,cosh u) if typ='o'.
*/
h_type  MathuYn(char typ, int n, h_type Bm[], h_type Jms[], h_type Ym[], int CDIM)
{
  int  m,ni,p;
  h_type Sum,MathuY;

    p=n%2;
    ni=(n-p)/2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"MathuYn: Type has not been specified correctly:STOP";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"MathuYn: n can not be zero for odd functions :STOP";
     exit(0);
     }

    if ( typ=='e' )
    { if ( p==0 )
	{   Sum=signf(ni)*Bm[0]*Ym[0]*Jms[0];
		for (m=1;m<CDIM;m++)
		    {Sum=Sum+signf(ni-m)*Bm[2*m]*Ym[m]*Jms[m];}
	}

     if ( p==1 )
	{     Sum=signf(ni)*Bm[1]*(Jms[0]*Ym[1]+Jms[1]*Ym[0]);
		for (m=1;m<CDIM-1;m++)
			{Sum=Sum+signf(ni-m)*Bm[2*m+1]*(Jms[m]*Ym[m+1]+Jms[m+1]*Ym[m]);}
	}

	  MathuY=pow(pi/2.,.5)*Sum/(Bm[p]+1e-30);
     }

     if ( typ=='o' )
	{
	   Sum=signf(ni)*Bm[p]*(Jms[1-p]*Ym[1]-Jms[1]*Ym[1-p]);
	   for (m=1;m<CDIM-1;m++)
		 {Sum=Sum+signf(ni-m)*Bm[2*m+p]*(Jms[m-1+p]*Ym[m+1]-Jms[m+1]*Ym[m-1+p]);}

	  MathuY=pow(pi/2.,.5)*Sum/(Bm[2-p]+1e-30);
	 }
	   return MathuY;
}
/*-----------------------------------------------------------------------*/
/*- Differential of the Radial Second Kind Even and Odd Mathieu Function -*/
/*
  This Function takes the coefficients Bm[], the Bessel functions
  Jms[]=Jm(0.5h exp(-u)) and JDms[]. Ym[]=Ym(0.5h exp(u)) and  YDm[].
  Sums the series giving  Ye'n(h,cosh u) if typ='e'
  or Yo'n(h, cosh u) if typ='o'
  Note that differentiation is wrt u.
*/
h_type MathuYDn(char typ, int n, h_type h, double u, h_type Bm[], h_type Jms[], h_type JDms[], h_type Ym[], h_type YDm[], int CDIM)
{
 int  m,ni,p;
 h_type Sum,MathuYD;
 double up,un;

 un=exp(-u);
 up=exp(u);

    p=n%2;
    ni=(n-p)/2;
    if ( (typ!='e') && (typ!='o'))
    { std::cout <<"\n"<<"MathuYDn: Type has not been specified correctly: type=" << typ <<"\n";
      exit(0);
    }
    if ( n==0 && typ=='o' )
     { std::cout <<"\n"<<"MathuYDn: n can not be zero for odd functions :STOP"<<"\n";
     exit(0);
     }


    Sum=0;
    if ( typ=='e' )
     {
	if ( p==0 )
	{
	   for ( m=0;m<CDIM;m++)
		{Sum=Sum+signf(ni-m)*Bm[2*m]*(up*YDm[m]*Jms[m]-un*Ym[m]*JDms[m]);}
	}

	if ( p==1 )
	{
	   for (m=0;m<CDIM-1;m++)
	    {Sum=Sum+signf(ni-m)*Bm[2*m+1]*(up*YDm[m+1]*Jms[m]-un*Ym[m+1]*JDms[m]+up*YDm[m]*Jms[m+1]-un*Ym[m]*JDms[m+1]);}
	}

	MathuYD=pow(pi/2.,.5)*Sum*.5*h/(Bm[p]+1e-30);
      }

    if ( typ=='o' )
     {
	Sum=signf(ni)*Bm[p]*(up*YDm[1]*Jms[0]-un*JDms[0]*Ym[1]+un*JDms[1]*Ym[0]-up*Jms[1]*YDm[0]);
	for ( m=1 ; m<CDIM-1; m++)
	  { Sum=Sum+signf(ni-m)*Bm[2*m+p]*(up*YDm[m+1]*Jms[m-1+p]-un*Ym[m+1]*JDms[m-1+p]+un*JDms[m+1]*Ym[m-1+p]-up*Jms[m+1]*YDm[m-1+p]);}

	MathuYD=pow(pi/2.,.5)*Sum*.5*h/(Bm[2-p]+1e-30);
      }

     return MathuYD;
}
/*------------------------------------------------------------------------*/
/*====================================================================*/
/*---------- Integrals of Sen and Son  -------------------------------*/
/*
     This uses Romberg integration Routine
     The routine is general only the function call is chaged
*/
h_type ISn(char typ,int ni, double vi, double ti, double ui, h_type Bm[], int CDIM)
{
  int K,M,N,N1,N2,TMP;
  double a,b,ACC,DX,R,X;
  h_type T[2][11],INTG,Z;

  a=vi-ti;b=vi+ti;    // Limits: a lower limit, b upper limit
  INTG =MathuSn(typ,ni,a,Bm,CDIM)*pow(pow(sinh(ui),2)+pow(sin(a),2),.5);
  INTG+=MathuSn(typ,ni,b,Bm,CDIM)*pow(pow(sinh(ui),2)+pow(sin(b),2),.5);      //Function CALL
  T[0][0]=0.5*(b-a)*INTG;

  N1=0;N2=1;K=0;ACC=1;

  while((ACC>1e-6 && K<10) || K<4)    // Max NUMBER OF ITERATIONS
   {
    K++;
    N=(int)pow(2,K); DX=(b-a)/N;Z=0;

   for (R=1;R<=N;R+=2)
	{
	 X=a+(R*DX);
	 Z +=MathuSn(typ,ni,X,Bm,CDIM)*pow(pow(sinh(ui),2)+pow(sin(X),2),.5);      // Function CALL
	}

   T[N2][0]=.5*T[N1][0]+DX*Z;   // VALUE OF T(K,0)

   for ( M=1 ; M<=K; M++)
    { T[N2][M]=(pow(4,M)*T[N2][M-1]-T[N1][M-1])/(pow(4,M)-1);}

   ACC=abs(T[N2][K]-T[N1][K-1]);
   TMP=N1;N1=N2;N2=TMP; //SWITCH ARRAY
  }

  INTG=T[N1][K];
  if (ACC>0.05) std::cout <<"INTERGRATION DID NOT CONVERGE";


   return INTG;

}

/*========================================================================*/
/*--------------------------------------------------------------------*/
#undef h_type
