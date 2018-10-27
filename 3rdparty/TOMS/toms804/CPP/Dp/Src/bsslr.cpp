#include   <iostream>
#include   <stdlib.h>
#include   <cmath>

/*
    Standard Bessel functions of first and second kind.
    Computation is based on routines from Numerical Recipes in C with
    some modifications.
    In most cases accuracy is more than 14 decimal places.
    The routines are specifically designed to compute arrays of
    Jn,J'n,Yn and Y'n

*/

inline int IMAX(int a,int b)
	{if (a>b) return a; else return b;}
inline double SIGN(double a, double b)
	      { if (b>=0.0) return fabs(a); else return -fabs(a);}
inline double abs(double z) {return fabs(z);}
const double pi=3.141592653589793;

       /* Routines in this file*/
       double lnGamma(double xx);
       double Gamma(double n);
       double Gamma(int  n);

       double BesselJn(double order, double z);
       double BesselJDn(double order, double z);
       double BesselJ(double z, int M, double bJn[]);
       double BesselJ(double z, int M, double bJn[], double JDn[]);

	double BesselYn(double order, double z);
	double BesselYDn(double order, double z);
	double BesselY(double z, int M, double bYn[]);
	double BesselY(double z, int M, double bYn[], double bYDn[]);

	void bessjy(double x, double xnu, double *Jv, double *Yv, double *JDv, double *YDv); /* fractional order as well as integer */

	double bessj(int n, double x);
	double bessjD(int n, double x);
	double bessy(int n, double x);
	double bessyD(int n, double x);

	/* Arrays of fractional order */
	double BesselJ(double z,double v, int M, double pJn[],double mJn[]);

/*------------- Gamma Function  ---------------------*/
double lnGamma(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double Gamma(double n)  /*   G(n)=(n-1)!  */
{
  int i,nend=20;
  double z,apx,gn;

  if (n==1) return 1;
  if (n+int(abs(n))==0) return 1e300;
  nend=nend+int(abs(n));
  z=n+nend;
  apx=exp((z-.5)*log(z)-z+.5*log(2*pi)+1/(12*z)-pow(z,-3)/360+pow(z,-5)/1260-pow(z,-7)/1680);
  gn=apx;
  for (i=nend;i>=1; i--)
  { z=z-1;
    gn=gn/z;
   }
   return gn;
}

double  Gamma(int n)  /*   G(n)=(n-1)!  */
{
  double gn,i;

  if (n<=0) return 1e300;
  gn=1;
  for (i=1;i<n; i++)   gn=gn*i;

   return gn;
}
/*===============================================================*/
/*--------------- Bessel functions Jn Real----------------------*/
double BesselJn(double order , double z)
{
   double Jn,Yn,JDn,YDn,acc;

   if (z==0)   { if (order==0) return 1; else return 0;}

	bessjy(z,order,&Jn,&Yn,&JDn,&YDn);
	acc=fabs((Jn*YDn-JDn*Yn)*0.5*z*pi-1.0);
    if(acc>1.0e-9) {std::cout <<"\n BesselJn: Compuation is not accurate: acc="<<acc;}

  return  Jn;

}
/*------------------------------------------------*/
double BesselJDn(double order , double z)
{
   double Jn,Yn,JDn,YDn,acc;

   if (z==0)
   { if (order==1) return 0.5; else return 0;}

	bessjy(z,order,&Jn,&Yn,&JDn,&YDn);
	acc=fabs((Jn*YDn-JDn*Yn)*0.5*z*pi-1.0);
    if(acc>1.0e-9) {std::cout <<"\n BesselJDn: Compuation is not accurate: acc="<<acc;}

   return  JDn;

}
/*------------------ Arrays of Bessel Real -----------------*/
double BesselJ(double z, int M, double bJn[])
{       /* Array size is M, i.e bJn[0],...,bJn[M-1]  */
  int i;

  if (z==0)
   { bJn[0]=1;
     for (i=1;i<M;i++)
      bJn[i]=0;
      return 1e-20;
   }
	int j,SJ;
	double az,bj,bjm,bjp,toz,fact;
	double ACR=72.0,sign,Sum;
	double J0,Y0,JD0,YD0,acc=0,J1;

	az=fabs(z);
	toz=2.0/az;

	if (az > (double) M) {
		bessjy(az,0,&J0,&Y0,&JD0,&YD0);
		acc=fabs((J0*YD0-JD0*Y0)*0.5*az*pi-1.0);
		J1=-JD0;
		bJn[0]=J0;
		bJn[1]=J1;
		for (j=1;j<M-1;j++) {
			bJn[j+1]=j*toz*bJn[j]-bJn[j-1];
		}
	} else {
		for(i=0;i<M;i++)   bJn[i]=0.0;
		Sum=0.0;
		bjp=0.0;
		bj=1.0;
		SJ=int(2*((M+int(sqrt(ACR*M))/2.0)));
		for (j=SJ;j>0;j--) {
			bjm=j*toz*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj) > 1.0e10) {
			     for(i=0;i<M;i++)   bJn[i] *=1.0e-10;
				Sum *=1.0e-10;
				bj *= 1.0e-10;
				bjp *= 1.0e-10;
			}
			if(j%2) Sum +=bj;
			if (j<M) bJn[j]=bjp;
		}
	     if(az<1)
		   {
		     bJn[0]=bj;
		     Sum=2.0*Sum-bJn[0];
		     J0=bJn[0]/Sum;
		   }
		else
		  {
		    bessjy(az,0,&J0,&Y0,&JD0,&YD0);
		    acc=fabs((J0*YD0-JD0*Y0)*0.5*az*pi-1.0);
		  }

		fact=J0/bj;
		bJn[0]=J0;
		sign=1.0;
		for(j=1;j<M;j++)
		{  if(z<0.0) sign=-sign;
		   bJn[j]=sign*bJn[j]*fact;
		}
	}
   if(acc>1.0e-9) {std::cout <<"\n BesselJ: Compuation is not accurate: acc="<<acc;}
    return acc;
}
/*-------------------------------------------------------*/
double BesselJ(double z, int M, double bJn[], double JDn[])
{
   int n;
   double acc;

   acc=BesselJ(z,M,bJn);
   JDn[0]=-bJn[1];
   for(n=1;n<M-1;n++)
    { JDn[n]=0.5*(bJn[n-1]-bJn[n+1]);}
   if(z==0) JDn[M-1]=0; else JDn[M-1]=bJn[M-2]-(M-1)*bJn[M-1]/z;

	return acc;
}

/*====================  Bessel Y   =================================*/
/*------------- Bessel Y real -------------------------------------*/
double BesselYn(double order, double z)
{
   double Jn,Yn,JDn,YDn,acc;
  if (abs(z)==0) return 1e200;

	bessjy(z,order,&Jn,&Yn,&JDn,&YDn);
	acc=fabs((Jn*YDn-JDn*Yn)*0.5*z*pi-1.0);
   if(acc>1.0e-9) {std::cout <<"\n BesselYn: Compuation is not accurate: acc="<<acc;}
    return Yn;

}
double BesselYDn(double order, double z)
{
   double Jn,Yn,JDn,YDn,acc;
  if (abs(z)==0) return 1e200;

	bessjy(z,order,&Jn,&Yn,&JDn,&YDn);
	acc=fabs((Jn*YDn-JDn*Yn)*0.5*z*pi-1.0);
    if(acc>1.0e-9) {std::cout <<"\n BesselYDn: Compuation is not accurate: acc="<<acc;}

    return YDn;

}
/*------------------ Bessel Arrays Y real -------------------------*/
double BesselY(double z, int M, double bYn[])
{

 if (z==0)
     {
      std::cout <<"BesselY d: zero argument i.e infinit result can not containue";
      return 1e9;
      }

	double by,bym,byp,toz,acc;
	double Jn,Yn,JDn,YDn;
	int j,i;

	toz=2.0/z;
	bessjy(z,0,&Jn,&Yn,&JDn,&YDn);
	acc=fabs((Jn*YDn-JDn*Yn)*0.5*z*pi-1.0);
	by=-YDn; // for n=1
	bym=Yn; // for n=0
	bYn[0]=bym;
	bYn[1]=by;
	for (j=1;j<M-1;j++) {
		byp=j*toz*by-bym;
		bym=by;
		by=byp;
		bYn[j+1]=by;
	  if(j>4 && fabs(by)>1e260) {break;}
//	       {std::cout<<"BesselY: Yn array is larger than largest machine number z="<<z;exit(0);}
	}
	if(j<M-1)
	      for(i=j+2;i<M;i++) bYn[i]=1e260;

	if(acc>1.0e-9) {std::cout <<"\n BesselY: Compuation is not accurate: acc="<<acc;}
	return acc;

}
double BesselY(double z, int M, double bYn[], double bYDn[])
{
   int n;
   double acc;

    if (z==0)
     {
      std::cout <<"BesselY d: zero argument i.e infinit result can not containue";
      return 1e9;
      }

   acc=BesselY(z,M,bYn);
   bYDn[0]=-bYn[1];
   for(n=1;n<M-1;n++)
    { bYDn[n]=0.5*(bYn[n-1]-bYn[n+1]);}
  bYDn[M-1]=bYn[M-2]-(M-1.0)*bYn[M-1]/z;

	return acc;
}
/*================================================================*/
/*-----------------------------------------------------------*/

double bessj(int n, double z)
{
	int j,jsum,m;
	double az,bj,bjm,bjp,sum,toz,ans,sign=1.0;
	double ACR=400.0;
	double J0,Y0,JD0,YD0,acc,J1;

	if (z == 0.0)	if(n==0) return 1; else return 0.0;
	az=fabs(z);
	if(z<0) sign=pow(-1,n);
	else if (az > (double) n) {
		bessjy(az,0,&J0,&Y0,&JD0,&YD0);
		acc=fabs((J0*YD0-JD0*Y0)*0.5*z*pi-1.0);
		J1=-JD0;
		if(n==0) return J0;
		if(n==1) return sign*J1;
		toz=2.0/az;
		bjm=J0;
		bj=J1;
		for (j=1;j<n;j++) {
			bjp=j*toz*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=sign*bj;
	} else {
		toz=2.0/az;
		m=int(2*((n+(int) sqrt(ACR*n))/2.0));
		jsum=0;
		bjp=ans=sum=0.0;
		bj=1.0;
		for (j=m;j>0;j--) {
			bjm=j*toz*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj) > 1.0e10) {
				bj *= 1.0e-10;
				bjp *= 1.0e-10;
				ans *= 1.0e-10;
				sum *= 1.0e-10;
			}
			if (jsum) sum += bj;
			jsum=!jsum;
			if (j == n) ans=bjp;
		}
		sum=2.0*sum-bj;
		ans /= sum;
		ans=ans*sign;
	}
	if(acc>1.0e-9) {std::cout <<"\n bessj: Compuation is not accurate: acc="<<acc;}

	return ans;
}

/*====================================================*/
double bessy(int n, double z)
{
	int j;
	double by,bym,byp,toz;
	double J0,Y0,JD0,YD0,acc,Y1;

	bessjy(z,0,&J0,&Y0,&JD0,&YD0);
	acc=fabs((J0*YD0-JD0*Y0)*0.5*z*pi-1.0);
	Y1=-YD0;
	if(n==0) return Y0;
	if(n==1) return Y1;

	toz=2.0/z;
	by=Y1;
	bym=Y0;
	for (j=1;j<n;j++) {
		byp=j*toz*by-bym;
		bym=by;
		by=byp;
	}
	if(acc>1.0e-9) {std::cout <<"\n bessy: Compuation is not accurate: acc="<<acc;}
	return by;
}
double bessyD(int n, double z)
{
	int j;
	double by,bym,byp,toz;
	double J0,Y0,JD0,YD0,acc,Y1;

	bessjy(z,0,&J0,&Y0,&JD0,&YD0);
	acc=fabs((J0*YD0-JD0*Y0)*0.5*z*pi-1.0);
	Y1=-YD0;
	if(n==0) return -YD0;

	toz=2.0/z;
	by=Y1;
	bym=Y0;
	for (j=1;j<=n;j++) {
		byp=j*toz*by-bym;
		bym=by;
		by=byp;
	}
	if(acc>1.0e-9) {std::cout <<"\n bessyD: Compuation is not accurate: acc="<<acc;}

	return bym*n/z-byp;
}
/*-------------------------------------------------------------------*/
/*         Fractional Order Bessel Function                          */
/*-------------------------------------------------------------------*/
void bessjy(double x, double xnu, double *Jv, double *Yv, double *JDv, double *YDv)
{
	void beschB(double x, double *gam1, double *gam2, double *gampl,
		double *gammi);

	int i,isign,l,nl;
	double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2,
		fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,Jvl,
		Jvl1,Jvmu,JDv1,JDvl,Jvtemp,Yv1,Yvmu,Yvmup,Yvtemp,sum,sum1,
		temp,w,x2,xi,xi2,xmu,xmu2;
	double EPS=1.0e-16,FPMIN=1.0e-30;
	int MAXIT=10000,XMIN=2;

	if (x <= 0.0 || xnu < 0.0) {std::cout<<"bad arguments in bessjy";exit(0);}
	nl=(x < XMIN ? (int)(xnu+0.5) : IMAX(0,(int)(xnu-x+1.5)));
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	w=xi2/pi;
	isign=1;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) {
		b += xi2;
		d=b-d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b-1.0/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=c*d;
		h=del*h;
		if (d < 0.0) isign = -isign;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) {std::cout<<"x too large in bessjy; tYv asymptotic expansion";exit(0);}
	Jvl=isign*FPMIN;
	JDvl=h*Jvl;
	Jvl1=Jvl;
	JDv1=JDvl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--) {
		Jvtemp=fact*Jvl+JDvl;
		fact -= xi;
		JDvl=fact*Jvtemp-Jvl;
		Jvl=Jvtemp;
	}
	if (Jvl == 0.0) Jvl=EPS;
	f=JDvl/Jvl;
	if (x < XMIN) {
		x2=0.5*x;
		pimu=pi*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschB(xmu,&gam1,&gam2,&gampl,&gammi);
		ff=2.0/pi*fact*(gam1*cosh(e)+gam2*fact2*d);
		e=exp(e);
		p=e/(gampl*pi);
		q=1.0/(e*pi*gammi);
		pimu2=0.5*pimu;
		fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
		r=pi*pimu2*fact3*fact3;
		c=1.0;
		d = -x2*x2;
		sum=ff+r*q;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*(ff+r*q);
			sum += del;
			del1=c*p-i*del;
			sum1 += del1;
			if (fabs(del) < (1.0+fabs(sum))*EPS) break;
		}
		if (i > MAXIT) {std::cout<<"bessy series failed to converge";exit(0);}
		Yvmu = -sum;
		Yv1 = -sum1*xi2;
		Yvmup=xmu*xi*Yvmu-Yv1;
		Jvmu=w/(Yvmup-f*Yvmu);
	} else {
		a=0.25-xmu2;
		p = -0.5*xi;
		q=1.0;
		br=2.0*x;
		bi=2.0;
		fact=a*xi/(p*p+q*q);
		cr=br+q*fact;
		ci=bi+p*fact;
		den=br*br+bi*bi;
		dr=br/den;
		di = -bi/den;
		dlr=cr*dr-ci*di;
		dli=cr*di+ci*dr;
		temp=p*dlr-q*dli;
		q=p*dli+q*dlr;
		p=temp;
		for (i=2;i<=MAXIT;i++) {
			a += 2*(i-1);
			bi += 2.0;
			dr=a*dr+br;
			di=a*di+bi;
			if (fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
			fact=a/(cr*cr+ci*ci);
			cr=br+cr*fact;
			ci=bi-ci*fact;
			if (fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
			den=dr*dr+di*di;
			dr /= den;
			di /= -den;
			dlr=cr*dr-ci*di;
			dli=cr*di+ci*dr;
			temp=p*dlr-q*dli;
			q=p*dli+q*dlr;
			p=temp;
			if (fabs(dlr-1.0)+fabs(dli) < EPS) break;
		}
		if (i > MAXIT)
		     { std::cout <<"cf2 failed in bessjy"; exit(0); }
		gam=(p-f)/q;
		Jvmu=sqrt(w/((p-f)*gam+q));
		Jvmu=SIGN(Jvmu,Jvl);
		Yvmu=Jvmu*gam;
		Yvmup=Yvmu*(p+q/gam);
		Yv1=xmu*xi*Yvmu-Yvmup;
	}
	fact=Jvmu/Jvl;
	*Jv=Jvl1*fact;
	*JDv=JDv1*fact;
	for (i=1;i<=nl;i++) {
		Yvtemp=(xmu+i)*xi2*Yv1-Yvmu;
		Yvmu=Yv1;
		Yv1=Yvtemp;
	}
	*Yv=Yvmu;
	*YDv=xnu*xi*Yvmu-Yv1;
}
void beschB(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
	double Chebev(double a, double b, double c[], int m, double x);
	double xx;
	int NUSE1=7,NUSE2=8;
	static double c1[] = {
		-1.142022680371172e0,6.516511267076e-3,
		3.08709017308e-4,-3.470626964e-6,6.943764e-9,
		3.6780e-11,-1.36e-13};
	static double c2[] = {
		1.843740587300906e0,-0.076852840844786e0,
		1.271927136655e-3,-4.971736704e-6,-3.3126120e-8,
		2.42310e-10,-1.70e-13,-1.0e-15};

	xx=8.0*x*x-1.0;
	*gam1=Chebev(-1.0,1.0,c1,NUSE1,xx);
	*gam2=Chebev(-1.0,1.0,c2,NUSE2,xx);
	*gampl= *gam2-x*(*gam1);
	*gammi= *gam2+x*(*gam1);
}
double Chebev(double a, double b, double c[], int m, double x)
{
	double d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0) {std::cout<<"x not in range in routine Chebev";exit(0);}
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>=1;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}
/*-------------------------------------------------------------------------*/
/*  Takes care of negative orders  */
void nbessjy(double x, double xnu, double *Jv, double *Yv, double *JDv, double *YDv)
{

	double	Jvl1,Jvtemp,Yv1,Yvtemp,xmu;
	int sign=1;

	xmu=xnu;
	if (xmu<0)  {xmu=-xmu;sign=-1;}
	bessjy(x,xmu,Jv,Yv,JDv,YDv);

	if(sign<0)  /* use reflection formula */
	{
	 Jvl1= *Jv*cos(xmu*pi) - *Yv*sin(xmu*pi);
	 Yv1= *Jv*sin(xmu*pi) + *Yv*cos(xmu*pi);
	 Jvtemp= *JDv*cos(xmu*pi) - *YDv*sin(xmu*pi);
	 Yvtemp= *JDv*sin(xmu*pi) + *YDv*cos(xmu*pi);
	 *Jv=Jvl1;
	 *Yv=Yv1;
	 *JDv=Jvtemp;
	 *YDv=Yvtemp;
	}
}
/*===================================================================*/
