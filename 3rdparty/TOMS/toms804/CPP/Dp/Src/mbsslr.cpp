#include <cmath>
#include <iostream>
#include <cstdlib>

/*
   Implementation of Modified Bessel functions based on numerical recipes.
   The routins give accuracy of over 14 decimal places for most cases.
   Accuracy can be checked using the Wronskians

*/

	  double BesselIn(int n, double x);
	  double BesselKn(int n, double x);
	  void bessik(double x, double xnu, double *Iv, double *Kv, double *IDv, double *KDv);

	  void BesselI(double z, int M, double bIn[]);
	  void BesselI(double z, int M, double bIn[], double IDn[]);
	  void BesselK(double z, int M, double bKn[]);
	  void BesselK(double z, int M, double bKn[], double KDn[]);

/*-----------------------------------------------------------------*/
double BesselIn(int n, double x)
{
	int j,jsum=0,m;
	double bi,bim,bip,tox,ans,Sum=0;
	double ACR=40.0,sign=1;
	double I0,K0,ID0,KD0,I1;

	if (x == 0.0)	{if (n==0) return 1; else return 0.0;}

	bessik(x,0,&I0,&K0,&ID0,&KD0);
	I1=ID0;
	if (n==0) return I0;
	if (n==1) return I1;

		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;
		m=(n+(int) sqrt(ACR*n));
		if(m%2==0) sign=-1;
		for (j=2*m;j>0;j--) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if (fabs(bi) > 1.0e10) {
				ans *= 1.0e-10;
				bi *= 1.0e-10;
				bip *= 1.0e-10;
				Sum *=1.0e-10;
			}
			if (jsum) {sign=-sign;Sum+=sign*bi;}
			jsum=!jsum;
			if (j == n) ans=bip;
		}
		Sum=2*Sum-bi;
		ans *= I0/bi;
		return x < 0.0 && (n & 1) ? -ans : ans;

}
/*---------------------------------------------------------------*/
double BesselKn(int n, double x)
{
	int j;
	double bk,bkm,bkp,tox;
	double I0,K0,ID0,KD0,K1;
	bessik(x,0,&I0,&K0,&ID0,&KD0);
	K1=-KD0;

	if(n==0) return K0;
	if(n==1) return K1;

	tox=2.0/x;
	bkm=K0;
	bk=K1;
	for (j=1;j<n;j++) {
		bkp=bkm+j*tox*bk;
		bkm=bk;
		bk=bkp;
	}
	return bk;
}
/*----------------------------------------------------------*/
/*		Arrays of Modified Bessel Functions         */
/*==========================================================*/

void BesselI(double z, int M, double bIn[])
{
	int j,m,N;
	double bi,bim,bip,toz,fact,sign;
	double ACR=72.0,Sum;
	double I0,K0,ID0,KD0;

	if (z == 0.0)
	     {
	      bIn[0]=1.0;
	      for(j=1;j<M;j++) bIn[j]=0.0;
	      return;
	     }


		toz=2.0/fabs(z);
		bip=0.0;
		for(m=0;m<M;m++)   bIn[m]=0.0;
		bi=1.0;Sum=0.0;
		N=(M+(int) sqrt(ACR*M));
		if(N%2) sign=1; else sign=-1;
		for (j=2*N;j>0;j--) {
			bim=bip+j*toz*bi;
			bip=bi;
			bi=bim;
			if (fabs(bi) > 1.0e10) {
				for(m=0;m<M;m++) bIn[m]=bIn[m]*1.0e-10;
				bi *= 1.0e-10;
				bip *= 1.0e-10;
				Sum *= 1.0e-10;
			}
			if (j<M) bIn[j]=bip;
			if (j%2==0) {if(j%4==0) Sum+=bip; else Sum -=bip;}
		}
		Sum=2*Sum+bi;
		I0=bi/Sum;
		if(fabs(z)>1) bessik(z,0,&I0,&K0,&ID0,&KD0);
		bIn[0]=I0;
		fact=bIn[0]/bi;
		sign=1.0;
		for(j=1;j<M;j++)
		{  if(z<0.0) sign=-sign;
		   bIn[j]=sign*bIn[j]*fact;
		}

}

void BesselI(double z, int M, double bIn[], double IDn[])
{
   int n;

   BesselI(z,M,bIn);
   IDn[0]=bIn[1];
   for(n=1;n<M-1;n++)
    { IDn[n]=0.5*(bIn[n-1]+bIn[n+1]);}
   if(z==0) IDn[M-1]=0; else IDn[M-1]=bIn[M-2]-(M-1.0)*bIn[M-1]/z;

	return;
}

/*---------------------------------------------------------*/
void BesselK(double z,int M,double bKn[])
{
	int j;
	double bk,bkm,bkp,toz;
	double I0,K0,ID0,KD0,K1;

	bessik(z,0,&I0,&K0,&ID0,&KD0);
	for (j=0;j<M;j++) bKn[j]=0;

	K1=-KD0;
	toz=2.0/z;
	bkm=K0;
	bk=K1;
	bKn[0]=bkm;
	bKn[1]=bk;
	for (j=1;j<M-1;j++) {
		bkp=bkm+j*toz*bk;
		bkm=bk;
		bk=bkp;
		bKn[j+1]=bkp;
//		if(bkp>1e250) j=M;
	}
	return;
}
void BesselK(double z, int M, double bKn[], double KDn[])
{
   int n;

   BesselK(z,M,bKn);
   KDn[0]=-bKn[1];
   for(n=1;n<M-1;n++)
    {KDn[n]=-0.5*(bKn[n-1]+bKn[n+1]);}

   KDn[M-1]=-bKn[M-2]-(M-1.0)*bKn[M-1]/z;

	return;
}
/*--------------------------------------------------------------*/
/*---------   Fractional Order Modified Bessel Functions -------*/
/*--------------------------------------------------------------*/
void bessik(double x, double xnu, double *Iv, double *Kv, double *IDv, double *KDv)
{
	void beschb(double x, double *gam1, double *gam2, double *gampl,
		double *gammi);
	void nrerror(char error_text[]);
	int i,l,nl;
	double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
		gammi,gampl,h,p,pimu,q,q1,q2,qnew,Ivl,Ivl1,Ivmu,IDv1,IDvl,
		Ivtemp,Kv1,Kvmu,Kvmup,Kvtemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;

	double EPS=1.0e-16,FPMIN=1.0e-30,PI=3.141592653589793;
	int MAXIT=10000,XMIN=2;


	if (x <= 0.0 || xnu < 0.0) std::cout<<"\n bad arguments in bessik";
	nl=(int)(xnu+0.5);
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) {
		b += xi2;
		d=1.0/(b+d);
		c=b+1.0/c;
		del=c*d;
		h=del*h;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) {std::cout<<"x too large in bessik; try asymptotic expansion";exit(0);}
	Ivl=FPMIN;
	IDvl=h*Ivl;
	Ivl1=Ivl;
	IDv1=IDvl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--) {
		Ivtemp=fact*Ivl+IDvl;
		fact -= xi;
		IDvl=fact*Ivtemp+Ivl;
		Ivl=Ivtemp;
	}
	f=IDvl/Ivl;
	if (x < XMIN) {
		x2=0.5*x;
		pimu=PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,&gam1,&gam2,&gampl,&gammi);
		ff=fact*(gam1*cosh(e)+gam2*fact2*d);
		sum=ff;
		e=exp(e);
		p=0.5*e/gampl;
		q=0.5/(e*gammi);
		c=1.0;
		d=x2*x2;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*ff;
			sum += del;
			del1=c*(p-i*ff);
			sum1 += del1;
			if (fabs(del) < fabs(sum)*EPS) break;
		}
		if (i > MAXIT) {std::cout<<"bessk seIves failed to converge";exit(0);}
		Kvmu=sum;
		Kv1=sum1*xi2;
	} else {
		b=2.0*(1.0+x);
		d=1.0/b;
		h=delh=d;
		q1=0.0;
		q2=1.0;
		a1=0.25-xmu2;
		q=c=a1;
		a = -a1;
		s=1.0+q*delh;
		for (i=2;i<=MAXIT;i++) {
			a -= 2*(i-1);
			c = -a*c/i;
			qnew=(q1-b*q2)/a;
			q1=q2;
			q2=qnew;
			q += c*qnew;
			b += 2.0;
			d=1.0/(b+a*d);
			delh=(b*d-1.0)*delh;
			h += delh;
			dels=q*delh;
			s += dels;
			if (fabs(dels/s) < EPS) break;
		}
		if (i > MAXIT) {std::cout<<"bessik: failure to converge in cf2";exit(0);}
		h=a1*h;
		Kvmu=sqrt(PI/(2.0*x))*exp(-x)/s;
		Kv1=Kvmu*(xmu+x+0.5-h)*xi;
	}
	Kvmup=xmu*xi*Kvmu-Kv1;
	Ivmu=xi/(f*Kvmu-Kvmup);
	*Iv=(Ivmu*Ivl1)/Ivl;
	*IDv=(Ivmu*IDv1)/Ivl;
	for (i=1;i<=nl;i++) {
		Kvtemp=(xmu+i)*xi2*Kv1+Kvmu;
		Kvmu=Kv1;
		Kv1=Kvtemp;
	}
	*Kv=Kvmu;
	*KDv=xnu*xi*Kvmu-Kv1;
}

void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
	double chebev(double a, double b, double c[], int m, double x);
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
	*gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
	*gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
	*gampl= *gam2-x*(*gam1);
	*gammi= *gam2+x*(*gam1);
}
double chebev(double a, double b, double c[], int m, double x)
{
	double d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0) {std::cout<<"x not in range in routine chebev";exit(0);}
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>=1;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}
