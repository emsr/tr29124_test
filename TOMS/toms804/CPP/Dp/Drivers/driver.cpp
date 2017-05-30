/*            Standard include libraries     */
#include   <cstdio>
#include   <cstdlib>
#include   <iostream>
#include   <fstream>
#include   <cmath>

/*       Include libraries  part of Mathieu functions Library  */
#include   "mathur.h"
#include   "bsslr.h"
#include   "mmathur.h"
#include   "mbsslr.h"
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

  Last modified 20/5/1999

*/

struct point {
   double x;
   double y;
};
/*
  This a driver for the Mathieu functions library
  The following files must be compiled and linked:
  MATHUR.CPP     Mathieu functions
  MMATHUR.CPP    Modified Mathieu functions
  MCNR.CPP       Mathieu charactristic Numbers
  BSSLR.CPP      Bessel functions
  MBSSLR.CPP     Modified Bessel functions
*/
    double Radial(char typ,int n,double h,double u,int kind, double mdf);
    double Circumf(char typ,int n,double h,double v,int kind, double mdf);

/*----------------- Main Porgram --------------------------*/
int  main (void)
{
  point sn1[6];
  int i;
  int n,M;
  double h,v;
  double t,tm;
  double pi;
  char typ,r,mdc,kd;
  short kind=1,md=0;
  char num[]="0123456789";
  char ordc[]="00";
  int n1,n2;
  pi=4*atan(1.0);

  std::cout <<"\n typ=";
  std::cin >>typ;
  std::cout <<" order n=";
  std::cin >>n;
  std::cout <<" Circumferential(c) or Radial(r): ";
  std::cin >>r;
  std::cout <<" kind: First(f) Second(d): ";
  std::cin >>kd;
  std::cout <<" Standard(s) or Modified(m): ";
  std::cin >>mdc;

  if(n<10) ordc[1]=num[n];
	   else
	   {
	    n1=n/10;
	    n2=n-n1*10;
	    ordc[0]=num[n1];
	    ordc[1]=num[n2];
	   }

  if (kd=='d') kind=2; else kind=1;
  if (mdc=='m') md=1; else md=0;

   M=5;

 printf("\n\n%4d \n",n);
 if(r=='r')
 {
  M=4;
  tm=0.0125;
  if(kind==2) tm=10*0.0125;
  for(t=tm;t<=2.5001+tm;t +=0.0125)
  {
    h=t*n;
    for(i=1;i<=M;i++) sn1[i].x=t;

    sn1[1].y=Radial(typ,n,h,0.25,kind,md);  // u=2.5
    sn1[2].y=Radial(typ,n,h,0.5,kind,md);
    sn1[3].y=Radial(typ,n,h,1,kind,md);
    sn1[4].y=Radial(typ,n,h,2,kind,md);

    printf("\n%3.3f   ",sn1[1].x);
    for (i=1;i<=M;i++)
     {
       printf("%+1.6f    ",sn1[i].y);
     }

     }
  }
 else
 {
  M=4;
  for(v=0;v<=2*pi;v +=0.0314)
  {
    for(i=1;i<=M;i++) sn1[i].x=v*180/pi;

    sn1[1].y=Circumf(typ,n,0.5*n,v,kind,md);  // h=0.5*n
    sn1[2].y=Circumf(typ,n,1.0*n,v,kind,md);
    sn1[3].y=Circumf(typ,n,2.0*n,v,kind,md);
    sn1[4].y=Circumf(typ,n,3.0*n,v,kind,md);

    printf("\n%3.3f   ",sn1[1].x);
    for (i=1;i<=M;i++)
     {
       printf("%+1.6f     ",sn1[i].y);
     }

     }

 }


  printf("\n");
  exit(0);
}  // End Main

double Radial(char typ,int n,double h,double u,int kind, double mdf)
{

  if(mdf==1) return MathuMZn(typ,n,h,u,kind);
  return MathuZn(typ,n,h,u,kind);
}
double Circumf(char typ,int n,double h,double v,int kind, double mdf)
{

  if(mdf==1) return MathuQn(typ,n,h,v,kind);
  return MathuSn(typ,n,h,v,kind);
}
