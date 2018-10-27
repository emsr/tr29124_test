#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "wright.h"

/* Test driver to evaluate the wright omega function along the */
/* boundaries of the different regions of the approximations.  */

int
main(void)
{
  int i;
  double complex z,w,e,r,cond;
  double pi = M_PI;
  int n=100;

  double td;
  double x[2],y[2];
  double exp_num=160.0;
  FILE *fp=fopen("Results.dat","w");

  //region 1;
  //x=(-2.0,1.0] ,y=2*pi
  x[0]=nextafter(-2.0,1.0);
  x[1]=1.0;
  y[0] = nextafter(2.0*pi,-1.0);
  td= (x[1]-x[0])/(double)n;
  for(i=0;i<n;i++)
    {
      z = x[0]+td*(double)i +I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 2;
  //x=1.0 ,y=(1.0,2*pi)
  y[0]=nextafter(1.0,2.0);
  y[1]=nextafter(2.0*pi,1.0);;
  td= -(y[1]-y[0])/(double)n;
  x[0] = 1.0;
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[1]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
    
  //region 3;
  //x=(-2.0,1.0] ,y=1.0
  x[0]=nextafter(-2.0,1.0);
  x[1]=1.0;
  td= -(x[1]-x[0])/(double)n;
  y[0] = nextafter(1.0,2.0);
  for(i=0;i<n;i++)
    {
      z = x[1]+td*(double)i +y[0]*I;
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
 
  //region 4;
  //x=-2.0 ,y=(1.0,2*pi)
  y[0]=nextafter(1.0,2.0);
  y[1]=nextafter(2.0*pi,-1.0);
  td= (y[1]-y[0])/(double)n;
  x[0] = nextafter(-2.0,1.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[0]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }

  //region 5;
  //x=(-2.0,1.0] ,y=-2*pi
  x[0]=nextafter(-2.0,1.0);
  x[1]=1.0;
  y[0] = nextafter(-2.0*pi,1.0);
  td= (x[1]-x[0])/(double)n;
  for(i=0;i<n;i++)
    {
      z = x[0]+td*(double)i +I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
 
  //region 6;
  //x=1.0 ,y=(-2*pi,-1.0)
  y[0]=nextafter(-2.0*pi,1.0);
  y[1]=nextafter(-1.0,-2.0);
  td= (y[1]-y[0])/(double)n;
  x[0] = 1.0;
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[0]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }

  //region 7;
  //x=(-2.0,1.0] ,y=-1.0
  x[0]=nextafter(-2.0,1.0);
  x[1]=1.0;
  td= (x[0]-x[1])/(double)n;
  y[0] = nextafter(-1.0,-2.0);
  for(i=0;i<n;i++)
    {
      z = x[1]+td*(double)i +y[0]*I;
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 8;
  //x=-2.0 ,y=(-2*pi,-1.0)
  y[0]=nextafter(-2.0*pi,1.0);
  y[1]=nextafter(-1.0,-2.0);
  td= (y[0]-y[1])/(double)n;
  x[0] = nextafter(-2.0,1.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[1]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 9;
  // x=-2.0 y=[-1.0,1.0]
  y[0]=-1.0;
  y[1]=1.0;
  td= (y[1]-y[0])/(double)n;
  x[0] = nextafter(-2.0,1.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[0]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }

  //region 10;
  // x=(-2.0,1.0] y=1.0
  x[0]=nextafter(-2.0,1.0);
  x[1]=1.0;
  td= (x[1]-x[0])/(double)n;
  y[0] = 1.0;
  for(i=0;i<n;i++)
    {
      z = x[0]+td*(double)i +y[0]*I;
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  // region 11
  // x=1.0 y=[1.0,pi]
  y[0]=1.0;
  y[1]=pi;
  td= (y[1]-y[0])/(double)n;
  x[0] = nextafter(1.0,2.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[0]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }

  //region 12
  // (x-0.1e1)*(x-0.1e1)+y*y=pi*pi)
  //        (on inside)
  td = pi/(double)n;
  x[0] = pi/2.0;
  for(i=0;i<n;i++)
    {
      z = nextafter(pi,-1.0)*(cos(x[0]-td*(double)i)+I*sin(x[0]-td*(double)i))+nextafter(1.0,-1.0);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  // region 13
  // x=1.0 y=[-pi,-1.0]
  y[0]=-pi;
  y[1]=-1.0;
  td= (y[1]-y[0])/(double)n;
  x[0] = nextafter(1.0,2.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[0]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
 
  // region 14
  // x=(-2.0,1.0] y=-1.0
  x[0]=nextafter(-2.0,1.0);
  x[1]=1.0;
  td= (x[1]-x[0])/(double)n;
  y[0] = -1.0;
  for(i=0;i<n;i++)
    {
      z = x[0]+td*(double)i +y[0]*I;
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
 
  //region 15
  // x=(-inf,-2) y=pi^+
  for(i=0;i<n;i++)
    {
      x[0] = nextafter(-1.0-exp((double)(n-1-i)/exp_num),HUGE_VAL);
      y[0] = nextafter(pi-0.75*(x[0]+1.0),HUGE_VAL);
      z = x[0]  + I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  // Region 16
    y[0]=0.75+pi;
  y[1]=2.0*pi;
  td= (y[1]-y[0])/(double)n;
  x[0] = nextafter(-2.0,-3.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[0]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 17
  //x=(-2.0,1.0] ,y=2*pi
  x[0]=nextafter(-2.0,1.0);
  x[1]=1.0;
  y[0] =2.0*pi;
  td= (x[1]-x[0])/(double)n;
  for(i=0;i<n;i++)
    {
      z = x[0]+td*(double)i +I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 18
  //x=1.0 ,y=(pi,2*pi)
  y[0]=nextafter(pi,6.0);
  y[1]=nextafter(2.0*pi,1.0);
  td= -(y[1]-y[0])/(double)n;
  x[0] = nextafter(1.0,2.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[1]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 19
  // (x-0.1e1)*(x-0.1e1)+y*y=pi*pi)
  //        (on outside)
  td = pi/(double)(n-1);
  y[0] = pi/2.0;
  for(i=0;i<n;i++)
    {
      y[1]=pi*sin(y[0]-td*i);
      x[0]=sqrt(pi*pi-y[1]*y[1])+1.0;
      if(y[1]<0)
        z=nextafter(x[0],HUGE_VAL)+I*nextafter(y[1],-HUGE_VAL);
      else
        z=nextafter(x[0],HUGE_VAL)+I*nextafter(y[1],HUGE_VAL);
      
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 20;
  //x=1.0 ,y=(-2*pi,-pi)
  y[0]=nextafter(-2.0*pi,1.0);
  y[1]=nextafter(-pi,-6.0);
  td= -(y[1]-y[0])/(double)n;
  x[0] = nextafter(1.0,2.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[1]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 21;
  //x=(-2.0,1.0] ,y=-2*pi
  x[0]=nextafter(-2.0,1.0);
  x[1]=1.0;
  y[0] = nextafter(-2.0*pi,-7.0);
  td= -(x[1]-x[0])/(double)n;
  for(i=0;i<n;i++)
    {
      z = x[1]+td*(double)i +I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  // Region 22
  y[0]=-0.75-pi;
  y[1]=-2.0*pi;
  td= -(y[1]-y[0])/(double)n;
  x[0] = nextafter(-2.0,-3.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[1]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }

  // Region 23
  // x=(-inf,-2) y=pi^+
  for(i=0;i<n;i++)
    {
      x[0] = nextafter(-1.0-exp((double)(i)/exp_num),HUGE_VAL);
      y[0] = nextafter(-pi+0.75*(x[0]+1.0),-HUGE_VAL);
      z = x[0]  + I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 24
  // x=(-inf,-2) y=pi^+
  for(i=0;i<n;i++)
    {
      x[0] = nextafter(-1.0-exp((double)(n-1-i)/exp_num),-HUGE_VAL);
      y[0] = nextafter(-pi+0.75*(x[0]+1.0),HUGE_VAL);
      z = x[0]  + I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }

  
  // Region 25
  y[0]=-pi;
  y[1]=-0.75-pi;
  td= -(y[1]-y[0])/(double)n;
  x[0] = nextafter(-2.0,-3.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[1]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 26
  // x=(-inf,-2) y=pi^+
  y[0] = nextafter(-pi,-7.0);
  for(i=0;i<n;i++)
    {
      x[0] = -1.0-exp((double)(i)/exp_num);
      
      z = x[0]  + I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 27
  // x=(-inf,-2) y=pi^+
  y[0] = nextafter(-pi,1.0);
  for(i=0;i<n;i++)
    {
      x[0] = -1.0-exp((double)(n-1-i)/exp_num);
      z = x[0]  + I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }

  // Region 28
  y[0]=nextafter(-pi,1.0);
  y[1]=pi;
  td= (y[1]-y[0])/(double)n;
  x[0] = nextafter(-2.0,-3.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[0]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 29
  // x=(-inf,-2) y=pi^+
  y[0] = nextafter(pi,1.0);
  for(i=0;i<n;i++)
    {
      x[0] = -1.0-exp((double)(i)/exp_num);
      
      z = x[0]  + I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }

  //region 30
  // x=(-inf,-2) y=pi^+
  y[0] = nextafter(pi,7.0);
  for(i=0;i<n;i++)
    {
      x[0] = -1.0-exp((double)(n-1-i)/exp_num);
      z = x[0]  + I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  // Region 31
  y[0]=nextafter(pi,7.0);
  y[1]=0.75+pi;
  td= (y[1]-y[0])/(double)n;
  x[0] = nextafter(-2.0,-3.0);
  for(i=0;i<n;i++)
    {
      z = x[0]+I*(y[0]+td*(double)i);
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  //region 32
  // x=(-inf,-2) y=pi^+
  for(i=0;i<n;i++)
    {
      x[0] = -1.0-exp((double)(n-1-i)/exp_num);
      y[0] = nextafter(pi-0.75*(x[0]+1.0),0.1);
      z = x[0]  + I*y[0];
      wrightomega_ext(z,&w,&e,&r,&cond);
      fprintf(fp,"%f %f %f %f\n",creal(z),cimag(z),creal(w),cimag(w));
    }
  
  fclose(fp);
  
  return 0;
}



