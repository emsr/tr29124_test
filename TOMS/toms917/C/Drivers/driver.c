#include<stdio.h>
#include "wright.h"


int
main(void)
{
  double x,y;
  double complex z,w,e,r,condest,r_ult;

  printf("Simple driver usage\nEnter value to evaluate in form: Re(z),Im(z)\n");
  scanf("%lf,%lf",&x,&y);
  
  z=x+I*y;

  w=wrightomega(z);

  printf("\nomega(%f %+f I)=%f %+f I\n\n",creal(z),cimag(z),creal(w),cimag(w));

  
  printf("Extended driver usage\n\n");
  
  wrightomega_ext(z,&w,&e,&r,&condest);
  
  /* Calculate ultimate residual */
  r_ult=(2.0*w*w-8.0*w-1.0)/cpow(1.0+w,6.0)*r*r*r*r;

  printf("omega(%f %+f I)=%f %+f I\n",creal(z),cimag(z),creal(w),cimag(w));
  printf("last update step: %g %+g I\n",creal(e),cimag(e));
  printf("Penultimate residual: %g %+g I\n",creal(r),cimag(r));
  printf("Ultimate residual: %g %+g I\n",creal(r_ult),cimag(r_ult));
  printf("Condition number estimate: %g %+g I\n",creal(condest),cimag(condest));


  

  return 0;
}
