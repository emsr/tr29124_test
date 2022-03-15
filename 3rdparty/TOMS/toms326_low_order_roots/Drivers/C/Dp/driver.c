#include <stdio.h>
/*
Driver program to demponstrate QUADROOTS, CUBICROOTS, BIQUADROOTS
*/
void main()
{
  int i,k,n;
  double p[5],r[3][5];
  int QUADROOTS(double[],double[][]);
  int CUBICROOTS(double[],double[][]);
  int BIQUADROOTS(double[],double[][]);

  for (n=2; n<=4; n++){
/* Various test cases */
  if(n==4)
  { p[0]=1; p[1]=-10; p[2]=35; p[3]=-50; p[4]=24; /* 1,2,3,4 */
    p[0]=1; p[1]=0; p[2]=4; p[3]=0; p[4]=4; /* +-i sqrt(2) */
  }
  if(n==3) { p[0]=1; p[1]=-6; p[2]=11; p[3]=-6; } /* 1,2,3 */
  if(n==2) { p[0]=1; p[1]=-2; p[2]=2;} /* 1+- i */

  if(n==2){
    printf("Testing quadratic\n");
    i=QUADROOTS(p,r);
  }
  else if(n==3){
    printf("\nTesting cubic\n");
    i=CUBICROOTS(p,r);
  }
  else if(n==4){
    printf("\nTesting biquadratic\n");
    i=BIQUADROOTS(p,r);
  }

  for(k=1;k<=n;k++)printf(" Re=%14.10f Im=%14.10f \n",r[1][k],r[2][k]);
 }
}
