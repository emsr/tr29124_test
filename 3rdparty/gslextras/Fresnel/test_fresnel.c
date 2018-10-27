#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#pragma hdrstop
#include "fresnel.h"
/* Results from Abramovitz & Stegun */
static double fc01 = 0.0999975; /* =C(0.1) */
static double fs01 = 0.0005236; /* =S(0.1) */

int main(void)
{
/*  FILE *out;
*/
  double x;
  time_t t_beg, t_end;

  printf("fresnel_c(0.1):\n");
  printf("Exact: %2.7f \t Calculated: %2.7f \n ", fc01, fresnel_c(0.1));
  printf("fresnel_s(0.1):\n");
  printf("Exact: %2.7f \t Calculated: %2.7f \n ", fs01, fresnel_s(0.1));

  printf("Test for speed of calculations...\n");
  time(&t_beg);
  for(x=0.0; x<=100.; x+=0.0001) /* 1000000 values*/
  {
   fresnel_c(x);
  }
  time(&t_end);
  printf("Time of calculation in seconds of 1000000 values of fresnel_c(x) = %d\n",
               t_end-t_beg);

  time(&t_beg);
  for(x=0.0; x<=100.; x+=0.0001) /* 1000000 values*/
  {
   fresnel_s(x);
  }
  time(&t_end);
  printf("Time of calculation in seconds of 1000000 values of fresnel_s(x) = %d\n",
               t_end-t_beg);

  /* Table of calculated values*/
/*
  if ((out = fopen("fresnel_c.dat", "wt")) == NULL)
   {
      fprintf(stderr, "Cannot open output file.\n");
      return 1;
   }
  for (x=0.; x<=10.; x+=0.1)
  {
      fprintf(out,"%2.1f \t %2.15f\n", x, fresnel_c(x));
  }
  fclose(out);
  if ((out = fopen("fresnel_s.dat", "wt")) == NULL)
   {
      fprintf(stderr, "Cannot open output file.\n");
      return 1;
   }
  for (x=0.; x<=10.; x+=0.1)
  {
      fprintf(out,"%2.1f \t %2.15f\n", x, fresnel_s(x));
  }
  fclose(out);
*/  
  getchar();
  return 0;
}
//---------------------------------------------------------------------------
