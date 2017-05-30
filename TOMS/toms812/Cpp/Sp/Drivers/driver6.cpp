#include"Bernstein.h"

// Testing for the Bernstein library: 
//   gcd
//
// This function finds the gcd of two arbitrary
// polynomials.  The gcd found is then used to
// divide the original two polynomials, the norms
// of the two remainders are examined to see if
// they are small enough.
//
// by Evan Yi-Feng Tsai, 2001

void main(void)
{
  int ii;
  double aa[2], bb[2], cc[2], dd[2], ee[2];
  double rf, rg, epsilon = 1e-5;
  Bernstein AA, BB, CC, DD, EE, F, G, gcd, truegcd;
  FILE *out;

  out = fopen("res6","w");

  aa[0] = -0.53;   aa[1] = 1.0 + aa[0];
  bb[0] = -0.81;   bb[1] = 1.0 + bb[0];
  cc[0] = -0.19;   cc[1] = 1.0 + cc[0];
  dd[0] = -0.24;   dd[1] = 1.0 + dd[0];
  ee[0] = -0.66;   ee[1] = 1.0 + ee[0];

  AA = Bernstein(aa,1);
  BB = Bernstein(bb,1);
  CC = Bernstein(cc,1);
  DD = Bernstein(dd,1);
  EE = Bernstein(ee,1);

  F = (AA^4) * (BB^4) * (CC^6);
  G = (AA^4) * (DD^3) * (EE^4);

  fprintf(out,"Test case -- Greatest common divisor:\n\n");
  fprintf(out,"f = (cc^6) (aa^4) (bb^4),\n");
  fprintf(out,"g = (dd^3) (aa^4) (ee^4),\n");

  fprintf(out,"where aa = (x - 0.53), bb = (x - 0.81), cc = (x - 0.19),\n");
  fprintf(out,"and dd = (x - 0.24), ee = (x - 0.66).\n\n");

  gcd = GCD(F,G,epsilon);

  rf = Norm(rem((Normalize(F) / gcd),(gcd.dgr - 1)));
  rg = Norm(rem((Normalize(G) / gcd),(gcd.dgr - 1)));

  fprintf(out,"Norms of the remainders:\n");
  fprintf(out,"rf = %1.18f\n",rf);
  fprintf(out,"rg = %1.18f\n",rg);

  gcd = Normalize(gcd);
  truegcd = Normalize(AA^4);


  fprintf(out,"\nBernstein coefficients of the gcd (Normalized):\n");
  fprintf(out,"Computed              True                  Error\n");
  
  for ( ii = 0; ii <= truegcd.dgr; ii++ )
    fprintf(out,"%1.16f    %1.16f    %1.16f\n",gcd.cf[ii],truegcd.cf[ii],gcd.cf[ii]-truegcd.cf[ii]);

  fclose(out);

}
