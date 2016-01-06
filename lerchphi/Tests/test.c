/*

---------------------------
Test program for lerchphi() 
---------------------------

This program is copyright by

Sergej V. Aksenov (http://www.geocities.com/saksenov) and 
Ulrich D. Jentschura (jentschura@physik.tu-dresden.de), 2002.

Version 1.00 (May 1, 2002)

*/

#include <math.h>
#include <stdio.h>

int main(void)
    {
  
        double z, s, v, res, eps;
        int stat, it;

        printf("\n--------------------------- \nTest program for lerchphi()\nv1.0 May 1, 2002\n--------------------------- \n\n");
        
        eps = 1.e-14;
        printf(" Accuracy eps=%-2.16e\n\n\n",eps);
	
	z=-1.0;
	s=2.0;
	v=1.0;
	printf(" Run 1: \n ------ \n 1<=|z| \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
	
	z=0.99999;
	s=2.0;
	v=-1.0;
	printf("\n\n Run 2: \n ------ \n v<0 int \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
	
	z=0.99999;
	s=2.3;
	v=-1.5;
	printf("\n\n Run 3: \n ------ \n v<0 not int, s not int \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
	
	z=0.99999999;
	s=1.0;
	v=1.0;
	printf("\n\n Run 4: \n ------ \n overflow in a_j \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
	
	z=0.99999;
	s=2.0;
	v=1.0;
	printf("\n\n Run 5: \n ------ \n regular case (CNCT) \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

	z=-0.99999;
	s=2.0;
	v=1.0;
	printf("\n\n Run 6: \n ------ \n regular case (delta) \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

	z=0.99999;
	s=2.0;
	v=1.0e3;
	printf("\n\n Run 7: \n ------ \n regular case (CNCT) \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
	
        z=0.3;
	s=2.0;
	v=-4.5;
	printf("\n\n Run 8: \n ------ \n regular case (pow ser) \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
        
	z=0.00001;
	s=2.0;
	v=1.0;
	printf("\n\n Run 9: \n ------ \n regular case (pow ser) \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
        
        z=-0.000063;
	s=2.0;
	v=1.0;
	printf("\n\n Run 10: \n ------ \n regular case (pow ser) \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
        
        z=3.4709929976435479e-06;
	s=1.0;
	v=1.5172413793103448e+00;
	printf("\n\n Run 11: \n ------ \n regular case (pow ser) \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);
        
        z=0.0003;
	s=2.0;
	v=-3.00000000000001;
	printf("\n\n Run 12: \n ------ \n regular case (pow ser) \n ------ \n");
	printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

        return ;
	
    }
