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

/* LINUX: Uncomment the following line when compiling on a Linux system 
   with Intel processor.
#include <x86.h>
*/


/* Include headers to interface properly with the 
multiprecision software library libqd by Bailey. */

#include "qd.h"
#include <iostream>

/* Redefine type double as quad-double */

#define double qd_real

int main(void)
    {
  
        int lerchphi(double *z, double *s, double *v, double *acc, 
             double *result, int *iter);
            
        double z, s, v, res, eps;
        int stat, it;
        
/* LINUX: Uncomment the following two lines (and the corresponding one at 
   the bottom) when compiling on a Linux system with Intel processor.
        unsigned short old_cw;
        x86_fix_start(&old_cw);
  */

        cout<<endl<<"--------------------------- "<<endl<<"Test program for lerchphi()"<<endl<<"v1.0 May 1, 2002"<<endl<<"--------------------------- "<<endl<<endl;
        
        eps = "1.e-32";
        cout<<" Accuracy eps="<<eps<<endl<<endl<<endl;
	
	z="0.99999";
	s="2.0";
	v="1.0";
	cout<<endl<<endl<<" Run 5: "<<endl<<" ------ "<<endl<<" regular case (CNCT) "<<endl<<" ------ "<<endl;
	cout<<" z="<<z<<endl<<" s="<<s<<endl<<" v="<<v<<endl;	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	cout<<" Flag="<<stat<<endl<<" Phi="<<res<<endl<<" Iterations="<<it<<endl;

	z="-0.99999";
	s="2.0";
	v="1.0";
        cout<<endl<<endl<<" Run 6: "<<endl<<" ------ "<<endl<<" regular case (delta) "<<endl<<" ------ "<<endl;
	cout<<" z="<<z<<endl<<" s="<<s<<endl<<" v="<<v<<endl;	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	cout<<" Flag="<<stat<<endl<<" Phi="<<res<<endl<<" Iterations="<<it<<endl;

	z="0.99999";
	s="2.0";
	v="1.0e3";
	cout<<endl<<endl<<" Run 7: "<<endl<<" ------ "<<endl<<" regular case (CNCT) "<<endl<<" ------ "<<endl;
	cout<<" z="<<z<<endl<<" s="<<s<<endl<<" v="<<v<<endl;	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	cout<<" Flag="<<stat<<endl<<" Phi="<<res<<endl<<" Iterations="<<it<<endl;
	
        z="0.3";
	s="2.0";
	v="-4.5";
        cout<<endl<<endl<<" Run 8: "<<endl<<" ------ "<<endl<<" regular case (pow ser) "<<endl<<" ------ "<<endl;
	cout<<" z="<<z<<endl<<" s="<<s<<endl<<" v="<<v<<endl;	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	cout<<" Flag="<<stat<<endl<<" Phi="<<res<<endl<<" Iterations="<<it<<endl;
        
	z="0.00001";
	s="2.0";
	v="1.0";
	cout<<endl<<endl<<" Run 9: "<<endl<<" ------ "<<endl<<" regular case (pow ser) "<<endl<<" ------ "<<endl;
	cout<<" z="<<z<<endl<<" s="<<s<<endl<<" v="<<v<<endl;	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	cout<<" Flag="<<stat<<endl<<" Phi="<<res<<endl<<" Iterations="<<it<<endl;
        
        z="-0.000063";
	s="2.0";
	v="1.0";
	cout<<endl<<endl<<" Run 10: "<<endl<<" ------ "<<endl<<" regular case (pow ser) "<<endl<<" ------ "<<endl;
	cout<<" z="<<z<<endl<<" s="<<s<<endl<<" v="<<v<<endl;	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	cout<<" Flag="<<stat<<endl<<" Phi="<<res<<endl<<" Iterations="<<it<<endl;
                
        z="3.4709929976435479e-06";
	s="1.0";
	v="1.5172413793103448e+00";
	cout<<endl<<endl<<" Run 11: "<<endl<<" ------ "<<endl<<" regular case (pow ser) "<<endl<<" ------ "<<endl;
	cout<<" z="<<z<<endl<<" s="<<s<<endl<<" v="<<v<<endl;	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	cout<<" Flag="<<stat<<endl<<" Phi="<<res<<endl<<" Iterations="<<it<<endl;
        
        z="0.0003";
	s="2.0";
	v="-3.00000000000001";
	cout<<endl<<endl<<" Run 12: "<<endl<<" ------ "<<endl<<" regular case (pow ser) "<<endl<<" ------ "<<endl;
	cout<<" z="<<z<<endl<<" s="<<s<<endl<<" v="<<v<<endl;	
        stat = lerchphi(&z, &s, &v, &eps, &res, &it);
	cout<<" Flag="<<stat<<endl<<" Phi="<<res<<endl<<" Iterations="<<it<<endl;

/* LINUX: Uncomment the following line (and the corresponding ones on top)  
   when compiling on a Linux system with Intel processor.
    x86_fix_end(&old_cw);
*/

        return 0;
	
    }
    
#undef double
