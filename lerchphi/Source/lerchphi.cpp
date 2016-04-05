/*

-------------------------------
Lerch's transcendent Phi(z,s,v)
-------------------------------

This program is copyright by

Sergej V. Aksenov (http://www.geocities.com/saksenov) and 
Ulrich D. Jentschura (jentschura@physik.tu-dresden.de), 2002.

Version 1.00 (May 1, 2002)

Calling sequence:

int lerchphi(double *z, double *s, double *v, double *acc,
             double *result, int *iter)

calculates Lerch's Phi transcendent Phi(z,s,v) with *result to a specified 
accuracy *acc after *iter iterations. Double precision is used throughout the calculation. 
The program uses direct summation of the defining series for |z| <= 0.5
and CNCT for 0.5 < |z| < 1.0.
The integer return code has to be interpreted as follows.

-------------
Return codes:
-------------

0 - Normal termination.
1 - Lerch Phi diverges for 1 <= |z|.
2 - Lerch Phi is not defined for integer v <= 0.
3 - pow() is not defined for v < 0 and s not integer.
4 - Long integer overflow in aj().
5 - Underflow in remainder estimate omega in lerchphi().
6 - No convergence within the maximum number of iterations.


Implementation note:

In subroutine aj(), defining variables ind and two2k as type double
instead of long int might eliminate overflow error which occurs for
high indices (error code 4).

*/ 

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <limits>

/* Function that computes van Wijngaarden's A_j for a given j. */

static int aj(double *z, double *s, double *v, int j, double *acc, double *res)
    {
  
	double sum, bjk, z2ind;
	int k, flag;
	unsigned long int ind, two2k;
	const double machmin = std::numeric_limits<double>::min();
	const double macheps = std::numeric_limits<double>::epsilon();
	
	sum = bjk = 0.0;
	k = -1;
	two2k = 1;
	flag = 0;

	/* Sum b^j_k's over k. */

	for (;;)
	    {
	  
		k++;

		/* Index for the term of the original series. */

		if (k > 0) two2k *= 2;
		ind = two2k * (j + 1) - 1;
        
		/* If long integer overflow occurs, variables become zero. 
		Not relevant in v1.0 because two2k and ind are double type. */
        
		if (k > 0 && (two2k == 0 || ind == 0)) 
		    {
			flag = 4;
			break;
		    }

		/* Increment the sum. */

		z2ind = pow(*z, ind);
		bjk = two2k * z2ind / pow(*v + ind, *s);
		sum += bjk;

		/* Stop summation if either sum is zero or
		    |term/sum| is below requested accuracy. */

		if (fabs(sum) <= machmin || fabs(bjk/sum) < 1.0e-2 * (*acc)) break;
        
	    }
    
        *res = sum;
	return flag;
	
    }

/* Function that computes approximation to Lerch Phi as a 
   converging sequence of CNC transforms S^n_k. */

int lerchphi(double *z, double *s, double *v, double *acc, 
             double *result, int *iter)
    {
  
	const double machmin = std::numeric_limits<double>::min();
	const double macheps = std::numeric_limits<double>::epsilon();
	const unsigned short int beta = 1, n = 0, imax = 100;
	unsigned short int j, m; 
	int i, sign, flag;
	double v1, sn, eps0, eps, skn, skn0, omega, *num, *den, *StoreAj,
	    factor, factor1, x, est, iom, sum1, cacc;
	
	/* Local copy of v. */

	v1 = *v;

	/* Special cases. */

	/* 1 <= |z|.  
        (Return error, Lerch Phi diverges.) */

	if (1.0 <= fabs(*z))
	    {
		*result = 1.0;
		*iter = 0;
		return 1;
	    }
    
	/* v <= 0 is integer. 
        (Return error, Lerch Phi is not defined.) */

	if (fabs(floor(*v) - *v) <= macheps*fabs(*v) && *v <= 0.0)
	    {
		*result = 1.0;
		*iter = 0;
		return 2;
	    }

	/* v < 0 is not integer or zero and z != 0 (z == 0 considered below) ... */

	if (*v < 0.0 && fabs(*z) > machmin)
	    {

		/* s is not an integer.
		(Return error because pow() is not defined.) */

		if (fabs(floor(*s) - *s) > macheps*fabs(*s))
		    {
			*result = 1.0;
			*iter = 0;
			return 3;
		    }

		/* s is an integer.
		(Transform v to positive). */

		else
		    {
			m = - (int) floor(*v);
			v1 += m;
			sum1 = 0.0;
			if ((int) *s % 2 == 0) sign = 1;
			else sign = -1;
			for (i = 0; i <= m-1; i++)
			    {
				if ((i > 0) && (*z < 0)) sign = -sign;
				sum1 += sign*pow(fabs(*z),i)/pow(fabs(*v+i),*s);
			    }
		    }
	    }

	/* z = 0 and ... */
 
	if (fabs(*z) <= machmin)
	    {

		/* ... v < 0 is not integer or zero and ... */

		if (*v < 0)
		    {
            
			/* s is not an integer.
			(Return error because pow() is not defined.) */

			if (fabs(floor(*s) - *s) > macheps*fabs(*s))
			    {
				*result = 1.0;
				*iter = 0;
				return 3;
			    }

			/* s is an integer. 
			(Return first term of series.)*/

			else
			    {
				if ((int) *s % 2 == 0) sign = 1;
				else sign = -1;
				*result = sign * 1.0 / pow(fabs(*v), *s);
			    }
		    }

		/* ... v > 0. 
		(Return first term of series.) */

		else
		    {
			*result = 1.0 / pow(*v, *s);
			*iter = 1;
			return 0;
		    }
	    }

	/* General case. */

	/* Some initializations. */

	/* sn denotes current partial sum of defining series:
            z > 0.5: sn is partial sum S_n of the van 
		Wijngaarden transformed series.
	    z <= 0.5: sn is the partial sum of the
		power series defining LerchPhi.
        skn0 and skn denote successive partial sums S^k_n
        that are same as sn in case of direct summation and
        delta-transformed in case of CNCT.
        eps0 and eps denote successive differences between 
        partial sums S^k_n. */

	eps0 = skn = skn0 = sn = 0.0;

	/* omega is next term of a partial sum (of defining power series 
        for direct summation, of van Wijngaarden transformed series 
        for CNCT) and also becomes a remainder estimate in the delta
        transformation in CNCT). */
	
        /* For z <= 0.5 van Wijngaarden transformation is not used
        [hence no calls to aj()]. */
        
	/* Direct summation and CNCT (z < -0.5) case. */

	if (*z <= 0.5) 
	omega = 1.0 / pow(v1, *s);
    
	/* CNCT (z > 0.5) case. */
    
	else
	    {
		flag = aj(z, s, &v1, 0, acc, &omega); 
		if (flag) 
		    {
			*result = 1.0;
			*iter = 0;
			return flag;
		    }
	    }
	  
	/* Allocate memory for working arrays. */

	num = (double *) malloc(imax * sizeof(double));
	den = (double *) malloc(imax * sizeof(double));
	/* StoreAj is used only in CNCT */
	if (*z > 0.5) StoreAj = (double *) malloc(imax * sizeof(double)); 

	flag = 0;
	i = -1;
	sign = -1;

	/* Main loop: iterations for S^k_n. */

	for (;;)
	    {
	  
		/* i points to current iterate. */
	    
		i++;

		/* Increment the sum. */
	    
		sign = -sign;
		sn += omega;

		/* Next term: omega. */
            
		if (*z < 0.0) /* Direct summation and CNCT (z < -0.5) case. */
                    /* Recurrence for power series. */
                    omega = (*z) * pow((v1+i)/(v1+i+1), *s) * omega;
		else /* z > 0 */ 
		    {
			if (*z <= 0.5) /* "Direct summation". */ 
			omega = (*z) * pow((v1+i)/(v1+i+1), *s) * omega;
			else /* CNCT (z > 0.5) case. */
			    {
				*(StoreAj+i) = sign * omega;
				if (i % 2 == 0) 
                                    /* Recurrence for odd pointer i. */
				    {omega = -sign * 0.5 * (*(StoreAj+i/2) - pow(*z, i/2) /
					pow(v1+i/2, *s));}
				else 
				    {
					flag = aj(z, s, &v1, i+1, acc, &omega);
					if (flag) break;
                                        else  omega = -sign * omega;
				    }
			    }
		    }
	
		/* Direct summation case: store current sum and remainder estimate. */
	    
		if (fabs(*z) <= 0.5)
		    {
			skn = sn;
			est = 2.0 * pow(fabs(*z), (i+1)) / pow(v1+i+1, *s);
		    }
	    
		/* CNCT case. */
	    
		else
		    {	
			
			/* Make sure omega is representable machine number. */
		
			if (fabs(omega) <= machmin)
			    {
				flag = 5;
				break;
			    }
			else iom = 1.0 / omega;
   
			/* Last terms in sums of numerator and denominator of
			i-th partial sum. */

			*(num+i) = sn * iom;
			*(den+i) = iom;

			/* Recurrence computation of numerator and
			denominator of a S_k^n. */

			if (i > 0)
			    {
				factor = 1.0;
				*(num+i-1) = *(num+i) - factor * (*(num+i-1));
				*(den+i-1) = *(den+i) - factor * (*(den+i-1));
			    }

			factor1 = (double) (beta+n+i-1) * (beta+n+i-2);
			for(j = 2; j <= i; j++)
			    {
				factor = factor1 / (beta+n+i+j-2) / (beta+n+i+j-3);
				*(num+i-j) = *(num+i-j+1) - factor * (*(num+i-j));
				*(den+i-j) = *(den+i-j+1) - factor * (*(den+i-j));
			    }

			/* Current approximation of the sum S_k^n. */

			skn = *num / *den;
	
		    } /* else CNCT case. */
	
		eps = fabs(skn - skn0);
	
		/* Check the three termination criteria. */

		/* |est/skn| is less than the requested accuracy
		(est is a remainder estimate). */

		if (i > 0 && eps < eps0)
		    {
			if (fabs(*z) > 0.5)
			    {
				x = eps/eps0;
				est = 2.0/x/(1.0-x)*eps;
			    }
                        cacc = fabs(est/skn);    
			if (cacc < (*acc)) break;
		    }
	  
		/* Successive iterates skn are the same. */

		if (eps <= 0.0) break;

		/* Maximum number of iterations is exceeded. */

		if (i > imax-2)
		    {
			flag = 6;
			break;
		    }
	  
	  	/* Go on to the next iteration. */

		skn0 = skn;
		eps0 = eps;
		
	    } /* for */

	/* Store the resulting sum. */

	if (*v < 0)
	    {
		sign = 1;
		if ((*z < 0) && (m % 2 != 0)) sign = -1;
		*result = sum1 + skn * sign * pow(fabs(*z),m);
	    }
	else *result = skn;

	/* Store the number of iterations. */

	*iter = i + 1;

	/* Clean up. */

	free(num);
	free(den);
	if (*z > 0.5) free(StoreAj);

	return flag;
	
    }
