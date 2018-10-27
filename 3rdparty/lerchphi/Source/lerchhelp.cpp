/*
$HOME/bin/bin/g++ -o lerchhelp lerchhelp.cpp
./lerchhelp > lerchhelp.txt
kdiff3 lerchhelp.txt ../Tests/test.out
*/

#include <cmath>
#include <vector>
#include <limits>
#include <cmath>
#include <cstdio>


template<typename _Tp>
  int
  aj(_Tp z, _Tp s, _Tp a, int j, _Tp acc, _Tp& res)
  {
    _Tp sum, bjk, z2ind;
    long long int k, flag;
    unsigned long long int ind, two2k;
    const _Tp machmin = std::numeric_limits<_Tp>::min();

    sum = bjk = _Tp{0};
    k = -1;
    two2k = 1;
    flag = 0;

    // Sum b^j_k's over k.
    while (true)
      {
	++k;

	// Index for the term of the original series.
	if (k > 0)
	  two2k *= 2;
	ind = two2k * (j + 1) - 1;

	/* If long integer overflow occurs, variables become zero. 
	Not relevant in v1.0 because two2k and ind are _Tp type. */

	if (k > 0 && (two2k == 0 || ind == 0)) 
	  {
	    flag = 4;
	    break;
	  }

	// Increment the sum.
	z2ind = pow(z, ind);
	bjk = two2k * z2ind / pow(a + ind, s);
	sum += bjk;

	/* Stop summation if either sum is zero or
	    |term/sum| is below requested accuracy. */
	if (fabs(sum) <= machmin || fabs(bjk/sum) < _Tp{1.0e-2} * acc)
	  break;
      }

    res = sum;
    return flag;
  }

template<typename _Tp>
  int
  lerch_phi(_Tp z, _Tp s, _Tp a, _Tp acc, 
            _Tp& result, int& iter)
  {
    const _Tp machmin = std::numeric_limits<_Tp>::min();
    const _Tp macheps = std::numeric_limits<_Tp>::epsilon();

    const unsigned int beta = 1, n = 0, imax = 100;
    _Tp omega, x, est, iom, sum1, cacc;

    bool aneg = a < _Tp{0};
    bool ztiny = std::abs(z) <= machmin;
    bool zsmall = !ztiny && z <= _Tp{0.5};
    bool sint = std::abs(std::floor(s) - s) <= macheps * std::abs(s);

    int m = 0;
    _Tp a1 = a;

    /* Special cases. */

    if (std::abs(z) >= _Tp{1})
      { // |z| >= 1. (Return error, Lerch Phi diverges.)
	result = _Tp{1};
	iter = 0;
	return 1;
      }

    if (std::abs(std::floor(a) - a) <= macheps * std::abs(a) && a <= _Tp{0})
      { // a <= 0 is integer. (Return error, Lerch Phi is not defined.)
	result = _Tp{1};
	iter = 0;
	return 2;
      }

    if (aneg && !ztiny)
      { // a < 0 is not integer or zero and z != 0 (z == 0 considered below) ...

	if (!sint)
	  { // ... s is not an integer. (Return error because pow() is not defined.)
	    result = _Tp{1};
	    iter = 0;
	    return 3;
	  }
	else
	  { // ... s is an integer. (Transform a to positive).
	    m = -int(std::floor(a));
	    a1 += m;
	    sum1 = _Tp{0};
	    int sign = (int(s) % 2 == 0) ? 1 : -1;
	    for (int i = 0; i <= m - 1; ++i)
	      {
		sum1 += sign * std::pow(std::abs(z), i)
			     / std::pow(std::abs(a + i), s);
		if (z < _Tp{0})
		  sign = -sign;
	      }
	  }
      }

    if (ztiny)
      { // z = 0 and ...
	if (aneg)
	  { // ... a < 0 is not integer or zero and ...
	    if (!sint)
	      { // ... s is not an integer. (Return error because pow() is not defined.)
		result = _Tp{1};
		iter = 0;
		return 3;
	      }
	    else
	      { // ... s is an integer. (Return first term of series.)
		int sign = (int(s) % 2 == 0) ? 1 : -1;
		result = sign * _Tp{1} / std::pow(std::abs(a), s);
	      }
	  }
	else
	  { // ... a > 0. Return first term of series.)
	    result = _Tp{1} / std::pow(a, s);
	    iter = 1;
	    return 0;
	  }
      }

    /* General case. */

    /* Some initializations. */

    /* sn denotes current partial sum of defining series:
	z > 0.5: sn is partial sum S_n of the van 
	    Wijngaarden transformed series.
	z <= 0.5: sn is the partial sum of the
	    power series defining Lerch Phi.
    skn0 and skn denote successive partial sums S^k_n
    that are same as sn in case of direct summation and
    delta-transformed in case of CNCT.
    eps0 and eps denote successive differences between 
    partial sums S^k_n. */

    auto eps0 = _Tp{0};
    auto skn = _Tp{0};
    auto skn0 = _Tp{0};
    auto sn = _Tp{0};

    /* omega is next term of a partial sum (of defining power series 
    for direct summation, of van Wijngaarden transformed series 
    for CNCT) and also becomes a remainder estimate in the delta
    transformation in CNCT). */

    /* For z <= 0.5 van Wijngaarden transformation is not used [hence no calls to aj()]. */

    if (zsmall) // Direct summation and CNCT (z < -0.5) case.
      omega = _Tp{1} / std::pow(a1, s);
    else
      { // CNCT (z > 0.5) case.
	auto flag = aj(z, s, a1, 0, acc, omega); 
	if (flag) 
	  {
	    result = _Tp{1};
	    iter = 0;
	    return flag;
	  }
      }

    std::vector<_Tp> num(imax);
    std::vector<_Tp> den(imax);
    // StoreAj is used only in CNCT.
    std::vector<_Tp> StoreAj;
    if (!zsmall)
      StoreAj.resize(imax); 

    int flag = 0;
    int i = -1;
    int sign = -1;

    // Main loop: iterations for S^k_n.

    while (true)
      {
	// i points to current iterate.
	++i;

	// Increment the sum.
	sign = -sign;
	sn += omega;

	// Next term: omega.

	if (z < _Tp{0}) // Direct summation and CNCT (z < -0.5) case.
          /* Recurrence for power series. */
          omega = z * omega * std::pow((a1 + i) / (a1 + i + 1), s);
	else // z >= 0 */ 
	  {
	    if (zsmall) // Direct summation.
	      omega = z * omega * std::pow((a1 + i) / (a1 + i + 1), s);
	    else
	      { // CNCT (z > 0.5) case.
		StoreAj[i] = sign * omega;
		if (i % 2 == 0) 
		  { // Recurrence for odd pointer i.
		    omega = -sign * _Tp{0.5} * (StoreAj[i / 2]
		          - std::pow(z, i / 2) / std::pow(a1 + i / 2, s));
		  }
		else 
		  {
		    flag = aj(z, s, a1, i + 1, acc, omega);
		    if (flag)
		      break;
                    else
		      omega *= -sign;
		  }
	      }
	  }

	if (std::abs(z) <= _Tp{0.5})
	  { // Direct summation case: store current sum and remainder estimate.
	    skn = sn;
	    est = _Tp{2} * std::pow(std::abs(z), i + 1) / std::pow(a1 + i + 1, s);
	  }
	else 
	  { // CNCT case.
	    if (std::abs(omega) <= machmin)
	      { // Make sure omega is representable machine number.
		flag = 5;
		break;
	      }
	    else
	      iom = _Tp{1} / omega;

	    // Last terms in sums of numerator and denominator of i-th partial sum.

	    num[i] = sn * iom;
	    den[i] = iom;

	    // Recurrence computation of numerator and denominator of a S_k^n.

	    if (i > 0)
	      {
		_Tp factor = _Tp{1};
		num[i-1] = num[i] - factor * num[i-1];
		den[i-1] = den[i] - factor * den[i-1];
	      }

	    _Tp factor1 = _Tp((beta + n + i - 1) * (beta + n + i - 2));
	    for(int j = 2; j <= i; ++j)
	      {
		_Tp factor = factor1 / (beta + n + i + j - 2) / (beta + n + i + j - 3);
		num[i-j] = num[i - j + 1] - factor * num[i - j];
		den[i-j] = den[i - j + 1] - factor * den[i - j];
	      }

	    // Current approximation of the sum S_k^n.
	    skn = num[0] / den[0];
	  }

	auto eps = std::abs(skn - skn0);

	/* Check the three termination criteria. */

	/* |est/skn| is less than the requested accuracy
	(est is a remainder estimate). */

	if (i > 0 && eps < eps0)
	  {
	    if (std::abs(z) > _Tp{0.5})
	      {
		x = eps / eps0;
		est = _Tp{2} / x / (_Tp{1} - x) * eps;
	      }
            cacc = std::abs(est/skn);    
	    if (cacc < acc)
	      break;
	  }

	if (eps <= _Tp{0})
	  break;

	if (i > imax - 2)
	  {
	    flag = 6;
	    break;
	  }

	skn0 = skn;
	eps0 = eps;
      }

    if (aneg)
      {
	sign = 1;
	if ((z < 0) && (m % 2 != 0))
	  sign = -1;
	result = sum1 + skn * sign * std::pow(std::abs(z), m);
      }
    else
      result = skn;

    iter = i + 1;

    return flag;
  }

/*
---------------------------
Test program for lerch_phi() 
---------------------------

This program is copyright by

Sergej V. Aksenov (http://www.geocities.com/saksenov) and 
Ulrich D. Jentschura (jentschura@physik.tu-dresden.de), 2002.

Version 1.00 (May 1, 2002)
*/

int
main()
{
  double z, s, v, res, eps;
  int stat, it;

  printf("\n--------------------------- \nTest program for lerch_phi()\nv1.0 May 1, 2002\n--------------------------- \n\n");

  eps = 1.e-14;
  printf(" Accuracy eps=%-2.16e\n\n\n",eps);

  z=-1.0;
  s=2.0;
  v=1.0;
  printf(" Run 1: \n ------ \n 1<=|z| \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=0.99999;
  s=2.0;
  v=-1.0;
  printf("\n\n Run 2: \n ------ \n v<0 int \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=0.99999;
  s=2.3;
  v=-1.5;
  printf("\n\n Run 3: \n ------ \n v<0 not int, s not int \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=0.99999999;
  s=1.0;
  v=1.0;
  printf("\n\n Run 4: \n ------ \n overflow in a_j \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=0.99999;
  s=2.0;
  v=1.0;
  printf("\n\n Run 5: \n ------ \n regular case (CNCT) \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=-0.99999;
  s=2.0;
  v=1.0;
  printf("\n\n Run 6: \n ------ \n regular case (delta) \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=0.99999;
  s=2.0;
  v=1.0e3;
  printf("\n\n Run 7: \n ------ \n regular case (CNCT) \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=0.3;
  s=2.0;
  v=-4.5;
  printf("\n\n Run 8: \n ------ \n regular case (pow ser) \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=0.00001;
  s=2.0;
  v=1.0;
  printf("\n\n Run 9: \n ------ \n regular case (pow ser) \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=-0.000063;
  s=2.0;
  v=1.0;
  printf("\n\n Run 10: \n ------ \n regular case (pow ser) \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=3.4709929976435479e-06;
  s=1.0;
  v=1.5172413793103448e+00;
  printf("\n\n Run 11: \n ------ \n regular case (pow ser) \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  z=0.0003;
  s=2.0;
  v=-3.00000000000001;
  printf("\n\n Run 12: \n ------ \n regular case (pow ser) \n ------ \n");
  printf(" z=%-2.16e\n s=%-2.16e\n v=%-2.16e\n",z,s,v);	
  stat = lerch_phi(z, s, v, eps, res, it);
  printf(" Flag=%d\n Phi=%-2.16e\n Iterations=%d\n",stat,res,it);

  return 0;
}
