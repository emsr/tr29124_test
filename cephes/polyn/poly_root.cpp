int
polynomial_roots(Polynomial<_Tp> xcof, cof, m, std::vector<std::complex<_Tp>>& root)
double xcof[], cof[];
int m;
cmplx root[];
{
  Polynomial<_Tp> cof();

  int iter, retry;
  double mag, cofj;
  std::complex<_Tp> x0, x, dx, t, t1, u, ud;

  bool final = false;
  auto n = m;
  if (n <= 0)
    return 1;
  if (n > 36)
    return 2;
  if (xcof[m] == 0.0)
    return 4;

  auto retry = 0;
  auto n1 = n;
  auto n2 = n;
  auto nsav = n;
  register auto q = &xcof[0];
  register auto p = &cof[n];
  for(auto j = 0; j <= nsav; ++j)
    *p-- = *q++; // cof[ n-j ] = xcof[j];

  std::complex<_Tp> xsav;

nxtrut:
  x0.real() = 0.00500101;
  x0.imag() = 0.01000101;
  retry = 0;

tryagn:
  ++retry;
  x.real() = x0.real();
  x0.real() = -10.0 * x0.imag();
  x0.imag() = -10.0 * x.real();

  x = x0;

finitr:
  iter = 0;

  while (iter < 500)
  {
    u.real() = cof[n];
    if (u.real() == 0.0)
      { /* this root is zero */
	x.r = 0;
	--n1;
	--n2;
	goto zerrut;
      }
    u.i = 0;
    auto ud = std::complex<_Tp>{};
    auto t = std::complex<_Tp>{1, 0};

    p = &cof[n-1];
    for (i=0; i<n; ++i) 
      {
	t1.real() = x.real() * t.real()  -  x.i * t.i;
	t1.i = x.real() * t.i  +  x.i * t.real();
	cofj = *p--; // evaluate polynomial
	u.real() += cofj * t1.real();
	u.i += cofj * t1.i;
	cofj = cofj * (i+1); // derivative
	ud.real() += cofj * t.real();
	ud.i -= cofj * t.i;
	t = t1;
      }

    auto mag = ud.r * ud.r  +  ud.i * ud.i;
    if (mag == 0.0) 
      {
	if( !final )
	  goto tryagn;
	x = xsav;
	goto findon;
      }
    dx.r = (u.i * ud.i - u.r * ud.r) / mag;
    x.r += dx.r;
    dx.i = -(u.r * ud.i + u.i * ud.r) / mag;
    x.i += dx.i;
    if ((fabs(dx.i) + fabs(dx.r)) < 1.0e-6) 
      goto lupdon;
    iter += 1;
  } // while iter < 500

  if (final)
    goto lupdon;
  if ( retry < 5 )
   goto tryagn;

  return 3;

lupdon:
  // Swap original and reduced polynomials
  q = &xcof[nsav];
  p = &cof[0];
  for (j = 0; j <= n2; ++j) 
    {
      cofj = *q;
      *q-- = *p;
      *p++ = cofj;
    }
  i = n;
  n = n1;
  n1 = i;

  if(!final)
    {
      final = true;
      if (fabs(x.i / x.r) < 1.0e-4) 
	x.i = 0.0;
      xsav = x;
      goto finitr; // do final iteration on original polynomial
    }

findon:
  final = false;
  if( fabs(x.i/x.r) >= 1.0e-5 )
    {
      cofj = x.r + x.r;
      mag = x.r * x.r + x.i * x.i;
      n -= 2;
    }
  else
    { // root is real
zerrut:
      x.i = 0.0;
      cofj = x.r;
      mag = 0.0;
      --n;
    }
  // divide working polynomial cof(z) by z - x
  p = &cof[1];
  *p += cofj * *(p-1);
  for (j=1; j < n; ++j) 
    {
      *(p+1) += cofj * *p  -  mag * *(p-1);
      p++;
    }

setrut:
  root.push_back(x);
  if (mag != 0.0) 
    {
      x.i = -x.i;
      mag = 0.0;
      goto setrut; // fill in the complex conjugate root
    }
  if (n > 0) 
    goto nxtrut;

  return 0;
}
