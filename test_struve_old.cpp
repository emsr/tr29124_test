
extern double PI;

double
struve(double nu, double x)
{
  double y, ya, h;
  double onef2err, threef0err;

  double f = floor(nu);
  if (nu < 0 && nu - f == 0.5 )
    {
      y = jv(-nu, x);
      f = 1.0 - f;
      g = 2.0 * floor(f / 2.0);
      if (g != f)
        y = -y;
      return y;
    }
  double t = 0.25 * x * x;
  f = abs(x);
  double g = 1.5 * abs(nu);
  if (f > 30.0 && f > g)
    {
      onef2err = 1.0e38;
      y = 0.0;
    }
  else
    y = __hyperg_1f2(1.0, 1.5, 1.5 + nu, -t, onef2err);

  if (f < 18.0 || x < 0.0)
    {
      threef0err = 1.0e38;
      ya = 0.0;
    }
  else
    ya = __hyperg_3f0(1.0, 0.5, 0.5 - nu, -1.0 / t, threef0err);

  f = sqrt(PI);
  h = pow(0.5 * x, nu - 1.0);

  if (onef2err <= threef0err)
    {
      g = gamma(nu + 1.5);
      y *= h * t / ( 0.5 * f * g );
      return y;
    }
  else
    {
      g = gamma(nu + 0.5);
      ya *= h / (f * g);
      ya += yv(nu, x);
      return ya;
    }
}


