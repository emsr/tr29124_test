
extern double PI;

double
struve(double v, double x)
{
  double y, ya, h;
  double onef2err, threef0err;

  double f = floor(v);
  if (v < 0 && v - f == 0.5 )
    {
      y = jv(-v, x);
      f = 1.0 - f;
      g = 2.0 * floor(f / 2.0);
      if (g != f)
        y = -y;
      return y;
    }
  double t = 0.25 * x * x;
  f = abs(x);
  double g = 1.5 * abs(v);
  if (f > 30.0 && f > g)
    {
      onef2err = 1.0e38;
      y = 0.0;
    }
  else
    y = __hyperg_1f2(1.0, 1.5, 1.5 + v, -t, onef2err);

  if (f < 18.0 || x < 0.0)
    {
      threef0err = 1.0e38;
      ya = 0.0;
    }
  else
    ya = __hyperg_3f0(1.0, 0.5, 0.5-v, -1.0 / t, threef0err);

  f = sqrt(PI);
  h = pow(0.5 * x, v - 1.0);

  if (onef2err <= threef0err)
    {
      g = gamma(v + 1.5);
      y *= h * t / ( 0.5 * f * g );
      return y;
    }
  else
    {
      g = gamma(v + 0.5);
      ya *= h / (f * g);
      ya += yv(v, x);
      return ya;
    }
}


