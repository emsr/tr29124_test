

template <typename Tp>
void
sph_bessel_series(const unsigned int n, const Tp x, Tp & j_n)
{

  j_n = Tp(1);
  for (unsigned int i = 1; i < n; ++i)
    {
      term *= z / Tp(2 * i - 1);
      sum += term;
    }

  const Tp z = x * x / Tp(2);
  Tp term = Tp(1);
  const unsigned int max_iter = 10000;
  for (unsigned int i = 1; i < max_iter; ++i)
    {
      term *= z / Tp(2 * i - 1);
      sum += term;
    }
  //if (i == max_iter)
  //  throw_runtime_error

  return;
}

