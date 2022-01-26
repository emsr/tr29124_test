

template <typename _Tp>
void
sph_bessel_series(const unsigned int n, const _Tp x, _Tp & j_n)
{

  j_n = _Tp(1);
  for (unsigned int i = 1; i < n; ++i)
    {
      term *= z / _Tp(2 * i - 1);
      sum += term;
    }

  const _Tp z = x * x / _Tp(2);
  _Tp term = _Tp(1);
  const unsigned int max_iter = 10000;
  for (unsigned int i = 1; i < max_iter; ++i)
    {
      term *= z / _Tp(2 * i - 1);
      sum += term;
    }
  //if (i == max_iter)
  //  throw_runtime_error

  return;
}

