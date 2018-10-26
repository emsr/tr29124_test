

template <typename _Tp>
void
__sph_bessel_series(const unsigned int __n, const _Tp __x, _Tp & __j_n)
{

  __j_n = _Tp(1);
  for (unsigned int __i = 1; __i < __n; ++__i)
    {
      __term *= __z / _Tp(2 * __i - 1);
      __sum += __term;
    }

  const _Tp __z = __x * __x / _Tp(2);
  _Tp __term = _Tp(1);
  const unsigned int __max_iter = 10000;
  for (unsigned int __i = 1; __i < __max_iter; ++__i)
    {
      __term *= __z / _Tp(2 * __i - 1);
      __sum += __term;
    }
  //if (__i == __max_iter)
  //  __throw_runtime_error

  return;
}

