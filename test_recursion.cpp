#include <iostream>

  template<typename _Tp>
    _Tp
    __gegenbauer_poly(unsigned int __n, _Tp __alpha, _Tp __x)
    {
      auto _C0 = _Tp(1);
      if (__n == 0)
        return _C0;

      auto _C1 = _Tp(2) * __alpha * __x;
      if (__n == 1)
        return _C1;

      auto _Cn = _Tp(0);
      for (unsigned int __nn = 2; __nn <= __n; ++__nn)
        {
          _Cn = (_Tp(2) * (_Tp(__nn) - _Tp(1) + __alpha) * __x * _C1
              - (_Tp(__nn) - _Tp(2) + _Tp(2) * __alpha) * _C0)
              / _Tp(__nn);
          _C0 = _C1;
          _C1 = _Cn;
        }
      return _Cn;
    }

int
main()
{
  std::cout.precision(18);

  std::cout <<  __gegenbauer_poly(100, 5.0, 10.000000) << '\n';
  std::cout <<  __gegenbauer_poly(100, 5.0, 10.000001) << '\n';
}

