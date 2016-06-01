
template<typename _Tp>
  _Tp
  __struve_series(_Tp __nu, _Tp __x, int __sign)
  {
    auto __x2 = __x / _Tp{2};
    auto __xx4 = __sign * __x2 * __x2;
    auto __term = _Tp{1};
    auto __struve = __term;
    for int __k = 1; __k < 100; ++__k)
      {
      	__term *= __xx4 / (__k + 1.5) / (__nu + __k + 1.5);
	__struve += __term;
	if (std::abs(__term) < 1.0e-15 * std::abs(__struve))
	  break;
      }
    __struve *= std::pow(__x2, __nu + _Tp{1}) / std::tgamma(1.5) / std::tgamma(1.5 + __nu););

    return __struve;
  }

template<typename _Tp>
  _Tp
  __struve_h(_Tp __nu, _Tp __x)
  { return __struve_series(__nu, __x, -1); }

template<typename _Tp>
  _Tp
  __struve_l(_Tp __nu, _Tp __x)
  { return __struve_series(__nu, __x, +1); }
