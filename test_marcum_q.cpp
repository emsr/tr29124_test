
//
// @brief The Marcum Q function.
//
template<typename _Tp>
  _Tp
  marcum_q_series(unsigned int __m, _Tp __a, _Tp __b)
  {
    if ( __b == _Tp(0))
      return _Tp(1);

    _Tp __q1 = _Tp(0);
    const _Tp __arg = __a * __b;
    const _Tp __ab = __a / __b;
    _Tp __temp = _Tp(1);
    unsigned int __k = 0;
    while (true)
      {
	_Tp __term = __temp * std::cyl_bessel_i(__k, __arg);
	__q1 += __term;
	if (std::abs(__term) < std::numeric_limits<_Tp>::epsilon())
          break;
	__temp *= __ab;
	++__k;
      }

    _Tp __q = __q1;

    const _Tp __ba = __b / __a;
    __temp = _Tp(1);
    for (unsigned int __k = 1; __k < __m; ++__k)
      {
	__temp *= __ba;
	__q += __temp * std::cyl_bessel_i(__k, __arg);
      }

    __q *= std::exp(-(__a * __a + __b * __b) / 2.0);

    return __q;
  }

//
// Recent software developments for special functions
// in the Santander-Amsterdam project
//
template<typename _Tp>
  std::pair<_Tp, _Tp>
  marcum_q_gamma_series(_Tp __mu, _Tp __a, _Tp __b)
  {
    auto __fact = _Tp{1};
    for (unsigned int __k = 1; __k < __S_max_iter; ++__k)
      {
	auto _PQ = std::__detail::__gamma(__k + __mu, __b);
	__fact *= __a / __k;
      }
    auto __exp = std::exp(-__a);
    return std::make_pair(__exp * , __exp * );
  }

//
// @brief The Marcum Q function.
//
// Recent software developments for special functions
// in the Santander-Amsterdam project
//
template<typename _Tp>
  _Tp
  marcum_q_integral(unsigned int __m, _Tp __a, _Tp __b)
  {
    auto __rho = [](_Tp __theta, _Tp __xi)
		 -> _Tp
		 { return ; };
    
    return __q;
  }
