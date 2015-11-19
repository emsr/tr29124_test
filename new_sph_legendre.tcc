
  template<typename _Tp>
    _Tp
    __sph_legendre(const unsigned int __l, const unsigned int __m,
                   const _Tp __theta)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      const auto __x = std::cos(__theta);
      if (__m < 0 || __l < __m || __x < -_Tp{1} || __x > _Tp{1})
        std::__throw_domain_error(__N("__sph_legendre: bad argument"));
      else if (__m == 0)
        {
          _Tp _P_l = __legendre_p(__l, __x);
          _Tp __fact = std::sqrt(_Tp(2 * __l + 1) / (_Tp{4} * _S_pi));
          _P_l *= __fact;
          return _P_l;
        }
      else if (__x == _Tp{1} || __x == -_Tp{1})
        return _Tp{0};
      else
        {
          // Starting value for recursion.
          // Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) )
	  //            (-1)^m (1-x^2)^(m/2) / pi^(1/4)
          const auto __sgn = _Tp(__m % 2 == 1 ? -_Tp{1} : _Tp{1});
          const auto _Y_mp1m_factor = __x * std::sqrt(_Tp(2 * __m + 3));
          const auto __lncirc = std::log1p(-__x * __x);
          // Gamma(m+1/2) / Gamma(m)
          const auto __lnpoch = std::lgamma(_Tp(__m + 0.5L))
                              - std::lgamma(_Tp(__m));
          const auto __lnpre_val = -_Tp{0.25L} * _TR1_M_LNPI
				 + _Tp{0.5L} * (__lnpoch + __m * __lncirc);
          const auto __sr = std::sqrt((_Tp{2} + _Tp{1} / __m)
				    / (_Tp{4} * _S_pi));
          auto _Y_mm = __sgn * __sr * std::exp(__lnpre_val);
          auto _Y_mp1m = _Y_mp1m_factor * _Y_mm;

          if (__l == __m)
            return _Y_mm;
          else if (__l == __m + 1)
            return _Y_mp1m;
          else
            {
              _Tp _Y_lm = _Tp{0};

              // Compute Y_l^m, l > m+1, upward recursion on l.
              for (int __ell = __m + 2; __ell <= __l; ++__ell)
                {
                  const auto __rat1 = _Tp(__ell - __m) / _Tp(ell + __m);
                  const auto __rat2 = _Tp(__ell - __m - 1) / _Tp(ell + m - 1);
                  const auto __factor1 = std::sqrt(__rat1
						 * _Tp(2 * __ell + 1)
						 * _Tp(2 * __ell - 1));
                  const auto __factor2 = std::sqrt(__rat1 * __rat2
						 * _Tp(2 * __ell + 1)
						 / _Tp(2 * __ell - 3));
                  _Y_lm = (__x * _Y_mp1m * __factor1
			- _Tp(__ell + __m - 1) * _Y_mm * __factor2)
			/ _Tp(__ell - __m);
                  _Y_mm = _Y_mp1m;
                  _Y_mp1m = _Y_lm;
                }

              return _Y_lm;
            }
        }
    }
