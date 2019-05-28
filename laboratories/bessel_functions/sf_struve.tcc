
namespace std
{
namespace __detail
{

  /**
   * An enum to dispatch Struve function summation.
   */
  enum _StruveType
  : int
  {
    _StruveH,
    _StruveK,
    _StruveL,
    _StruveM
  };

  /**
   * Return either the Struve function of the first kind
   * @f$ \boldmath{H}_\nu(x) @f$ or the modified Struve function
   * of the first kind @f$ \boldmath{L}_\nu(x) @f$
   * depending on whether @c sign is -1 or +1 respectively.
   */
  template<_StruveType _Type, typename _Tp>
    _Tp
    __struve_series(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;

      using _BasicSum = __gnu_cxx::_BasicSum<_Val>;
      using _WenigerBasSum = __gnu_cxx::_WenigerDeltaSum<_BasicSum>;
      using _WijnSum = __gnu_cxx::_VanWijngaardenSum<_Val>;
      using _WenigerWijnSum = __gnu_cxx::_WenigerDeltaSum<_WijnSum>;
      using _WenigerSum = std::conditional_t<_Type == _StruveH,
					     _WenigerWijnSum, _WenigerBasSum>;
      int __sign = (_Type == _StruveH ? -1 : _Type == _StruveL ? +1 : 0);
      assert(__sign != 0);

      constexpr int _S_max_iter = 1000;
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__x));
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(std::real(__x));

      auto __x2 = __x / _Val{2};
      auto __xx4 = _Tp(__sign) * __x2 * __x2;
      auto __term = _Val{1};
      auto __struve = _WenigerSum(_Val{1});
      __struve += __term;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
      	  __term *= __xx4 / _Val(__k + 0.5L) / (__nu + _Val(__k + 0.5L));
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      auto __factor = _Val{2} * std::pow(__x2, __nu + _Val{1})
		    / std::__detail::__gamma(__nu + _Val{1.5L}) / _S_sqrt_pi;

      return __factor * __struve();
    }

  /**
   * Return either the Struve function of the second kind
   * @f$ \boldmath{K}_\nu(x) @f$ or the modified Struve function
   * of the second kind @f$ \boldmath{M}_\nu(x) @f$
   * depending on whether @c sign is +1 or -1 respectively.
   */
  template<_StruveType _Type, typename _Tp>
    _Tp
    __struve_asymp(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;

      using _BasicSum = __gnu_cxx::_BasicSum<_Val>;
      using _WenigerBasSum = __gnu_cxx::_WenigerDeltaSum<_BasicSum>;
      using _WijnSum = __gnu_cxx::_VanWijngaardenSum<_Val>;
      using _WenigerWijnSum = __gnu_cxx::_WenigerDeltaSum<_WijnSum>;
      using _WenigerSum = std::conditional_t<_Type == _StruveM,
					     _WenigerWijnSum, _WenigerBasSum>;

      int __sign = (_Type == _StruveK ? +1 : _Type == _StruveM ? -1 : 0);
      assert(__sign != 0);

      constexpr int _S_max_iter = 1000;
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__x));
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(std::real(__x));

      auto __x2 = __x / _Val{2};
      auto __xx4 = _Val(__sign) * __x2 * __x2;
      auto __term = _Val{1};
      auto __struve = _WenigerSum(_Val{1});
      __struve += __term;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  const auto __term_prev = __term;
      	  __term *= _Val(__k - 0.5L) / (_Val(-__k - 0.5L) + __nu) / __xx4;
	  if (std::abs(__term) > std::abs(__term_prev))
	    break;
	  __struve += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__struve))
	    break;
	}
      auto __fact = _Val(__sign) * std::pow(__x2, __nu - _Val{1})
		  / std::__detail::__gamma(__nu + _Val{0.5L}) / _S_sqrt_pi;

      return __fact * __struve();
    }

  /**
   * Return the Struve function of the first kind
   * @f$ \boldmath{H}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    __struve_h(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_nan = __gnu_cxx::__quiet_NaN(std::real(__x));
      const auto _S_max = __gnu_cxx::__digits10(std::real(__x));

      if (std::real(__x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	std::__throw_domain_error(__N("__struve_h: bad argument"));
      else if (std::isnan(__nu) || std::isnan(__x))
	return _S_nan;
      else if (std::abs(__x) < _S_max)
	return __struve_series<_StruveH>(__nu, __x);
      else
	{
	  auto _Nnu = __cyl_neumann_n(__nu, __x);
	  auto _Knu = __struve_asymp<_StruveK>(__nu, __x);
	  return _Knu + _Nnu;
	}
    }

  /**
   * Return the Struve function of the second kind
   * @f$ \boldmath{K}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    __struve_k(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_nan = __gnu_cxx::__quiet_NaN(std::real(__x));
      const auto _S_max = __gnu_cxx::__digits10(std::real(__x));

      if (std::real(__x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	std::__throw_domain_error(__N("__struve_k: bad argument"));
      else if (std::isnan(__nu) || std::isnan(__x))
	return _S_nan;
      else if (std::abs(__x) >= _S_max)
	return __struve_asymp<_StruveK>(__nu, __x);
      else
	{
	  auto _Nnu = __cyl_neumann_n(__nu, __x);
	  auto _Hnu = __struve_series<_StruveH>(__nu, __x);
	  return _Hnu - _Nnu;
	}
    }

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    __struve_l(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_nan = __gnu_cxx::__quiet_NaN(std::real(__x));
      const auto _S_max = __gnu_cxx::__digits10(std::real(__x));

      if (std::real(__x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	std::__throw_domain_error(__N("__struve_l: bad argument"));
      else if (std::isnan(__nu) || std::isnan(__x))
	return _S_nan;
      else if (std::abs(__x) < _S_max)
	return __struve_series<_StruveL>(__nu, __x);
      else
	{
	  auto _Inu = __cyl_bessel_i(__nu, __x);
	  auto _Mnu = __struve_asymp<_StruveM>(__nu, __x);
	  return _Mnu + _Inu;
	}
    }

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    __struve_m(_Tp __nu, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_nan = __gnu_cxx::__quiet_NaN(std::real(__x));
      const auto _S_max = __gnu_cxx::__digits10(std::real(__x));

      if (std::real(__x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	std::__throw_domain_error(__N("__struve_k: bad argument"));
      else if (std::isnan(__nu) || std::isnan(__x))
	return _S_nan;
      else if (std::abs(__x) >= _S_max)
	return __struve_asymp<_StruveM>(__nu, __x);
      else
	{
	  auto _Inu = __cyl_bessel_i(__nu, __x);
	  auto _Lnu = __struve_series<_StruveL>(__nu, __x);
	  return _Lnu - _Inu;
	}
    }

} // namespace __detail
} // namespace std
