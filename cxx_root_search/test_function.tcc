#ifndef CXX_ROOT_SEARCH_TEST_FUNCTION_TCC
#define CXX_ROOT_SEARCH_TEST_FUNCTION_TCC 1

  // f(x) = sin(x)
  // f'(x) = cos(x)
  // f''(x) = -sin(x)
  template<typename _Tp>
    _Tp
    sin_f(_Tp x)
    { return std::sin(x); }

  template<typename _Tp>
    _Tp
    sin_df(_Tp x)
    { return std::cos(x); }

  template<typename _Tp>
    _Tp
    sin_df2(_Tp x)
    { return -std::sin(x); }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    sin_fdf(_Tp x)
    { return {std::sin(x), std::cos(x)}; }

  // f(x) = cos(x)
  // f'(x) = -sin(x)
  // f''(x) = -cos(x)
  template<typename _Tp>
    _Tp
    cos_f(_Tp x)
    { return std::cos(x); }

  template<typename _Tp>
    _Tp
    cos_df(_Tp x)
    { return -std::sin(x); }

  template<typename _Tp>
    _Tp
    cos_df2(_Tp x)
    { return -std::cos(x); }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    cos_fdf(_Tp x)
    { return {std::cos(x), -std::sin(x)}; }

  // f(x) = x^{20} - 1
  // f'(x) = 20 x^{19}
  // f''(x) = 380 x^{18}
  // zero at x = 1 or -1
  template<typename _Tp>
    _Tp
    func1(_Tp x)
    { return std::pow(x, _Tp{20}) - _Tp{1}; }

  template<typename _Tp>
    _Tp
    func1_df(_Tp x)
    { return _Tp{20} * std::pow(x, _Tp{19}); }

  template<typename _Tp>
    _Tp
    func1_df2(_Tp x)
    { return _Tp{20} * _Tp{19} * std::pow(x, _Tp{18}); }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    func1_fdf(_Tp x)
    { return {func1(x), func1_df(x)}; }

  // f(x) = signum(x) sqrt(std::abs(x))
  // f'(x) = 1 / sqrt(std::abs(x))
  // f''(x) = -signum(x) / pow(abs(x), 3/2) / 2
  // zero at x = 0
  template<typename _Tp>
    _Tp
    func2(_Tp x)
    {
      _Tp delta;

      if (x > _Tp{0})
	delta = _Tp{1};
      else if (x < _Tp{0})
	delta = _Tp{-1};
      else
	delta = _Tp{0};

      return std::sqrt(std::abs(x)) * delta;
    }

  template<typename _Tp>
    _Tp
    func2_df(_Tp x)
    { return _Tp{1} / std::sqrt(std::abs(x)); }

  template<typename _Tp>
    _Tp
    func2_df2(_Tp x)
    {
      _Tp delta;

      if (x == _Tp{0})
	return _Tp{0};
      else
	return _Tp(x < 0 ? +1 : -1)
		/ std::sqrt(std::abs(x)) / std::abs(x) / _Tp{2};
    }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    func2_fdf(_Tp x)
    { return {func2(x), func2_df(x)}; }

  // f(x) = x^2 - 1e-8
  // f'(x) = 2x
  // f''(x) = 2
  // zero at x = sqrt(1e-8) or -sqrt(1e-8)
  template<typename _Tp>
    _Tp
    func3(_Tp x)
    { return std::pow(x, _Tp{2}) - _Tp{1e-8}; }

  template<typename _Tp>
    _Tp
    func3_df(_Tp x)
    { return _Tp{2} * x; }

  template<typename _Tp>
    _Tp
    func3_df2(_Tp)
    { return _Tp{2}; }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    func3_fdf(_Tp x)
    { return {func3(x), func3_df(x)}; }

  // f(x) = x exp(-x)
  // f'(x) = (1 - x) exp(-x)
  // f''(x) = (x - 2) exp(-x)
  // zero at x = 0
  template<typename _Tp>
    _Tp
    func4(_Tp x)
    { return x * std::exp(-x); }

  template<typename _Tp>
    _Tp
    func4_df(_Tp x)
    { return (_Tp{1} - x) * std::exp(-x); }

  template<typename _Tp>
    _Tp
    func4_df2(_Tp x)
    { return  (x - _Tp{2}) * std::exp(-x); }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    func4_fdf(_Tp x)
    { return {func4(x), func4_df(x)}; }

  // f(x) = 1 / (1 + exp(x))
  // f'(x) = -exp(x) / (1 + exp(x))^2
  // f''(x) = -exp(x) / (1 + exp(x))^2 + 2 exp(2x) / (1 + exp(x))^3
  // no roots!
  template<typename _Tp>
    _Tp
    func5(_Tp x)
    { return _Tp{1} / (_Tp{1} + std::exp(x)); }

  template<typename _Tp>
    _Tp
    func5_df(_Tp x)
    { return -std::exp(x) / std::pow(_Tp{1} + std::exp(x), _Tp{2}); }

  template<typename _Tp>
    _Tp
    func5_df2(_Tp x)
    { return func5_df(x) * (_Tp{1} - _Tp{2} * std::exp(x) / (_Tp{1} + std::exp(x))); }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    func5_fdf(_Tp x)
    { return {func5(x), func5_df(x)}; }

  // f(x) = (x - 1)^7
  // f'(x) = 7 (x - 1)^6
  // f''(x) = 42 (x - 1)^5
  // zero at x = 1
  template<typename _Tp>
    _Tp
    func6(_Tp x)
    { return std::pow(x - _Tp{1}, _Tp{7}); }

  template<typename _Tp>
    _Tp
    func6_df(_Tp x)
    { return _Tp{7} * std::pow(x - _Tp{1}, _Tp{6}); }

  template<typename _Tp>
    _Tp
    func6_df2(_Tp x)
    { return _Tp{42} * std::pow(x - _Tp{1}, _Tp{5}); }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    func6_fdf(_Tp x)
    { return {func6(x), func6_df(x)}; }

  template<typename _Tp>
    std::tuple<_Tp, _Tp, _Tp>
    func6_fdf2(_Tp x)
    { return {func6(x), func6_df(x), func6_df(x)}; }

  // Linear function to test that solvers exit correctly 
  // when entered with an exact root
  // f(x) = -pi x + e
  // f'(x) = -pi
  // f''(x) = 0
  template<typename _Tp>
    _Tp
    func7(_Tp x)
    { return -pi<_Tp> * x + e<_Tp>; }

  template<typename _Tp>
    _Tp
    func7_df(_Tp)
    { return -pi<_Tp>; }

  template<typename _Tp>
    _Tp
    func7_df2(_Tp)
    { return _Tp{0}; }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    func7_fdf(_Tp x)
    { return {func7(x), func7_df(x)}; }

#endif // CXX_ROOT_SEARCH_TEST_FUNCTION_TCC
