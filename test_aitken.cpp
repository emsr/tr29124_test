// $HOME/bin_specfun/bin/g++ -std=gnu++1z -o test_aitken test_aitken.cpp

// ./test_aitken > test_aitken.txt

// g++ -std=c++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_aitken test_aitken.cpp

// ./test_aitken > test_aitken.txt

#include <cmath>
#include <cstdlib>
#include <vector>
#include <limits>

#include <iostream>
#include <iomanip>

//
//  Sum adapters...
//

  /**
   * The Aitken's delta-squared summation process.
   */
  template<typename _Sum>
    class _AitkenDeltaSqaredSum
    {
    public:

      using value_type = typename _Sum::value_type;

      ///  Default constructor.
      _AitkenDeltaSqaredSum()
      : _M_part_sum{_Sum{}}, _M_a{}, _M_sum{}, _M_converged{false}
      { }

      /// Add a new term to the sum.
      _AitkenDeltaSqaredSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    this->_M_part_sum += __term;
	    this->_M_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _AitkenDeltaSqaredSum&
      operator-=(value_type __term)
      { return operator+=(-__term); }

      /// Return true if the sum converged.
      operator
      bool() const
      { return !this->_M_converged; }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms; }

      ///  Reset the sum to it's initial state.
      _AitkenDeltaSqaredSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_a.clear();
	this->_M_sum = value_type{0};
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _AitkenDeltaSqaredSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      void _M_update();

      _Sum _M_part_sum;
      std::vector<value_type> _M_a;
      value_type _M_sum;
      bool _M_converged;
    };

  template<typename _Sum>
    void
    _AitkenDeltaSqaredSum<_Sum>::_M_update()
    {
      using _Tp = value_type;
      constexpr auto _S_huge = __gnu_cxx::__root_max(_Tp{5}); // 1.0e+60
      constexpr auto _S_tiny = __gnu_cxx::__root_min(_Tp{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      auto& __a = this->_M_a;

      __a.push_back(__s_n);
      if (__n < 2)
	this->_M_sum = __s_n;
      else
	{
	  auto __lowmax = __n / 2;
	  for (auto __j = 1; __j <= __lowmax; ++__j)
	    {
	      auto __m = __n - 2 * __j;
	      auto __denom = (__a[__m + 2] - __a[__m + 1])
			   - (__a[__m + 1] - __a[__m]);
	      if (std::abs(__denom) < _S_tiny)
		__a[__m] = _S_huge;
	      else
		{
		  auto __del = __a[__m] - __a[__m + 1];
		  __a[__m] -= __del * __del / __denom;
		}
	    }
	  this->_M_sum = __a[__n % 2];
	}
    }

  /**
   * The Winn's epsilon summation process.
   */
  template<typename _Sum>
    class _WinnEpsilonSum
    {
    public:

      using value_type = typename _Sum::value_type;

      ///  Default constructor.
      _WinnEpsilonSum()
      : _M_part_sum{_Sum{}}, _M_e{}, _M_sum{}, _M_converged{false}
      { }

      /// Add a new term to the sum.
      _WinnEpsilonSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    this->_M_part_sum += __term;
	    this->_M_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _WinnEpsilonSum&
      operator-=(value_type __term)
      { return operator+=(-__term); }

      /// Return true if the sum converged.
      operator
      bool() const
      { return !this->_M_converged; }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms; }

      ///  Reset the sum to it's initial state.
      _WinnEpsilonSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_e.clear();
	this->_M_sum = value_type{0};
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _WinnEpsilonSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      void _M_update();

      _Sum _M_part_sum;
      std::vector<value_type> _M_e;
      value_type _M_sum;
      bool _M_converged;
    };

  /**
   * Single step of Winn's epsilon transformation.
   */
  template<typename _Sum>
    void
    _WinnEpsilonSum<_Sum>::_M_update()
    {
      using _Tp = value_type;
      constexpr auto _S_huge = __gnu_cxx::__root_max(_Tp{5}); // 1.0e+60
      constexpr auto _S_tiny = __gnu_cxx::__root_min(_Tp{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      auto& __e = this->_M_e;

      __e.push_back(__s_n);
      if (__n == 0)
        this->_M_sum = __s_n;
      else
        {
          auto __aux2 = _Tp{0};
          for (auto __j = __n; __j >= 1; --__j)
            {
	      auto __aux1 = __aux2;
	      __aux2 = __e[__j - 1];
	      auto __diff = __e[__j] - __aux2;
	      if (std::abs(__diff) < _S_tiny)
		__e[__j - 1] = _S_huge;
	      else
		__e[__j - 1] = __aux1 + _Tp{1} / __diff;
	    }
	  this->_M_sum = __e[__n % 2];
        }
      return;
    }

  /**
   * The Brezinski theta summation process.
   */
  template<typename _Sum>
    class _BrezinskiThetaSum
    {
    public:

      using value_type = typename _Sum::value_type;

      ///  Default constructor.
      _BrezinskiThetaSum()
      : _M_part_sum{_Sum{}}, _M_arj{}, _M_sum{}, _M_converged{false}
      { }

      /// Add a new term to the sum.
      _BrezinskiThetaSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    this->_M_part_sum += __term;
	    this->_M_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _BrezinskiThetaSum&
      operator-=(value_type __term)
      { return operator+=(-__term); }

      /// Return true if the sum converged.
      operator
      bool() const
      { return !this->_M_converged; }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms; }

      ///  Reset the sum to it's initial state.
      _BrezinskiThetaSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_arj.clear();
	this->_M_sum = value_type{0};
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _BrezinskiThetaSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      void _M_update();

      _Sum _M_part_sum;
      std::vector<value_type> _M_arj;
      value_type _M_sum;
      bool _M_converged;
    };

  /**
   * Single step of Brezinski's Theta transformation.
   */
  template<typename _Sum>
    void
    _BrezinskiThetaSum<_Sum>::_M_update()
    {
      using _Tp = value_type;
      constexpr auto _S_huge = __gnu_cxx::__root_max(_Tp{5}); // 1.0e+60
      constexpr auto _S_tiny = __gnu_cxx::__root_min(_Tp{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      auto& __arj = this->_M_arj;

      __arj.push_back(__s_n);
      if (__n < 3)
	this->_M_sum = __s_n;
      else
	{
	  auto __lmax = __n / 3;
	  auto __m = __n;
	  for (auto __l = 1; __l <= __lmax; ++__l)
	    {
	      __m -= 3;
	      auto __diff0 = __arj[__m + 1] - __arj[__m];
	      auto __diff1 = __arj[__m + 2] - __arj[__m + 1];
	      auto __diff2 = __arj[__m + 3] - __arj[__m + 2];
	      auto __denom = __diff2 * (__diff1 - __diff0)
			   - __diff0 * (__diff2 - __diff1);
	      if (std::abs(__denom) < _S_tiny)
		__arj[__m] = _S_huge;
	      else
		__arj[__m] = __arj[__m + 1]
			   - __diff0 * __diff1 * (__diff2 - __diff1) / __denom;
	    }
	  this->_M_sum = __arj[__n % 3];
	}
    }

  /**
   * The Levin's summation process.
   */
  template<typename _Sum>
    class _LevinSum
    {
    public:

      using value_type = typename _Sum::value_type;

      ///  Default constructor.
      _LevinSum(value_type __beta = value_type{1})
      : _M_part_sum{_Sum{}}, _M_e{}, _M_beta{__beta},
	_M_sum{}, _M_converged{false}
      { }

      /// Get the beta parameter.
      value_type
      beta() const
      { return this->_M_beta; }

      /// Set the beta parameter.
      _LevinSum&
      beta(value_type __beta)
      {
	this->_M_beta = __beta;
	return *this;
      }

      /// Add a new term to the sum.
      _LevinSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    this->_M_part_sum += __term;
	    this->_M_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _LevinSum&
      operator-=(value_type __term)
      { return operator+=(-__term); }

      /// Return true if the sum converged.
      operator
      bool() const
      { return !this->_M_converged; }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms; }

      ///  Reset the sum to it's initial state.
      _LevinSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_e.clear();
	this->_M_sum = value_type{0};
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _LevinSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      void _M_update(value_type __r_n);

      _Sum _M_part_sum;
      std::vector<value_type> _M_e;
      value_type _M_beta;
      value_type _M_sum;
      bool _M_converged;
    };

  /**
   * One step of Levin's summation process.
   */
  template<typename _Sum>
    void
    _LevinSum<_Sum>::_M_update(value_type __r_n)
    {
      using _Tp = value_type;
      constexpr auto _S_huge = __gnu_cxx::__root_max(_Tp{5}); // 1.0e+60
      constexpr auto _S_tiny = __gnu_cxx::__root_min(_Tp{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      const auto __beta = this->_M_beta;
      auto& __anum = this->_M_anum;
      auto& __aden = this->_M_aden;

      __anum.push_back(__s_n / __r_n);
      __aden.push_back(value_type{1} / __r_n);
      if (__n > 0)
	{
	  __anum[__n - 1] = __anum[__n] - __anum[__n - 1];
	  __aden[__n - 1] = __aden[__n] - __aden[__n - 1];
	  if (__n > 1)
	    {
	      auto __bn1 = __beta + _Tp(__n - 1);
	      auto __bn2 = __beta + _Tp(__n);
	      auto __coef = __bn1 / __bn2;
	      auto __coefp = _Tp{1};
	      for (auto __j = 2; __j <= __n; ++__j, __coefp *= __coef)
		{
		  auto __fact = (__beta + _Tp(__n - __j)) * __coefp / __bn2;
		  __anum[__n - __j] = __anum[__n - __j + 1]
				    - __fact * __anum[__n - __j];
		  __aden[__n - __j] = __aden[__n - __j + 1]
				    - __fact * __aden[__n - __j];
		}
	    }
	}
      if (std::abs(__aden[0]) < _S_tiny)
	this->_M_sum = _S_huge;
      else
	this->_M_sum = __anum[0] / __aden[0];
    }

  /**
   * The Levin's summation process.
   */
  template<typename _Sum>
    class _WenigerSum
    {
    public:

      using value_type = typename _Sum::value_type;

      ///  Default constructor.
      _WenigerSum(value_type __beta = value_type{1})
      : _M_part_sum{_Sum{}}, _M_e{}, _M_beta{__beta},
	_M_sum{}, _M_converged{false}
      { }

      /// Get the beta parameter.
      value_type
      beta() const
      { return this->_M_beta; }

      /// Set the beta parameter.
      _WenigerSum&
      beta(value_type __beta)
      {
	this->_M_beta = __beta;
	return *this;
      }

      /// Add a new term to the sum.
      _WenigerSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    this->_M_part_sum += __term;
	    this->_M_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _WenigerSum&
      operator-=(value_type __term)
      { return operator+=(-__term); }

      /// Return true if the sum converged.
      operator
      bool() const
      { return !this->_M_converged; }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms; }

      ///  Reset the sum to it's initial state.
      _WenigerSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_e.clear();
	this->_M_sum = value_type{0};
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _WenigerSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      void _M_update(value_type __r_n);

      _Sum _M_part_sum;
      std::vector<value_type> _M_e;
      value_type _M_beta;
      value_type _M_sum;
      bool _M_converged;
    };

  /**
   * One step of Levin's summation process.
   */
  template<typename _Sum>
    void
    _WenigerSum<_Sum>::_M_update(value_type __r_n)
    {
      using _Tp = value_type;
      constexpr auto _S_huge = __gnu_cxx::__root_max(_Tp{5}); // 1.0e+60
      constexpr auto _S_tiny = __gnu_cxx::__root_min(_Tp{5}); // 1.0e-60;

      const auto __n = this->_M_part_sum.num_terms() - 1;
      const auto __s_n = this->_M_part_sum();
      const auto __beta = this->_M_beta;
      auto& __anum = this->_M_anum;
      auto& __aden = this->_M_aden;

      __anum.push_back(__s_n / __r_n);
      __aden.push_back(value_type{1} / __r_n);
      if (__n > 0)
	{
	  __anum[__n - 1] = __anum[__n] - __anum[__n - 1];
	  __aden[__n - 1] = __aden[__n] - __aden[__n - 1];
	  if (__n > 1)
	    {
	      auto __bn1 = __beta + _Tp(__n - 1);
	      auto __bn2 = __beta + _Tp(__n);
	      auto __coef = __bn1 / __bn2;
	      auto __coefp = _Tp{1}
	      for (auto __j = 2; __j <= __n; ++__j, __coefp *= __coef)
		{
		  auto __fact = (__beta + _Tp(__n - __j))
			      * __coefp / __bn2;
		  __anum[__n - __j] = __anum[__n - __j + 1]
				    - __fact * __anum[__n - __j];
		  __aden[__n - __j] = __aden[__n - __j + 1]
				    - __fact * __aden[__n - __j];
		}
	    }
	}
      if (std::abs(__aden[0]) < _S_tiny)
	this->_M_sum = _S_huge;
      else
	this->_M_sum = __anum[0] / __aden[0];
    }

template<typename Tp>
  void
  test()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();

    std::__detail::_BasicSum<Tp> BS;
    _AitkenDeltaSqaredSum<std::__detail::_BasicSum<Tp>> ABS;
    _AitkenDeltaSqaredSum<std::__detail::_KahanSum<Tp>> AKS;
    _WinnEpsilonSum<std::__detail::_BasicSum<Tp>> WBS;
    _WinnEpsilonSum<std::__detail::_KahanSum<Tp>> WKS;
    _BrezinskiThetaSum<std::__detail::_BasicSum<Tp>> BBS;
    //_LevinSum<std::__detail::_BasicSum<Tp>> LS;

    auto s = Tp{1.2};
    auto zeta = Tp{5.591582441177750776536563193423143277642L};
    std::cout << "\n\nzeta(1.2) = 5.59158244117775077653\n";
    for (auto k = 0; k < 100; ++k)
      {
	auto term = std::pow(k + 1, -s);
	BS += term;
	ABS += term;
	AKS += term;
	WBS += term;
	WKS += term;
	BBS += term;
	//LS += term;
	std::cout << std::setw(w) << k
		  << std::setw(w) << BS()
		  << std::setw(w) << ABS()
		  << std::setw(w) << AKS()
		  << std::setw(w) << WBS()
		  << std::setw(w) << WKS()
		  << std::setw(w) << BBS()
	  	  //<< std::setw(w) << LS()
		  << '\n';
      }

    // 2F0(1,1;;-1/3)
    std::cout << "\n\n2F0(1,1;;-1/3) = 0.78625122076596\n";
    auto a = Tp{1};
    auto b = Tp{1};
    auto z = Tp{-1} / Tp{3};
    auto term = Tp{1};
    BS.reset(term);
    ABS.reset(term);
    AKS.reset(term);
    WBS.reset(term);
    WKS.reset(term);
    BBS.reset(term);
    //LS.reset(term);
    for (auto k = 1; k < 100; ++k)
      {
	std::cout << std::setw(w) << (k - 1)
		  << std::setw(w) << BS()
		  << std::setw(w) << ABS()
		  << std::setw(w) << AKS()
		  << std::setw(w) << WBS()
		  << std::setw(w) << WKS()
		  << std::setw(w) << BBS()
		  //<< std::setw(w) << LS()
		  << '\n';
	term *= (a + k - 1) * (b + k - 1) * z / k;
	BS += term;
	ABS += term;
	AKS += term;
	WBS += term;
	WKS += term;
	BBS += term;
	//LS += term;
      }
  }

int
main()
{
  //std::cout << "\nfloat\n=====\n\n";
  //test<float>();

  std::cout << "\ndouble\n======\n";
  test<double>();

  std::cout << "\nlong double\n===========\n";
  test<long double>();

  //std::cout << "\n__float128\n===========\n";
  //test<__float128>();
}
