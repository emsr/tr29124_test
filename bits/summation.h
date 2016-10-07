// math special functions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/summation.h
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SUMMATION_H
#define _GLIBCXX_BITS_SUMMATION_H 1

#pragma GCC visibility push(default)

#include <complex>
#include <vector>
#include <array>
#include <bits/c++config.h>
#include <bits/complex_util.h> // for complex __isnan, __isinf

#pragma GCC system_header

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A Sum is default constructable
   * It has operator+=(Tp)
   * It has operator-=(Tp)
   * It has operator() const -> _Tp
   * It has converged() const -> Tp
   * It has operator bool() const // !converged() i.e. still needs work.
   * It has num_terms() const -> std::size_t
   * It has term() const -> _Tp // last term
   * 
   */

  /**
   * This is a basic naive sum.
   */
  template<typename _Tp>
    class _BasicSum
    {
    public:

      using value_type = _Tp;

      ///  Default constructor.
      _BasicSum()
      : _M_sum{}, _M_term{}, _M_num_terms{0}, _M_converged{false}
      { }

      /// Add a new term to the sum.
      _BasicSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    if (std::__detail::__isnan(__term))
	      std::__throw_runtime_error(__N("_BasicSum: bad term"));
	    if (std::__detail::__isinf(__term))
	      std::__throw_runtime_error(__N("_BasicSum: infinite term"));
	    ++this->_M_num_terms;
	    this->_M_term = __term;
	    this->_M_sum += this->_M_term;
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _BasicSum&
      operator-=(value_type __term)
      { return this->operator+=(-__term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->_M_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_num_terms; }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->_M_term; }

      ///  Reset the sum to it's initial state.
      _BasicSum&
      reset()
      {
	this->_M_sum = value_type{};
	this->_M_term = value_type{};
	this->_M_num_terms = 0;
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _BasicSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      value_type _M_sum;
      value_type _M_term;
      std::size_t _M_num_terms;
      bool _M_converged;
    };

  /**
   * This is a Kahan sum which tries to account for roundoff error.
   */
  template<typename _Tp>
    class _KahanSum
    {
    public:

      using value_type = _Tp;

      ///  Default constructor.
      _KahanSum()
      : _M_sum{}, _M_term{}, _M_temp{}, _M_num_terms{0}, _M_converged{false}
      { }

      /// Add a new term to the sum.
      _KahanSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    if (std::__detail::__isnan(__term))
	      std::__throw_runtime_error(__N("_KahanSum: bad term"));
	    if (std::__detail::__isinf(__term))
	      std::__throw_runtime_error(__N("_KahanSum: infinite term"));
	    ++this->_M_num_terms;
	    this->_M_term = __term - this->_M_temp;
	    this->_M_temp = this->_M_sum;
	    this->_M_sum += this->_M_term;
	    this->_M_temp = this->_M_term - (this->_M_sum - this->_M_temp);
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _KahanSum&
      operator-=(value_type __term)
      { return this->operator+=(-__term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->_M_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_num_terms; }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->_M_term; }

      ///  Reset the sum to it's initial state.
      _KahanSum&
      reset()
      {
	this->_M_sum = value_type{};
	this->_M_term = value_type{};
	this->_M_temp = value_type{};
	this->_M_num_terms = 0;
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _KahanSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      value_type _M_sum;
      value_type _M_term;
      value_type _M_temp;
      std::size_t _M_num_terms;
      bool _M_converged;
    };

  /**
   * 
   */
  template<typename _Tp>
    class _VanWijngaardenSum
    {
    public:

      using value_type = _Tp;

      ///  Default constructor.
      _VanWijngaardenSum()
      : _M_sum{}, _M_term{}, _M_delta{}, _M_num_terms{0}, _M_converged{false}
      { }

      /// Add a new term to the sum.
      _VanWijngaardenSum&
      operator+=(value_type __term);

      /// Subtract a new term from the sum.
      _VanWijngaardenSum&
      operator-=(value_type __term)
      { return this->operator+=(-__term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->_M_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_num_terms; }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->_M_term; }

      ///  Reset the sum to it's initial state.
      _VanWijngaardenSum&
      reset()
      {
	this->_M_sum = value_type{};
	this->_M_term = value_type{};
	this->_M_delta.clear();
	this->_M_num_terms = 0;
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _VanWijngaardenSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      value_type _M_sum;
      value_type _M_term;
      std::vector<value_type> _M_delta;
      std::size_t _M_num_terms;
      bool _M_converged;
    };

  /**
   *  This performs a series compression on a monotone series - converting
   *  it to an alternating series - for the regular van Wijngaarden sum.
   *  ADL for ctors anyone?  I'd like to put a lambda in here*
   */
  template<typename _TermFn>
    class _VanWijngaardenCompressor
    {
    public:

      _VanWijngaardenCompressor(_TermFn __term_fn)
      : _M_term_fn{__term_fn}
      { }

      auto
      operator[](std::size_t __j) const;

    private:

      _TermFn _M_term_fn;
    };

  /**
   * The Aitken's delta-squared summation process.
   */
  template<typename _Sum>
    class _AitkenDeltaSquaredSum
    {
    public:

      using value_type = typename _Sum::value_type;

      ///  Default constructor.
      _AitkenDeltaSquaredSum()
      : _M_part_sum{_Sum{}}, _M_a{}, _M_sum{}, _M_converged{false}
      { }

      /// Add a new term to the sum.
      _AitkenDeltaSquaredSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    if (std::__detail::__isnan(__term))
	      std::__throw_runtime_error(__N("_AitkenDeltaSquaredSum: bad term"));
	    if (std::__detail::__isinf(__term))
	      std::__throw_runtime_error(__N("_AitkenDeltaSquaredSum: infinite term"));
	    this->_M_part_sum += __term;
	    this->_M_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _AitkenDeltaSquaredSum&
      operator-=(value_type __term)
      { return this->operator+=(-__term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->_M_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->_M_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      _AitkenDeltaSquaredSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_a.clear();
	this->_M_sum = value_type{};
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _AitkenDeltaSquaredSum&
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

  /**
   * The Winn epsilon summation process.
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
	    if (std::__detail::__isnan(__term))
	      std::__throw_runtime_error(__N("_WinnEpsilonSum: bad term"));
	    if (std::__detail::__isinf(__term))
	      std::__throw_runtime_error(__N("_WinnEpsilonSum: infinite term"));
	    this->_M_part_sum += __term;
	    this->_M_update();
	  }
	return *this;
      }
  
      /// Subtract a new term from the sum.
      _WinnEpsilonSum&
      operator-=(value_type __term)
      { return this->operator+=(-__term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->_M_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->_converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->_M_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      _WinnEpsilonSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_e.clear();
	this->_M_sum = value_type{};
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
	    if (std::__detail::__isnan(__term))
	      std::__throw_runtime_error(__N("_BrezinskiThetaSum: bad term"));
	    if (std::__detail::__isinf(__term))
	      std::__throw_runtime_error(__N("_BrezinskiThetaSum: infinite term"));
	    this->_M_part_sum += __term;
	    this->_M_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _BrezinskiThetaSum&
      operator-=(value_type __term)
      { return this->operator+=(-__term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->_M_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->_converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->_M_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      _BrezinskiThetaSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_arj.clear();
	this->_M_sum = value_type{};
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


  // These sequence transformations depend on the provision of remainder
  // estimates. The update methods do not depend on the remainder model
  // and could be provided in CRTP derived classes.


  /**
   * 
   */
  template<typename _Tp>
    struct _RemainderTerm
    {
      using value_type = _Tp;

      value_type term = value_type{};
      value_type remainder = value_type{};
    };


  /**
   * This class implements an explicit remainder model where
   * the caller supplies the corresponding ramainder estimate after the
   * corresponding term.
   */
  template<typename _Tp>
    class _ExplicitRemainderModel
    {
    public:

      using value_type = _Tp;

      constexpr _ExplicitRemainderModel()
      : _M_n{0}
      { }

      void
      operator<<(value_type __term)
      {
	if (this->_M_n < 2)
	  {
	    this->_M_term[this->_M_n] = __term;
	    ++this->_M_n;
	  }
	else
	  /* error */;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return this->_M_n == 2; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      _RemainderTerm<value_type>
      operator()()
      {
	this->_M_n = 0;
	return _RemainderTerm<value_type>{this->_M_term[0], this->_M_term[1]};
      }

      _ExplicitRemainderModel&
      reset()
      {
	this->_M_n = 0;
	return *this;
      }

    private:

      int _M_n;
      std::array<value_type, 2> _M_term;
    };


  /**
   * This class implements the Levin U remainder model.
   * The remainder for term @f$ s_n @f$ is @f$ (n + \beta)a_n @f$
   * or in terms of the backward difference of the partial sums
   * @f$ (n + \beta)(s_n - s_{n-1}) @f$.
   */
  template<typename _Tp>
    class _URemainderModel
    {
    public:

      using value_type = _Tp;

      constexpr _URemainderModel(value_type __beta)
      : _M_n{0}, _M_term{}, _M_beta{__beta}
      { }

      constexpr void
      operator<<(value_type __term)
      {
	this->_M_term = __term;
	++this->_M_n;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return this->_M_n >= 2; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      _RemainderTerm<value_type>
      constexpr operator()()
      {
	this->_M_ok = false;
	return _RemainderTerm<value_type>{this->_M_term,
				(this->_M_n + this->_M_beta) * this->_M_term};
      }

      _URemainderModel&
      reset()
      {
	this->_M_n = 0;
	return *this;
      }

    private:

      int _M_n;
      value_type _M_term;
      value_type _M_beta;
    };


  /**
   * This class implements the Levin T remainder model.
   * The remainder for term @f$ s_n @f$ is the current term @f$ a_n @f$
   * or the backward difference of the partial sums @f$ s_n - s_{n-1} @f$.
   */
  template<typename _Tp>
    class _TRemainderModel
    {
    public:

      using value_type = _Tp;

      constexpr _TRemainderModel()
      : _M_ok{false}
      { }

      constexpr void
      operator<<(value_type __term)
      {
	if (!this->_M_ok)
	  {
	    this->_M_term = __term;
	    this->_M_ok = true;
	  }
	else
	  /* error */;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return this->_M_ok; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      _RemainderTerm<value_type>
      constexpr operator()()
      {
	this->_M_ok = false;
	return _RemainderTerm<value_type>{this->_M_term, this->_M_term};
      }

      _TRemainderModel&
      reset()
      {
	this->_M_ok = false;
	return *this;
      }

    private:

      bool _M_ok;
      value_type _M_term;
    };


  /**
   * This class implements the Levin D remainder model.
   * The remainder for term @f$ s_n @f$ is simply the next term @f$ a_{n+1} @f$
   * or the forward difference of the partial sums @f$ s_{n+1} - s_n @f$.
   */
  template<typename _Tp>
    class _DRemainderModel
    {
    public:

      using value_type = _Tp;

      constexpr _DRemainderModel()
      : _M_n{0}
      { }

      void
      operator<<(value_type __term)
      {
	this->_M_term[(this->_M_n) % 2] = __term;
	++this->_M_n;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return this->_M_n >= 2; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      _RemainderTerm<value_type>
      operator()()
      {
	return _RemainderTerm<value_type>{this->_M_term[(this->_M_n) % 2],
					  this->_M_term[(this->_M_n + 1) % 2]};
      }

      _DRemainderModel&
      reset()
      {
	this->_M_n = 0;
	return *this;
      }

    private:

      int _M_n;
      std::array<value_type, 2> _M_term;
    };


  /**
   * This class implements the Levin V remainder model.
   * The remainder for term @f$ s_n @f$ is
   * @f[
   *   r_n = \frac{a_n a_{n+1}}{a_n - a_{n+1}}
   * @f]
   */
  template<typename _Tp>
    class _VRemainderModel
    {
    public:

      using value_type = _Tp;

      constexpr _VRemainderModel()
      : _M_n{0}
      { }

      void
      operator<<(value_type __term)
      {
	this->_M_term[(this->_M_n) % 2] = __term;
	++this->_M_n;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return this->_M_n >= 2; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      _RemainderTerm<value_type>
      operator()()
      {
	auto __anp1 = this->_M_term[(this->_M_n + 1) % 2];
	auto __an = this->_M_term[(this->_M_n) % 2];
	return _RemainderTerm<value_type>{__an,
					  __an * __anp1 / (__an - __anp1)};
      }

      _VRemainderModel&
      reset()
      {
	this->_M_n = 0;
	return *this;
      }

    private:

      int _M_n;
      std::array<value_type, 2> _M_term;
    };


  /**
   * The Levin summation process.
   */
  template<typename _Sum,
	   typename _RemainderModel
		      = _ExplicitRemainderModel<typename _Sum::value_type>>
    class _LevinSum
    {
    public:

      using value_type = typename _Sum::value_type;

      ///  Default constructor.
      _LevinSum(value_type __beta = value_type{1})
      : _M_part_sum{_Sum{}}, _M_num{}, _M_den{},
	_M_beta{__beta},
	_M_sum{}, _M_converged{false}, _M_rem_mdl{}
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
	    if (std::__detail::__isnan(__term))
	      std::__throw_runtime_error(__N("_LevinSum: bad term"));
	    if (std::__detail::__isinf(__term))
	      std::__throw_runtime_error(__N("_LevinSum: infinite term"));
	    this->_M_rem_mdl << __term;
	    if (this->_M_rem_mdl.ready())
	      {
		auto __thing = this->_M_rem_mdl();
		this->_M_part_sum += __thing.term;
		this->_M_update(__thing.remainder);
	      }
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _LevinSum&
      operator-=(value_type __term)
      { return this->operator+=(-__term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->_M_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->_M_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      ///  The beta parameter is unchanged.
      _LevinSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_num.clear();
	this->_M_den.clear();
	this->_M_sum = value_type{};
	this->_M_converged = false;
	this->_M_rem_mdl.reset();
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

      const _RemainderModel&
      _M_self() const
      { return static_cast<const _RemainderModel&>(*this); }

      _RemainderModel&
      _M_self()
      { return static_cast<_RemainderModel&>(*this); }

      void _M_update(value_type __r_n);

      _Sum _M_part_sum;
      std::vector<value_type> _M_num;
      std::vector<value_type> _M_den;
      value_type _M_beta;
      value_type _M_sum;
      bool _M_converged;
      _RemainderModel _M_rem_mdl;
    };

  /**
   * The Weniger summation process.
   */
  template<typename _Sum,
	   typename _RemainderModel
		      = _ExplicitRemainderModel<typename _Sum::value_type>>
    class _WenigerSum
    {
    public:

      using value_type = typename _Sum::value_type;

      ///  Default constructor.
      _WenigerSum(value_type __beta = value_type{1})
      : _M_part_sum{_Sum{}}, _M_num{}, _M_den{},
	_M_beta{__beta},
	_M_sum{}, _M_converged{false}, _M_rem_mdl{}
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
	    if (std::__detail::__isnan(__term))
	      std::__throw_runtime_error(__N("_WenigerSum: bad term"));
	    if (std::__detail::__isinf(__term))
	      std::__throw_runtime_error(__N("_WenigerSum: infinite term"));
	    this->_M_rem_mdl << __term;
	    if (this->_M_rem_mdl.ready())
	      {
		auto __thing = this->_M_rem_mdl();
		this->_M_part_sum += __thing.term;
		this->_M_update(__thing.remainder);
	      }
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _WenigerSum&
      operator-=(value_type __term)
      { return this->operator+=(-__term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->_M_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->_M_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      ///  The beta parameter is unchanged.
      _WenigerSum&
      reset()
      {
	this->_M_part_sum.reset();
	this->_M_num.clear();
	this->_M_den.clear();
	this->_M_sum = value_type{};
	this->_M_converged = false;
	this->_M_rem_mdl.reset();
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
      std::vector<value_type> _M_num;
      std::vector<value_type> _M_den;
      value_type _M_beta;
      value_type _M_sum;
      bool _M_converged;
      _RemainderModel _M_rem_mdl;
    };

  // Specializations for specific remainder models.

  /**
   * The Levin T summation process.
   */
  template<typename _Sum>
    using _LevinTSum
      = _LevinSum<_Sum, _TRemainderModel<typename _Sum::value_type>>;

  /**
   * The Levin D summation process.
   */
  template<typename _Sum>
    using _LevinDSum
      = _LevinSum<_Sum, _DRemainderModel<typename _Sum::value_type>>;

  /**
   * The Levin V summation process.
   */
  template<typename _Sum>
    using _LevinUSum
      = _LevinSum<_Sum, _VRemainderModel<typename _Sum::value_type>>;

  /**
   * The Weniger Tau summation process.
   */
  template<typename _Sum>
    using _WenigerTauSum
      = _WenigerSum<_Sum, _TRemainderModel<typename _Sum::value_type>>;

  /**
   * The Weniger Delta summation process.
   */
  template<typename _Sum>
    using _WenigerDeltaSum
      = _WenigerSum<_Sum, _DRemainderModel<typename _Sum::value_type>>;

  /**
   * The Weniger Phi summation process.
   */
  template<typename _Sum>
    using _WenigerPhiSum
      = _WenigerSum<_Sum, _VRemainderModel<typename _Sum::value_type>>;

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#pragma GCC visibility pop

#include <bits/summation.tcc>

#endif // _GLIBCXX_BITS_SUMMATION_H
