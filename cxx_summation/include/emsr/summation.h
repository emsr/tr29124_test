
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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

#ifndef SUMMATION_H
#define SUMMATION_H 1

#include <complex>
#include <vector>
#include <array>
#include <bits/c++config.h>
#include <emsr/complex_util.h> // for complex std::isnan, std::isinf

namespace emsr
{

  /**
   * A Sum is default constructable
   * It has operator+=(Tp)
   * It has operator-=(Tp)
   * It has operator() const -> Tp
   * It has converged() const -> Tp
   * It has operator bool() const // !converged() i.e. still needs work.
   * It has num_terms() const -> std::size_t
   * It has term() const -> Tp // last term
   * 
   */

  /**
   * This is a basic naive sum.
   */
  template<typename Tp>
    class BasicSum
    {
    public:

      using value_type = Tp;

      ///  Default constructor.
      BasicSum()
      : m_sum{}, m_term{}, m_num_terms{0}, m_converged{false}
      { }

      /// Add a new term to the sum.
      BasicSum&
      operator+=(value_type term)
      {
	if (!this->m_converged)
	  {
	    if (std::isnan(term))
	      throw std::runtime_error("BasicSum: bad term");
	    if (std::isinf(term))
	      throw std::runtime_error("BasicSum: infinite term");
	    ++this->m_num_terms;
	    this->m_term = term;
	    this->m_sum += this->m_term;
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      BasicSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_num_terms; }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_term; }

      ///  Reset the sum to it's initial state.
      BasicSum&
      reset()
      {
	this->m_sum = value_type{};
	this->m_term = value_type{};
	this->m_num_terms = 0;
	this->m_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      BasicSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    private:

      value_type m_sum;
      value_type m_term;
      std::size_t m_num_terms;
      bool m_converged;
    };

  /**
   * This is a Kahan sum which tries to account for roundoff error.
   */
  template<typename Tp>
    class KahanSum
    {
    public:

      using value_type = Tp;

      ///  Default constructor.
      KahanSum()
      : m_sum{}, m_term{}, m_temp{}, m_num_terms{0}, m_converged{false}
      { }

      /// Add a new term to the sum.
      KahanSum&
      operator+=(value_type term)
      {
	if (!this->m_converged)
	  {
	    if (std::isnan(term))
	      throw std::runtime_error("KahanSum: bad term");
	    if (std::isinf(term))
	      throw std::runtime_error("KahanSum: infinite term");
	    ++this->m_num_terms;
	    this->m_term = term - this->m_temp;
	    this->m_temp = this->m_sum;
	    this->m_sum += this->m_term;
	    this->m_temp = this->m_term - (this->m_sum - this->m_temp);
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      KahanSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_num_terms; }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_term; }

      ///  Reset the sum to it's initial state.
      KahanSum&
      reset()
      {
	this->m_sum = value_type{};
	this->m_term = value_type{};
	this->m_temp = value_type{};
	this->m_num_terms = 0;
	this->m_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      KahanSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    private:

      value_type m_sum;
      value_type m_term;
      value_type m_temp;
      std::size_t m_num_terms;
      bool m_converged;
    };

  /**
   * 
   */
  template<typename Tp>
    class VanWijngaardenSum
    {
    public:

      using value_type = Tp;

      ///  Default constructor.
      VanWijngaardenSum(std::size_t start_term = 0u)
      : m_sum{}, m_term{}, m_delta{}, m_num_terms{0},
	m_start_term{start_term},
	m_converged{false}
      { }

      /// Add a new term to the sum.
      VanWijngaardenSum&
      operator+=(value_type term);

      /// Subtract a new term from the sum.
      VanWijngaardenSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_num_terms; }

      /// Return the number of initial terms to add to the sum before
      /// switching to the vanWijngaarden algorithm.
      std::size_t
      start_term() const
      { return this->m_start_term; }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_term; }

      ///  Reset the sum to it's initial state.
      VanWijngaardenSum&
      reset()
      {
	this->m_sum = value_type{};
	this->m_term = value_type{};
	this->m_delta.clear();
	this->m_num_terms = 0;
	this->m_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      VanWijngaardenSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    private:

      value_type m_sum;
      value_type m_term;
      std::vector<value_type> m_delta;
      std::size_t m_num_terms;
      std::size_t m_start_term;
      bool m_converged;
    };

  /**
   *  This performs a series compression on a monotone series - converting
   *  it to an alternating series - for the regular van Wijngaarden sum.
   *  ADL for ctors anyone?  I'd like to put a lambda in here*
   */
  template<typename TermFn>
    class VanWijngaardenCompressor
    {
    private:

      TermFn m_term_fn;

    public:

      using return_t = decltype(m_term_fn.operator()(std::size_t{}));

      VanWijngaardenCompressor(TermFn term_fn)
      : m_term_fn{term_fn}
      { }

      return_t
      operator[](std::size_t j) const;
    };

  /**
   * The Aitken's delta-squared summation process.
   */
  template<typename Sum>
    class AitkenDeltaSquaredSum
    {
    public:

      using value_type = typename Sum::value_type;

      ///  Default constructor.
      AitkenDeltaSquaredSum()
      : m_part_sum{Sum{}}, m_a{}, m_sum{}, m_converged{false}
      { }

      /// Add a new term to the sum.
      AitkenDeltaSquaredSum&
      operator+=(value_type term)
      {
	if (!this->m_converged)
	  {
	    if (std::isnan(term))
	      throw std::runtime_error("AitkenDeltaSquaredSum: bad term");
	    if (std::isinf(term))
	      throw std::runtime_error("AitkenDeltaSquaredSum: infinite term");
	    this->m_part_sum += term;
	    this->m_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      AitkenDeltaSquaredSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      AitkenDeltaSquaredSum&
      reset()
      {
	this->m_part_sum.reset();
	this->m_a.clear();
	this->m_sum = value_type{};
	this->m_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      AitkenDeltaSquaredSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    private:

      void m_update();

      Sum m_part_sum;
      std::vector<value_type> m_a;
      value_type m_sum;
      bool m_converged;
    };

  /**
   * The Winn epsilon summation process.
   */
  template<typename Sum>
    class WinnEpsilonSum
    {
    public:

      using value_type = typename Sum::value_type;

      ///  Default constructor.
      WinnEpsilonSum()
      : m_part_sum{Sum{}}, m_e{}, m_sum{}, m_converged{false}
      { }

      /// Add a new term to the sum.
      WinnEpsilonSum&
      operator+=(value_type term)
      {
	if (!this->m_converged)
	  {
	    if (std::isnan(term))
	      throw std::runtime_error("WinnEpsilonSum: bad term");
	    if (std::isinf(term))
	      throw std::runtime_error("WinnEpsilonSum: infinite term");
	    this->m_part_sum += term;
	    this->m_update();
	  }
	return *this;
      }
  
      /// Subtract a new term from the sum.
      WinnEpsilonSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->_converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      WinnEpsilonSum&
      reset()
      {
	this->m_part_sum.reset();
	this->m_e.clear();
	this->m_sum = value_type{};
	this->m_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      WinnEpsilonSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    private:

      void m_update();

      Sum m_part_sum;
      std::vector<value_type> m_e;
      value_type m_sum;
      bool m_converged;
    };

  /**
   * The Brezinski theta summation process.
   */
  template<typename Sum>
    class BrezinskiThetaSum
    {
    public:

      using value_type = typename Sum::value_type;

      ///  Default constructor.
      BrezinskiThetaSum()
      : m_part_sum{Sum{}}, m_arj{}, m_sum{}, m_converged{false}
      { }

      /// Add a new term to the sum.
      BrezinskiThetaSum&
      operator+=(value_type term)
      {
	if (!this->m_converged)
	  {
	    if (std::isnan(term))
	      throw std::runtime_error("BrezinskiThetaSum: bad term");
	    if (std::isinf(term))
	      throw std::runtime_error("BrezinskiThetaSum: infinite term");
	    this->m_part_sum += term;
	    this->m_update();
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      BrezinskiThetaSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->_converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      BrezinskiThetaSum&
      reset()
      {
	this->m_part_sum.reset();
	this->m_arj.clear();
	this->m_sum = value_type{};
	this->m_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      BrezinskiThetaSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    private:

      void m_update();

      Sum m_part_sum;
      std::vector<value_type> m_arj;
      value_type m_sum;
      bool m_converged;
    };


  // These sequence transformations depend on the provision of remainder
  // estimates. The update methods do not depend on the remainder model
  // and could be provided in CRTP derived classes.


  /**
   * 
   */
  template<typename Tp>
    struct RemainderTerm
    {
      using value_type = Tp;

      explicit RemainderTerm(Tp val, Tp rem)
      : term(val), remainder(rem)
      { }

      value_type term = value_type{};
      value_type remainder = value_type{};
    };


  /**
   * This class implements an explicit remainder model where
   * the caller supplies the corresponding ramainder estimate after the
   * corresponding term.
   */
  template<typename Tp>
    class ExplicitRemainderModel
    {
    public:

      using value_type = Tp;

      constexpr ExplicitRemainderModel()
      : m_n{0}, m_term{{Tp{}, Tp{}}}
      { }

      void
      operator<<(value_type term)
      {
	if (this->m_n < 2)
	  {
	    this->m_term[this->m_n] = term;
	    ++this->m_n;
	  }
	else
	  {
            /* error */;
          }
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return this->m_n == 2; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      RemainderTerm<value_type>
      operator()()
      {
	this->m_n = 0;
	return RemainderTerm<value_type>{this->m_term[0], this->m_term[1]};
      }

      ExplicitRemainderModel&
      reset()
      {
	this->m_n = 0;
	return *this;
      }

    private:

      int m_n;
      std::array<value_type, 2> m_term;
    };


  /**
   * This class implements the Levin U remainder model.
   * The remainder for term @f$ s_n @f$ is @f$ (n + \beta)a_n @f$
   * or in terms of the backward difference of the partial sums
   * @f$ (n + \beta)(s_n - s_{n-1}) @f$.
   */
  template<typename Tp>
    class URemainderModel
    {
    public:

      using value_type = Tp;

      constexpr URemainderModel(value_type beta)
      : m_n{0}, m_term{}, m_beta{beta}
      { }

      constexpr void
      operator<<(value_type term)
      {
	this->m_term = term;
	++this->m_n;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return true; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      RemainderTerm<value_type>
      constexpr operator()()
      {
	return RemainderTerm<value_type>{this->m_term,
				(this->m_n + this->m_beta) * this->m_term};
      }

      URemainderModel&
      reset()
      {
	this->m_n = 0;
	return *this;
      }

    private:

      int m_n;
      value_type m_term;
      value_type m_beta;
    };


  /**
   * This class implements the Levin T remainder model.
   * The remainder for term @f$ s_n @f$ is the current term @f$ a_n @f$
   * or the backward difference of the partial sums @f$ s_n - s_{n-1} @f$.
   */
  template<typename Tp>
    class TRemainderModel
    {
    public:

      using value_type = Tp;

      constexpr TRemainderModel()
      : m_term{}
      { }

      constexpr void
      operator<<(value_type term)
      {
	this->m_term = term;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return true; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      RemainderTerm<value_type>
      constexpr operator()()
      { return RemainderTerm<value_type>{this->m_term, this->m_term}; }

      TRemainderModel&
      reset()
      {
	return *this;
      }

    private:

      value_type m_term;
    };


  /**
   * This class implements the Levin D remainder model.
   * The remainder for term @f$ s_n @f$ is simply the next term @f$ a_{n+1} @f$
   * or the forward difference of the partial sums @f$ s_{n+1} - s_n @f$.
   */
  template<typename Tp>
    class DRemainderModel
    {
    public:

      using value_type = Tp;

      constexpr DRemainderModel()
      : m_n{0}, m_term{{Tp{}, Tp{}}}
      { }

      void
      operator<<(value_type term)
      {
	this->m_term[(this->m_n) % 2] = term;
	++this->m_n;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return this->m_n >= 2; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      RemainderTerm<value_type>
      operator()()
      {
	return RemainderTerm<value_type>{this->m_term[(this->m_n) % 2],
					  this->m_term[(this->m_n + 1) % 2]};
      }

      DRemainderModel&
      reset()
      {
	this->m_n = 0;
	return *this;
      }

    private:

      int m_n;
      std::array<value_type, 2> m_term;
    };


  /**
   * This class implements the Levin V remainder model.
   * The remainder for term @f$ s_n @f$ is
   * @f[
   *   r_n = \frac{a_n a_{n+1}}{a_n - a_{n+1}}
   * @f]
   */
  template<typename Tp>
    class VRemainderModel
    {
    public:

      using value_type = Tp;

      constexpr VRemainderModel()
      : m_n{0}, m_term{{Tp{}, Tp{}}}
      { }

      void
      operator<<(value_type term)
      {
	this->m_term[(this->m_n) % 2] = term;
	++this->m_n;
      }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr bool
      ready() const
      { return this->m_n >= 2; }

      // Return true if the remainder model has accumulated enough terms
      // to start work on the sum.
      constexpr operator
      bool() const
      { return this->ready(); }

      RemainderTerm<value_type>
      operator()()
      {
	auto anp1 = this->m_term[(this->m_n + 1) % 2];
	auto an = this->m_term[(this->m_n) % 2];
	return RemainderTerm<value_type>{an,
					  an * anp1 / (an - anp1)};
      }

      VRemainderModel&
      reset()
      {
	this->m_n = 0;
	return *this;
      }

    private:

      int m_n;
      std::array<value_type, 2> m_term;
    };


  /**
   * The Levin summation process.
   */
  template<typename Sum,
	   typename RemainderModel
		      = ExplicitRemainderModel<typename Sum::value_type>>
    class LevinSum
    {
    public:

      using value_type = typename Sum::value_type;

      ///  Default constructor.
      LevinSum(value_type beta = value_type{1})
      : m_part_sum{Sum{}}, m_num{}, m_den{},
	m_beta{beta},
	m_sum{}, m_converged{false}, m_rem_mdl{}
      { }

      /// Get the beta parameter.
      value_type
      beta() const
      { return this->m_beta; }

      /// Set the beta parameter.
      LevinSum&
      beta(value_type beta)
      {
	this->m_beta = beta;
	return *this;
      }

      /// Add a new term to the sum.
      LevinSum&
      operator+=(value_type term)
      {
	if (!this->m_converged)
	  {
	    if (std::isnan(term))
	      throw std::runtime_error("LevinSum: bad term");
	    if (std::isinf(term))
	      throw std::runtime_error("LevinSum: infinite term");
	    this->m_rem_mdl << term;
	    if (this->m_rem_mdl.ready())
	      {
		auto thing = this->m_rem_mdl();
		this->m_part_sum += thing.term;
		this->m_update(thing.remainder);
	      }
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      LevinSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      ///  The beta parameter is unchanged.
      LevinSum&
      reset()
      {
	this->m_part_sum.reset();
	this->m_num.clear();
	this->m_den.clear();
	this->m_sum = value_type{};
	this->m_converged = false;
	this->m_rem_mdl.reset();
	return *this;
      }

      ///  Restart the sum with the first new term.
      LevinSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    private:

      const RemainderModel&
      m_self() const
      { return static_cast<const RemainderModel&>(*this); }

      RemainderModel&
      m_self()
      { return static_cast<RemainderModel&>(*this); }

      void m_update(value_type r_n);

      Sum m_part_sum;
      std::vector<value_type> m_num;
      std::vector<value_type> m_den;
      value_type m_beta;
      value_type m_sum;
      bool m_converged;
      RemainderModel m_rem_mdl;
    };

  // Specializations for specific remainder models.

  /**
   * The Levin U summation process.
  template<typename Sum>
    using LevinUSum
      = LevinSum<Sum, URemainderModel<typename Sum::value_type>>;
   */

  /**
   * The Levin T summation process.
   */
  template<typename Sum>
    using LevinTSum
      = LevinSum<Sum, TRemainderModel<typename Sum::value_type>>;

  /**
   * The Levin D summation process.
   */
  template<typename Sum>
    using LevinDSum
      = LevinSum<Sum, DRemainderModel<typename Sum::value_type>>;

  /**
   * The Levin V summation process.
   */
  template<typename Sum>
    using LevinVSum
      = LevinSum<Sum, VRemainderModel<typename Sum::value_type>>;

  /**
   * The Weniger summation process.
   */
  template<typename Sum,
	   typename RemainderModel
		      = ExplicitRemainderModel<typename Sum::value_type>>
    class WenigerSum
    {
    public:

      using value_type = typename Sum::value_type;

      ///  Default constructor.
      WenigerSum(value_type beta = value_type{1})
      : m_part_sum{Sum{}}, m_num{}, m_den{},
	m_beta{beta},
	m_sum{}, m_converged{false}, m_rem_mdl{}
      { }

      /// Get the beta parameter.
      value_type
      beta() const
      { return this->m_beta; }

      /// Set the beta parameter.
      WenigerSum&
      beta(value_type beta)
      {
	this->m_beta = beta;
	return *this;
      }

      /// Add a new term to the sum.
      WenigerSum&
      operator+=(value_type term)
      {
	if (!this->m_converged)
	  {
	    if (std::isnan(term))
	      throw std::runtime_error("WenigerSum: bad term");
	    if (std::isinf(term))
	      throw std::runtime_error("WenigerSum: infinite term");
	    this->m_rem_mdl << term;
	    if (this->m_rem_mdl.ready())
	      {
		auto thing = this->m_rem_mdl();
		this->m_part_sum += thing.term;
		this->m_update(thing.remainder);
	      }
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      WenigerSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_part_sum.num_terms(); }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_part_sum.term(); }

      ///  Reset the sum to it's initial state.
      ///  The beta parameter is unchanged.
      WenigerSum&
      reset()
      {
	this->m_part_sum.reset();
	this->m_num.clear();
	this->m_den.clear();
	this->m_sum = value_type{};
	this->m_converged = false;
	this->m_rem_mdl.reset();
	return *this;
      }

      ///  Restart the sum with the first new term.
      WenigerSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    private:

      void m_update(value_type r_n);

      Sum m_part_sum;
      std::vector<value_type> m_num;
      std::vector<value_type> m_den;
      value_type m_beta;
      value_type m_sum;
      bool m_converged;
      RemainderModel m_rem_mdl;
    };

  // Specializations for specific remainder models.

  /**
   * The Weniger Upsilon summation process.
  template<typename Sum>
    using WenigerUpsilonSum
      = WenigerSum<Sum, URemainderModel<typename Sum::value_type>>;
   */

  /**
   * The Weniger Tau summation process.
   */
  template<typename Sum>
    using WenigerTauSum
      = WenigerSum<Sum, TRemainderModel<typename Sum::value_type>>;

  /**
   * The Weniger Delta summation process.
   */
  template<typename Sum>
    using WenigerDeltaSum
      = WenigerSum<Sum, DRemainderModel<typename Sum::value_type>>;

  /**
   * The Weniger Phi summation process.
   */
  template<typename Sum>
    using WenigerPhiSum
      = WenigerSum<Sum, VRemainderModel<typename Sum::value_type>>;

} // namespace emsr

#include <emsr/summation.tcc>

#endif // SUMMATION_H
