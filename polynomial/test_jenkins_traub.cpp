/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_jenkins_traub test_jenkins_traub.cpp -lquadmath
*/

#include <vector>
#include <limits>
#include <cmath>

#include "bits/specfun_state.h"
#include "polynomial/polynomial.h"
#include "solver_low_degree.h"

namespace __gnu_cxx
{

template<typename _Real>
  class _JenkinsTraubSolver
  {
  public:

    _JenkinsTraubSolver(const std::vector<_Real>& op);
    _JenkinsTraubSolver(std::vector<_Real>&& op);

    std::vector<solution_t<_Real>> solve();

  private:

    void quadratic(_Real a, _Real b, _Real c,
		   solution_t<_Real> &z_small, solution_t<_Real> &z_large);
    int fxshfr(int l2);
    int iter_quadratic(_Real uu, _Real vv);
    int iter_real(_Real sss, int& iflag);
    int calcsc();
    void next_k_poly(int type);
    std::pair<_Real, _Real> quadratic_coefficients(int type);
    void remquo_quadratic(int n, _Real u, _Real v,
			  std::vector<_Real>& poly, std::vector<_Real>& quot,
			  _Real& a, _Real& b);

    static constexpr auto eta = std::numeric_limits<_Real>::epsilon();
    static constexpr auto base = _Real{std::numeric_limits<_Real>::radix};
    static constexpr auto infty = 1.0e+50;//std::numeric_limits<_Real>::max();
    static constexpr auto small = 1.0e-50;//std::numeric_limits<_Real>::min();

    _Real min_log_deriv = _Real{0.005L};
    int max_iter_real = 10;
    // Epsilon parameter.
    _Real are = this->eta;
    // Epsilon parameter.
    _Real mre = this->eta;

    std::vector<_Real> P;
    std::vector<_Real> P_quot;
    std::vector<_Real> H, H_quot, H_save;
    _Real sr, si;
    _Real u, v;
    _Real a;
    _Real b;
    _Real c;
    _Real d;
    _Real e;
    _Real f;
    _Real g;
    _Real h;
    _Real a1;
    _Real a2;
    _Real a3;
    _Real a6;
    _Real a7;
    solution_t<_Real> z_small;
    solution_t<_Real> z_large;
    int order;
    bool zerok;
    int num_iters = 0;
  };

/**
 * Constructor from input polygon.
 */
template<typename _Real>
  _JenkinsTraubSolver<_Real>::_JenkinsTraubSolver(const std::vector<_Real>& op)
  : P(op)
  {
    if (this->P.size() == 0)
      std::__throw_domain_error("Polynomial degree must be at least 1.");

    // Algorithm fails of the leading coefficient is zero.
    // We could erase leading-order zero coefficients.
    if (this->P[0] == _Real{0})
      std::__throw_domain_error("Leading coefficient must be nonzero.");

    const auto degree = this->P.size() - 1;
    this->order = degree;
    this->P_quot.resize(degree + 1);
    this->H.resize(degree + 1);
    this->H_quot.resize(degree + 1);
    this->H_save.resize(degree + 1);
  }

/**
 * 
 */
template<typename _Real>
  std::vector<solution_t<_Real>>
  _JenkinsTraubSolver<_Real>::solve()
  {
    auto lo = this->small / this->eta;
    // Initialization of constants for shift rotation.
    auto xx = std::sqrt(_Real{0.5});
    auto yy = -xx;
    const auto rot = _Real{94} * _Real{0.017453293};
    const auto cosr = std::cos(rot);
    const auto sinr = std::sin(rot);

    std::vector<solution_t<_Real>> zero;
    zero.reserve(this->P.size());

    // Remove the zeros at the origin, if any.
    while (this->P[this->order] == _Real{0})
      {
	zero.push_back(_Real{0});
	--this->order;
      }
    if (this->order < 1)
      return zero;

    std::vector<_Real> pt(this->order + 1);
    std::vector<_Real> H_temp(this->order + 1);

    while (true)
      {
	// Start the algorithm for one zero.
	this->num_iters = 0;
	if (this->order == 1)
	  {
	    zero.push_back(-this->P[1] / this->P[0]);
	    --this->order;
	    return zero;
	  }
	// Calculate the final zero or pair of zeros.
	if (this->order == 2)
	  {
	    solution_t<_Real> z_small, z_large;
	    this->quadratic(this->P[0], this->P[1], this->P[2], z_small, z_large);
	    if (z_small.index() != 0)
	      {
		zero.push_back(z_small);
		--this->order;
	      }
	    if (z_large.index() != 0)
	      {
		zero.push_back(z_large);
		--this->order;
	      }
	    return zero;
	  }
	// Find largest and smallest moduli of coefficients.
	auto max = _Real{0};
	auto min = this->infty;
	for (int i = 0; i <= this->order; ++i)
	  {
	    auto x = std::abs(this->P[i]);
	    if (x > max)
	      max = x;
	    if (x != _Real{0} && x < min)
	      min = x;
	  }
	// Scale if there are large or very small coefficients.
	// Computes a scale factor to multiply the coefficients
	// of the polynomial. The scaling is done to avoid overflow
	// and to avoid undetected underflow interfering
	// with the convergence criterion.
	// The factor is a power of the base.
	auto scale = lo / min;
	bool rescale = true;
	if (scale > _Real{1} && this->infty / scale < max)
	  rescale = false;
	if (scale <= _Real{1} && max < _Real{10})
	  rescale = false;

	if (rescale)
	  {
	    // Scale polynomial.
	    if (scale == _Real{0})
	      scale = this->small;
	    auto l = int(std::log(scale) / std::log(base) + _Real{0.5});
	    auto factor = std::pow(base, l);
	    if (factor != _Real{1})
	      for (int i = 0; i <= this->order; ++i)
		this->P[i] *= factor;
	  }

	// Compute lower bound on moduli of roots.
	for (int i = 0; i <= this->order; ++i)
	  pt[i] = std::abs(this->P[i]);
	pt[this->order] = -pt[this->order];
	// Compute upper estimate of bound.
	auto x = std::exp((std::log(-pt[this->order])
			 - std::log(pt[0])) / _Real(this->order));
	// If Newton step at the origin is better, use it.	
	if (pt[this->order - 1] != _Real{0})
	  {
	    auto xm = -pt[this->order] / pt[this->order - 1];
	    if (xm < x)
	      x = xm;
	  }
	// Chop the interval (0,x) until ff <= 0.
	while (true)
	  {
	    auto xm = x * _Real{0.1L};
	    auto ff = pt[0];
	    for (int i = 1; i <= this->order; ++i)
	      ff = ff * xm + pt[i];
	    if (ff <= _Real{0})
	      break;
	    x = xm;
	  }
	// Do Newton interation until x converges to two decimal places.
	auto dx = x;
	while (std::abs(dx / x) > this->min_log_deriv)
	  {
	    auto ff = pt[0];
	    auto df = ff;
	    for (int i = 1; i < this->order; ++i)
	      {
		ff = ff * x + pt[i];
		df = df * x + ff;
	      }
	    ff = ff * x + pt[this->order];
	    dx = ff / df;
	    x -= dx;
	    ++this->num_iters;
	  }
	auto bound = x;
	// Compute the derivative as the initial H polynomial
	// and do 5 steps with no shift.
	auto nm1 = this->order - 1;
	for (int i = 1; i < this->order; ++i)
	  this->H[i] = _Real(this->order - i) * this->P[i] / _Real(this->order);
	this->H[0] = this->P[0];
	auto aa = this->P[this->order];
	auto bb = this->P[this->order - 1];
	this->zerok = (this->H[this->order - 1] == _Real{0});
	for(int jj = 0; jj < 5; ++jj)
	  {
	    ++this->num_iters;
	    auto cc = this->H[this->order - 1];
	    if (!this->zerok)
	      {
		// Use a scaled form of recurrence if value of H at 0 is nonzero.	
		auto t = -aa / cc;
		for (int i = 0; i < nm1; ++i)
		  {
		    const auto j = this->order - i - 1;
		    this->H[j] = t * this->H[j - 1] + this->P[j];
		  }
		this->H[0] = this->P[0];
		this->zerok = (std::abs(this->H[this->order - 1])
			    <= _Real{10} * this->eta * std::abs(bb));
	    }
	    else
	      {
		// Use unscaled form of recurrence.
		for (int i = 0; i < nm1; ++i)
		  {
		    const auto j = this->order - i - 1;
		    this->H[j] = this->H[j - 1];
		  }
		this->H[0] = _Real{0};
		this->zerok = (this->H[this->order - 1] == _Real{0});
	      }
	}
	// Save H for restarts with new shifts.
	H_temp = this->H;

	// Loop to select the quadratic corresponding to each new shift.
	for (int count = 0; count < 20; ++count)
	  {
	    /*  Quadratic corresponds to a _Real shift to a	
	     *  non-real point and its complex conjugate. The point
	     *  has modulus bound and amplitude rotated by 94 degrees
	     *  from the previous shift.
	     */
	    auto xxx = cosr * xx - sinr * yy;
	    yy = sinr * xx + cosr * yy;
	    auto xx = xxx;
	    this->sr = bound * xx;
	    this->si = bound * yy;
	    this->u = -_Real{2} * this->sr;
	    this->v = bound;
	    auto num_zeros = this->fxshfr(20 * (count + 1));
	    bool cycle = false;
	    if (num_zeros != 0)
	      {
	      /*  The second stage jumps directly to one of the third
	       *  stage iterations and returns here if successful.
	       *  Deflate the polynomial, store the zero or zeros
	       *  and return to the main algorithm.
	       */
		zero.push_back(this->z_small);
		this->order -= num_zeros;
		this->P = this->P_quot;
		if (num_zeros != 1)
		  zero.push_back(this->z_large);
		cycle = true;
		break;
	      }
	    if (cycle)
	      continue;

	    // If the iteration is unsuccessful another quadratic
	    // is chosen after restoring H.
	    this->H = H_temp;
	 }
      }
  }


/**
 * Computes up to L2 fixed shift H-polynomials, testing for convergence
 * in the linear or quadratic case.
 * Initiates one of the variable shift iterations and returns
 * the number of zeros found.
 */
template<typename _Real>
  int
  _JenkinsTraubSolver<_Real>::fxshfr(int l2)
  {
    _Real ots, otv;
    int iflag;

    int num_zeros = 0;

    auto betav = _Real{0.25};
    auto betas = _Real{0.25};
    auto oss = this->sr;
    auto ovv = this->v;
    // Evaluate polynomial by synthetic division.
    this->remquo_quadratic(this->order, this->u, this->v,
			   this->P, this->P_quot,
			   this->a, this->b);
    auto type = this->calcsc();
    for (int j = 0; j < l2; ++j)
      {
	// Calculate next H polynomial and estimate v.
	this->next_k_poly(type);
	type = this->calcsc();
	auto [ui, vi] = this->quadratic_coefficients(type);
	auto vv = vi;
	// Estimate s.
	auto ss = _Real{0};
	if (this->H[this->order - 1] != _Real{0})
	  ss = -this->P[this->order] / this->H[this->order - 1];
	auto tv = _Real{1};
	auto ts = _Real{1};
	if (j == 0 || type == 3)
	  {
	    ovv = vv;
	    oss = ss;
	    otv = tv;
	    ots = ts;
	    continue;
	  }
	// Compute relative measures of convergence of s and v sequences.
	if (vv != _Real{0})
	  tv = std::abs((vv - ovv) / vv);
	if (ss != _Real{0})
	  ts = std::abs((ss - oss) / ss);
	// If decreasing, multiply two most recent convergence measures.
	auto tvv = _Real{1};
	if (tv < otv)
	  tvv = tv * otv;
	auto tss = _Real{1};
	if (ts < ots)
	  tss = ts * ots;
	// Compare with convergence criteria.
	auto vpass = tvv < betav;
	auto spass = tss < betas;
	if (!(spass || vpass))
	  {
	    ovv = vv;
	    oss = ss;
	    otv = tv;
	    ots = ts;
	    continue;
	  }
	// At least one sequence has passed the convergence test.
	// Store variables before iterating.
	auto u_save = this->u;
	auto v_save = this->v;
	this->H_save = this->H;
	auto s = ss;
	// Choose iteration according to the fastest converging sequence.
	auto vtry = false;
	auto stry = false;
	if ((spass && !vpass) || tss < tvv)
	  goto _40;
  _20:	
	num_zeros = this->iter_quadratic(ui, vi);
	if (num_zeros > 0)
	  return num_zeros;
	// Quadratic iteration has failed. Flag that it has
	// been tried and decrease the convergence criterion.
	vtry = true;
	betav *= _Real{0.25};
	// Try linear iteration if it has not been tried and
	// the S sequence is converging.
	if (stry || !spass)
	  goto _50;
	this->H = this->H_save;

  _40:
	num_zeros = this->iter_real(s, iflag);
	if (num_zeros > 0)
	  return num_zeros;
	// Linear iteration has failed. Flag that it has been
	// tried and decrease the convergence criterion.
	stry = true;
	betas *= _Real{0.25};
	if (iflag == 0)
	  goto _50;
	// If linear iteration signals an almost _Real real
	// zero attempt quadratic iteration.
	ui = -_Real{2} * s;
	vi = s * s;
	goto _20;
  _50:
	// Restore variables.
	this->u = u_save;
	this->v = v_save;
	this->H = this->H_save;

	// Try quadratic iteration if it has not been tried
	// and the V sequence is convergin.
	if (vpass && !vtry)
	  goto _20;
	// Recompute polynomial quotient and remainder
        // to continue the second stage.
	this->remquo_quadratic(this->order, this->u, this->v,
			       this->P, this->P_quot,
			       this->a, this->b);
	type = this->calcsc();
      }
    return num_zeros;
  }


/**
 * Variable-shift H-polynomial iteration for a quadratic factor
 * converges only if the zeros are equimodular or nearly so.
 * @param uu The linear coefficient of the starting quadratic equation.
 * @param vv The constant  coefficient of the starting quadratic equation.
 * @return The number of zeros found.
 */
template<typename _Real>
  int
  _JenkinsTraubSolver<_Real>::iter_quadratic(_Real uu, _Real vv)
  {
    _Real mp, omp, ee, relstp, t, zm;
    int type;

    int num_zeros = 0;
    int tried = 0;
    this->u = uu;
    this->v = vv;
    int j = 0;

    while (true)
      {
	++this->num_iters;
	this->quadratic(_Real{1}, this->u, this->v,
			this->z_small, this->z_large);
	// Return if roots of the quadratic are real and not
	// close to multiple or nearly equal and of opposite sign.
	if (std::abs(std::abs(real(this->z_small))
		   - std::abs(real(this->z_large)))
	       > _Real{0.01L} * std::abs(real(this->z_large)))
	  return num_zeros;
	// Evaluate polynomial by quadratic synthetic division.
	this->remquo_quadratic(this->order, this->u, this->v,
			       this->P, this->P_quot, this->a, this->b);
	mp = std::abs(this->a - real(this->z_small) * this->b)
	   + std::abs(imag(this->z_small) * this->b);
	// Compute a rigorous bound on the rounding error in evaluating P.
	zm = std::sqrt(std::abs(this->v));
	ee = _Real{2} * std::abs(this->P_quot[0]);
	t = -real(this->z_small) * this->b;
	for (int i = 1; i < this->order; ++i)
	  ee = ee * zm + std::abs(this->P_quot[i]);
	ee = ee * zm + std::abs(this->a + t);
	ee *= (_Real{5} * this->mre + _Real{4} * this->are);
	ee -= (_Real{5} * this->mre + _Real{2} * this->are)
	    * (std::abs(this->a + t) + std::abs(this->b) * zm);
	ee += _Real{2} * this->are * std::abs(t);
	// Iteration has converged sufficiently if the
	// polynomial value is less than 20 times this bound.
	if (mp <= _Real{20} * ee)
	  {
	    num_zeros = 2;
	    return num_zeros;
	  }
	j++;
	// Stop iteration after 20 steps.
	if (j > 20)
	  return num_zeros;
	if (j < 2 || relstp > _Real{0.01L} || mp < omp || tried)
	  {
	    omp = mp;
	    // Calculate next H polynomial and new u and v.
	    type = this->calcsc();
	    this->next_k_poly(type);
	    type = this->calcsc();
	    auto [ui, vi] = this->quadratic_coefficients(type);
	    // If vi is zero the iteration is not converging.
	    if (vi == _Real{0})
	      return num_zeros;
	    relstp = std::abs((vi - this->v) / vi);
	    this->u = ui;
	    this->v = vi;
	    continue;
	  }
	// A cluster appears to be stalling the convergence.
	// Five fixed shift steps are taken with a u, v close to the cluster.
	if (relstp < this->eta)
	  relstp = this->eta;
	relstp = std::sqrt(relstp);
	this->u -= this->u * relstp;
	this->v += this->v * relstp;
	this->remquo_quadratic(this->order, this->u, this->v,
			       this->P, this->P_quot,
			       this->a, this->b);
	for (int i = 0; i < 5; ++i)
	  {
	    type = this->calcsc();
	    this->next_k_poly(type);
	  }
	tried = 1;
	j = 0;
      }
  }


/**
 * Variable-shift H polynomial iteration for a real zero.
 * @param sss Starting iterate.
 * @param iflag Flag to indicate a pair of zeros near real axis.
 * @return The number of zeros found.
 */
template<typename _Real>
  int
  _JenkinsTraubSolver<_Real>::iter_real(_Real sss, int& iflag)
  {
    _Real t = _Real{0};
    _Real omp;

    int num_zeros = 0;
    _Real s = sss;
    iflag = 0;
    int i_real = 0;

    while (true)
      {
	++this->num_iters;
	auto pval = this->P[0];
	// Evaluate P at s.
	this->P_quot[0] = pval;
	for (int i = 1; i <= this->order; ++i)
	  {
	    pval = pval * s + this->P[i];
	    this->P_quot[i] = pval;
	  }
	auto mp = std::abs(pval);
	// Compute a rigorous bound on the error in evaluating P.
	auto ms = std::abs(s);
	auto ee = (this->mre / (this->are + this->mre)) * std::abs(this->P_quot[0]);
	for (int i = 1; i <= this->order; ++i)
	  ee = ee * ms + std::abs(this->P_quot[i]);
	// Iteration has converged sufficiently if the polynomial
	// value is less than 20 times this bound.
	if (mp <= _Real{20} * ((this->are + this->mre) * ee - this->mre * mp))
	  {
	    num_zeros = 1;
	    this->z_small = s;
	    return num_zeros;
	  }
	++i_real;
	// Stop iteration after max_iter_real steps.
	if (i_real > this->max_iter_real)
	  return num_zeros;
	else if (i_real < 2
		 || std::abs(t) > _Real{0.001L} * std::abs(s - t) || mp < omp)
	  {
	    // Return if the polynomial value has increased significantly.
	    omp = mp;

	    // Compute t, the next polynomial, and the new iterate.
	    auto kval = this->H[0];
	    this->H_quot[0] = kval;
	    for (int i = 1; i < this->order; ++i)
	      {
		kval = kval * s + this->H[i];
		this->H_quot[i] = kval;
	      }

	    if (std::abs(kval)
		 <= std::abs(this->H[this->order - 1]) * _Real{10} * this->eta)
	      { // HVE n -> n-1
		// Use unscaled form.
		this->H[0] = _Real{0};
		for (int i = 1; i < this->order; ++i)
		  this->H[i] = this->H_quot[i-1];
	      }
	    else
	      {
		// Use the scaled form of the recurrence if the value
		// of H at s is nonzero.
		t = -pval / kval;
		this->H[0] = this->P_quot[0];
		for (int i = 1; i < this->order; ++i)
		  this->H[i] = t * this->H_quot[i - 1] + this->P_quot[i];
	      }

	    kval = this->H[0];
	    for (int i = 1; i < this->order; ++i)
	      kval = kval * s + this->H[i];
	    auto t = _Real{0};
	    if (std::abs(kval)
		 > std::abs(this->H[this->order - 1] * _Real{10} * this->eta))
	      t = -pval / kval;
	    s += t;
	  }
	else
	  {
	    // A cluster of zeros near the real axis has been encountered.
	    // Return with iflag set to initiate a quadratic iteration.
	    iflag = 1;
	    sss = s; // HVE sss = s added
	    return num_zeros;
	  }
      }
  }

/**
 * This routine calculates scalar quantities used to compute
 * the next H-polynomial and new estimates of the quadratic coefficients.
 *
 * @return Flag indicating how the calculations are normalized
 *         to avoid overflow.
 */
template<typename _Real>
  int
  _JenkinsTraubSolver<_Real>::calcsc()
  {
    const auto eps = _Real{100} * this->eta;
    // Synthetic division of H by the quadratic 1, u, v
    int type = 0;
    this->remquo_quadratic(this->order - 1, this->u, this->v,
			   this->H, this->H_quot, this->c, this->d);
    if (std::abs(this->c) > eps * std::abs(this->H[this->order - 1])
     || std::abs(this->d) > eps * std::abs(this->H[this->order - 2]))
      {
	if (std::abs(this->d) < std::abs(this->c))
	  {
	    // Type = 1 indicates that all formulas are divided by c.
	    type = 1;
	    this->e = this->a / this->c;
	    this->f = this->d / this->c;
	    this->g = this->u * this->e;
	    this->h = this->v * this->b;
	    this->a3 = this->a * this->e
		     + (this->h / this->c + this->g) * this->b;
	    this->a1 = this->b - this->a * (this->d / this->c);
	    this->a7 = this->a + this->g * this->d + this->h * this->f;
	    return type;
	  }
	else
	  {
	    // Type = 2 indicates that all formulas are divided by d.
	    type = 2;
	    this->e = this->a / this->d;
	    this->f = this->c / this->d;
	    this->g = this->u * this->b;
	    this->h = this->v * this->b;
	    this->a3 = (this->a + this->g) * this->e
		     + this->h * (this->b / this->d);
	    this->a1 = this->b * this->f - this->a;
	    this->a7 = (this->f + this->u) * this->a + this->h;
	    return type;
	  }
      }
    else
      {
	// Type == 3 indicates the quadratic is almost a factor of H.
	type = 3;
	return type;
      }
  }


/**
 * Computes the next H polynomials using scalars computed in calcsc.
 */
template<typename _Real>
  void
  _JenkinsTraubSolver<_Real>::next_k_poly(int type)
  {
    if (type == 3)
      {
	// Use unscaled form of the recurrence if type is 3.
	this->H[0] = _Real{0};
	this->H[1] = _Real{0};
	for (int i = 2; i < this->order; ++i)
	  this->H[i] = this->H_quot[i-2];
	return;
      }
    auto ab_temp = this->a;
    if (type == 1)
      ab_temp = this->b;
    if (std::abs(this->a1) <= std::abs(ab_temp) * this->eta * _Real{10})
      {
	// If a1 is nearly zero then use a special form of the recurrence.
	this->H[0] = _Real{0};
	this->H[1] = -this->a7 * this->P_quot[0];
	for(int i = 2; i < this->order; ++i)
	  this->H[i] = this->a3 * this->H_quot[i-2] - this->a7 * this->P_quot[i-1];
	return; // HVE return added
      }
    else
      {
	// Use scaled form of the recurrence.
	this->a7 /= this->a1;
	this->a3 /= this->a1;
	this->H[0] = this->P_quot[0];
	this->H[1] = this->P_quot[1] - this->a7 * this->P_quot[0];
	for (int i = 2; i < this->order; ++i)
	  this->H[i] = this->a3 * this->H_quot[i-2]
		     - this->a7 * this->P_quot[i-1] + this->P_quot[i];
      }
  }


/**
 * Compute new estimates of the quadratic coefficients
 * using the scalars computed in calcsc.
 */
template<typename _Real>
  std::pair<_Real, _Real>
  _JenkinsTraubSolver<_Real>::quadratic_coefficients(int type)
  {
    if (type == 3)
      // If type=3 the quadratic is zeroed.
      return std::make_pair(_Real{0}, _Real{0});

    _Real a4, a5;
    if (type == 2)
      {
	a4 = (this->a + this->g) * this->f + this->h;
	a5 = (this->f + this->u) * this->c + this->v * this->d;
      }
    else
      {
	a4 = this->a + this->u * this->b + this->h * this->f;
	a5 = this->c + (this->u + this->v * this->f) * this->d;
      }

    // Evaluate new quadratic coefficients.
    const auto n = this->order;
    auto b1 = -this->H[n - 1] / this->P[n];
    auto b2 = -(this->H[n - 2] + b1 * this->P[n - 1])
	    / this->P[n];
    auto c1 = this->v * b2 * this->a1;
    auto c2 = b1 * this->a7;
    auto c3 = b1 * b1 * this->a3;
    auto c4 = c1 - c2 - c3;
    if (auto temp = a5 + b1 * a4 - c4; temp == _Real{0})
      return std::make_pair(_Real{0}, _Real{0});
    else
      {
	auto uu = this->u - (this->u * (c3 + c2)
		+ this->v * (b1 * this->a1 + b2 * this->a7)) / temp;
	auto vv = this->v * (_Real{1} + c4 / temp);
	return std::make_pair(uu, vv);
      }
  }

/**
 * Divides the polynomial P by the quadratic 1x^2 + ux + v
 * placing the quotient in q and the remainder in a, b.
 */
template<typename _Real>
  void
  _JenkinsTraubSolver<_Real>::remquo_quadratic(int nn, _Real u, _Real v,
					       std::vector<_Real>& poly,
					       std::vector<_Real>& quot,
					       _Real& a, _Real& b)
  {
    b = poly[0];
    quot[0] = b;
    a = poly[1] - b * u;
    quot[1] = a;
    for (int i = 2; i <= nn; ++i)
      {
	auto c = poly[i] - a * u - b * v;
	quot[i] = c;
	b = a;
	a = c;
      }	
  }


/**
 * Calculate the zeros of the quadratic az^2 + bz + c.
 * The quadratic formula, modified to avoid overflow, is used
 * to find the larger zero if the zeros are real and both
 * are complex. The smaller real zero is found directly from
 * the product of the zeros c/a.
 */
template<typename _Real>
  void
  _JenkinsTraubSolver<_Real>::quadratic(_Real a, _Real b, _Real c,
					solution_t<_Real>& z_small,
					solution_t<_Real>& z_large)
  {
    z_small = {};
    z_large = {};
    if (a == _Real{0})
      { // Less than two roots.
	if (b != _Real{0})
	  z_small = -c / b;
	return;
      }

    if (c == _Real{0})
      { // one real root, one zero root.
	z_small = _Real{0};
	z_large = -b / a;
	return;
      }

    // Compute discriminant avoiding overflow.
    auto b2 = b / _Real{2};

    _Real d, e;
    if (std::abs(b2) < std::abs(c))
      {
	e = std::copysign(a, c);
	e = b2 * (b2 / std::abs(c)) - e;
	d = std::sqrt(std::abs(e)) * std::sqrt(std::abs(c));
      }
    else
      {
	e = _Real{1} - (a / b2) * (c / b2);
	d = std::sqrt(std::abs(e)) * std::abs(b2);
      }

    if (e < _Real{0})
      { // complex conjugate zeros.
	z_small = std::complex<_Real>(-b2 / a, +std::abs(d / a));
	z_large = std::complex<_Real>(-b2 / a, -std::abs(d / a));
      }
    else
      {
	if (b2 >= _Real{0})
	  d = -d; // Real zeros.
	z_large = (-b2 + d) / a;
	z_small = _Real{0};
	if (z_large != _Real{0})
	  z_small = (c / z_large) / a;
      }
  }

} // namespace __gnu_cxx


template<typename _Real>
  void
  test_jenkins_traub()
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);

    int order = 0;
    int MAX_TERMS = 100;
    while ((order < 2) || (order > MAX_TERMS - 1))
      {
	std::cout << "Polynomial order (2 - " << MAX_TERMS - 1 << "): ";
	std::cin >> order;
      }
    std::vector<_Real> a(order + 1);

    std::cout << "Enter coefficients, high order to low order.\n";
    for (int i = 0; i <= order; ++i)
      {
	std::cout << "a[" << i << "] = ";
	std::cin >> a[i];
      }

    __gnu_cxx::_JenkinsTraubSolver jenkins_traub(a);
    const auto zeros = jenkins_traub.solve();
    std::cout << "\nThe zeros are:\n";
    for (const auto& z : zeros)
      std::cout << z << '\n';
/*
    const auto eq = jenkins_traub.equations();
    std::cout << "\nThe quadratic factors are:\n";
    for (int p = 0; p < eq.size() / 2; ++p)
      std::cout << "t^2 + " << eq[2 * p + 1] << " t + " << eq[2 * p] << '\n';
    if ((eq.size() % 2) == 1)
      std::cout << "The linear term is: \nt - " << eq.back() << '\n';
*/
    std::cout << "\nSolution tests:\n";
    __gnu_cxx::_Polynomial<_Real> poly(a.begin(), a.end());
    for (const auto& z : zeros)
      {
	const auto idx = z.index();
	std::cout << "P(" << z << ") = ";
	if (idx == 1)
	  std::cout << poly(std::get<1>(z));
	else if (idx == 2)
	  std::cout << poly(std::get<2>(z));
	std::cout << '\n';
      }
  }

int
main()
{
  std::cout << "\ndouble\n\n======\n" << std::flush;
  test_jenkins_traub<double>();

  //std::cout << "\nlong double\n\n===========\n" << std::flush;
  //test_jenkins_traub<long double>();
}

