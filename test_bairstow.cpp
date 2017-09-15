/*
$HOME/bin/bin/g++ -std=c++17 -g -I. -o test_bairstow test_bairstow.cpp
./test_bairstow < test_bairstow.in > test_bairstow.txt
*/

#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>

#include "bits/specfun_state.h"
#include "solver.h"

namespace __gnu_cxx
{

///
/// @todo: If you *don't* reverse the input array you solve for 1/z.
///
template<typename _Real>
  class _BairstowSolver
  {
  public:

    _BairstowSolver(const std::vector<_Real>& __coeff)
    : _M_coeff{__coeff.rbegin(), __coeff.rend()},
      _M_b(__coeff.size()), _M_c(__coeff.size()),
      _M_order(__coeff.size() - 1)
    {
      if (this->_M_coeff.size() == 0)
	std::__throw_domain_error("_BairstowSolver: "
				  "Coefficient size must be nonzero.");

      if (this->_M_coeff[0] == _Real{0})
	std::__throw_domain_error("_BairstowSolver: "
				  "Leading-order coefficient must be nonzero.");

      const auto __scale = this->_M_coeff[0];
      for (int i = 0; i <= this->_M_order; ++i)
	this->_M_coeff[i] /= __scale;

      this->_M_zero.reserve(this->_M_coeff.size());
    }

    std::vector<solution_t<_Real>> solve();

  private:

    void _M_iterate();

    void
    _M_add_zero(solution_t<_Real> __z)
    {
      this->_M_zero.push_back(__z);
      --this->_M_order;
    }

    static constexpr int _S_max_rand_iter = 200;
    static constexpr int _S_max_error_iter = 500;
    static constexpr auto _S_eps
      = _Real{10} * std::numeric_limits<_Real>::epsilon();

    std::vector<_Real> _M_coeff;
    std::vector<_Real> _M_b;
    std::vector<_Real> _M_c;
    std::vector<solution_t<_Real>> _M_zero;
    _Real _M_eps = _S_eps;
    int _M_order;
    bool _M_precision_error = false;
  };

template<typename _Real>
  void
  _BairstowSolver<_Real>::_M_iterate()
  {
    auto __r = _Real{0};
    auto __s = _Real{0};
    auto __dr = _Real{1};
    auto __ds = _Real{0};

    int __iter = 1;
    this->_M_precision_error = false;

    while ((std::abs(__dr) + std::abs(__ds)) > this->_M_eps)
      {
	if (__iter % _S_max_rand_iter == 0)
	  {
	    std::random_device rd;
	    std::mt19937 urng(rd());
	    std::uniform_real_distribution<_Real> pdf(0, 2);
	    __r = pdf(urng);
	  }
	if (__iter % _S_max_error_iter == 0)
	  {
	    this->_M_eps *= _Real{10};
	    this->_M_precision_error = true;
	    std::cout << "Loss of precision\n";
	  }
	this->_M_b[1] = this->_M_coeff[1] - __r;
	this->_M_c[1] = this->_M_b[1] - __r;

	for (int i = 2; i <= this->_M_order; ++i)
	  {
	    this->_M_b[i] = this->_M_coeff[i]
			  - __r * this->_M_b[i-1] - __s * this->_M_b[i-2];
	    this->_M_c[i] = this->_M_b[i]
			  - __r * this->_M_c[i-1] - __s * this->_M_c[i-2];
	  }

	auto __dn = this->_M_c[this->_M_order-1] * this->_M_c[this->_M_order-3]
		  - this->_M_c[this->_M_order-2] * this->_M_c[this->_M_order-2];
	if (std::abs(__dn) < std::numeric_limits<_Real>::epsilon())
	  {
	    __dr = _Real{1};
	    __ds = _Real{1};
	  }
	else
	  {
	    auto __drn = this->_M_b[this->_M_order]
		       * this->_M_c[this->_M_order-3]
		       - this->_M_b[this->_M_order-1]
		       * this->_M_c[this->_M_order-2];
	    auto __dsn = this->_M_b[this->_M_order-1]
		       * this->_M_c[this->_M_order-1]
		       - this->_M_b[this->_M_order]
		       * this->_M_c[this->_M_order-2];
	    __dr = __drn / __dn;
	    __ds = __dsn / __dn;
	  }

	__r += __dr;
	__s += __ds;
	++__iter;
      }
    for (int i = 0; i < this->_M_order - 1; ++i) 
      this->_M_coeff[i] = this->_M_b[i];

    std::array<solution_t<_Real>, 2> __z2 = __quadratic(__s, __r, _Real{1});

    this->_M_add_zero(__z2[0]);
    this->_M_add_zero(__z2[1]);
  }

template<typename _Real>
  std::vector<solution_t<_Real>>
  _BairstowSolver<_Real>::solve()
  {
    this->_M_eps = _S_eps;

    this->_M_b[0] = _Real{1};
    this->_M_c[0] = _Real{1};

    while (this->_M_order > 2)
      this->_M_iterate();
    this->_M_add_zero(this->_M_coeff[1]);

    return this->_M_zero;
  }

} // namespace __gnu_cxx

template<typename _Real>
  void
  test_bairstow()
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

    __gnu_cxx::_BairstowSolver bairstow(a);
    for (const auto& z : bairstow.solve())
      std::cout << z << '\n';
    //std::cout << "\nThe quadratic factors are:\n";

    //for (int i = order; i >= 2; i -= 2)
    //  std::cout << "t^2 + " << a[i - 1] << " t + " << a[i] << '\n';
    //if ((order % 2) == 1)
    //  std::cout << "The linear term is: \nt - " << a[1] << '\n';
  }

int
main()
{
  test_bairstow<double>();
}
