#ifndef _CHEBYSHEV_H
#define _CHEBYSHEV_H 1

#include <initializer_list>
#include <vector>
#include <iosfwd>

#include "polynomial.h"

namespace __gnu_cxx
{

  /**
   *  @brief  _Chebyshev represents a Chebyshev fit of a function.
   *  Given a function @c func, the lower limit @c a, the upper limit @c b,
   *  and the number of points @c n _Chebyshev contains the coefficients
   *  of a Chebyshev polynomial expansion such that
   *    func(x) =~ kSum_0^{n-1} c_k T_k(x) - c_0/2.
   *  Use these constructors with moderately large n of 30 or 50.
   *  Then truncate the Chebyshev fit to a smaller number of terms to satisfy
   *  the accuracy requirements.
   */
  template<typename _Tp>
    class _Chebyshev
    {
    public:

      using value_type = _Tp;

      _Chebyshev()
      : _M_lower{-1}, _M_upper{+1},
	_M_coef{}
      { }

      _Chebyshev(_Tp __a, _Tp __b, unsigned __n, _Tp __func(_Tp));

      _Chebyshev(_Tp __a, _Tp __b, std::initializer_list<_Tp> __il)
      : _M_lower{__a}, _M_upper{__b},
	_M_coef{__il}
      { }

      _Chebyshev(_Tp __a, _Tp __b, const std::vector<_Tp>& __coeff)
      : _M_lower{__a}, _M_upper{__b},
	_M_coef{__coeff}
      { }

      _Chebyshev(_Tp __a, _Tp __b, std::vector<_Tp>&& __coeff)
      : _M_lower{__a}, _M_upper{__b},
	_M_coef{std::move(__coeff)}
      { }

      template<typename _Up>
        _Chebyshev(_Up __a, _Up __b, const _Polynomial<_Up>& __poly);

      _Chebyshev derivative() const;

      _Chebyshev integral() const;

      std::size_t
      order() const
      { return this->_M_coef.size(); }

      value_type
      lower() const
      { return this->_M_lower; }

      value_type
      upper() const
      { return this->_M_upper; }

      _Polynomial<_Tp>
      to_polynomial() const;

      value_type operator()(value_type __x) const;

      template<typename _Up>
        void truncate(_Up __eps = std::numeric_limits<_Up>::epsilon());

      template<typename _Tp1, typename _CharT, typename _Traits>
	friend std::basic_ostream<_CharT, _Traits>&
	operator<<(std::basic_ostream<_CharT, _Traits>& __os,
		   const _Chebyshev<_Tp1>& __cheb);

    private:
      _Tp _M_lower;
      _Tp _M_upper;
      std::vector<_Tp> _M_coef;
    };

} // namespace __gnu_cxx

#include "chebyshev.tcc"

#endif // _CHEBYSHEV_H
