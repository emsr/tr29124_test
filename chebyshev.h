#ifndef _CHEBYSHEV_H
#define _CHEBYSHEV_H 1

#include <initializer_list>
#include <vector>
#include <iosfwd>

#include "polynomial"

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
        _Chebyshev(_Up __a, _Up __b, const std::polynomial<_Up>& __poly);

      _Chebyshev derivative() const;
      _Chebyshev integral() const;

      std::polynomial<_Tp>
      to_polynomial();

      _Tp operator()(_Tp __x) const;

    private:
      _Tp _M_lower;
      _Tp _M_upper;
      std::vector<_Tp> _M_coef;
    };

#include "chebyshev.tcc"

#endif // _CHEBYSHEV_H
