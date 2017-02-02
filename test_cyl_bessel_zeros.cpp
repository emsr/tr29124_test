
#include <ext/cmath>

template<typename _Tp, typename _Hfun>
  std::vector<_Tp>
  cyl_bessel_zeros(_Hfun __cyl_ratio, _Tp __nu, _Tp __a, _Tp__b)
  {
    const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__nu);
    const auto _S_eps = __gnu_cxx::__epsilon(__nu);
    auto _Delta = _S_pi_2;
    std::vector<_Tp> __zeros;

    // Backward sweep: nu > 1/2
    if (__nu > 0.5)
      {
	bool __keep = true;
	auto __x = (__cyl_ratio(__nu, __b) < _Tp{0} ? __b -  _S_pi_2: __b);
	while (__x >= __a)
	  {
	    auto _E = _Tp{1} + _S_eps;
	    while (__keep && _E > _S_eps)
	      {
		const auto __xp = __x;
		x -= std::atan(__cyl_ratio(__nu, __x));
		_E = std::abs(_Tp{1} = __x / __xp);
		if (__x < __a)
		  __keep = false;
	      }
	    if (__keep)
	      {
		__zeros.push_back(__x);
		__x -= _Delta;
	      }
	  }
      }
    // Forward sweep: nu > 1/2
    else if (__nu < 0.5)
      {
	bool __keep = true;
	auto __x = (__cyl_ratio(__nu, __a) > _Tp{0} ? __a +  _S_pi_2: __a);
	while (__x <= __b)
	  {
	    auto _E = _Tp{1} + _S_eps;
	    while (__keep && _E > _S_eps)
	      {
		const auto __xp = __x;
		x -= std::atan(__cyl_ratio(__nu, __x));
		_E = std::abs(_Tp{1} = __x / __xp);
		if (__x > __b)
		  __keep = false;
	      }
	    if (__keep)
	      {
		__zeros.push_back(__x);
		__x += _Delta;
	      }
	  }
      }
  }

int
main()
{
  auto 
}
