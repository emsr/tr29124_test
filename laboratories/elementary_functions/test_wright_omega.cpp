/**
 *
 */

#include <cmath>
#include <cfenv>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  std::complex<_Tp>
  __wright_omega(const std::complex<_Tp>& __z,
		 std::complex<_Tp>& __err, std::complex<_Tp>& __res,
		 std::complex<_Tp>& __condn)
  {
    using _Cmplx = std::complex<_Tp>;
    const auto _S_NaN = __gnu_cxx::__quiet_NaN(__z.real());
    const auto _S_eps = __gnu_cxx::__epsilon(__z.real());
    const auto _S_pi = emsr::pi_v<_Tp>;
    const auto _S_i = _Cmplx(_Tp{0}, _Tp{1});
    const auto _S_0 = _Cmplx{};
    auto [__x, __y] = reinterpret_cast<const _Tp(&)[2]>(__z);
    auto __ympi = __y - _S_pi;
    auto __yppi = __y + _S_pi;
    const auto _S_near = _Tp{0.1};
    __err = _S_0;
    __res = _S_0;

    if (std::isnan(__x) || std::isnan(__y))
      return _Cmplx(_S_NaN, _S_NaN);
    else if (std::isinf(__x) && (__x < _Tp{0})
	  && (-_S_pi < __y) && (__y <= _S_pi))
      { // Signed zeros between branches.
	_Cmplx __w;
	if (std::abs(__y) <= _S_pi / _Tp{2})
	  __w = +_Tp{0};
	else
	  __w = -_Tp{0};

	if (__y < _Tp{0})
	  __w += -_S_i * _Tp{0};

	return __w;
      }
    else if (std::isinf(__x) || std::isinf(__y))
      return _Cmplx(__x, __y);
    else if (__x == _Tp{-1} && std::abs(__y) == _S_pi)
      return _Cmplx(_Tp{-1}, _Tp{0});
    else
      {
	//  Choose approximation based on region.

	_Cmplx __w;
	if ((_Tp{-2} < __x && __x <= _Tp{1}
	  && _Tp{1} < __y && __y < _Tp{2} * _S_pi))
	  {
	    // Region 1: upper branch point.
	    // Series about z = -1 + i*pi.
	    const auto __dz = __z + _Tp{1} - _S_i * _S_pi;
	    const auto __pz = std::conj(std::sqrt(std::conj(_Tp{2} * __dz)));

	    __w = _Tp{-1}
		+ (_S_i
		+ (_Tp{1} / _Tp{3}
		+ (_Tp{-1} / _Tp{36} * _S_i
		+ (_Tp{1} / _Tp{270} + _Tp{1} / _Tp{4320} * _S_i * __pz)
		* __pz) * __pz) * __pz) * __pz;
	  }
	else if ((_Tp{-2} < __x && __x <= _Tp{1}
	       && _Tp{-2} * _S_pi < __y && __y < _Tp{-1}))
	  {
	    // Region 2: lower branch point.
	    // Series about z = -1 - i*pi.
	    const auto __dz = __z + _Tp{1} + _S_i * _S_pi;
	    const auto __pz = std::conj(std::sqrt(std::conj(_Tp{2} * __dz)));

	    __w = _Tp{-1}
		+ (-_S_i
		+ (_Tp{1} / _Tp{3}
		+ (_Tp{1} / _Tp{36} * _S_i
		+ (_Tp{1} / _Tp{270} - _Tp{1} / _Tp{4320} * _S_i * __pz)
		* __pz) * __pz) * __pz) * __pz;
	  }
	else if (__x <= _Tp{-2} && -_S_pi < __y && __y <= _S_pi)
	  {
	    // Region 3: between branch cuts.
	    // Series: About -infinity.
	    const auto __pz = std::exp(__z);
	    __w = (_Tp{1}
		+ (_Tp{-1}
		+ (_Tp{3} / _Tp{2}
		+ (_Tp{-8} / _Tp{3}
		+ _Tp{125} / _Tp{24} * __pz) * __pz) * __pz) * __pz) * __pz;
	  }
	else if (((_Tp{-2} < __x) && (__x <= _Tp{1})
	       && (_Tp{-1} <= __y) && (__y <= _Tp{1}))
		|| ((_Tp{-2} < __x)
		 && (__x - _Tp{1}) * (__x - _Tp{1}) + __y * __y
			 <= _S_pi * _S_pi))
	  {
	    // Region 4: Mushroom.
	    // Series about z = 1.
	    const auto __pz = __z - _Tp{1};
	    __w = _Tp{1} / _Tp{2} + _Tp{1} / _Tp{2} * __z
		+ (_Tp{1} / _Tp{16}
		+ (_Tp{-1} / _Tp{192}
		+ (_Tp{-1} / _Tp{3072}
		+ _Tp{13} / _Tp{61440} * __pz) * __pz) * __pz) * __pz * __pz;
	  }
	else if (__x <= -_Tp{3} / _Tp{2}
		 && _S_pi < __y
		 && __y - _S_pi <= _Tp{-3} / _Tp{4} * (__x + _Tp{1}))
	  {
	    // Region 5: Top wing.
	    // Negative log series.
	    const auto __t = __z - _S_i * _S_pi;
	    const auto __pz = std::log(-__t);
	    __w = ((_Tp{1}
		  + (-_Tp{3} / _Tp{2}
		  + _Tp{1} / _Tp{3} * __pz) * __pz) * __pz
		 + ((_Tp{-1}
		  + _Tp{1} / _Tp{2} * __pz) * __pz
		  + (__pz + (-__pz + __t) * __t) * __t) * __t)
		/ (__t * __t * __t);
	  }
	else if (__x <= -_Tp{3} / _Tp{2}
		 && _Tp{3} / _Tp{4} * (__x + _Tp{1}) < __y + _S_pi
				    && __y + _S_pi <= _Tp{0})
	  {
	    // Region 6: Bottom wing.
	    // Negative log series.
	    const auto __t = __z + _S_i * _S_pi;
	    const auto __pz = std::log(-__t);
	    __w = ((_Tp{1}
		 + (_Tp{-3} / _Tp{2}
		 + _Tp{1} / _Tp{3} * __pz) * __pz) * __pz
		+ ((_Tp{-1}
		 + _Tp{1} / _Tp{2} * __pz) * __pz
		 + (__pz + (-__pz + __t) * __t) * __t) * __t)
		/ (__t * __t * __t);
	  }
	else
	  {
	    // Region 7: Everywhere else.
	    // Series solution about infinity.
	    const auto __pz = std::log(__z);
	    __w = ((_Tp{1}
		 + (_Tp{-3} / _Tp{2}
		 + _Tp{1} / _Tp{3} * __pz) * __pz) * __pz
		+ ((_Tp{-1}
		 + _Tp{1} / _Tp{2} * __pz) * __pz
		 + (__pz + (-__pz + __z) * __z) * __z) * __z)
		/ (__z * __z * __z);
	  }

	auto __sgn = _Tp{0};
	auto __zr = __z;
	if (__x <= _Tp{-1} + _S_near
	    && (std::abs(__ympi) <= _S_near || std::abs(__yppi) <= _S_near))
	  {
	    __sgn = _Tp{-1};
	    // Regularize if near branch cuts.
	    if (std::abs(__ympi) <= _S_near)
	      {
		// Recompute ympi with directed rounding.
		fesetround(FE_UPWARD);

		__ympi = __y - _S_pi;

		if (__ympi <= _Tp{0})
		  {
        	    fesetround(FE_DOWNWARD);
        	    __ympi = __y - _S_pi;
		  }

		__zr = _Cmplx(__x, __ympi);

		// Return rounding to default.
		fesetround(FE_TONEAREST);
	      }
	    else
	      {
		// Recompute yppi with directed rounding.
		fesetround(FE_UPWARD);

		__yppi = __y + _S_pi;

		if (__yppi <= _Tp{0})
		  {
        	    fesetround(FE_DOWNWARD);
        	    __yppi = __y + _S_pi;
		  }

		__zr = _Cmplx(__x, __yppi);

		//  Return rounding to default.
		fesetround(FE_TONEAREST);
	      }
	  }
	else
	  __sgn = _Tp{+1};

	__w *= __sgn;

	while (true)
	  {
	    const auto __res = __zr - __sgn * __w - std::log(__w);
	    const auto __wp1 = __sgn * __w + _Tp{1};
	    const auto __yy = _Tp{2} * __wp1 * (__wp1 + _Tp{2} / _Tp{3} * __res);
	    const auto __err = __res / __wp1 * (__yy - __res)
				/ (__yy - _Tp{2} * __res);
	    __w *= _Tp{1} + __err;
	    const auto __res4 = std::pow(std::abs(__res), _Tp{4});
	    const auto __wpol = _Tp{-1} + __w * (_Tp{-8} + _Tp{2} * __w);
	    const auto __test = std::abs(__res4 * __wpol);
	    if (__test < _S_eps * _Tp{72} * std::pow(std::abs(__wp1), _Tp{6}))
	      break;
	  }

	// Undo regularization.
	__w *= __sgn;

	// Provide condition number estimate.
	__condn = __zr / (_Tp{1} + __w);

	return __w;
      }
  }

/**
 * Return the Wright omega function for complex argument.
 */
template<typename _Tp>
  std::complex<_Tp>
  wright_omega(const std::complex<_Tp>& __z)
  {
    std::complex<_Tp> __err, __res, __condn;
    return __wright_omega(__z, __err, __res, __condn);
  }

/**
 * Return the Wright omega function for real argument.
 */
template<typename _Tp>
  _Tp
  wright_omega(_Tp __x)
  {
    using _Cmplx = std::complex<_Tp>;
    _Cmplx __z(__x), __err, __res, __condn;
    return std::real(__wright_omega(__z, __err, __res, __condn));
  }

/**
 * Return the 0-branch of the Lambert W function of real argument.
 * This is denoted Wp in the DLMF.
 */
template<typename _Tp>
  _Tp
  lambert_wp(_Tp __x)
  {
    using _Cmplx = std::complex<_Tp>;
    const auto _S_1de = emsr::inv_e_v<_Tp>;
    if (__x < _S_1de)
      std::__throw_domain_error("lambert_wp: Argument out of range.");
    else if (__x == _S_1de)
      return _Tp{-1};
    else
      {
	_Cmplx __z(std::log(_Cmplx(__x))),
	       __err, __res, __condn;
	return std::real(__wright_omega(__z, __err, __res, __condn));
      }
  }

/**
 * Return the -1-branch of the Lambert W function of real argument.
 * This is denoted Wm in the DLMF.
 */
template<typename _Tp>
  _Tp
  lambert_wm(_Tp __x)
  {
    using _Cmplx = std::complex<_Tp>;
    const auto _S_2pi = emsr::tau_v<_Tp>;
    const auto _S_i = _Cmplx(_Tp{0}, _Tp{1});
    const auto _S_1de = emsr::inv_e_v<_Tp>;
    if (__x < _S_1de || __x > _Tp{0})
      std::__throw_domain_error("lambert_wm: Argument out of range.");
    else if (__x == _S_1de)
      return _Tp{-1};
    else if (__x == _Tp{0})
      return -std::numeric_limits<_Tp>::infinity();
    else
      {
	_Cmplx __z(std::log(_Cmplx(__x)) - _S_i * _S_2pi),
	       __err, __res, __condn;
	return std::real(__wright_omega(__z, __err, __res, __condn));
      }
  }

/*
 * Burkhardt driver.
 */
void
driver(std::complex <double> z)
{
  std::complex <double> condn;
  std::complex <double> err;
  std::complex <double> res;
  std::complex <double> res_ult;
  std::complex <double> w;

  std::cout << "\n";
  std::cout << "DRIVER:\n";
  std::cout << "  Demonstrate simple and extended Wright Omega evaluators.\n";

  // Simple evaluator.
  w = wright_omega(z);

  std::cout << "\n";
  std::cout << "  Calling:\n";
  std::cout << "    w = wrightomega(z);\n";
  std::cout << "  returns:\n";
  std::cout << "    w = omega(" << std::real ( z ) 
       << ", " << std::imag ( z ) 
       << ") =  ( " << std::real ( w )
       << ", " << std::imag ( w ) << ")\n";

  // Extended evaluator.
  w = __wright_omega(z, err, res, condn);

  std::cout << "\n";
  std::cout << "  Calling:\n";
  std::cout << "    wrightomega_ext ( z, w, e, r, condest );\n";
  std::cout << "  returns:\n";
  std::cout << "    w = omega(" << std::real ( z ) 
       << ", " << std::imag ( z ) 
       << ") =  ( " << std::real ( w ) 
       << ", " << std::imag ( w ) << ")\n";
  std::cout << "  e = last update step = ( " << std::real ( err ) 
       << ", " << std::imag ( err ) << ")\n";
  std::cout << "  r = penultimate residual = ( " << std::real ( res ) 
       << ", " << std::imag ( res ) << ")\n";
  std::cout << "  condest = condition number estimate = ( " << std::real ( condn ) 
       << ", " << std::imag ( condn ) << ")\n";

  // Calculate and print ultimate residual.
  res_ult = ( 2.0 * w * w - 8.0 * w - 1.0 ) 
    / std::pow ( 1.0 + w, 6.0 ) * std::pow ( res, 4.0 );
  std::cout << "\n";
  std::cout << "  ultimate residual = ( " << std::real ( res_ult )
       << ", " << std::imag ( res_ult ) << ")\n";

  return;
}

/*
 * Make sure we didn't blow anything.
 */
void
burkhadt_test()
{
  double a;
  double b;
  double pi = M_PI;
  std::complex <double> z;

  //timestamp ( );
  std::cout << "\n";
  std::cout << "TOMS917_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the TOMS917 library.\n";

  a = 0.0;
  b = 0.0;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 1.0;
  b = 0.0;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 1.0 + exp ( 1.0 );
  b = 0.0;
  z = std::complex<double>(a, b);
  driver ( z );

  a = - 1.0;
  b = pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = - 1.0;
  b = - pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = - 2.0 + log ( 2.0 );
  b = pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = - 2.0 + log ( 2.0 );
  b = - pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 0.0;
  b = 1.0 + pi / 2.0;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 0.0;
  b = pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 1.0;
  b = 1.0;
  z = std::complex<double>(a, b);
  driver ( z );
};

void
test_burkhardt_boundary()
{
  std::complex<double> cond;
  std::complex<double> e;
  double exp_num = 160.0;
  std::string filename = "results.txt";
  std::ofstream fp;
  int i;
  int n = 100;
  double pi = M_PI;
  std::complex<double> r;
  double td;
  std::complex<double> w;
  double x[2];
  double y[2];
  std::complex<double> z;

  std::cout << "\n";
  std::cout << "TEST_BOUNDARY:\n";
  std::cout << "  Test wrightomega_ext() near approximation region boundaries.\n";
  std::cout << "  Store results in a file for comparison with benchmark data.\n";

  fp.open ( filename.c_str ( ) );

  //  We want trailing zeros, to make comparison with benchmark easier.
  fp << std::fixed;

  //  Region 1;
  //  x=(-2.0,1.0] ,y=2*pi
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  y[0] = std::nextafter ( 2.0 * pi, -1.0 );
  td = ( x[1] - x[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 2;
  //  x=1.0 ,y=(1.0,2*pi)
  x[0] = 1.0;
  y[0] = std::nextafter ( 1.0, 2.0 );
  y[1] = std::nextafter ( 2.0 * pi, 1.0 );
  td = - ( y[1] - y[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 3;
  //  x=(-2.0,1.0] ,y=1.0
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  td = - ( x[1] - x[0] ) / double(n);
  y[0] = std::nextafter ( 1.0, 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[1] + td *double(i), y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 4;
  //  x=-2.0 ,y=(1.0,2*pi)
  y[0] = std::nextafter ( 1.0, 2.0 );
  y[1] = std::nextafter ( 2.0 * pi, - 1.0 );
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 5;
  //  x=(-2.0,1.0] ,y=-2*pi
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  y[0] = std::nextafter ( - 2.0 * pi, 1.0 );
  td = ( x[1] - x[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 6;
  //  x=1.0, y=(-2*pi,-1.0)
  y[0] = std::nextafter ( - 2.0 * pi, 1.0 );
  y[1] = std::nextafter ( - 1.0, - 2.0 );
  td = ( y[1] - y[0] ) / double(n);
  x[0] = 1.0;

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 7;
  //  x=(-2.0,1.0] ,y=-1.0
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  td = ( x[0] - x[1] ) / double(n);
  y[0] = std::nextafter ( - 1.0, - 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[1] + td *double(i), y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
    }

  //  Region 8;
  //  x=-2.0 ,y=(-2*pi,-1.0)
  y[0] = std::nextafter ( - 2.0 * pi, 1.0 );
  y[1] = std::nextafter ( - 1.0, -2.0 );
  td = ( y[0] - y[1] ) / double(n);
  x[0] = std::nextafter ( - 2.0, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 9;
  //  x=-2.0 y=[-1.0,1.0]
  y[0] = - 1.0;
  y[1] = 1.0;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 10;
  //  x=(-2.0,1.0] y=1.0
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  td = ( x[1] - x[0] ) / double(n);
  y[0] = 1.0;

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 11
  //  x=1.0 y=[1.0,pi]
  y[0] = 1.0;
  y[1] = pi;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( 1.0, 2.0 );
  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 12
  //  (x-0.1e1)*(x-0.1e1)+y*y=pi*pi)
  //  (on inside)
  td = pi / double(n);
  x[0] = pi / 2.0;

  for ( i = 0; i < n; i++ )
  {
    double a = std::nextafter ( pi, -1.0 ) * cos ( x[0] - td *double(i) ) 
             + std::nextafter ( 1.0, - 1.0 );
    double b = std::nextafter ( pi, -1.0 ) * sin ( x[0] - td *double(i) );
    z = std::complex<double> ( a, b );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 13
  //  x=1.0 y=[-pi,-1.0]
  y[0] = - pi;
  y[1] = - 1.0;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( 1.0, 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 14
  //  x=(-2.0,1.0] y=-1.0
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  td = ( x[1] - x[0] ) / double(n);
  y[0] = -1.0;

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 15
  //  x=(-inf,-2) y=pi^+
  for ( i = 0; i < n; i++ )
  {
    x[0] = std::nextafter ( - 1.0 - exp (double( n - 1 - i ) / exp_num ), HUGE_VAL );
    y[0] = std::nextafter ( pi - 0.75 * ( x[0] + 1.0 ), HUGE_VAL );
    z = std::complex<double> ( x[0], y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 16
  y[0] = 0.75 + pi;
  y[1] = 2.0 * pi;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, - 3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 17
  //  x=(-2.0,1.0] ,y=2*pi
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  y[0] = 2.0 * pi;
  td = ( x[1] - x[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 18
  //  x=1.0 ,y=(pi,2*pi)
  y[0] = std::nextafter ( pi, 6.0 );
  y[1] = std::nextafter ( 2.0 * pi, 1.0 );
  td = - ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( 1.0, 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 19
  //  (x-0.1e1)*(x-0.1e1)+y*y=pi*pi)
  //  (on outside)
  td = pi / double( n - 1 );
  y[0] = pi / 2.0;

  for ( i = 0; i < n; i++ )
  {
    y[1] = pi * sin ( y[0] - td * i );
    x[0] = sqrt ( pi * pi - y[1] * y[1] ) + 1.0;
    if ( y[1] < 0 )
    {
      z = std::complex<double> ( std::nextafter ( x[0], HUGE_VAL ), std::nextafter ( y[1], - HUGE_VAL ) );
    }
    else
    {
      z = std::complex<double> ( std::nextafter ( x[0], HUGE_VAL ), std::nextafter ( y[1], HUGE_VAL ) );
    } 
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 20;
  //  x=1.0 ,y=(-2*pi,-pi)
  y[0] = std::nextafter ( - 2.0 * pi, 1.0 );
  y[1] = std::nextafter ( - pi, - 6.0 );
  td = - ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( 1.0, 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 21;
  //  x=(-2.0,1.0] ,y=-2*pi
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  y[0] = std::nextafter ( - 2.0 * pi, - 7.0 );
  td = - ( x[1] - x[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[1] + td *double(i), y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 22
  y[0] = - 0.75 - pi;
  y[1] = - 2.0 * pi;
  td = - ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, - 3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 23
  //  x=(-inf,-2) y=pi^+
  for ( i = 0; i < n; i++ )
  {
    x[0] = std::nextafter ( - 1.0 - exp (double( i ) / exp_num ), HUGE_VAL );
    y[0] = std::nextafter ( - pi + 0.75 * ( x[0] + 1.0 ), - HUGE_VAL );
    z = std::complex<double> ( x[0], y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 24
  //  x=(-inf,-2) y=pi^+
  for ( i = 0; i < n; i++ )
  {
    x[0] = std::nextafter ( - 1.0 - exp (double( n - 1 - i ) / exp_num ), - HUGE_VAL );
    y[0] = std::nextafter ( - pi + 0.75 * ( x[0] + 1.0 ), HUGE_VAL );
    z = std::complex<double> ( x[0], y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 25
  y[0] = - pi;
  y[1] = - 0.75 - pi;
  td = - ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, - 3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td * double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 26
  //  x=(-inf,-2) y=pi^+
  y[0] = std::nextafter ( - pi, -7.0 );

  for ( i = 0; i < n; i++ )
  {
    x[0] = - 1.0 - exp (double( i ) / exp_num );
    z = std::complex<double> ( x[0], y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 27
  //  x=(-inf,-2) y=pi^+
  y[0] = std::nextafter ( - pi, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    x[0] = - 1.0 - exp (double( n - 1 - i ) / exp_num );
    z = std::complex<double> ( x[0], y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 28
  y[0] = std::nextafter ( - pi, 1.0 );
  y[1] = pi;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, -3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 29
  //  x=(-inf,-2) y=pi^+
  y[0] = std::nextafter ( pi, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    x[0] = - 1.0 - exp (double( i ) / exp_num );
    z = std::complex<double> ( x[0], y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 30
  //  x=(-inf,-2) y=pi^+
  y[0] = std::nextafter ( pi, 7.0 );

  for ( i = 0; i < n; i++ )
  {
    x[0] = - 1.0 - exp (double( n - 1 - i ) / exp_num );
    z = std::complex<double> ( x[0], y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 31
  y[0] = std::nextafter ( pi, 7.0 );
  y[1] = 0.75 + pi;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, - 3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 32
  //  x=(-inf,-2) y=pi^+
  for ( i = 0; i < n; i++ )
  {
    x[0] = -1.0 - exp (double( n - 1 - i ) / exp_num );
    y[0] = std::nextafter ( pi - 0.75 * ( x[0] + 1.0 ), 0.1 );
    z = std::complex<double> ( x[0], y[0] );
    w = __wright_omega( z,  e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Terminate.
  //
  fp.close ( );

  std::cout << "\n";
  std::cout << "TEST_BOUNDARY:\n";
  std::cout << "  Results saved in file '" << filename << "'\n";

  return;
}

/**
 *
 */
template<typename _Tp>
  void
  test_wright_omega(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << '\n';
    for (int i = -40; i <= 100; ++i)
      {
	auto x = i * _Tp{0.1L};
	auto y = wright_omega(x);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << y
		  << '\n';
      }
  }

int
main()
{
  burkhadt_test();
  test_burkhardt_boundary();

  test_wright_omega<double>();
}
