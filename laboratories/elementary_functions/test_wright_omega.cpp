/**
 *
 */

#include <cmath>
#include <cfenv>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

template<typename Tp>
  std::complex<Tp>
  wright_omega(const std::complex<Tp>& z,
		 std::complex<Tp>& err, std::complex<Tp>& res,
		 std::complex<Tp>& condn)
  {
    using _Cmplx = std::complex<Tp>;
    const auto _S_NaN = emsr::quiet_NaN(z.real());
    const auto _S_eps = emsr::epsilon(z.real());
    const auto _S_pi = emsr::pi_v<Tp>;
    const auto _S_i = _Cmplx(Tp{0}, Tp{1});
    const auto _S_0 = _Cmplx{};
    auto [x, y] = reinterpret_cast<const Tp(&)[2]>(z);
    auto ympi = y - _S_pi;
    auto yppi = y + _S_pi;
    const auto _S_near = Tp{0.1};
    err = _S_0;
    res = _S_0;

    if (std::isnan(x) || std::isnan(y))
      return _Cmplx(_S_NaN, _S_NaN);
    else if (std::isinf(x) && (x < Tp{0})
	  && (-_S_pi < y) && (y <= _S_pi))
      { // Signed zeros between branches.
	_Cmplx w;
	if (std::abs(y) <= _S_pi / Tp{2})
	  w = +Tp{0};
	else
	  w = -Tp{0};

	if (y < Tp{0})
	  w += -_S_i * Tp{0};

	return w;
      }
    else if (std::isinf(x) || std::isinf(y))
      return _Cmplx(x, y);
    else if (x == Tp{-1} && std::abs(y) == _S_pi)
      return _Cmplx(Tp{-1}, Tp{0});
    else
      {
	//  Choose approximation based on region.

	_Cmplx w;
	if ((Tp{-2} < x && x <= Tp{1}
	  && Tp{1} < y && y < Tp{2} * _S_pi))
	  {
	    // Region 1: upper branch point.
	    // Series about z = -1 + i*pi.
	    const auto dz = z + Tp{1} - _S_i * _S_pi;
	    const auto pz = std::conj(std::sqrt(std::conj(Tp{2} * dz)));

	    w = Tp{-1}
		+ (_S_i
		+ (Tp{1} / Tp{3}
		+ (Tp{-1} / Tp{36} * _S_i
		+ (Tp{1} / Tp{270} + Tp{1} / Tp{4320} * _S_i * pz)
		* pz) * pz) * pz) * pz;
	  }
	else if ((Tp{-2} < x && x <= Tp{1}
	       && Tp{-2} * _S_pi < y && y < Tp{-1}))
	  {
	    // Region 2: lower branch point.
	    // Series about z = -1 - i*pi.
	    const auto dz = z + Tp{1} + _S_i * _S_pi;
	    const auto pz = std::conj(std::sqrt(std::conj(Tp{2} * dz)));

	    w = Tp{-1}
		+ (-_S_i
		+ (Tp{1} / Tp{3}
		+ (Tp{1} / Tp{36} * _S_i
		+ (Tp{1} / Tp{270} - Tp{1} / Tp{4320} * _S_i * pz)
		* pz) * pz) * pz) * pz;
	  }
	else if (x <= Tp{-2} && -_S_pi < y && y <= _S_pi)
	  {
	    // Region 3: between branch cuts.
	    // Series: About -infinity.
	    const auto pz = std::exp(z);
	    w = (Tp{1}
		+ (Tp{-1}
		+ (Tp{3} / Tp{2}
		+ (Tp{-8} / Tp{3}
		+ Tp{125} / Tp{24} * pz) * pz) * pz) * pz) * pz;
	  }
	else if (((Tp{-2} < x) && (x <= Tp{1})
	       && (Tp{-1} <= y) && (y <= Tp{1}))
		|| ((Tp{-2} < x)
		 && (x - Tp{1}) * (x - Tp{1}) + y * y
			 <= _S_pi * _S_pi))
	  {
	    // Region 4: Mushroom.
	    // Series about z = 1.
	    const auto pz = z - Tp{1};
	    w = Tp{1} / Tp{2} + Tp{1} / Tp{2} * z
		+ (Tp{1} / Tp{16}
		+ (Tp{-1} / Tp{192}
		+ (Tp{-1} / Tp{3072}
		+ Tp{13} / Tp{61440} * pz) * pz) * pz) * pz * pz;
	  }
	else if (x <= -Tp{3} / Tp{2}
		 && _S_pi < y
		 && y - _S_pi <= Tp{-3} / Tp{4} * (x + Tp{1}))
	  {
	    // Region 5: Top wing.
	    // Negative log series.
	    const auto t = z - _S_i * _S_pi;
	    const auto pz = std::log(-t);
	    w = ((Tp{1}
		  + (-Tp{3} / Tp{2}
		  + Tp{1} / Tp{3} * pz) * pz) * pz
		 + ((Tp{-1}
		  + Tp{1} / Tp{2} * pz) * pz
		  + (pz + (-pz + t) * t) * t) * t)
		/ (t * t * t);
	  }
	else if (x <= -Tp{3} / Tp{2}
		 && Tp{3} / Tp{4} * (x + Tp{1}) < y + _S_pi
				    && y + _S_pi <= Tp{0})
	  {
	    // Region 6: Bottom wing.
	    // Negative log series.
	    const auto t = z + _S_i * _S_pi;
	    const auto pz = std::log(-t);
	    w = ((Tp{1}
		 + (Tp{-3} / Tp{2}
		 + Tp{1} / Tp{3} * pz) * pz) * pz
		+ ((Tp{-1}
		 + Tp{1} / Tp{2} * pz) * pz
		 + (pz + (-pz + t) * t) * t) * t)
		/ (t * t * t);
	  }
	else
	  {
	    // Region 7: Everywhere else.
	    // Series solution about infinity.
	    const auto pz = std::log(z);
	    w = ((Tp{1}
		 + (Tp{-3} / Tp{2}
		 + Tp{1} / Tp{3} * pz) * pz) * pz
		+ ((Tp{-1}
		 + Tp{1} / Tp{2} * pz) * pz
		 + (pz + (-pz + z) * z) * z) * z)
		/ (z * z * z);
	  }

	auto sgn = Tp{0};
	auto zr = z;
	if (x <= Tp{-1} + _S_near
	    && (std::abs(ympi) <= _S_near || std::abs(yppi) <= _S_near))
	  {
	    sgn = Tp{-1};
	    // Regularize if near branch cuts.
	    if (std::abs(ympi) <= _S_near)
	      {
		// Recompute ympi with directed rounding.
		fesetround(FE_UPWARD);

		ympi = y - _S_pi;

		if (ympi <= Tp{0})
		  {
        	    fesetround(FE_DOWNWARD);
        	    ympi = y - _S_pi;
		  }

		zr = _Cmplx(x, ympi);

		// Return rounding to default.
		fesetround(FE_TONEAREST);
	      }
	    else
	      {
		// Recompute yppi with directed rounding.
		fesetround(FE_UPWARD);

		yppi = y + _S_pi;

		if (yppi <= Tp{0})
		  {
        	    fesetround(FE_DOWNWARD);
        	    yppi = y + _S_pi;
		  }

		zr = _Cmplx(x, yppi);

		//  Return rounding to default.
		fesetround(FE_TONEAREST);
	      }
	  }
	else
	  sgn = Tp{+1};

	w *= sgn;

	while (true)
	  {
	    const auto res = zr - sgn * w - std::log(w);
	    const auto wp1 = sgn * w + Tp{1};
	    const auto yy = Tp{2} * wp1 * (wp1 + Tp{2} / Tp{3} * res);
	    const auto err = res / wp1 * (yy - res)
				/ (yy - Tp{2} * res);
	    w *= Tp{1} + err;
	    const auto res4 = std::pow(std::abs(res), Tp{4});
	    const auto wpol = Tp{-1} + w * (Tp{-8} + Tp{2} * w);
	    const auto test = std::abs(res4 * wpol);
	    if (test < _S_eps * Tp{72} * std::pow(std::abs(wp1), Tp{6}))
	      break;
	  }

	// Undo regularization.
	w *= sgn;

	// Provide condition number estimate.
	condn = zr / (Tp{1} + w);

	return w;
      }
  }

/**
 * Return the Wright omega function for complex argument.
 */
template<typename Tp>
  std::complex<Tp>
  wright_omega(const std::complex<Tp>& z)
  {
    std::complex<Tp> err, res, condn;
    return wright_omega(z, err, res, condn);
  }

/**
 * Return the Wright omega function for real argument.
 */
template<typename Tp>
  Tp
  wright_omega(Tp x)
  {
    using _Cmplx = std::complex<Tp>;
    _Cmplx z(x), err, res, condn;
    return std::real(wright_omega(z, err, res, condn));
  }

/**
 * Return the 0-branch of the Lambert W function of real argument.
 * This is denoted Wp in the DLMF.
 */
template<typename Tp>
  Tp
  lambert_wp(Tp x)
  {
    using _Cmplx = std::complex<Tp>;
    const auto _S_1de = emsr::inv_e_v<Tp>;
    if (x < _S_1de)
      throw std::domain_error("lambert_wp: Argument out of range.");
    else if (x == _S_1de)
      return Tp{-1};
    else
      {
	_Cmplx z(std::log(_Cmplx(x))),
	       err, res, condn;
	return std::real(wright_omega(z, err, res, condn));
      }
  }

/**
 * Return the -1-branch of the Lambert W function of real argument.
 * This is denoted Wm in the DLMF.
 */
template<typename Tp>
  Tp
  lambert_wm(Tp x)
  {
    using _Cmplx = std::complex<Tp>;
    const auto _S_2pi = emsr::tau_v<Tp>;
    const auto _S_i = _Cmplx(Tp{0}, Tp{1});
    const auto _S_1de = emsr::inv_e_v<Tp>;
    if (x < _S_1de || x > Tp{0})
      throw std::domain_error("lambert_wm: Argument out of range.");
    else if (x == _S_1de)
      return Tp{-1};
    else if (x == Tp{0})
      return -std::numeric_limits<Tp>::infinity();
    else
      {
	_Cmplx z(std::log(_Cmplx(x)) - _S_i * _S_2pi),
	       err, res, condn;
	return std::real(wright_omega(z, err, res, condn));
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
  w = wright_omega(z, err, res, condn);

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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
    w = wright_omega( z,  e, r, cond );
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
template<typename Tp>
  void
  test_wright_omega(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << '\n';
    for (int i = -40; i <= 100; ++i)
      {
	auto x = i * Tp{0.1L};
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
