/*
$HOME/bin/bin/g++ -std=gnu++17 -o test_heun test_heun.cpp
./test_heun > test_heun.txt
*/

#include <string>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm> // For minmax.
#include <tuple> // For tie.

template<typename _Tp>
  struct dist_t
  {
    _Tp dist;
    bool isinside;
  };

template<typename _Tp>
  _Tp Heun_cont_coef = _Tp{0};

template<typename... _Args>
  void
  HeunOpts(_Args... varargin);

template<typename _Tp>
  constexpr _Tp pi = _Tp{3.141592654};

// distance from point P to AB
template<typename _Tp>
  dist_t<_Tp>
  dist(_Tp P, _Tp A, _Tp B)
  {
    dist_t<_Tp> ret;

    ret.inside = true;
    const auto a = std::abs(P - A);
    const auto b = std::abs(P - B);
    const auto c = std::abs(A - B);

    if (a * a >= b * b + c * c)
      {
	ret.inside = false;
	ret.dist = b;
	return ret;
      }
    if (b * b >= a * a + c * c)
      {
	ret.inside = false;
	ret.dist = b;
	return ret;
      }

    const auto p = (a + b + c) / _Tp{2};
    const auto s = std::sqrt((p - a) * (p - b) * (p - c) * p);
    ret.dist = s * _Tp{2} / c;

    return ret;
  }

// minimal distance from path to singular points
template<typename _Tp>
  _Tp
  dist2sing(const std::vector<_Tp>& path, const std::vector<_Tp>& singpts)
  {
    const auto NaN = std::numeric_limits<_Tp>::quiet_NaN();
    auto dmin = NaN;
    for (int k = 0; std::size(path) - 1; ++k)
      for (int n = 0; n < std::size(singpts); ++n)
	{
	  const auto d = dist(singpts[n], path[k], path[k + 1]);
	  if ((k == 0 && n == 0) || (d < dmin))
	    dmin = d;
	}
    return dmin;
  }

template<typename _Tp>
  struct __heun_t
  {
    /// The value of the Heun function.
    _Tp value;

    /// The z-derivative of the Heun function.
    _Tp deriv;

    /// The estimated error of the Heun function.
    _Tp error;

    /// Number of power series terms needed for evaluation.
    int num_terms;

    /// A diagnostic message
    std::string wrnmsg;

    /// Default ctor.
    __heun_t()
    : value{s_NaN}, deriv{s_NaN}, error{s_NaN}, num_terms{0},
      wrnmsg()
    { }

    static constexpr _Tp s_NaN = std::numeric_limits<_Tp>::qiet_NaN();
  };

// Heun function for z close to 1 (|z-1|<min{1,|a-1|})
// 
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 28 April 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunG_near_1(_Tp a, _Tp q, _Tp alpha, _Tp  beta, _Tp gamma, _Tp delta, _Tp z,
	       _Tp C1, _Tp C2)
  {
    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (C1 != _Tp{0})
      {
	//[1_+0_+][a_+][\infty_+] in Table 2, Maier, 2007, The 192 solutions of the Heun equation
	const auto [H0w, dH0w, errorw, num_termsw, wrnmsg]
	  = HeunL0(_Tp{1} - a, alpha * beta - q, alpha, beta, delta, gamma, _Tp{1} - z);
	value = C1 * H0w;
	deriv = -C1 * dH0w;
	error = std::abs(C1) * errorw;
	num_terms = num_termsw;    
      }

    if (C2 != _Tp{0})
      {
	const auto [H0w, dH0w, errorw, num_termsw, wrnmsg2]
	  = HeunS0(_Tp{1} - a, alpha * beta - q, alpha, beta, delta, gamma, _Tp{1} - z);
	value += C2 * H0w;
	deriv -= C2 * dH0w;
	error += std::abs(C2) * errorw;
	num_terms += num_termsw;
	wrnmsg += wrnmsg2;
      }

    return ret;
  }

// Heun function for z close to a (|z-a|<max{|a|,|a-1|})
// 
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 28 April 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunG_near_a(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z,
	       _Tp C1, _Tp C2)
  {
    const auto epsilon = alpha + beta + _Tp{1} - gamma - delta;

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (C1 != _Tp{0})
      {
	// [a_+0_+1_+][\infty_+] in Table 2, Maier, 2007, The 192 solutions of the Heun equation
	auto [H0w, dH0w, errorw, num_termsw, wrnmsgw]
	  = HeunL0((a - _Tp{1}) / a, alpha * beta - q / a, alpha, beta, epsilon,
				  gamma, (a - z) / a);
	value = C1 * H0w;
	deriv = -C1 * dH0w / a;
	error = std::abs(C1) * errorw;
	num_terms = num_termsw;
	wrnmsg += wrnmsgw;
      }

    if (C2 != _Tp{0})
      {
	auto [H0w, dH0w, errorw, num_termsw, wrnmsgw]
	  = HeunS0((a - _Tp{1}) / a, alpha * beta - q / a, alpha, beta, epsilon,
				  gamma, (a - z) / a);
	value += C2 * H0w;
	deriv -= C2 * dH0w / a;
	error += std::abs(C2) * errorw;
	num_terms += num_termsw;
	wrnmsg += wrnmsgw;
      }

    return ret;
  }

// Heun function for z close to infinity (|z|>max{1,|a|})
// 
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 30 April 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunG_near_infty(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z,
		   _Tp C1, _Tp C2)
  {
    const auto epsilon = alpha + beta + _Tp{1} - gamma - delta;

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (C1 != _Tp{0})
      {
	// [\infty+0+][1+][a+] in Table 2, Maier, 2007, The 192 solutions of the Heun equation
	auto [H0w, dH0w, errorw, num_termsw, wrnmsgw]
	  = HeunL0(_Tp{1} / a, (q + alpha * (delta - beta)) / a + alpha * (epsilon - beta),
		  alpha, alpha - gamma + _Tp{1}, alpha - beta + _Tp{1}, delta, _Tp{1} / z);
	value = C1 * std::pow(z, -alpha) * H0w;
	deriv = -C1 * std::pow(z, -alpha - _Tp{1}) * (dH0w / z + alpha * H0w);
	error = std::abs(C1 * std::pow(z,-alpha)) * errorw;
	num_terms = num_termsw;
	wrnmsg += wrnmsgw;
      }

    if (C2 != _Tp{0})
      {
	auto [H0w, dH0w, errorw, num_termsw, wrnmsgw]
	  = HeunS0(_Tp{1} / a, (q + alpha * (delta - beta)) / a + alpha * (epsilon - beta),
		  alpha, alpha - gamma + _Tp{1}, alpha - beta + _Tp{1}, delta, _Tp{1} / z);
	value += C2 * std::pow(z, -alpha) * H0w;
	deriv -= C2 * std::pow(z, -alpha - _Tp{1}) * (dH0w / z + alpha * H0w);
	error += std::abs(C2 * std::pow(z, -alpha)) * errorw;
	num_terms += num_termsw;
	wrnmsg += wrnmsgw;
      }

    return ret;
  }

// Heun function, by power series at Z0 for the given values H(Z0)=H0, H'(Z0)=dH0 
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 04 June 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunGfromZ0(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z,
	      _Tp Z0, _Tp H0, _Tp dH0)
  {
    int Heun_klimit = 200;
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto NaN = std::numeric_limits<_Tp>::quiet_NaN();

    //if (Heun_klimit == 0)
    //  HeunOpts();

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    auto R = std::min({std::abs(Z0), std::abs(Z0 - _Tp{1}), std::abs(Z0 - a)});

    auto epsilon = alpha + beta + _Tp{1} - gamma - delta;

    if (std::abs(z - Z0) >= R)
      {
	wrnmsg = "HeunGfromZ0: z is out of the convergence radius = " << R;
	return ret;
      }
    else if (std::min(std::abs(z - _Cmplx(_Tp{1}, a))) < eps
	  || std::min(std::abs(Z0 - _Cmplx(_Tp{1}, a))) < eps)
      {
	wrnmsg = "HeunGfromZ0: z or Z0 is too close to the singular points"; 
	return ret;
      }
    else if (z == Z0)
      {
	value = H0;
	deriv = dH0; 
	error = _Tp{0};
	num_terms = 0;
	return ret;
      }
    else
      {
	auto zeta = z - Z0;

	auto recur = [alpha, beta, gamma, delta, epsilon, a, q, zeta, Z0](int k, _Tp ckm1, _Tp ckm2, _Tp ckm3)
	  { 
	    return -(ckm1 * zeta * (_Tp{1} - _Tp{1} / k) * ((gamma + delta + epsilon + 3 * (k - 2)) * Z0 * Z0
	     - ((a + _Tp{1}) * (gamma + 2 * k - 4) + epsilon + a * delta) * Z0 + a * (gamma + k - 2))
	    + ckm2 * zeta * zeta * ((2 * (_Tp{1} - 2 / k) * (gamma + delta + epsilon + 3 / 2 * (k - 3)) + alpha * beta / k) * Z0
	    - q / k - (_Tp{1} - 2 / k) * ((a + _Tp{1}) * (gamma + k - 3) + epsilon + a * delta))
	    + ckm3 * zeta*zeta*zeta * ((_Tp{1} - 3 / k) * (gamma + epsilon + delta + k - 4) + alpha * beta / k))
	    / ((k - _Tp{1}) * Z0 * (Z0 - _Tp{1}) * (Z0 - a));
	  };

	auto ckm3 = H0;
	auto ckm2 = dH0 * zeta;
	auto ckm1 = recur(2, ckm2, ckm3, 0);

	value = ckm3 + ckm2 + ckm1;
	auto vm1 = value;
	auto vm2 = NaN;
	auto dm2 = dH0;
	auto dm1 = dH0 + _Tp{2} * ckm1 / zeta;
	auto deriv = dm1;
	auto dderiv = _Tp{2} * ckm1 / zeta / zeta; 

	auto k = 3;
	auto ckm0 = _Tp{1};

	while ((k <= Heun_klimit) && ((vm2 != vm1) || (dm2 != dm1) || (std::abs(ckm0) > eps)))
	  {
	    auto ckm0 = recur(k, ckm1, ckm2, ckm3);
	    value += ckm0;
	    deriv = dm1 + _Tp(k) * ckm0 / zeta;
	    dderiv += k * _Tp(k - 1) * ckm0 / zeta / zeta;
	    ckm3 = ckm2;
	    ckm2 = ckm1;
	    ckm1 = ckm0;
	    vm2 = vm1;
	    vm1 = value;
	    dm2 = dm1;
	    dm1 = deriv;
	    ++k;
	  }

	num_terms = k - 1;

	if (isinf(value) || isinf(deriv) || isnan(value) || isnan(deriv))
	  {
	    wrnmsg = "failed convergence of recurrence and summation"; 
	    value = NaN;
	    deriv = NaN;
	    error = NaN; 
 	    return ret;
	  }
	else
	  {
	    auto value2 = (z * (z - _Tp{1}) * (z - a) * dderiv
		   + (gamma * (z - _Tp{1}) * (z - a) + delta * z * (z - a) + epsilon * z * (z - _Tp{1})) * deriv)
		   / (q - alpha * beta *z);
	    auto error1 = std::abs(value - value2);

	    if (std::abs(q - alpha * beta * z) < 0.01)
	      {
		auto error2 = std::abs(ckm0) * sqrt(num_terms) + std::abs(value) * eps * num_terms;
		error =  std::min(error1, error2);
	      }
	    else
	      error = error1;

	    return ret;
	  }
      }  
  }

// (local) Heun function, based on HeunL0 (analytic continuation from z=0)
// but with improvements near singular points z = 1, a, infinity 
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 15 June 2015
//
template<typename _Tp, typename... _Args>
  __heun_t<_Tp>
  HeunL(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z, _Args... varargin)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    bool gamma_nonpos = std::abs(std::ceil(gamma - 5 * eps) + std::abs(gamma)) < 5 * eps;

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (gamma_nonpos)
      gamma = std::ceil(gamma - 5 * eps);

    if ((std::real(z / a) >= _Tp{1}) && (angle(z / a) == _Tp{0}) || (std::real(z) >= _Tp{1}) && (angle(z) == _Tp{0}) ||
	(gamma_nonpos && (angle(z) == _Tp{0}) && (sign(std::real(z)) < _Tp{0})))
      {
	wrnmsg = "HeunL: z belongs to a possible branch cut; "; 
	return ret;
      }
    else if (std::abs(angle(z)) == pi<_Tp>)
      return HeunL0(a, q, alpha, beta, gamma, delta, z);
    else
      {
	if (sizeof...(varargin) > 0)
	  return HeunLS(_Tp{1}, a, q, alpha, beta, gamma, delta, z, varargin...);
	else
	  return HeunLS(_Tp{1}, a, q, alpha, beta, gamma, delta, z);
      }
    return ret;
  }


// local Heun function, by analytic continuation from z=0 to z
// using a consequence of power expansions
//
// it is assumed that z is not equal to 1 or a;
//
// The function uses a parameter Heun_cont_coef which can be changed by HeunOpts: 
// for each power expansion of the analytic continuation procedure Heun_cont_coef
// is the relative distance from the centre to the calculated point 
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 04 June 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunL0(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z)
  {
    using namespace std::complex_literals;
    using namespace std::string_literals;

    const auto eps = std::numeric_limits<_Tp>::epsilon();
    if (Heun_cont_coef<_Tp> == _Tp{0})
      HeunOpts();

    bool gamma_nonpos = std::abs(std::ceil(gamma - 5 * eps) + std::abs(gamma)) < 5 * eps;

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (gamma_nonpos)
      gamma = std::ceil(gamma - 5 * eps);

    if (std::real(z / a) >= _Tp{1}
	&& std::imag(z / a) == _Tp{0}
	|| (std::real(z) >= _Tp{1}) && (std::imag(z) == _Tp{0})
	|| (gamma_nonpos && std::imag(z) == _Tp{0} && sign(z) < _Tp{0}))
      {
	wrnmsg = "HeunL0: z belongs to a possible branch cut; "; 
	return ret;
      }
    else if (std::min(std::abs(z - std::complex<_Tp>(_Tp{1}, a))) < eps)
      {
	wrnmsg = "HeunL0: z is too close to one of the singular points; "; 
	return ret;
      }
    else
      {
	const auto R0 = std::min(_Tp{1}, std::abs(a));

	if (std::abs(z) <= R0 * (Heun_cont_coef<_Tp> + 0.01))
	  {
	    if (gamma_nonpos)
	      return HeunL00log(a, q, alpha, beta, gamma, delta, z);
	    else
	      return HeunL00(a, q, alpha, beta, gamma, delta, z);
	  }
	else
	  {
	    auto [P1, P2] = std::minmax(a, _Tp{1});

	    const auto R1 = std::min(std::abs(P1), std::abs(a - _Tp{1}));
	    const auto R2 = std::min(std::abs(P2), std::abs(a - _Tp{1}));

	    std::vector<_Tp> zz;
	    zz.reserve(4);

	    zz.push_back(_Tp{0});

	    auto [d, isinside] = dist(P1, _Tp{0}, z);
	    if (isinside && (d < R1 / 2))
	      zz.push_back(P1 + std::exp(1i * (pi<_Tp> / 2 + angle(z)))
				 * std::min(R1 / 2, std::abs(z - P1))
				 * sign(std::imag(z / P1)));

	    std::tie(d, isinside) = dist(P2, _Tp{0}, z);
	    if (isinside && (d < R2 / 2))
	      zz.push_back(P2 + std::exp(1i * (pi<_Tp> / 2 + angle(z)))
				 * std::min(R2 / 2, std::abs(z - P2))
				 * sign(std::imag(z / P2)));

	    zz.push_back(z);

	    auto sum_error = _Tp{0};
	    auto sum_terms = _Tp{0};

	    auto failure = false;

	    for (int k = 1; k < std::size(zz); ++k)
	      {
		auto z0 = zz[k];
		auto theta = angle(zz[k + 1] - z0);
		auto insearch = true;

		while (insearch && !failure)
		  {
		    const auto R = z0 == _Tp{0}
				 ? R0
				 : std::min({std::abs(z0), std::abs(z0 - 1),
					     std::abs(z0 - a)});

		    decltype(zz[0]) z1;
		    if (std::abs(zz[k + 1] - z0) <= R * Heun_cont_coef<_Tp>)
		      {
			z1 = zz[k + 1];
			insearch = false;
		      }
		    else
		      z1 = z0 + Heun_cont_coef<_Tp> * R * std::exp(1i * theta);

		    if (z0 == _Tp{0})
		      {
			std::tie(value, deriv, error, num_terms, wrnmsg)
			  = gamma_nonpos
			  ? HeunL00log(a, q, alpha, beta, gamma, delta, z1)
			  : HeunL00(a, q, alpha, beta, gamma, delta, z1);

			if (wrnmsg.size() != 0)
			  {
			    wrnmsg = "HeunL0: "
				     "problem invoking HeunL00";
			    if (gamma_nonpos)
			      wrnmsg << "log";
			    wrnmsg
			      << '(' << a << ','
			      << q << ',' << alpha << ',' << beta << ','
			      << gamma << ','<< delta << ',' << z1
			      << "); warning: " << wrnmsg << "; ";
			    failure = true;
			  }
		      }
		    else
		      {
			std::tie(value, deriv, error, num_terms, wrnmsg)
			   = HeunGfromZ0(a, q, alpha, beta, gamma, delta,
					 z1, z0, value, deriv);
			if (wrnmsg.size() != 0)
			  {
			    wrnmsg = "HeunL0: "
				     "problem invoking HeunGfromZ0("
				   << a << ',' << q
				   << ',' << alpha << ',' << beta
				   << ',' << gamma << ',' << delta
				   << ',' << z0 << ',' << value << ',' << deriv
				   << "); warning: "<< wrnmsg << "; ";
			    failure = true;
			  }
		      }

		    sum_error += error;
		    sum_terms += num_terms;

		    z0 = z1;
		  }

		if (failure)
		  break;

	      }

	    num_terms = sum_terms;
	    error = sum_error;
	    return ret;
	  }
      }
  }

// local Heun function, by power series at z=0 for Hl(0)=1, Hl'(0)=q/(a*gamma)
// the case when gamma is not equal to 0, -1, -2, 
//
// |z| should not exceed the convergency radius min{1,|a|}
//
// The function uses a parameter Heun_klimit which can be changed by HeunOpts: 
// Heun_klimit is the maximum number of series' terms
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 28 April 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunL00(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto NaN = std::numeric_limits<_Tp>::quiet_NaN();
    int Heun_klimit = 200;

    //if (isempty(Heun_klimit))
    //  HeunOpts();

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    auto R = std::min(1, std::abs(a));

    auto epsilon = alpha + beta + 1 - gamma - delta;

    auto gamma_nonpos = std::abs(std::ceil(gamma - 5 * eps) + std::abs(gamma)) < 5 * eps;

    if (gamma_nonpos)
      {
	wrnmsg = "HeunL00: gamma is a non-positive integer, use HeunL00log; ";
	return ret;
      }
    else if (std::abs(z) >= R)
      {
	wrnmsg = "HeunL00: z is out of the convergence radius = " << R << "; ";
	return ret;
      }
    else if (z == _Tp{0})
      {
	value = 1;
	deriv = q / (a * gamma);
	error = _Tp{0};
	num_terms = 1;
	return ret;
      }
    else
      {
	auto recur = [a, q, alpha, beta, gamma, delta, epsilon, z](int k, _Tp ckm1, _Tp ckm2)
	  {
	    return (ckm1 * z * (q / k + (1 - 1 / k) * ((a + 1) * gamma + epsilon + a * delta) + (k + 2 / k - 3) * (a + 1)) -
	   ckm2 * z*z * ((1 - 2 / k) * (gamma + epsilon + delta) + 1 / k * alpha * beta + k + 6 / k - 5)) / (a * gamma + (k - 1) * a);
	  };

	auto ckm2 = _Tp{1};
	auto ckm1 = z * q / (a * gamma);

	value = ckm2 + ckm1;
	auto vm1 = value;
	auto vm2 = NaN;
	deriv = q / (a * gamma);
	auto dm1 = deriv;
	auto dm2 = NaN;
	auto dderiv = _Tp{0};

	int k = 2;
	auto ckm0 = _Tp{1};

	while ((k <= Heun_klimit) && ((vm2 != vm1) || (dm2 != dm1) || (std::abs(ckm0) > eps)))
	  {
	    ckm0 = recur(k, ckm1, ckm2);
	    value += ckm0;
	    deriv = dm1 + k * ckm0 / z;
	    dderiv += k * (k - 1) * ckm0 / z / z;
	    ckm2 = ckm1;
	    ckm1 = ckm0;
	    vm2 = vm1;
	    vm1 = value;
	    dm2 = dm1;
	    dm1 = deriv;
	    ++k;
	  }

	num_terms = k - 1;

	if (isinf(value) || isinf(deriv) || isnan(value) || isnan(deriv))
	  {
	    wrnmsg = "HeunL00: failed convergence of recurrence and summation; "; 
	    value = NaN;
	    deriv = NaN;
	    error = NaN;
	    return ret;
	  }
	else
	  {
	    auto value2 = (z * (z - 1) * (z - a) * dderiv + (gamma * (z - 1) * (z - a) + delta * z * (z - a) + epsilon * z * (z - 1)) * deriv) / (q - alpha * beta * z);
	    auto error1 = std::abs(value - value2);

	    if (std::abs(q - alpha * beta * z) < 0.01)
	      {
		auto error2 = std::abs(ckm0) * sqrt(num_terms) + eps * num_terms * std::abs(value);
		error = std::min(error1, error2);
	      }
	    else
	      error = error1;

	    return ret;
	  }
      }
  }

// local Heun function, by power series at z=0, for Hl(0)=1
// the case when gamma = 0, -1, -2, ...
//
// |z| should not exceed the convergency radius min{1,|a|}
//
// The function uses a parameter Heun_klimit which can be changed by HeunOpts: 
// Heun_klimit is the maximum number of series' terms
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 28 April 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunL00log(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto NaN = std::numeric_limits<_Tp>::quiet_NaN();
    int Heun_klimit = 200;

    //if (isempty(Heun_klimit))
    //  HeunOpts();

    auto R = std::min(_Tp{1}, std::abs(a));

    auto gamma_nonpos = std::abs(std::ceil(gamma - 5*eps) + std::abs(gamma)) < 5*eps;

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (!gamma_nonpos)
      std::__throw_domain_error("HeunL00log: gamma is not a non-positive integer");
    else if (std::abs(z) >= R)
      {
	wrnmsg = "HeunL00log: z is outside of the convergence radius = " << R;
	return ret;
      }
    else
      {
	const auto epsilon = alpha + beta + 1 - gamma - delta;

	auto N = 1 - gamma;

	auto recur0 = [a, q, alpha, beta, gamma, delta, epsilon, z](int k, _Tp ckm1, _Tp ckm2) 
	  {
	    return (ckm1 * z  * (q / k + (1 - 1 / k) * ((a + 1) * gamma + epsilon + a * delta) + (k + 2 / k - 3) * (a + 1)) - 
		  ckm2 * z * z * ((1 - 2 / k) * (gamma + epsilon + delta) + 1 / k * alpha * beta + k + 6 / k - 5))/(a * gamma +(k - 1) * a);
	  };

	auto recur1 = [recur0, a, q, alpha, beta, gamma, delta, epsilon, z](int k, _Tp ckm1, _Tp ckm2, _Tp dkm0, _Tp dkm1, _Tp dkm2)
	  {
	    return recur0(k, ckm1, ckm2) + 
		  (dkm0 * a * (1 - gamma - 2 * k) + dkm1 * z * (epsilon + a * delta + (a + 1) * (gamma + 2 * k - 3)) + 
	 	  dkm2 * z * z * (4 - 2 * k - alpha - beta)) / (a * k  * (gamma + k - 1));
	  };

	auto L1 = _Tp{1};
	auto dL1 = _Tp{0};
	auto ddL1 = _Tp{0};
	auto ckm0 = _Tp{1};
	auto ckm1 = _Tp{1};
	auto ckm2 = _Tp{0};

	for (int k = 1; k <= N - 1; ++k)
	  {  
	    ckm0 = recur0(k, ckm1, ckm2);
	    L1 += ckm0;
	    dL1 += k * ckm0 / z;
	    ddL1 += k * (k - 1) * ckm0 / z / z;
	    ckm2 = ckm1;
	    ckm1 = ckm0; 
	  }

	auto sN = (ckm1 * z * (q - gamma * (epsilon + a * delta - a - 1))
		- ckm2 * z * z * (alpha * beta - (gamma + 1) * (epsilon + delta - 2)))
			 / (a * (1 - gamma));

	auto L2 = _Tp{0};
	auto dL2 = _Tp{0};
	auto ddL2 = _Tp{0};
	auto dm1 = dL2;
	auto dm2 = NaN;
	ckm1 = _Tp{0};
	ckm2 = ckm0; 

	auto L3 = sN;
	auto skm2 = _Tp{0};
	auto skm1 = sN;
	auto dL3 = N * sN / z;
	auto ddL3 = N * (N - 1) * sN / z / z; 
	auto dsm1 = dL3;
	auto dsm2 = NaN;
	skm0 = NaN;

	int k = N + 1;

	while ((k <= Heun_klimit) && ((dsm2 != dsm1) || (std::abs(skm0) > eps) || 
	  (dm2 != dm1) || (std::abs(ckm0) > eps)))
	  {
	    skm0 = recur0(k, skm1, skm2);
	    ckm0 = recur1(k, ckm1, ckm2, skm0, skm1, skm2);

	    L2 += ckm0;
	    dL2 = dm1 + k * ckm0 / z;
	    ddL2 += k * (k - 1) * ckm0 / z / z;
	    ckm2 = ckm1;
	    ckm1 = ckm0;
	    dm2 = dm1;
	    dm1 = dL2;

	    L3 += skm0;
	    dL3 = dsm1 + k * skm0 / z;
	    ddL3 += k * (k - 1) * skm0 / z / z;
	    skm2 = skm1;
	    skm1 = skm0;
	    dsm2 = dsm1;
	    dsm1 = dL3;

	    ++k;
	  }

	num_terms = k - 1;

	value = L1 + L2 + std::log(z) * L3;
	deriv = dL1 + dL2 + std::log(z) * dL3 + L3 / z;
	auto dderiv = ddL1 + ddL2 - L3 / z / z + 2 * dL3 / z + std::log(z) * ddL3;

	if (isinf(value) || isinf(deriv) || isnan(value) || isnan(deriv))
	  {
	    wrnmsg = "HeunL00log: failed convergence of recurrence and summation; "; 
	    value = NaN;
	    deriv = NaN;
	    error = NaN;
	    return ret;
	  }
	else
	  {
	    auto value2 = (z * (z - 1) * (z - a) * dderiv + (gamma * (z - 1) * (z - a) + delta * z * (z - a) + epsilon * z * (z - 1)) * deriv) / 
		    (q - alpha * beta * z);
	    auto error1 = std::abs(value - value2);

	    if (std::abs(q - alpha * beta * z) < 0.01)
	      {
		auto error2 = std::abs(L1) * eps * N + std::abs(ckm0) * sqrt(num_terms - N + 1) + std::abs(L2) * eps * (num_terms - N + 1) + 
		     std::abs(std::log(z)) * (std::abs(skm0) * sqrt(num_terms - N + 1) + std::abs(L3) * eps * (num_terms - N + 1));
		error =  std::min(error1,error2);
	      }
	    else
	      error = error1;
	    return ret;
	  }
      }
  }

// local Heun function, based on HeunL0 or HeunS0
// by analytic continuation from z=0
// with improvements near singular points z = 1, a, infinity 
//
// if numfunc=1 HeunLS computes the value of the first local solution HeunL
// otherwise HeunLS computes the value of the second local solution HeunS
//
// the optional parameter memlimit (500 as default) is the maximum number of already
//   computed matching data which are kept in memory
//
// The function uses parameters Heun_proxco and Heun_proxcoinf which can be changed
// by HeunOpts: they specify relative proximity to singular point where we use the
// special representations
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 10 June 2015
//
template<typename _Tp, typename... _Args>
  __heun_t<_Tp>
  HeunLS(int numfunc, _Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z,
	 _Args... varargin)
  {
    if (sizeof...(varargin) > 0)
      memlimit = varargin{1};
    else
      memlimit = 500;

    global Heun_proxco1st Heun_proxcoinf1st Heun_proxco Heun_proxcoinf;

    if (isempty(Heun_proxco1st)||isempty(Heun_proxco)||isempty(Heun_proxcoinf1st)||isempty(Heun_proxcoinf))
      HeunOpts();

    Heun0 = HeunL0;

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (numfunc != 1)
      {
	numfunc = 2;
	Heun0 = HeunS0;
      }

    auto Ra = std::min(std::abs(a), std::abs(a-1));
    auto R1 = std::min(1, std::abs(a-1));
    auto Rinf = std::max(1, std::abs(a));

    std::string wrnmsg;
    bool failed = true;

    if (std::abs(z - a) < Heun_proxco * Ra)
      {
	singpt = "a";
	sdir = a * sign(std::imag(a)) / std::abs(a);

	if ((std::imag(a) == _Tp{0}) && ((a < _Tp{0}) || (a > 1)))
	  {
	    if (std::imag(z) > _Tp{0})
	      {
		singpt = "aU";
		sdir = 1;
	      }
	    else
	      {
		singpt = "aL";
		sdir = -1;
	      }
	  } 

	midpoint = a * 0.5 + 0.70710678i * sdir;
	impev = std::abs(z - a) < Heun_proxco1st * Ra;

	[value, deriv, error, num_terms, wrnmsg, failed]
	  = HeunLSnearsing(numfunc, Heun0, HeunG_near_a, singpt, midpoint, memlimit, impev,
			   a, q, alpha, beta, gamma, delta, z);

	return ret;
      }
    else if (std::abs(z - 1) < Heun_proxco * R1)
      {
	singpt = "1";
	sdir = -sign(std::imag(a));

	if ((std::imag(a) == _Tp{0}) && (a > _Tp{0}) && (a < 1))
	  {
	    if (std::imag(z) > _Tp{0})
	      {
		singpt = "1U";
		sdir = 1;
	      }
	    else
	      {
		singpt = "1L";
		sdir = -1;
	      }
	  } 

	midpoint = 0.5 + 0.70710678i * sdir;
	impev = std::abs(z - 1) < Heun_proxco1st * R1;

	[value, deriv, error, num_terms, wrnmsg, failed]
	  = HeunLSnearsing(numfunc, Heun0, HeunG_near_1, singpt, midpoint, memlimit, impev,
			   a, q, alpha, beta, gamma, delta, z);

	return ret;
      }
    else if (std::abs(z) > Heun_proxcoinf * Rinf)
      {
	ls = sort({-pi<_Tp>,_Tp{0},pi<_Tp>,angle(a)});
	idx = sum(angle(z)>ls);
	singpt = 'I' << idx;
	midarg = (ls(idx) +ls(idx+1))/2;

	midpoint = Heun_proxcoinf * Rinf * std::exp(1i*midarg);
	impev = std::abs(z)>Heun_proxcoinf1st*Rinf;

	[value, deriv, error, num_terms, wrnmsg, failed]
	  = HeunLSnearsing(numfunc, Heun0, HeunG_near_infty, singpt, midpoint, impev, memlimit,
			   a, q, alpha, beta, gamma, delta, z);

	return ret;
      }

    if (failed)
      return Heun0(a, q, alpha, beta, gamma, delta, z);
    else
      return ret; // Added EMSR
  }


  function [value, deriv, error, num_terms, wrnmsg, failed]
  HeunLSnearsing(numfunc, Heun0, HeunG_nearsing, singpt, midpoint,
		 memlimit, impev, _Tp a, q_Tp , _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z)
  {
    const auto NaN = std::numeric_limits<_Tp>::quiet_NaN();
    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    [A1, A2, errorco, consts_known] = extrdatfromsav(numfunc, a, q, alpha, beta, gamma, delta, singpt);

    if (!consts_known && impev)
    {
      auto [value1, deriv1, error1, num_terms1, wrnmsg1]
	= HeunG_nearsing(a, q, alpha, beta, gamma, delta, midpoint, 1, _Tp{0});
      auto [value2, deriv2, error2, num_terms2, wrnmsg2]
	= HeunG_nearsing(a, q, alpha, beta, gamma, delta, midpoint, _Tp{0}, 1);
      auto [value0, deriv0, error0, num_terms0, wrnmsg0]
	= Heun0(a, q, alpha, beta, gamma, delta, midpoint);
      auto M = [value1 value2;
       deriv1 deriv2];
      auto dcs = M \ [value0; deriv0];
      auto A1 = dcs(1);
      auto A2 = dcs(2);

      auto errorco = error0 + error1 + error2;
      auto num_terms0 = num_terms0 + num_terms1 + num_terms2;

      if (std::min(svd(M)) / std::max(svd(M)) < 10e-6) 
	{    
	  A1 = NaN;
	  A2 = NaN;
	}

      keepdattosav(numfunc, a, q, alpha, beta, gamma, delta, errorco, A1, A2, singpt, memlimit);
      wrnmsg = wrnmsg0+ wrnmsg1+ wrnmsg2;
    }

    failed = isnan(A1) || isnan(A2) || (!consts_known && !impev);

    if (!failed)
      {
      [value1, deriv1, error1, num_terms1, wrnmsg1]
	= HeunG_nearsing(a, q, alpha, beta, gamma, delta, z, 1, _Tp{0});
      if (!isempty(wrnmsg1))
	{
	  if (std::abs(A1) < errorco)
	    {
	       wrnmsg1 = "";
	       value1 = _Tp{0};
	       deriv1 = _Tp{0};
	       error1 = _Tp{0}; 
	    }
	}

      auto [value2, deriv2, error2, num_terms2, wrnmsg2]
	= HeunG_nearsing(a, q, alpha, beta, gamma, delta, z, _Tp{0}, 1);
      if (!isempty(wrnmsg2))
	{
	  if (std::abs(A2) < errorco)
	    {
	       wrnmsg2 = "";
	       value2 = _Tp{0};
	       deriv2 = _Tp{0};
	       error2 = _Tp{0}; 
	    }
	}
      wrnmsg += wrnmsg1 + wrnmsg2;

      value = A1 * value1 + A2 * value2;
      deriv = A1 * deriv1 + A2 * deriv2;
      error = std::abs(A1) * error1 + std::abs(A2) * error2 + errorco;
      num_terms = num_terms1 + num_terms2;

      if (!consts_known)
	num_terms += num_terms0;
    }

    failed ||= isnan(value);
  }
  
  function [A1, A2, errorco, consts_known]
  extrdatfromsav(numfunc, _Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, singpt)
  {
    global savdata;
    const auto NaN = std::numeric_limits<_Tp>::quiet_NaN();

    auto A1 = NaN;
    auto A2 = NaN;
    auto errorco = NaN;
    bool consts_known = false;
    if (std::size(savdata) != 0)
      for (k = 1 : std::size(savdata))
	if (savdata(k).numfunc==numfunc)&&strcmp(savdata(k).singpt,singpt)&&(savdata(k).a==a)&&
	   (savdata(k).q==q)&&(savdata(k).alpha==alpha)&&(savdata(k).beta==beta)&&
	   (savdata(k).gamma==gamma)&&(savdata(k).delta==delta)
	  {
	    A1 = savdata(k).A1;
	    A2 = savdata(k).A2;
	    errorco = savdata(k).errorco;
	    consts_known = true;
	    break;
	  }
  }

void
keepdattosav(numfunc, a, q, alpha, beta, gamma, delta, errorco, A1, A2, singpt, memlimit)
 {
  global savdata;

  if (std::size(savdata) == 0)
    savdata=struct("numfunc",numfunc,"singpt",singpt,"a",a,"q",q,"alpha",alpha,
      "beta",beta,"gamma",gamma,"delta",delta,"errorco",errorco,"A1",A1,"A2",A2);
  else
    {
      if (std::size(savdata) <= memlimit)
	{
	  savdata(end+1) = struct("numfunc",numfunc,"singpt",singpt,"a",a,"q",q,"alpha",alpha,
	  "beta",beta,"gamma",gamma,"delta",delta,"errorco",errorco,"A1",A1,"A2",A2);
	}
      else
	{
	  savdata(1).numfunc=numfunc; savdata(1).singpt=singpt; savdata(1).a=a;
	  savdata(1).q=q; savdata(1).alpha=alpha; savdata(1).beta=beta;
	  savdata(1).gamma=gamma; savdata(1).delta=delta;
	  savdata(1).errorco=errorco; savdata(1).A1=A1; savdata(1).A2=A2;
	  savdata = shift(savdata,-1);
	}  
    }
}

// local Heun function, by analytic continuation from z=0 to z
// using a consequence of power expansions
//
// path2z is a list of points from 0 to z
//
// The function uses a parameter Heun_cont_coef which can be changed by HeunOpts: 
// for each power expansion of the analytic continuation procedure Heun_cont_coef
// is the relative distance from the centre to the calculated point 
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 04 June 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunLmv(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta,
	  const std::vector<_Tp>& path2z)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();

    if (Heun_cont_coef<_Tp> == _Tp{0})
      HeunOpts();

    gamma_nonpos = std::abs(std::ceil(gamma - 5 * eps) + std::abs(gamma)) < 5 * eps;

    if (gamma_nonpos)
      gamma = std::ceil(gamma - 5 * eps);

    if (path2z[0] != _Tp{0})
      {
	wrnmsg = "HeunLmv: path2z should start from zero; "; 
	return ret;
      }
    else if (dist2sing(path2z,  {_Tp{1}, a}) < _Tp{10} * eps
	  || dist2sing(path2z(2 : end), {_Tp{0}}) < _Tp{10} * eps)
      {
	wrnmsg = "HeunLmv: path2z is too close to one of the singular points; "; 
	return ret;
      }
    else
      {
	R0 = std::min(1, std::abs(a));
	z = path2z.back();

	if (std::abs(z) <= R0 * (Heun_cont_coef<_Tp> + 0.01))
	  {
	    if (gamma_nonpos)
	      return HeunL00log(a, q, alpha, beta, gamma, delta, z);
	    else
	      return HeunL00(a, q, alpha, beta, gamma, delta, z);
	  }
	else
	  {
	    auto sum_error = _Tp{0};
	    auto sum_terms = _Tp{0};

	    auto failure = false;

	    for (k = 1 : (std::size(path2z) - 1))
	      {
		z0 = path2z[k];
		theta = angle(path2z[k + 1] - z0);
		insearch = true;

		while (insearch && !failure)
		  {
		    if (z0 == _Tp{0})
		      R = R0;
		    else
		      R = std::min({std::abs(z0), std::abs(z0 - 1), std::abs(z0 - a)});

		    if (std::abs(path2z[k + 1] - z0) <= R * Heun_cont_coef<_Tp>)
		      {
			z1 = path2z[k + 1];
			insearch = false;
		      }
		    else
		      z1 = z0 + Heun_cont_coef<_Tp> * R * std::exp(1i*theta);

		    if (z0 == _Tp{0})
		      {
			if (gamma_nonpos)
			  {
			    [H0, dH0, error, num_terms, wrnmsg]
			       = HeunL00log(a, q, alpha, beta, gamma, delta, z1);
			    st = "log";
			  }
			else
			  {
			    [H0, dH0, error, num_terms, wrnmsg] = HeunL00(a, q, alpha, beta, gamma, delta, z1);
			    st = "";
			  }

			if (wrnmsg.length() != 0)
			  {
			    wrnmsg = "HeunLmv: problem invoking HeunL00"<< st<< "("<< a << ','
			      << q << ',' << alpha << ',' << beta << ','
			      << gamma << ',' << delta << ',' << z1
			      << "); warning: " << wrnmsg << "; ";
			    failure = true;
			  }
		      }
		      else
		      {
			[H0, dH0, error, num_terms, wrnmsg]
			   = HeunGfromZ0(a,q,alpha,beta,gamma,delta,z1,z0,H0,dH0);
			if (wrnmsg.length() != 0)
			  {
			    wrnmsg = "HeunLmv: problem invoking HeunGfromZ0("<< a << ','
			      << q << ',' << alpha << ',' << beta << ','
			      << gamma << ',' << delta << ',' << z1 << ','
			      << z0 << ',' << H0 << ',' << dH0
			      << "); warning: " << wrnmsg << "; ";
			    failure = true;
			  }
		      }

		    sum_error += error;
		    sum_terms += num_terms;

		    z0 = z1;
		  }

		if failure
		  break;	
	      }

	    value = H0;
	    deriv = dH0;
	    num_terms = sum_terms;
	    error = sum_error;
	    return ret;
	  }
      }  
  }


// Change internal parameters of HeunL, HeunL00, HeunL00log, HeunGfromZ0
//
// Usage: HeunOpts("cont_coef",Heun_cont_coef,"klimit",Heun_klimit,
//		 "proxco1st",Heun_proxco_1st,"proxcoinf1st",Heun_proxcoinf_1st)
//		 "proxco",Heun_proxco,"proxcoinf",Heun_proxcoinf)
//
// parameters are optional
//
// if a pair of parameters is omitted, but the corresponding global variable
// is empty then it is set to the default value
//
// Use [] as the value to reset coefficient to its default
// e.g. HeunOpts("cont_coef",[])
//
// Heun_cont_coef is used in HeunL0 and HeunS0gamma1; for each power expansion
// of the analytic continuation procedure Heun_cont_coef is the relative
// (to the radius of convergence) distance from the centre to the calculated point.
// By default Heun_cont_coef = 0.38 (golden ratio)
//
// Heun_klimit is used in HeunL00, HeunL00log, HeunGfromZ0, HeunS0gamma1, 
// HeunS00gamma1; it is the maximum number of power series' terms.
// Default value is 1000 
//
// Heun_proxco1st, Heun_proxcoinf1st, Heun_proxco, Heun_proxcoinf are used in HeunLS;
// they specify relative proximity to singular point where special representation is used
// "1st" for the case when the matching coefficients are not known 
// By default Heun_proxco1st = 0.05; Heun_proxcoinf1st = 10; Heun_proxco = 0.25; Heun_proxcoinf = 2 
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 05 June 2015
//
template<typename... _Args>
  void
  HeunOpts(_Args... varargin)
  {
    int Heun_klimit = 200;
    poop Heun_proxco Heun_proxcoinf Heun_proxco1st Heun_proxcoinf1st;

    [reg, props] = parseparams(varargin...);
    opts = cell2struct(props(2:2:end),props(1:2:end),2);

    if (isfield(opts,"cont_coef"))
      {
	if (isempty(opts.cont_coef))
	  Heun_cont_coef = 0.38;
	else
	  Heun_cont_coef = opts.cont_coef;
      }
    else if (Heun_cont_coef<_Tp> == _Tp{0})
      Heun_cont_coef = 0.38;

    if (isfield(opts,"klimit"))
      {
	if (isempty(opts.klimit))
	  Heun_klimit = 1000;
	else
	  Heun_klimit = opts.klimit;
      }
    else if (isempty(Heun_klimit))
      Heun_klimit = 1000;

    if (isfield(opts,"proxco1st"))
      {
	if (isempty(opts.proxco1st))
	  Heun_proxco1st = 0.05;
	else
	  Heun_proxco1st = opts.proxco1st;
      }
    else if (isempty(Heun_proxco1st))
      Heun_proxco1st = 0.05;

    if (isfield(opts,"proxco"))
      {
	if (isempty(opts.proxco))
	  Heun_proxco = 0.25;
	else
	  Heun_proxco = opts.proxco;
      }
    else if (isempty(Heun_proxco))
      Heun_proxco = 0.25;

    if isfield(opts,"proxcoinf1st")
      {
	if (isempty(opts.proxcoinf1st))
	  Heun_proxcoinf1st = 5;
	else
	  Heun_proxcoinf1st = opts.proxcoinf1st;
      }
    else if (isempty(Heun_proxcoinf1st))
      Heun_proxcoinf1st = 5;

    if (isfield(opts,"proxcoinf"))
      {
	if (isempty(opts.proxcoinf))
	  Heun_proxcoinf = 2;
	else
	  Heun_proxcoinf = opts.proxcoinf;
      }
    else if (isempty(Heun_proxcoinf))
      Heun_proxcoinf = 2;  
  }

// second local Heun function, based on HeunS0 (analytic continuation from z=0)
//
// with improvements near singular points z = 1, a, infinity 
//
// the optional parameter memlimit (500 as default) is the maximum number of already
//   computed matching data which are kept in memory
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 05 June 2015
//
template<typename _Tp, typename... _Args>
  __heun_t<_Tp>
  HeunS(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z,
	_Args... varargin)
  {
    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if ((std::real(z / a) >= 1 && std::imag(z / a) == _Tp{0}) || (std::real(z) >= 1 && std::imag(z) == _Tp{0}))
      {
	wrnmsg = "HeunS: z belongs to a possible branch cut; "; 
	return ret;
      }
    else
      {
	if (sizeof...(varargin) > 0)
	  return HeunLS(2, a, q, alpha, beta, gamma, delta, z, varargin...);
	else
	  return HeunLS(2, a, q, alpha, beta, gamma, delta, z);
      }
    return ret;
  }

// the second local Heun function at z=0
//
// evaluation by analytic continuation from z=0
// to z using a consequence of power expansions
//
// it is assumed that z is not equal to 1 or a;
//
// The function uses a parameter Heun_cont_coef which can be changed by HeunOpts: 
// for each power expansion of the analytic continuation procedure Heun_cont_coef
// is the relative distance from the centre to the calculated point 
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 01 May 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunS0(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta, _Tp z)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();

    if (Heun_cont_coef<_Tp> == _Tp{0})
      HeunOpts();

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if ((std::real(z/a) >= 1 && std::imag(z/a) == _Tp{0})
     || (std::real(z) >= 1 && std::imag(z) == _Tp{0}))
      { 
	wrnmsg = "HeunS0: z belongs to a possible branch cut; "; 
	return ret;
      }
    else if (std::min(std::abs(z - [1 a])) < eps)
      {
	wrnmsg = "HeunS0: z is too close to one of the singular points; "; 
	return ret;
      }
    else
      {
	if (gamma == 1)
	  return HeunS0gamma1(a, q, alpha, beta, delta, z);
	else
	  {
	    epsilon = alpha + beta + 1 - gamma - delta;
	    auto [H0w, dH0w, errorw, num_terms, wrnmsg]
		   = HeunL0(a, q - (gamma - 1) * (epsilon + a * delta),
			    beta - gamma + 1, alpha - gamma + 1, 2 - gamma, delta, z);
	    value = std::pow(z, (1 - gamma) * H0w);
	    deriv = (1 - gamma) * std::pow(z, -gamma) * H0w + std::pow(z, 1 - gamma) * dH0w;
	    error = std::abs(std::pow(z, 1 - gamma) * errorw;
	    return ret;
	  }
      }
  }

// the second local Heun function at z=0
// the case when gamma = 1
//
// evaluation by power series
//
// |z| should not exceed the convergency radius min{1,|a|}
//
// The function uses a parameter Heun_klimit which can be changed by HeunOpts: 
// Heun_klimit is the maximum number of series" terms
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 28 April 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunS00gamma1(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp delta, _Tp z)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto NaN = std::numeric_limits<_Tp>::quiet_NaN();
    int Heun_klimit = 200;

    if (isempty(Heun_klimit))
      HeunOpts();

    auto R = std::min(_Tp{1}, std::abs(a));

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (std::abs(z) >= R)
      {
	wrnmsg = "HeunS00gamma1: z is out of the convergence radius = "<< R <<"; ";
	return ret;
      }
    else
      {
	epsilon = alpha + beta - delta;

	auto recur0 = [alpha, beta, gamma, delta, epsilon, a, z](int k, _Tp ckm1, _Tp ckm2) 
	  {
	    return (ckm1 * z * (q / k + (1 - 1 / k) * ((a + 1) + epsilon + a * delta) + (k + 2 / k - 3) * (a + 1)) - 
	   ckm2 * z * z * ((1 - 2 / k) * (1 + epsilon + delta) + 1 / k * alpha * beta + k + 6 / k - 5)) / (a * k);
	  };

	auto recur1 = [alpha, beta, epsilon, a, z](int k, _Tp ckm1, _Tp ckm2, _Tp dkm0, _Tp dkm1, _Tp dkm2)
	  {
	    return recur0(k, ckm1, ckm2) + 
		  (-dkm0 * a * 2 + dkm1 * z * ((epsilon + a * delta) / k + (a + 1) * 2 * (1 - 1 / k)) + 
		  dkm2 * z * z * ((4 - alpha - beta) / k - 2)) / (a * k);
	  };

	L1 = _Tp{0};
	dL1 = _Tp{0};
	ddL1 = _Tp{0};
	dm1 = _Tp{0};
	dm2 = NaN;
	ckm0 = NaN;
	ckm1 = _Tp{0};
	ckm2 = _Tp{0}; 

	L2 = 1;
	dL2 = _Tp{0};
	ddL2 = _Tp{0};
	skm2 = _Tp{0};
	skm1 = 1;

	dsm1 = _Tp{0};
	dsm2 = NaN;
	skm0 = NaN;

	k = 1;

	while ((k <= Heun_klimit) && ((dsm2!=dsm1) || (std::abs(skm0)>eps) || 
	  {
	    (dm2 != dm1) || (std::abs(ckm0) > eps)))

	    skm0 = recur0(k, skm1, skm2);
	    ckm0 = recur1(k, ckm1, ckm2, skm0, skm1, skm2);

	    L1 += ckm0;
	    dL1 = dm1 + k * ckm0 / z;
	    ddL1 += k * (k - 1) * ckm0 / z / z;
	    ckm2 = ckm1;
	    ckm1 = ckm0;
	    dm2 = dm1;
	    dm1 = dL1;

	    L2 += skm0;
	    dL2 = dsm1 + k * skm0 / z;
	    ddL2 += k * (k - 1) * skm0 / z / z;
	    skm2 = skm1;
	    skm1 = skm0;
	    dsm2 = dsm1;
	    dsm1 = dL2;

	    ++k;
	  }

	num_terms = k-1;

	value = L1 + std::log(z) * L2;
	deriv = dL1 + std::log(z) * dL2 + L2 / z;
	dderiv = ddL1 - L2 / z / z + 2 * dL2 / z + std::log(z) * ddL2;

	if (isinf(value) || isinf(deriv) || isnan(value) || isnan(deriv))
	  {
	    wrnmsg = "HeunS00gamma1: failed convergence of recurrence and summation; "; 
	    value = NaN;
	    deriv = NaN;
	    error = NaN;
	  }
	else
	  {
	    value2 = (z * (z - 1) * (z - a) * dderiv
		    + ((z - 1) * (z - a)
		    + delta * z * (z - a) + epsilon * z * (z - 1)) * deriv) / 
		    (q - alpha * beta * z);
	    error1 = std::abs(value - value2);

	    if (std::abs(q - alpha * beta * z) < 0.01)
	      {
		error2 = std::abs(ckm0) * sqrt(num_terms) + std::abs(L1) * eps * num_terms
		       + std::abs(std::log(z)) * (std::abs(skm0) *sqrt(num_terms) + std::abs(L2) * eps * num_terms);
		error =  std::min(error1, error2);
	      }
	    else
	      error = error1;
	  }
	return ret;
      }
  }

// the second local Heun function at z=0
// the case when gamma = 1
//
// evaluation by analytic continuation from z=0
// to z using consequence of power expansions
//
// it is assumed that z is not equal to 1 or a;
//
// The function uses a parameter Heun_cont_coef which can be changed by HeunOpts: 
// for each power expansion of the analytic continuation procedure Heun_cont_coef
// is the relative distance from the centre to the calculated point 
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 28 April 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunS0gamma1(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp delta, _Tp z)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();

    if (Heun_cont_coef<_Tp> == _Tp{0})
      HeunOpts();

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if ((std::real(z / a) >= _Tp{1} && std::imag(z / a) == _Tp{0})
     || (std::real(z) >= _Tp{1}   && std::imag(z) == _Tp{0})
     || (std::imag(z) == _Tp{0}   && sign(z) < _Tp{0}))
      {
	wrnmsg = "HeunS0gamma1: z belongs to a possible branch cut; "; 
	return ret;
      }
    else if (std::min(std::abs(z-[1 a])) < eps)
      {
	wrnmsg = "HeunS0gamma1: z is too close to one of the singular points; "; 
	return ret;
      }
    else
      {
	auto R0 = std::min(_Tp{1}, std::abs(a));

	if (std::abs(z) <= R0 * (Heun_cont_coef<_Tp> + 0.01))
	  return HeunS00gamma1(a, q, alpha, beta, delta, z);
	else
	  {
	    if (std::abs(a) < 1)
	      {
		P1 = a;
		P2 = 1;
	      }
	    else
	      {
		P1 = 1;
		P2 = a;
	      }
	    auto R1 = std::min(std::abs(P1), std::abs(a - 1));
	    auto R2 = std::min(std::abs(P2), std::abs(a - 1));

	    std::vector<_Tp> zz;
	    zz.reserve(4);

	    zz.push_back(_Tp{0});

	    auto [d, isinside] = dist(P1, _Tp{0}, z);
	    if (isinside && d < R1 / 2)
	      zz.push_back(P1 + std::exp(1i * (pi<_Tp> / 2 + angle(z)))
			* std::min(R1 / 2, std::abs(z - P1)) * sign(std::imag(z / P1)));

	    [d, isinside] = dist(P2, _Tp{0}, z);
	    if (isinside && d < R2 / 2)
	      zz.push_backP2 + std::exp(1i * (pi<_Tp> / 2 + angle(z)))
			* std::min(R2 / 2, std::abs(z - P2)) * sign(std::imag(z / P2)));

	    zz.push_back(z);

	    auto sum_error = _Tp{0};
	    auto sum_terms = _Tp{0};

	    auto failure = false;

	    _Tp H0, dH0;
	    for (int k = 0; k < std::size(zz) - 1; ++k)
	      {
		auto z0 = zz[k];
		auto theta = angle(zz[k + 1] - z0);
		auto insearch = true;

		while (insearch && !failure)
		  {
		    if (z0 == _Tp{0})
		      R = R0;
		    else
		      R = std::min(std::abs({z0, z0 - 1, z0 - a}));

		    decltype(zz[0]) z1;
		    if (std::abs(zz[k + 1] - z0) <= R * Heun_cont_coef)
		      {
			z1 = zz[k + 1];
			insearch = false;
		      }
		    else
		      z1 = z0 + Heun_cont_coef<_Tp> * R * std::exp(1i * theta);

		    if (z0 == _Tp{0})
		      {
			[H0, dH0, error, num_terms, wrnmsg]
			  = HeunS00gamma1(a, q, alpha, beta, delta, z1);
			if (wrnmsg.length() != 0)
			  {
			    wrnmsg = "HeunS0gamma1: "
				     "problem invoking HeunS00gamma1("
				   << a << ',' << q
				   << ',' << alpha << ',' << beta
				   << ',' << delta << ',' << z1
				   << "); warning: "<< wrnmsg << "; ";
			    failure = true;
			  }
		      }
		    else
		      {
			[H0, dH0, error, num_terms, wrnmsg]
			  = HeunGfromZ0(a, q, alpha, beta, _Tp{1}, delta, z1, z0, H0, dH0);
			if (wrnmsg.length() != 0)
			  {
			    wrnmsg = "HeunS0gamma1: "
				     "problem invoking HeunGfromZ0("
				   << a << ',' << q
				   << ',' << alpha << ',' << beta
				   << ',' << gamma << ',' << delta
				   << ',' << z1 << ',' << z0
				   << ',' << H0 << ',' << dH0
				   << "); warning: " << wrnmsg << "; ";
			    failure = true;
			  }
		      }

		    sum_error += error;
		    sum_terms += num_terms;

		    z0 = z1;
		  }

		if (failure)
		  break;
	      }

	    value = H0;
	    deriv = dH0;
	    num_terms = sum_terms;
	    error = sum_error;

	    return ret;
	  }
      }
  }

// The second local Heun function at z=0
//
// evaluation by analytic continuation from z=0
// to z using a consequence of power expansions
//
// path2z is a list of points from 0 to z
//
// The function uses a parameter Heun_cont_coef which can be changed by HeunOpts: 
// for each power expansion of the analytic continuation procedure Heun_cont_coef
// is the relative distance from the centre to the calculated point 
//
// if wrnmsg is not empty, then the function returns value, deriv = NaN
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 04 June 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunSmv(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp gamma, _Tp delta,
	  const std::vector<_Tp>& path2z)
  {
    using namespace std::complex_literals;
    const auto eps = std::numeric_limits<_Tp>::epsilon();

    if (Heun_cont_coef<_Tp> == _Tp{0})
      HeunOpts();

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (path2z[0] != _Tp{0})
      {
	wrnmsg = "HeunSmv: path2z should start from zero; "; 
	return ret;
      }
    else if (dist2sing(path2z, {_Tp{1}, a}) < _Tp{10} * eps
	  || dist2sing({std::begin(path2z) + 2, std::end(path2z)}, {_Tp{0}}) < _Tp{10} * eps)
      {
	wrnmsg = "HeunSmv: path2z is too close to one of the singular points; "; 
	return ret;
      }
    else
      {
	if (gamma == 1)
	  return HeunSmvgamma1(a, q, alpha, beta, delta, path2z);
	else
	  {
	    auto z = path2z.back();

	    auto epsilon = alpha + beta + 1 - gamma - delta;
	    auto [H0w, dH0w, errorw, num_terms, wrnmsg]
	      = HeunLmv(a, q - (gamma - 1) * (epsilon + a * delta),
			beta - gamma + 1,
			alpha - gamma + 1, 2 - gamma, delta, path2z);

	    auto zp0 = std::pow(path2z[1], -gamma);
	    // increment of argument on the line segment (z0,z1)
	    auto addarg = [](_Tp z0, _Tp z1) -> _Tp
			  { return std::asin((std::real(z0) * std::imag(z1)
					     -std::real(z1) * std::imag(z0))
					   / (std::abs(z0) * std::abs(z1))); };

	    for (int k = 3; k < std::size(path2z); ++k)
	      zp0 *= std::exp(-gamma *(std::log(std::abs(path2z[k]))
			    - std::log(std::abs(path2z[k-1]))
			    + 1i * addarg(path2z[k-1], path2z[k])));

	    value = z * zp0 * H0w;
	    deriv = (1 - gamma) * zp0 * H0w + z * zp0 * dH0w;
	    error = std::abs(z * zp0) * errorw;

	    return ret;
	  }
      }
  }


// the second local Heun function at z=0
// the case when gamma = 1
//
// evaluation by analytic continuation from z=0
// to z using consequence of power expansions
//
// path2z is a list of points from 0 to z
//
// The function uses a parameter Heun_cont_coef which can be changed by HeunOpts: 
// for each power expansion of the analytic continuation procedure Heun_cont_coef
// is the relative distance from the centre to the calculated point 
//
// Oleg V. Motygin, copyright 2015, license: GNU GPL v3
//
// 04 June 2015
//
template<typename _Tp>
  __heun_t<_Tp>
  HeunSmvgamma1(_Tp a, _Tp q, _Tp alpha, _Tp beta, _Tp delta,
		const std::vector<_Tp>& path2z)
  {
    using namespace std::complex_literals;
    const auto eps = std::numeric_limits<_Tp>::epsilon();

    if (Heun_cont_coef<_Tp> == _Tp{0})
      HeunOpts();

    __heun_t<_Tp> ret;
    auto& [value, deriv, error, num_terms, wrnmsg] = ret;

    if (path2z[0] != _Tp{0})
      {
	wrnmsg = "HeunSmvgamma1: path2z should start from zero; "; 
	return ret;
      }
    else if (dist2sing(path2z, {_Tp{1}, a}) < _Tp{10} * eps
	  || dist2sing(path2z(2 : end), {_Tp{0}}) < _Tp{10} * eps)
      {
	wrnmsg = "HeunSmvgamma1: path2z is too close to one of the singular points; "; 
	return ret;
      }
    else
      {
	auto R0 = std::min(_Tp{1}, std::abs(a));
	auto z = path2z.back();

	if (std::abs(z) <= R0 * (Heun_cont_coef<_Tp> + 0.01))
	  return HeunS00gamma1(a, q, alpha, beta, delta, z);
	else
	  {
	    auto sum_error = _Tp{0};
	    auto sum_terms = _Tp{0};

	    auto failure = false;

	    __heun_t<_Tp> heun0;
	    auto& [H0, dH0, error, num_terms, wrnmsg] = heun0;

	    for (int k = 0; k < std::size(path2z) - 1; ++k)
	      {
		auto z0 = path2z[k];
		auto theta = angle(path2z[k+1] - z0);
		auto insearch = true;

		while (insearch && !failure)
		  {
		    _Tp R;
		    if (z0 == _Tp{0})
		      R = R0;
		    else
		      R = std::min({std::abs(z0), std::abs(z0 - 1), std::abs(z0 - a)});

		    _Tp z1;
		    if (std::abs(path2z[k+1] - z0) <= R * Heun_cont_coef)
		      {
			z1 = path2z[k+1];
			insearch = false;
		      }
		    else
		      z1 = z0 + Heun_cont_coef<_Tp> * R * std::exp(1i * theta);

		    if (z0 == _Tp{0})
		      {
			[H0, dH0, error, num_terms, wrnmsg]
			  = HeunS00gamma1(a, q, alpha, beta, delta, z1);

			if (wrnmsg.length() != 0)
			  {
			    wrnmsg = "HeunSmvgamma1: "
				     "problem invoking HeunS00gamma1("
				   << a << ',' << q << ','
				   << alpha << ',' << beta << ',' << delta << ','
				   << z1
				   << "); warning: " << wrnmsg << "; ";
			    failure = true;
			  }
		      }
		    else
		      {
			[H0, dH0, error, num_terms, wrnmsg]
			  = HeunGfromZ0(a, q, alpha, beta, 1, delta, z1, z0, H0, dH0);
			if (wrnmsg.length() != 0)
			  {
			    wrnmsg = "HeunSmvgamma1: "
				     "problem invoking HeunGfromZ0("
				   << a << ',' << q << ','
				   << alpha << ',' << beta << ',' << gamma << ',' << delta << ','
				   << z1 << ',' << z0 << ','
				   << H0 << ',' << dH0
				   << "); warning: " << wrnmsg << "; ";
			    failure = true;
			  }
		      }

		    sum_error += sum_error;
		    sum_terms += sum_terms;

		    z0 = z1;
		  }

		if (failure)
		  break;
	      }
	    value = H0;
	    deriv = dH0;
	    num_terms = sum_terms;
	    error = sum_error;

	    return ret;
	  }
      }  
  }
