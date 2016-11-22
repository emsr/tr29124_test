//Provides test example of using integration suite with legendre polynomials
//Written by Jason Dick

#include <iostream>
#include <cmath>
#include <tr1/cmath>
#include <functional>
#include <stdexcept>
#include <sstream>
#include <string>

//#include "simple_integrate.h"
#include "factorial_table.h"
#include "integration.h"

using namespace __gnu_test;

//Function which should integrate to 1 for n1=n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_legendre(int l1, int l2, _Tp x)
  {
    return (_Tp(l1 + l2 + 1) / _Tp{2})
	 * std::tr1::legendre(l2,x) * std::tr1::legendre(l1,x);
  }

template<typename _Tp>
  void
  test_legendre()
  {
    // Neverending loop: runs until integration fails
    for (int n1 = 0; ; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    std::function<_Tp(_Tp)> func(std::bind(&normalized_legendre<_Tp>, n1, n2,
                                         std::placeholders::_1));
	    _Tp integ_precision = _Tp{1000} * std::numeric_limits<_Tp>::epsilon();
	    _Tp comp_precision = _Tp{10} * integ_precision;
	    _Tp integration_result, integration_error;

	    typedef std::pair<_Tp&,_Tp&> ret_type;
	    ret_type{integration_result, integration_error}
        	= integrate(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});

	    if (n1 == n2)
	      {
        	//integration_result should be 1, throw error if differs from one by
        	//more than integration precision (with an additional fudge factor in
        	//case integration isn't quite that accurate)
        	if (fabs(_Tp{1} - integration_result) > comp_precision)
        	  {
        	    std::stringstream ss;
        	    ss.precision(-int(log10(std::numeric_limits<_Tp>::epsilon())));
        	    ss << "Integration failed at n1=" << n1 << ", n2=" << n2
        	       << ", returning result " << integration_result
        	       << " instead of the expected 1" << std::endl;
        	    throw std::logic_error(ss.str());
        	  }
	      }
	    else
	      {
        	//integration_result should be 0, throw error if differs from zero by
        	//more than integration precision (with an additional fudge factor in
        	//case integration isn't quite that accurate)
        	if (fabs(integration_result) > comp_precision)
        	  {
        	    std::stringstream ss;
        	    ss.precision(14);
        	    ss << "Integration failed at n1=" << n1 << ", n2=" << n2
        	       << ", returning result " << integration_result
        	       << " instead of the expected 0" << std::endl;
        	    throw std::logic_error(ss.str());
        	  }
	      }
	  }
	std::cout << "Integration successful for legendre polynomials up to n=" << n1
             << std::endl;
      }
  }

int
main()
{
  try
    {
      test_legendre<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  try
    {
      test_legendre<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
