/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -I../../cxx_summation/include -I../../quadrature/include -o plot_airy plot_airy.cpp -L$HOME/bin/lib64 -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./plot_airy ../plot_data > ../output/plot_airy.txt
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <string>

/**
 * 
template<typename _Tp>
  void
  plot_airy(std::string filename)
  {
    using _Val = std::__detail::__num_traits_t<_Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Val>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(w) << "t"
	 << std::setw(w) << "Ai"
	 << std::setw(w) << "Aip"
	 << std::setw(w) << "Bi"
	 << std::setw(w) << "Bip"
	 << std::setw(w) << "Wronskian"
	 << '\n';
    data << "#"
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << '\n';
    for (int i = -2000; i <= +500; ++i)
      {
	auto t = _Tp(0.01Q * i);
	auto airy0 = __gnu_cxx::airy(t);
	data << std::setw(w) << std::real(airy0.z)
	     << std::setw(w) << std::real(airy0.Ai)
	     << std::setw(w) << std::real(airy0.Aip)
	     << std::setw(w) << std::real(airy0.Bi)
	     << std::setw(w) << std::real(airy0.Bip)
	     << std::setw(w) << std::real(airy0.Wronskian())
	     << std::setw(w) << std::real(airy0.true_Wronskian())
	     << '\n';
      }
    data << "\n\n";
  }
 */

/**
 * 
 */
template<typename _Tp>
  void
  splot_airy(std::string filename)
  {
    using _Val = std::__detail::__num_traits_t<_Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Val>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = __gnu_cxx::airy_ai(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::pow(std::abs(Ai), _Val{1} / _Val{6})
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = __gnu_cxx::airy_ai(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::real(Ai)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = __gnu_cxx::airy_ai(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::imag(Ai)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = __gnu_cxx::airy_ai(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::arg(Ai)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Bi = __gnu_cxx::airy_bi(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::pow(std::abs(Bi), _Val{1} / _Val{6})
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Bi = __gnu_cxx::airy_bi(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::real(Bi)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Bi = __gnu_cxx::airy_bi(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::imag(Bi)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Bi = __gnu_cxx::airy_bi(t);
	    data << std::setw(w) << std::real(t)
	         << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::arg(Bi)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';
/*
    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto airy0 = airy(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::abs(airy0.Wronskian())
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';
*/
  }

/**
 * 
template<typename _Tp>
  void
  plot_scorer(std::string filename)
  {
    using _Val = std::__detail::__num_traits_t<_Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Val>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    _Scorer<_Tp> scorer;

    data << "\n\n";
    data << "#"
	 << std::setw(w) << "t"
	 << std::setw(w) << "Gi"
	 << std::setw(w) << "Gip"
	 << std::setw(w) << "Hi"
	 << std::setw(w) << "Hip"
	 << std::setw(w) << "Wronskian"
	 << '\n';
    data << "#"
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << '\n';
    for (int i = -2000; i <= +500; ++i)
      {
	auto t = _Tp(0.01Q * i);
	auto scorer0 = scorer(t);
	data << std::setw(w) << std::real(scorer0.z)
	     << std::setw(w) << std::real(scorer0.Ai)
	     << std::setw(w) << std::real(scorer0.Aip)
	     << std::setw(w) << std::real(scorer0.Bi)
	     << std::setw(w) << std::real(scorer0.Bip)
	     << std::setw(w) << std::real(scorer0.Wronskian())
	     << std::setw(w) << std::real(scorer0.true_Wronskian())
	     << '\n';
      }
    data << "\n\n";
  }
 */

/**
 * 
template<typename _Tp>
  void
  plot_fgh(std::string filename)
  {
    using _Val = std::__detail::__num_traits_t<_Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Val>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(w) << "t"
	 << std::setw(w) << "fai"
	 << std::setw(w) << "fai'"
	 << std::setw(w) << "gai"
	 << std::setw(w) << "gai'"
	 << std::setw(w) << "hai"
	 << std::setw(w) << "hai'"
	 << '\n';
    data << "#"
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << std::setw(w) << "========="
	 << '\n';
    for (int i = -2000; i <= +500; ++i)
      {
	auto t = _Tp(0.01Q * i);
	auto fgh0 = _Airy_series<_Val>::_S_FGH(t);
	data << std::setw(w) << std::real(fgh0.z)
	     << std::setw(w) << std::real(fgh0.fai)
	     << std::setw(w) << std::real(fgh0.faip)
	     << std::setw(w) << std::real(fgh0.gai)
	     << std::setw(w) << std::real(fgh0.gaip)
	     << std::setw(w) << std::real(fgh0.hai)
	     << std::setw(w) << std::real(fgh0.haip)
	     << '\n';
      }
    data << "\n\n";
  }
 */

int
main(int n_app_args, char** arg)
{
  std::string plot_data_dir = ".";
  if (n_app_args > 1)
    plot_data_dir = arg[1];

  using fcmplx = std::complex<float>;
  using cmplx = std::complex<double>;
  using lcmplx = std::complex<long double>;

  //plot_airy<fcmplx>(plot_data_dir + '/' + "plot_airy_float.txt");
  //plot_airy<cmplx>(plot_data_dir + '/' + "plot_airy_double.txt");
  //plot_airy<lcmplx>(plot_data_dir + '/' + "plot_airy_long_double.txt");

  splot_airy<fcmplx>(plot_data_dir + '/' + "plot_airy_complex_float.txt");
  splot_airy<cmplx>(plot_data_dir + '/' + "plot_airy_complex_double.txt");
  splot_airy<lcmplx>(plot_data_dir + '/' + "plot_airy_complex_long_double.txt");

  //plot_scorer<cmplx>(plot_data_dir + '/' + "plot_scorer_double.txt");
  //plot_fgh<cmplx>(plot_data_dir + '/' + "plot_fgh_double.txt");
}
