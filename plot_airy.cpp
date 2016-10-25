/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -I. -g -Wall -Wextra -o plot_airy plot_airy.cpp -L$HOME/bin/lib64 -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./plot_airy > plot_airy.new

$HOME/bin_tr29124/bin/g++ -std=gnu++17  -I.-DOLD -g -Wall -Wextra -o plot_airy_old plot_airy.cpp -L$HOME/bin/lib64 -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./plot_airy_old > plot_airy.new.old

$HOME/bin_tr29124/bin/g++ -std=gnu++17 -I. -UOLD -g -Wall -Wextra -o plot_airy_new plot_airy.cpp -L$HOME/bin/lib64 -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./plot_airy_new > plot_airy.new.new
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
    auto width = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "t"
	 << std::setw(width) << "Ai"
	 << std::setw(width) << "Aip"
	 << std::setw(width) << "Bi"
	 << std::setw(width) << "Bip"
	 << std::setw(width) << "Wronskian"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';
    for (int i = -2000; i <= +500; ++i)
      {
	auto t = _Tp(0.01Q * i);
	auto airy0 = __gnu_cxx::airy(t);
	data << std::setw(width) << std::real(airy0.z)
	     << std::setw(width) << std::real(airy0.Ai)
	     << std::setw(width) << std::real(airy0.Aip)
	     << std::setw(width) << std::real(airy0.Bi)
	     << std::setw(width) << std::real(airy0.Bip)
	     << std::setw(width) << std::real(airy0.Wronskian())
	     << std::setw(width) << std::real(airy0.true_Wronskian())
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
    auto width = 8 + data.precision();

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = __gnu_cxx::airy_ai(t);
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::pow(std::abs(Ai), _Val{1} / _Val{6})
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
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::real(Ai)
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
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::imag(Ai)
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
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::arg(Ai)
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
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::pow(std::abs(Bi), _Val{1} / _Val{6})
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
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::real(Bi)
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
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::imag(Bi)
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
	    data << std::setw(width) << std::real(t)
	         << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::arg(Bi)
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
	    data << std::setw(width) << std::real(t)
		 << std::setw(width) << std::imag(t)
		 << std::setw(width) << std::abs(airy0.Wronskian())
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
    auto width = 8 + data.precision();

    _Scorer<_Tp> scorer;

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "t"
	 << std::setw(width) << "Gi"
	 << std::setw(width) << "Gip"
	 << std::setw(width) << "Hi"
	 << std::setw(width) << "Hip"
	 << std::setw(width) << "Wronskian"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';
    for (int i = -2000; i <= +500; ++i)
      {
	auto t = _Tp(0.01Q * i);
	auto scorer0 = scorer(t);
	data << std::setw(width) << std::real(scorer0.z)
	     << std::setw(width) << std::real(scorer0.Ai)
	     << std::setw(width) << std::real(scorer0.Aip)
	     << std::setw(width) << std::real(scorer0.Bi)
	     << std::setw(width) << std::real(scorer0.Bip)
	     << std::setw(width) << std::real(scorer0.Wronskian())
	     << std::setw(width) << std::real(scorer0.true_Wronskian())
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
    auto width = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "t"
	 << std::setw(width) << "fai"
	 << std::setw(width) << "fai'"
	 << std::setw(width) << "gai"
	 << std::setw(width) << "gai'"
	 << std::setw(width) << "hai"
	 << std::setw(width) << "hai'"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';
    for (int i = -2000; i <= +500; ++i)
      {
	auto t = _Tp(0.01Q * i);
	auto fgh0 = _Airy_series<_Val>::_S_FGH(t);
	data << std::setw(width) << std::real(fgh0.z)
	     << std::setw(width) << std::real(fgh0.fai)
	     << std::setw(width) << std::real(fgh0.faip)
	     << std::setw(width) << std::real(fgh0.gai)
	     << std::setw(width) << std::real(fgh0.gaip)
	     << std::setw(width) << std::real(fgh0.hai)
	     << std::setw(width) << std::real(fgh0.haip)
	     << '\n';
      }
    data << "\n\n";
  }
 */

int
main()
{
  using fcmplx = std::complex<float>;
  using cmplx = std::complex<double>;
  using lcmplx = std::complex<long double>;

  //plot_airy<fcmplx>("plot/airy_float_new.txt");
  //plot_airy<cmplx>("plot/airy_double_new.txt");
  //plot_airy<lcmplx>("plot/airy_long_double_new.txt");

  splot_airy<fcmplx>("plot/airy_complex_float_new.txt");
  splot_airy<cmplx>("plot/airy_complex_double_new.txt");
  splot_airy<lcmplx>("plot/airy_complex_long_double_new.txt");

  //plot_scorer<cmplx>("plot/scorer_double_new.txt");
  //plot_fgh<cmplx>("plot/fgh_double_new.txt");
}
