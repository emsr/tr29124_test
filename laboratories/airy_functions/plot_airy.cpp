/**
 *
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <string>

#include <emsr/sf_airy.h>

/**
 * 
template<typename Tp>
  void
  plot_airy(std::string filename)
  {
    using Val = emsr::num_traits_t<Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<Val>::digits10);
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
	auto t = Tp(0.01Q * i);
	auto airy0 = emsr::airy(t);
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
template<typename Tp>
  void
  splot_airy(std::string filename)
  {
    using Val = emsr::num_traits_t<Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<Val>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = emsr::airy_ai(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::pow(std::abs(Ai), Val{1} / Val{6})
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = emsr::airy_ai(t);
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
	    auto t = Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = emsr::airy_ai(t);
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
	    auto t = Tp(0.10Q * i, 0.10Q * j);
	    auto Ai = emsr::airy_ai(t);
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
	    auto t = Tp(0.10Q * i, 0.10Q * j);
	    auto Bi = emsr::airy_bi(t);
	    data << std::setw(w) << std::real(t)
		 << std::setw(w) << std::imag(t)
		 << std::setw(w) << std::pow(std::abs(Bi), Val{1} / Val{6})
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = Tp(0.10Q * i, 0.10Q * j);
	    auto Bi = emsr::airy_bi(t);
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
	    auto t = Tp(0.10Q * i, 0.10Q * j);
	    auto Bi = emsr::airy_bi(t);
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
	    auto t = Tp(0.10Q * i, 0.10Q * j);
	    auto Bi = emsr::airy_bi(t);
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
	    auto t = Tp(0.10Q * i, 0.10Q * j);
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
template<typename Tp>
  void
  plot_scorer(std::string filename)
  {
    using Val = emsr::num_traits_t<Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<Val>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    Scorer<Tp> scorer;

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
	auto t = Tp(0.01Q * i);
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
template<typename Tp>
  void
  plot_fgh(std::string filename)
  {
    using Val = emsr::num_traits_t<Tp>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<Val>::digits10);
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
	auto t = Tp(0.01Q * i);
	auto fgh0 = Airy_series<Val>::s_FGH(t);
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
