/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

#include <zeta/riemann_zeta.tcc>
#include <zeta/PolyLog.h>

int
main()
{
  std::cout.precision(emsr::digits10(1.0));
  std::cout << std::scientific;
  const auto w = 8 + std::cout.precision();

  //const auto s_pi = emsr::pi;
  const auto s_2pi = emsr::tau;

  auto li_ok = emsr::detail::polylog(2.0, std::complex<double>(-12.1800));
  auto li_bad = emsr::detail::polylog(2.0, std::complex<double>(-12.1850));

  //  Start with asymp. DONE
  //double s = -12.9;
  //double x = 0.0;
  //emsr::detail::polylog_exp_asymp(s, 100.0 * std::polar(1.0, s_2pi * x));
  //PolyLog_Exp_asym(s, 100.0 * std::polar(1.0, s_2pi * x));

  double rs = -4.0;
  std::cout << '\n' << '\n';
  for (int i = 0; i <= 360; ++i)
    {
      const auto theta = s_2pi * i / 360;
      const auto phase = std::polar(1.0, theta);
      const auto li_gnu = emsr::detail::polylog_exp(rs, phase);
      const auto li_zeta = PolyLog_Exp(rs, phase);
      const auto delta = li_gnu - li_zeta;
      const auto absdelta = std::abs(delta);
      std::cout << ' ' << std::setw(w) << theta
		//<< ' ' << phase
		<< ' ' << std::setw(w) << std::real(li_gnu)
		<< ' ' << std::setw(w) << std::imag(li_gnu)
		<< ' ' << std::setw(w) << std::real(li_zeta)
		<< ' ' << std::setw(w) << std::imag(li_zeta)
		<< ' ' << std::setw(w) << absdelta
		<< ' ' << std::setw(w) << std::real(delta)
		<< ' ' << std::setw(w) << std::imag(delta)
		<< '\n';
    }

  // Zoom in on the area around pi/2
  // I'm right!  My series is continuous across pi/2.
  // I bet the arctan in zeta has a wrong sign or something.
  std::cout << '\n' << '\n';
  for (int i = -20; i <= 20; ++i)
    {
      const auto theta = s_2pi / 4 + s_2pi * i / 3600;
      const auto phase = std::polar(1.0, theta);
      const auto li_gnu = emsr::detail::polylog_exp(rs, phase);
      const auto li_zeta = PolyLog_Exp(rs, phase);
      const auto delta = li_gnu - li_zeta;
      const auto absdelta = std::abs(delta);
      std::cout << ' ' << std::setw(w) << theta
		//<< ' ' << phase
		<< ' ' << std::setw(w) << std::real(li_gnu)
		<< ' ' << std::setw(w) << std::imag(li_gnu)
		<< ' ' << std::setw(w) << std::real(li_zeta)
		<< ' ' << std::setw(w) << std::imag(li_zeta)
		<< ' ' << std::setw(w) << 100 * absdelta
		<< ' ' << std::setw(w) << 100 * std::real(delta)
		<< ' ' << std::setw(w) << 100 * std::imag(delta)
		<< '\n';
    }

  // Zoom in on the area around -12.2 to -12.1
  // My series is discontinuous!
  std::cout << '\n' << '\n';
  for (int i = 0; i <= 100; ++i)
    {
      const auto x = double(-12.2L + i * 0.001L);
      const auto li_gnu = emsr::detail::polylog(2.0, std::complex<double>(x));
      const auto dilog = emsr::detail::dilog(x);
      const auto li_zeta = PolyLog(2.0, std::complex<double>(x));
      std::cout << ' ' << std::setw(w) << x
		<< ' ' << std::setw(w) << std::real(li_gnu)
		<< ' ' << std::setw(w) << std::imag(li_gnu)
		<< ' ' << std::setw(w) << dilog
		<< ' ' << std::setw(w) << std::real(li_zeta)
		<< ' ' << std::setw(w) << std::imag(li_zeta)
		<< ' ' << std::setw(w) << -1000 * std::real(li_gnu - dilog)
		<< ' ' << std::setw(w) << -1000 * std::real(li_gnu - li_zeta)
		<< '\n';
    }


  // OK, now on to Test series 1 [PolyLog_Exp_pos(k, exp(i2pix)]:
  //int s = 3;
  //double x = 0.05;
  //emsr::detail::polylog_exp_pos(s, std::polar(1.0, s_2pi * x));
  //PolyLog_Exp_pos(s, std::polar(1.0, s_2pi * x));

  // Same, [PolyLog_Exp_pos(k, exp(i2pix)] but we need within pi/3 of the negative axis it seems:
  int s = 4;
  double x = 0.60;
  emsr::detail::polylog_exp_pos(s, std::polar(1.0, s_2pi * x));
  PolyLog_Exp_pos(s, std::polar(1.0, s_2pi * x));


  //int m = 2;
  //double = 3.1;
  //emsr::detail::clausen_cl<type>(m, x);
  //Claussen_Cl(m, x);


  // OK, now on to Test series 6 [PolyLog_Exp_pos(s, exp(i2pix)]: DONE
  //s = 0.145;
  //x = 0.05;
  //emsr::detail::polylog_exp_pos(s, std::polar(1.0, s_2pi * x));
  //PolyLog_Exp_pos(s, std::polar(1.0, s_2pi * x));
  // Guess what? my sincos_pi (and actually builtin_sincos called sincos...
  // ... which was picked up by cephes/cmath/sincos.c.
  // Which assumes the angle is in degrees. LOL!

  // OK, now on to Test series 6 [PolyLog_Exp_pos(s, exp(i2pix)]:
  //int s = 3;
  //double x = 0.05;
  // Complex arg seems OK.. But why issue with Clausen?
  //emsr::detail::polylog_exp_pos(s, std::polar(1.0, s_2pi * x));
  //PolyLog_Exp_pos(s, std::polar(1.0, s_2pi * x));

  // Try float, complex...
  //double s = 2.1;
  //double x = 3.1;
  //emsr::detail::polylog_exp_pos(s, std::polar(1.0, s_2pi * x));
  //PolyLog_Exp_pos(s, std::polar(1.0, s_2pi * x));
  //Test function 3 [PolyLog_Exp_neg(s<0, exp(i2pik)]:
  //s = -7.9;
  //emsr::detail::polylog_exp_neg(s, std::polar(1.0, s_2pi * x));
  //PolyLog_Exp_neg(s, std::polar(1.0, s_2pi * x));
  // Is riemann_seta_m_1 broken for s < 0?

}
