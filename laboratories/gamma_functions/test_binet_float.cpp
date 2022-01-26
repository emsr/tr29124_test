/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>

#include <emsr/float128_io.h>
#include <emsr/polynomial.h>
#include <emsr/summation.h>
#include <emsr/fp_type_util.h>
#include <emsr/math_constants.h>
#include <emsr/sf_gamma.h>

namespace emsr
{
namespace detail
{

  /**
   * Computes the sequence of Bernoulli numbers @f$B_{2m}@f$ (m > 0)
   * using the Akiyama-Tanigawa algorithm.
   * This is unstable.
   */
  template<typename _Real>
    std::vector<_Real>
    bernoulli_a_t(std::size_t len)
    {
      auto n = 2 * len + 1;
      std::vector<_Real> t;
      std::vector<_Real> a;

      t.emplace_back(1LL);

      for (std::size_t m = 1; m < n; ++m)
	{
	  t.push_back(_Real{1} / _Real(m + 1));
	  for (int j = m; j > 0; --j)
            t[j - 1] = _Real(j) * (t[j - 1] - t[j]);

	  // Get all Bernoulli numbers by deleting the 'if' clause.
	  if ((m & 1) == 0)
	    a.push_back(t[0]);
	}

      return a;
    }

  /**
   * @see The Gamma function revisited.
   */
  template<typename _Tp>
    std::vector<_Tp>
    recursive_thing(_Tp c)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;

      const int _N = 100;

      std::vector<_Real> F;
      _Tp _Fprev{1}, _Gprev{0}, _Hprev{0};
      F.push_back(_Fprev);
      for (int n = _N; n > 0; --n)
	{
	  auto _Fcurr = _Fprev + _Hprev * c / (2 * n);
	  auto _Gcurr = ((2 * n) * _Gprev + c * _Fcurr) / (2 * n - 1);
	  auto _Hcurr = _Hprev + _Gcurr;
	  _Fprev = _Fcurr;
	  _Gprev = _Gcurr;
	  _Hprev = _Hcurr;
	  F.push_back(_Fprev);
	}
      return F;
    }

  /**
   * @see The Gamma function revisited.
   */
  template<typename _Tp>
    _Tp
    binet_recursive(_Tp z)
    {
      const auto s_2pi = emsr::tau_v<_Tp>;
      const auto c = s_2pi * z;
      for (int k = 1; k < 10; ++k)
	{
          std::vector<_Tp> F = recursive_thing(c * k);
	  //for (auto f : F)
	  //  
	}
    }


  /**
   * Computes the sequence of even Bernoulli numbers @f$ B_{2m} @f$ (m > 0)
   * using the old Riemann zeta function series.
   */
  template<typename _Real>
    std::vector<_Real>
    bernoulli_vec(std::size_t len)
    {
      std::vector<_Real> a;
      for (std::size_t m = 1; m <= len; ++m)
	a.push_back(bernoulli<_Real>(2 * m));
      return a;
    }

  /**
   * Scales the even Bernoulli numbers @f$ B_{2m+2} @f$ with weights
   * @f$ (-1)^m/((2m + 1)(2m + 2)) @f$, m >= 0.
   */
  template<typename _Real>
    std::vector<_Real>
    weights(std::vector<_Real> b)
    {
      int sgn = 1;
      for (std::size_t m = 0; m < b.size(); ++m)
	{
	  b[m] *= _Real(sgn) / _Real(2 * m + 1) / _Real(2 * m + 2);
	  sgn = -sgn;
	}

      return b;
    }

  /**
   * Computes the partial numerators for the Binet function
   * using Rutishauser's Quotient-Difference (QD) algorithm.
   */
  template<typename _Real>
    std::vector<_Real>
    quotient_difference(std::vector<_Real> s)
    {
      auto len = s.size();
      auto zero = _Real{0};
      std::vector<std::vector<_Real>> m;
      std::vector<_Real> r;

      for (std::size_t n = 0; n < len; ++n)
	{
	  m.push_back(std::vector<_Real>());
	  m.back().push_back(zero); // e[k+1,0] = 0, k >= 0
	}
      for (std::size_t n = 0; n < len - 1; ++n)
	m[n].push_back(s[n + 1] / s[n]); // q[k,1] = c[k+1]/c[k], k >= 0

      r.push_back(s[0]);
      r.push_back(m[0][1]);

      for (std::size_t k = 2; k < len; ++k)
	{
	  for (std::size_t n = 0; n < len - k; ++n)
            {
              auto a = m[n + 1][k - 2];
              auto b = m[n + 1][k - 1];
              auto c = m[n][k - 1];
              m[n].push_back((k & 1) == 0
				 ? a + b - c
				 : a * b / c);
            }
	  r.push_back(m[0][k]);
	}

      return r;
    }

  /**
   * For a series specified by coefficients @f$ c_k @f$:
   * @f[
   *  \Lambda(x) = \sum_{k=0}^{N} c_k x^k
   * @f]
   * form the coefficients @f$ d_k @f$ for the inverse
   * @f[
   *  \frac{1}{\Lambda(x)} = \sum_{k=0}^{N} d_k x^k
   * @f]
   */
  template<typename _Real>
    std::vector<_Real>
    series_reciprocal_old(std::vector<_Real> c)
    {
      if (c.size() == 0)
	return std::vector<_Real>{};
      else if (c[0] == _Real{0})
	throw std::domain_error("series_reciprocal: first (constant) coefficient is zero.");
      else if (c.size() == 1)
	return std::vector<_Real>{{_Real{1} / c[0]}};
      else
	{
	  auto n = c.size();
	  auto m = n + 1;
	  std::vector<_Real> lambda(n);
	  lambda[0] = _Real{0};
	  for (unsigned i = 1; i < n; ++i)
	    lambda[i] = c[i] / c[0];

	  std::vector<_Real> lambdak(m);
	  std::vector<_Real> d(m);
	  d[0] = lambdak[0] = _Real{1}; // k == 0
	  for (unsigned k = 1; k < m; ++k)
	    {
	      std::vector<_Real> work(m);
	      for (unsigned i = 1; i < n; ++i)
		for (unsigned j = k - 1; j < m; ++j)
		  if (i + j < m)
		    work[i + j] += lambda[i] * lambdak[j];
	      std::swap(work, lambdak);
	      auto sign = (k % 2 == 0 ? +1 : -1);
	      for (unsigned j = k; j < m; ++j)
		d[j] += sign * lambdak[j];
	    }
	  for (unsigned j = 0; j < m; ++j)
	    d[j] /= c[0];
	  return d;
	}
    }

  /**
   * For a series specified by coefficients @f$ c_k @f$:
   * @f[
   *  \Lambda(x) = \sum_{k=0}^{N} c_k x^k
   * @f]
   * form the coefficients @f$ d_k @f$ for the inverse
   * @f[
   *  \frac{1}{\Lambda(x)} = \sum_{k=0}^{N} d_k x^k
   * @f]
   */
  template<typename _Real>
    std::vector<_Real>
    series_reciprocal(std::vector<_Real> c)
    {
      if (c.size() == 0)
	return std::vector<_Real>{};
      else if (c[0] == _Real{0})
	throw std::domain_error("series_reciprocal: first (constant) coefficient is zero.");
      else if (c.size() == 1)
	return std::vector<_Real>{{_Real{1} / c[0]}};
      else
	{
	  auto n = c.size();
	  auto m = n + 1;
	  std::vector<_Real> d(m);
	  d[0] = _Real{1} / c[0];
	  for (auto k = 1u; k < m; ++k)
	    {
	      for (auto i = 0u; i < k; ++i)
	        if (k - i < n)
		  d[k] += c[k - i] * d[i];
	      d[k] *= -d[0];
	    }
	  return d;
	}
    }

  /**
   * For a series specified by coefficients @f$ c_k @f$:
   * @f[
   *  \Lambda(x) = \sum_{k=0}^{N} c_k x^k
   * @f]
   * form the coefficients @f$ d_k @f$ for the inverse
   * @f[
   *  \frac{1}{\Lambda(x)} = \sum_{k=0}^{N} d_k x^k
   * @f]
   */
  template<typename _Real>
    std::vector<_Real>
    series_reciprocal_vanWijn(std::vector<_Real> c)
    {
std::cerr << std::showpoint;
//      using _WijnSum = emsr::_VanWijngaardenSum<_Real>;
      if (c.size() == 0)
	return std::vector<_Real>{};
      else if (c[0] == _Real{0})
	throw std::domain_error("series_reciprocal_vanWijn: first (constant) coefficient is zero.");
      else if (c.size() == 1)
	return std::vector<_Real>{{_Real{1} / c[0]}};
      else
	{
	  auto n = c.size();
	  auto m = n + 1;
	  std::vector<_Real> d(m);
	  d[0] = _Real{1} / c[0];
	  for (auto k = 1u; k < m; ++k)
	    {
std::cerr << '\n';
//	      _WijnSum sum(8);
emsr::BasicSum<_Real> sum;
	      for (auto i = 0u; i < k; ++i)
	        if (k - i < n)
{
std::cerr << std::setw(12) << c[k - i] * d[i] << '\t' << std::setw(12) << c[k - i] << '\t' << std::setw(12) << d[i] << '\n';
		  sum += c[k - i] * d[i];
}
	      d[k] = -d[0] * sum();
	    }
	  return d;
	}
    }

  /**
   * Computes the partial numerators for the Binet function
   * using Rutishauser's Quotient-Difference (QD) algorithm.
   * Use the more numerically stable progressive version.
   */
  template<typename _Real>
    std::vector<_Real>
    quotient_difference_prog(std::vector<_Real> s)
    {
      auto n = s.size();
      auto zero = _Real{0};
      std::vector<_Real> r;

      std::vector<std::vector<_Real>> q;
      q.push_back(std::vector<_Real>{});
      q.back().push_back(-s[1] / s[0]);
      for (unsigned l = 1; l < n; ++l)
	q.back().push_back(zero);

      std::vector<std::vector<_Real>> e;
      e.push_back(std::vector<_Real>{});
      e.back().push_back(_Real{0});
      for (unsigned l = 1; l < n - 1; ++l)
	e.back().push_back(s[l + 1] / s[l]);

      r.push_back(e[0][0]);

      for (unsigned k = 1; k < n - 2; ++k)
	{
	  q.push_back(std::vector<_Real>{});
	  for (unsigned l = 0; l < n - 2 * k; ++l)
	    q.back().push_back(q[k - 1][l] + e[k - 1][l + 1] - e[k - 1][l]);

	  if (q[k].size() > 0)
	    r.push_back(q[k][0]);

	  e.push_back(std::vector<_Real>{});
	  for (unsigned l = 0; l < n - 2 * k - 1; ++l)
	    e.back().push_back(q[k][l + 1] * e[k - 1][l + 1] / q[k][l]);

	  if (e[k].size() > 0)
	    r.push_back(e[k][0]);
	  else
	    break;
	}

      return r;
    }

  // Computes the Stieltjes continued fraction for the
  // Gamma function using Rutishauser's QD-algorithm.
  template<typename _Real>
    std::vector<_Real>
    stieltjes_cont_frac_seq(int len)
    { return quotient_difference(weights(bernoulli_vec<_Real>(len))); }

  // Computes the Stieltjes continued fraction for the
  // Gamma function using the progressive QD-algorithm.
  template<typename _Real>
    std::vector<_Real>
    stieltjes_cont_frac_seq_prog(int len)
    { return quotient_difference_prog(series_reciprocal(weights(bernoulli_vec<_Real>(len)))); }

  /**
   * Compute the Binet function using the asymptotic series
   */
  template<typename _Tp>
    _Tp
    binet_asymp(_Tp z)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;
      constexpr auto s_eps = std::numeric_limits<_Real>::epsilon();

      // Weighted Bernoulli numbers: (-1)^k B_{2k + 2} / ((2k + 1)(2k + 2)), k >= 0.
      constexpr std::size_t s_n = 50;
      constexpr _Real
      s_b[s_n]
      {
	 0.0833333333333333333333333333333333Q,
	-0.00277777777777777777777777777777778Q,
	 0.000793650793650793650793650793650794Q,
	-0.000595238095238095238095238095238095Q,
	 0.000841750841750841750841750841750842Q,
	-0.00191752691752691752691752691752692Q,
	 0.00641025641025641025641025641025641Q,
	-0.0295506535947712418300653594771242Q,
	 0.179644372368830573164938490015889Q,
	-1.3924322169059011164274322169059Q,
	 13.4028640441683919944789510006901Q,
	-156.848284626002017246255102181262Q,
	 2193.10333333333333242281567137488Q,
	-36108.7712537249893410286881332229Q,
	 691472.268851313066777148250599689Q,
	-15238221.539407416184496902516819Q,
	 382900751.391414141206257412944758Q,
	-10882266035.7843910827594224065177Q,
	 347320283765.002252041501237032466Q,
	-12369602142269.2744463508996924125Q,
	 488788064793079.334748002432392346Q,
	-21320333960919373.8819953763850981Q,
	 1021775296525700076.81475571783662Q,
	-53575472173300203569.7635271795249Q,
	 3061578263704883412598.75652933154Q,
	-189999174263992040345172.024164544Q,
	 12763374033828834138229314.2866605Q,
	-925284717612041629895616935.647741Q,
	 72188225951856102911502977063.9351Q,
	-6045183405995856961951306804314.17Q,
	 542067047157009453982686049507708.0Q,
	-5.19295781531408194139317505902564e+34Q,
	 5.30365885511970059106530722178703e+36Q,
	-5.76332534816496400763640094632658e+38Q,
	 6.65115571484845393008203176382139e+40Q,
	-8.13737835813668052936054479212525e+42Q,
	 1.05369669533571417913038325599530e+45Q,
	-1.44181805999622062443077173359329e+47Q,
	 2.08173565220895654364963845771048e+49Q,
	-3.1670226634886661786956177503001e+51Q,
	 5.07000646121113733654063737462523e+53Q,
	-8.52997282030055187017934268544118e+55Q,
	 1.50641728093405985560080145938733e+58Q,
	-2.78934947038316368320924011776134e+60Q,
	 5.40935043528604149280237369227601e+62Q,
	-1.09753378215085198388932020854243e+65Q,
	 2.32748762026184791385428030490926e+67Q,
	-5.15392916206532138229395451444576e+69Q,
	 1.19062102308902264390529844167985e+72Q,
	-2.86689389602966736504472887935303e+74Q,
      };

      auto z2 = _Real{1} / (z * z);
      auto J = _Val{};
      auto zk = _Real{1} / z;
      for (auto b : s_b)
	{
	  auto term = b * zk;
	  J += term;
	  if (std::abs(term) < s_eps * std::abs(J))
	    break;
	  zk *= z2;
	}

      return J;
    }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    binet_cont_frac(_Tp z)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;

      // Stieltjes partial numerators.
      constexpr std::size_t s_n = 96;
      constexpr _Real
      s_a[s_n]
      {
	0.0833333333333333333333333333333333Q,
	0.0333333333333333333333333333333333Q,
	 0.252380952380952380952380952380952Q,
	 0.525606469002695417789757412398922Q,
	  1.01152306812684171174737212473062Q,
	  1.51747364915328739842849151949548Q,
	  2.26948897420495996090915067220988Q,
	  3.00991738325939817007314073420771Q,
	  4.02688719234390122616887595318146Q,
	  5.00276808075403005168850241227618Q,
	  6.28391137081578218007266315495159Q,
	  7.49591912238403380685495790732127Q,
	  9.04066023436772762998666609001229Q,
	  10.4893036545094782875289814299128Q,
	  12.2971936103862187466923692040639Q,
	  13.9828769539923961784331447385825Q,
	  16.0535514167050146361374915307607Q,
	  17.9766073998701127949281372367462Q,
	  20.3097620274419730230673739635349Q,
	  22.4704716399325574758783364097705Q,
	  25.0658465489469610173532914926089Q,
	  27.4644518250275190652737946955527Q,
	  30.3218212316756007623693149398914Q,
	  32.9585339299691045132634785379822Q,
	  36.0776989313050131199954841899512Q,
	  38.9527066823032420486221769350022Q,
	  42.3334900435887335309177330095052Q,
	  45.4469608500453554059239128674111Q,
	  49.0892031290347752081214826346055Q,
	  52.4412887513857265619233000597589Q,
	  56.3448453453809833753955797537205Q,
	  59.9356839071150333537518705178518Q,
	  64.1004227559858314780894631982923Q,
	  67.9301407879351470488016641155734Q,
	   72.355940555316437684665595839878Q,
	  76.4246546266993907322325843076967Q,
	     81.1114032374092435525653584692Q,
	  85.4192212762136787021289558722571Q,
	  90.3668147241044672462686630639823Q,
	  94.9138370997201052742490661620312Q,
	  100.122178464277395495903162701654Q,
	  104.908498885286170621478245464067Q,
	  110.377497511751573913393590634963Q,
	  115.403203777985798757433371085184Q,
	  121.132774587272518084243872533321Q,
	   126.39794922552925075229845755837Q,
	  132.388012128380142749569888206836Q,
	  137.892732934202376912269757821377Q,
	  144.143212329966869026943650253781Q,
	  149.887552832809305746006710750765Q,
	  156.398377177576973745096638699112Q,
	  162.382407042900426530471938775574Q,
	  169.153508474984932408059609929377Q,
	   175.37729385398857680658699587357Q,
	  182.408607867218065510716673661826Q,
	  188.872211702762410172665369254403Q,
	  196.163676859917004824417285070116Q,
	  202.867159155531461804697073884307Q,
	  210.418716835726303983539348533132Q,
	  217.362134893306299782673360590071Q,
	  225.173729068255076835933399830165Q,
	  232.357137699049610253773779798926Q,
	  240.428714734018620041024938840451Q,
	  247.852166446770354011995062455453Q,
	   256.18367492258298121720294924574Q,
	  263.847220092464803208411818347916Q,
	  272.438610644428059254443668928037Q,
	  280.342297668527701485815071219909Q,
	  289.193522832224603576535146833356Q,
	  297.337398292267149631717466586696Q,
	  306.448412310057104239926692820533Q,
	  314.832521248409818147155868398378Q,
	  324.203279591189949457394418186791Q,
	  332.827666467109303098738023006716Q,
	  342.458123766008487891255835335155Q,
	  351.322837085329931349530858685116Q,
	  361.212936622621040016587632599681Q,
	  370.318052897603932328668688760666Q,
	   380.46767185378546297365914347277Q,
	  389.813420780440596080148545734177Q,
	  400.222083801356338201887356456528Q,
	  409.809503100119305624896129735861Q,
	  420.474887464105522818668573985185Q,
	   430.30921928664565081854293143378Q,
	  441.219490742117007174187514884464Q,
	  451.327279780387873931254713293876Q,
	   462.42352388554253922581949114896Q,
	  472.933427852983052154394997964881Q,
	  483.940623946790019803987275719333Q,
	  495.423974552696641115707823565594Q,
	  505.200988444745669875055339184825Q,
	  519.819449137925693258054033672317Q,
	  524.613967559501453127692011747469Q,
	  547.989836985329038445605446265368Q,
	  542.239386997206019448762884936956Q,
	  570.030818671740310619828040336405Q,
      };

      // Backward recurrence.
      auto w = _Val{}; // The tail function.
      auto J = w;
      for (std::ptrdiff_t k = s_n - 1; k >= 0; --k)
        J = s_a[k] / (z + J);

      return J;
    }

  /**
   * Use the recurrance formula for the Binet function to evaluate the function
   * where asymptotic evaluations work well.
   *
   * @f[
   *    J(z) = J(z+1) + z + \frac{1}{2}ln(1+\frac{1}{z}) - 1
   * @f]
   */
  template<typename _Tp>
    _Tp
    binet_recur(_Tp z)
    {
      _Tp stuff{};
      while (std::real(z) < _Tp{10})
        {
	  stuff += z + _Tp{0.5} * std::log1p(_Tp{1} / z) - _Tp{1};
	  z += _Tp{1};
        }
      return stuff + binet_cont_frac(z);
    }

  /**
   * Compute the Binet function @f$ J(z) @f$ defined by
   * @f[
   *    J(z) = log\left(\Gamma(z)\right) + z
   *         - \left(z-\frac{1}{2}\right) log(z) - \frac{1}{2}log(2\pi)
   * @f]
   * or
   * @f[
   *    \Gamma(z) = \sqrt{2\pi}z^{z-\frac{1}{2}}e^{-z}e^{J(z)}
   * @f]
   * where @f$ \Gamma(z) @f$ is the gamma function.
   */
  template<typename _Tp>
    _Tp
    binet(_Tp z)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;

      /// @todo Find Binet function switch.
      constexpr auto s_switchover = _Real{10};

      if (std::isnan(z))
	return emsr::quiet_NaN<_Tp>();
      else if (z < s_switchover)
	return binet_cont_frac(z);
      else
	return binet_asymp(z);
    }

} // namespace detail
} // namespace emsr


namespace emsr
{

  /**
   * Return the Binet function @f$ J(z) @f$ or the scaled log gamma function
   * @f$ \Gamma^*(z) @f$ for @c float argument @f$ z @f$.
   *
   * @see lgamma_scaled for details.
  float
  lgamma_scaledf(float z)
  { return emsr::detail::binet<float>(z); }
   */

  /**
   * Return the Binet function @f$ J(z) @f$ or the scaled log gamma function
   * @f$ \Gamma^*(z) @f$ for <tt>long double</tt> argument @f$ z @f$.
   *
   * @see lgamma_scaled for details.
   */
  long double
  lgamma_scaledl(long double z)
  { return emsr::detail::binet<long double>(z); }

  /**
   * Return the Binet function @f$ J(z) @f$ or the scaled log gamma function
   * @f$ log(\Gamma^*(z)) @f$ defined by
   * @f[
   *    J(z) = log(\Gamma^*(z)) = log\left(\Gamma(z)\right) + z
   *         - \left(z-\frac{1}{2}\right) log(z) - log(2\pi)
   * @f]
   * or
   * @f[
   *    \Gamma(z) = \sqrt{2\pi}z^{z-\frac{1}{2}}e^{-z}e^{J(z)}
   * @f]
   * where @f$ \Gamma(z) @f$ is the gamma function.
   */
  template<typename _Tp>
    _Tp
    lgamma_scaled(_Tp z)
    { return emsr::detail::binet(z); }


  /**
   * Return the exponential of the Binet function or the scaled gamma function
   * @f$ \Gamma^*(z) @f$ for @c float argument @f$ z @f$.
   *
   * @see tgamma_scaled for details.
  float
  tgamma_scaledf(float z)
  { return std::exp(emsr::detail::binet<float>(z)); }
   */

  /**
   * Return the exponential of the Binet function or the scaled gamma function
   * @f$ \Gamma^*(z) @f$ for <tt>long double</tt> argument @f$ z @f$.
   *
   * @see tgamma_scaled for details.
   */
  long double
  tgamma_scaledl(long double z)
  { return std::exp(emsr::detail::binet<long double>(z)); }

  /**
   * Return the exponential of the Binet function or the scaled gamma function
   * @f$ \Gamma^*(z) @f$ defined by
   * @f[
   *    \Gamma^*(z) = \Gamma(z)/\sqrt{2\pi}z^{z-\frac{1}{2}}e^{-z} = e^{J(z)}
   * @f]
   * or
   * @f[
   *    \Gamma(z) = \sqrt{2\pi}z^{z-\frac{1}{2}}e^{-z}\Gamma^*(z)
   * @f]
   * where @f$ \Gamma(z) @f$ is the gamma function.
   */
  template<typename _Tp>
    _Tp
    tgamma_scaled(_Tp z)
    { return std::exp(emsr::detail::binet(z)); }

} // namespace emsr


template<typename _Tp>
  void
  test()
  {
    using _Real = _Tp;
    constexpr auto s_ln2pi = emsr::ln2_v<_Real> + emsr::lnpi_v<_Real>;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    auto width = std::cout.precision() + 8;

    std::cout << "\nBernoulli numbers\n";
    auto bern = emsr::detail::bernoulli_vec<_Real>(50);
    for (auto& b : bern)
      std::cout << ' ' << std::setw(width) << b << '\n';

    std::cout << "\nBernoulli numbers\n";
    auto bern_a_t = emsr::detail::bernoulli_a_t<_Real>(bern.size());
    std::cout << ' ' << std::setw(4) << "k"
	      << ' ' << std::setw(width) << "B"
	      << ' ' << std::setw(width) << "B_AT"
	      << '\n';
    for (std::size_t k = 0; k < bern.size(); ++k)
      {
        std::cout << ' ' << std::setw(4) << k
		  << ' ' << std::setw(width) << bern[k]
		  << ' ' << std::setw(width) << bern_a_t[k]
		  << '\n';
      }

    std::cout << "\nWeighted Bernoulli numbers\n";
    auto wts = emsr::detail::weights(bern);
    for (auto& w : wts)
      std::cout << ' ' << std::setw(width) << w << '\n';

    std::cout << "\nStieltjes partial numerators\n";
    int len = 100;
    auto cf = emsr::detail::stieltjes_cont_frac_seq<_Real>(len);
    for (std::size_t k = 0; k < cf.size(); ++k)
      std::cout << ' ' << std::setw(2) << k + 1 << ": "
		<< ' ' << std::setw(width) << cf[k]
		<< ' ' << std::setw(width) << cf[k] / ((k+1) * (k+1) / _Tp{16}) << '\n';

    std::cout << "\nStieltjes partial numerators (progressive)\n";
    auto cfp = emsr::detail::stieltjes_cont_frac_seq_prog<_Real>(len);
    for (std::size_t k = 0; k < cfp.size(); ++k)
      std::cout << ' ' << std::setw(2) << k + 1 << ": "
		<< ' ' << std::setw(width) << cfp[k]
		<< ' ' << std::setw(width) << cfp[k] / ((k+1) * (k+1) / _Tp{16}) << '\n';

    std::cout << "\nBinet asymptotic\n";
    std::cout << ' ' << std::setw(4) << "x"
	      << ' ' << std::setw(width) << "J_as"
	      << ' ' << std::setw(width) << "J_cf"
	      << ' ' << std::setw(width) << "J_fake"
	      << ' ' << std::setw(width) << "(J_as - J_cf) / J_cf"
	      << ' ' << std::setw(width) << "(J_as - J_f) / J_f"
	      << ' ' << std::setw(width) << "(J_cf - J_f) / J_f"
	      << ' ' << std::setw(width) << "(J_bern - J_f) / J_f"
	      << '\n';
    const auto half = _Real{1} / _Real{2};
    const auto del = _Real{1} / _Real{10};
    for (int k = 1; k <= 5000; ++k)
      {
	auto x = del * k;
	auto j_as = emsr::detail::binet_asymp(x);
	auto j_cf = emsr::detail::binet_cont_frac(x);
	auto j_fake = std::lgamma(x)
		    - (x - half) * std::log(x) + x - s_ln2pi / _Real{2};
	auto j_bern = emsr::detail::log_gamma_bernoulli(x)
		    - (x - half) * std::log(x) + x - s_ln2pi / _Real{2};
	std::cout << ' ' << std::setw(4) << x
		  << ' ' << std::setw(width) << j_as
		  << ' ' << std::setw(width) << j_cf
		  << ' ' << std::setw(width) << j_fake
		  << ' ' << std::setw(width) << (j_as - j_cf) / j_cf
		  << ' ' << std::setw(width) << (j_as - j_fake) / j_fake
		  << ' ' << std::setw(width) << (j_cf - j_fake) / j_fake
		  << ' ' << std::setw(width) << (j_bern - j_fake) / j_fake
		  << '\n';
      }
  }

// Test polynomial reciprocal...
template<typename _Tp>
  void
  test_exp()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto width = std::cout.precision() + 8;

    std::vector<_Tp> coeff;
    _Tp fact = 1;
    coeff.push_back(_Tp{1} / fact);
    for (int i = 1; i <= width; ++i) // Width as number of terms... Why not.
      coeff.push_back(_Tp{1} / (fact *= _Tp(i)));

    std::cout << "\n exp(x) oefficients:\n";
    for (auto cf : coeff)
      std::cout << std::setw(width) << cf << '\n';

    std::vector<_Tp> recip = emsr::detail::series_reciprocal(coeff);
    std::cout << "\n Reciprocal (hopefully exp(-x)) coefficients:\n";
    for (auto cf : recip)
      std::cout << std::setw(width) << cf << '\n';

/* Failed experiment
    std::cout << "\n Reciprocal (hopefully exp(-x)) coefficients using vanWijngaarden:\n";
    std::vector<_Tp> recip_vW = emsr::detail::series_reciprocal_vanWijn(coeff);
    for (auto k = 0u; k < recip_vW.size(); ++k)
      {
	std::cout << std::setw(width) << recip[k]
		  << std::setw(width) << recip_vW[k]
		  << std::setw(width) << recip_vW[k] - recip[k]
		  << '\n';
      }
*/

/* Try to build Pade approximants someday.
    std::cout << "\nTest exp(x)\n";
    emsr::Polynomial<_Tp> expoly(std::begin(coeff), std::end(coeff));
    emsr::Polynomial<_Tp> rat_numer(std::begin(coeff), std::begin(coeff) + 10);
    emsr::Polynomial<_Tp> rat_denom(std::begin(recip), std::begin(recip) + 10);
    const auto del = _Tp{1} / _Tp{10};
    for (int k = 0; k <= 500; ++k)
      {
	auto x = del * k;
	auto pexp = expoly(x);
	auto rexp = std::sqrt(rat_numer(x) / rat_denom(x));
	std::cout << ' ' << std::setw(4) << x
		  << ' ' << std::setw(width) << pexp
		  << ' ' << std::setw(width) << pexp - std::exp(x)
		  << ' ' << std::setw(width) << rexp
		  << ' ' << std::setw(width) << rexp - std::exp(x)
		  << '\n';
      }
*/
  }

int
main()
{
  std::cout << "\nTest polynomial inversion\n";
  test_exp<long double>();

  //std::cout << "\nTest polynomial inversion\n";
  //test_exp<__float128>();

  //std::cout << "\nfloat\n=====\n\n";
  //test<float>();

  std::cout << "\ndouble\n======\n";
  test<double>();

  std::cout << "\nlong double\n===========\n";
  test<long double>();

  //std::cout << "\nfloat128\n===========\n";
  //test<__float128>();
}
