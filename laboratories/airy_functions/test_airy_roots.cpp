/**
 *
 */

#include <limits>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

template<typename _Tp>
  void
  test_airy_roots()
  {
    const auto s_eps = _Tp{10} * std::numeric_limits<_Tp>::epsilon();
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    const auto w = 6 + std::cout.precision();
    const auto max_iter = 10000;

    // Roots of Ai(-x)
    std::vector<_Tp> zai
    {
       2.33811f,  4.08795f,  5.52056f,  6.78671f,  7.94413f,
       9.02265f, 10.04017f, 11.00852f, 11.93602f, 12.82878f,
      13.69149f, 14.52783f, 15.34076f, 16.13269f, 16.90563f,
      17.66130f, 18.40113f, 19.12638f, 19.83813f, 20.53733f,
      21.22483f, 21.90137f, 22.56761f, 23.22417f, 23.87156f,
      24.51030f, 25.14082f, 25.76353f, 26.37881f, 26.98699f,
      27.58839f, 28.18331f, 28.77201f, 29.35475f, 29.93176f,
      30.50327f, 31.06947f, 31.63056f, 32.18671f, 32.73810f,
      33.28488f, 33.82721f, 34.36523f, 34.89907f, 35.42886f,
      35.95471f, 36.47675f, 36.99507f, 37.50980f, 38.02101f
    };

    // Roots of Ai'(-x)
    std::vector<_Tp> zaip
    {
       1.01879f,  3.24820f,  4.82010f,  6.16331f,  7.37218f,
       8.48849f,  9.53545f, 10.52766f, 11.47506f, 12.38479f,
      13.26222f, 14.11150f, 14.93594f, 15.73820f, 16.52050f,
      17.28470f, 18.03234f, 18.76480f, 19.48322f, 20.18863f,
      20.88192f, 21.56389f, 22.23523f, 22.89659f, 23.54853f,
      24.19156f, 24.82616f, 25.45274f, 26.07171f, 26.68341f,
      27.28818f, 27.88632f, 28.47811f, 29.06381f, 29.64367f,
      30.21792f, 30.78676f, 31.35039f, 31.90899f, 32.46275f,
      33.01183f, 33.55638f, 34.09654f, 34.63246f, 35.16426f,
      35.69207f, 36.21601f, 36.73618f, 37.25270f, 37.76566f
    };

    std::cout << '\n'
	      << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "Ai'(x)"
	      << ' ' << std::setw(w) << "Ai(x)"
	      << '\n';
    for (unsigned i = 0; i < zai.size(); ++i)
      {
	auto xai = -zai[i];
	auto xai_prev = xai;
	auto aip_prev = _Tp(0);
	auto ai_prev = _Tp(0);
	auto iter = 0;
	do
	  {
	    xai_prev = xai;
	    auto [x, ai, aip, bi, bip] = emsr::detail::airy(xai);
	    xai -= ai / aip;
	    aip_prev = aip;
	    ai_prev = ai;
	    if (++iter > max_iter)
	      break;
	  }
	while (std::abs(xai - xai_prev) >= s_eps * std::abs(xai));
	zai[i] = -xai;
	std::cout << ' ' << std::setw(w) << xai
		  << ' ' << std::setw(w) << aip_prev
		  << ' ' << std::setw(w) << ai_prev
		  << '\n';
      }

    std::cout << '\n'
	      << ' ' << std::setw(w) << "x'"
	      << ' ' << std::setw(w) << "Ai(x')"
	      << ' ' << std::setw(w) << "Ai'(x')"
	      << '\n';
    for (unsigned i = 0; i < zai.size(); ++i)
      {
	auto xaip = -zaip[i];
	auto xaip_prev = xaip;
	auto ai_prev = _Tp(0);
	auto aip_prev = _Tp(0);
	//auto aipp_prev = _Tp(0);
	auto iter = 0;
	do
	  {
	    xaip_prev = xaip;
	    auto [x, ai, aip, bi, bip] = emsr::detail::airy(xaip);
	    auto aipp = xaip * ai; // Ai''(x) = xAi(x)
	    xaip -= aip / aipp;
	    ai_prev = ai;
	    aip_prev = aip;
	    //aipp_prev = aipp;
	    if (++iter > max_iter)
	      break;
	  }
	while (std::abs(xaip - xaip_prev) >= 2 * s_eps * std::abs(xaip));
	zaip[i] = -xaip;
	std::cout << ' ' << std::setw(w) << xaip
		  << ' ' << std::setw(w) << ai_prev
		  << ' ' << std::setw(w) << aip_prev
		  << '\n';
      }

    // The roots of Bi(x) interlace those of Ai(x)
    // and are close to the roots of Ai'(x).
    std::vector<_Tp> zbi;
    std::cout << '\n'
	      << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "Bi'(x)"
	      << ' ' << std::setw(w) << "Bi(x)"
	      << '\n';
    auto zai_lower = _Tp(0);
    for (unsigned i = 0; i < zai.size() - 1; ++i)
      {
	auto xbi = (-zai_lower - zai[i]) / 2;
	auto xbi_prev = xbi;
	auto bip_prev = _Tp(0);
	auto bi_prev = _Tp(0);
	auto iter = 0;
	do
	  {
	    xbi_prev = xbi;
	    auto [x, ai, aip, bi, bip] = emsr::detail::airy(xbi);
	    xbi -= bi / bip;
	    bip_prev = bip;
	    bi_prev = bi;
	    if (++iter > max_iter)
	      break;
	  }
	while (std::abs(xbi - xbi_prev) >= s_eps * std::abs(xbi));
	zbi.push_back(-xbi);
	std::cout << ' ' << std::setw(w) << xbi
		  << ' ' << std::setw(w) << bip_prev
		  << ' ' << std::setw(w) << bi_prev
		  << '\n';
	zai_lower = zai[i];
      }

    // The roots of Bi'(x) interlace those of Ai'(x)
    // and are close to the roots of Ai(x).
    std::vector<_Tp> zbip;
    std::cout << '\n'
	      << ' ' << std::setw(w) << "x'"
	      << ' ' << std::setw(w) << "Bi(x')"
	      << ' ' << std::setw(w) << "Bi'(x')"
	      << '\n';
    for (unsigned i = 0; i < zai.size() - 1; ++i)
      {
	auto xbip = -zai[i];
	auto xbip_prev = xbip;
	auto bi_prev = _Tp(0);
	auto bip_prev = _Tp(0);
	//auto bipp_prev = _Tp(0);
	auto iter = 0;
	do
	  {
	    xbip_prev = xbip;
	    auto [x, ai, aip, bi, bip] = emsr::detail::airy(xbip);
	    auto bipp = xbip * bi; // Bi''(x) = xBi(x)
	    xbip -= bip / bipp;
	    bi_prev = bi;
	    bip_prev = bip;
	    //bipp_prev = bipp;
	    if (++iter > max_iter)
	      break;
	  }
	while (std::abs(xbip - xbip_prev) >= s_eps * std::abs(xbip));
	zbip.push_back(-xbip);
	std::cout << ' ' << std::setw(w) << xbip
		  << ' ' << std::setw(w) << bi_prev
		  << ' ' << std::setw(w) << bip_prev
		  << '\n';
      }
  }


int
main()
{
  test_airy_roots<double>();
  test_airy_roots<long double>();
}
