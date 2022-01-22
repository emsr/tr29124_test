
#include <ext/float128_math.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <limits>

#include <iostream>
#include <iomanip>

#include <ext/float128_io.h>
#include <emsr/summation.h>

template<typename Tp>
  void
  test(Tp proto = Tp{})
  {
    const auto _S_max_log = emsr::log_max(proto);

    using ABS_t = emsr::AitkenDeltaSquaredSum<emsr::BasicSum<Tp>>;
    using ShankS_t = emsr::AitkenDeltaSquaredSum<ABS_t>;

    std::cout.precision(emsr::digits10(proto));
    auto w = 8 + std::cout.precision();
    emsr::BasicSum<Tp> BS;
    emsr::AitkenDeltaSquaredSum<emsr::BasicSum<Tp>> ABS;
    emsr::AitkenDeltaSquaredSum<emsr::KahanSum<Tp>> AKS;
    ShankS_t ShankS;
    emsr::WinnEpsilonSum<emsr::BasicSum<Tp>> WBS;
    emsr::WinnEpsilonSum<emsr::KahanSum<Tp>> WKS;
    emsr::BrezinskiThetaSum<emsr::BasicSum<Tp>> BTS;

    //emsr::LevinUSum<emsr::BasicSum<Tp>> LUS;
    emsr::LevinTSum<emsr::BasicSum<Tp>> LTS;
    emsr::LevinDSum<emsr::BasicSum<Tp>> LDS;
    emsr::LevinVSum<emsr::BasicSum<Tp>> LVS;

    emsr::WenigerTauSum<emsr::BasicSum<Tp>> WTS;
    emsr::WenigerDeltaSum<emsr::BasicSum<Tp>> WDS;
    emsr::WenigerPhiSum<emsr::BasicSum<Tp>> WPS;
    emsr::WenigerDeltaSum<emsr::VanWijngaardenSum<Tp>> WDvW; // Start term: (8)

    auto s = Tp{1.2};
    auto zetaterm = [s, _S_max_log](std::size_t k)
		    -> Tp
		    { return (s * std::log(Tp(k + 1)) < _S_max_log ? std::pow(Tp(k + 1), -s) : Tp{0}); };

    auto VwT = emsr::VanWijngaardenCompressor<decltype(zetaterm)>(zetaterm);

    auto zeta = Tp{5.591582441177750776536563193423143277642L};
    std::cout << "\n\nzeta(1.2) = " << std::setw(w) << zeta << '\n';
    std::cout << std::setw(w) << "k"
	      << std::setw(w) << "Basic"
	      << std::setw(w) << "Aitken-Basic"
	      << std::setw(w) << "Aitken-Kahan"
	      << std::setw(w) << "Winn-Basic"
	      << std::setw(w) << "Winn-Kahan"
	      << std::setw(w) << "BrezinskiT-Basic"
	      << std::setw(w) << "LevinT-Basic"
	      << std::setw(w) << "LevinD-Basic"
	      << std::setw(w) << "LevinV-Basic"
	      << std::setw(w) << "WenigerTau-Basic"
	      << std::setw(w) << "WenigerDelta-Basic"
	      << std::setw(w) << "WenigerPhi-Basic"
	      << std::setw(w) << "WenigerD-vW"
	      << std::setw(w) << "Shanks"
	      << '\n';
    std::cout << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << std::setw(w) << "------------------"
	      << '\n';
    for (auto k = 0; k < 100; ++k)
      {
	auto term = std::pow(Tp(k + 1), -s);
	BS += term;
	ABS += term;
	AKS += term;
	WBS += term;
	WKS += term;
	BTS += term;
	LTS += term;
	LDS += term;
	LVS += term;
	WTS += term;
	WDS += term;
	WPS += term;
	WDvW += VwT[k];
	ShankS += term;
	std::cout << std::setw(w) << k
		  << std::setw(w) << BS()
		  << std::setw(w) << ABS()
		  << std::setw(w) << AKS()
		  << std::setw(w) << WBS()
		  << std::setw(w) << WKS()
		  << std::setw(w) << BTS()
	  	  << std::setw(w) << LTS()
	  	  << std::setw(w) << LDS()
	  	  << std::setw(w) << LVS()
	  	  << std::setw(w) << WTS()
	  	  << std::setw(w) << WDS()
	  	  << std::setw(w) << WPS()
	  	  << std::setw(w) << WDvW()
		  << std::setw(w) << ShankS()
		  << '\n';
      }

    // 2F0(1,1;;-1/z) = z e^z E_1(z)
    const auto expint_scaled
      = [](Tp z)->Tp
        { return z * std::exp(z)* __gnu_cxx::expint(1,z); };

    // 2F0(1,1;;-1/z) = z e^zE_1(z)
    for (auto z : {Tp{3}, Tp{0.5L}})
      {
	std::cout << "\n\n";
	auto a = Tp{1};
	auto b = Tp{1};
	std::cout << "  2F0(1,1;;" << -1/z << ") = 0.78625122076596\n";
	std::cout << "  expint_scaled = " << expint_scaled(z) << '\n';
	auto term = Tp{1};
	BS.reset(term);
	ABS.reset(term);
	AKS.reset(term);
	WBS.reset(term);
	WKS.reset(term);
	BTS.reset(term);
	LTS.reset(term);
	LDS.reset(term);
	LVS.reset(term);
	WTS.reset(term);
	WDS.reset(term);
	WPS.reset(term);
	ShankS.reset(term);
	std::cout << std::setw(w) << "k"
		  << std::setw(w) << "Basic"
		  << std::setw(w) << "Aitken-Basic"
		  << std::setw(w) << "Aitken-Kahan"
		  << std::setw(w) << "Winn-Basic"
		  << std::setw(w) << "Winn-Kahan"
		  << std::setw(w) << "BrezinskiT-Basic"
		  << std::setw(w) << "LevinT-Basic"
		  << std::setw(w) << "LevinD-Basic"
		  << std::setw(w) << "LevinV-Basic"
		  << std::setw(w) << "WenigerTau-Basic"
		  << std::setw(w) << "WenigerDelta-Basic"
		  << std::setw(w) << "WenigerPhi-Basic"
		  << std::setw(w) << "Shanks"
		  << '\n';
	std::cout << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << '\n';
	for (auto k = 1; k < 100; ++k)
	  {
	    std::cout << std::setw(w) << (k - 1)
		      << std::setw(w) << BS()
		      << std::setw(w) << ABS()
		      << std::setw(w) << AKS()
		      << std::setw(w) << WBS()
		      << std::setw(w) << WKS()
		      << std::setw(w) << BTS()
		      << std::setw(w) << LTS()
		      << std::setw(w) << LDS()
		      << std::setw(w) << LVS()
		      << std::setw(w) << WTS()
		      << std::setw(w) << WDS()
		      << std::setw(w) << WPS()
		      << std::setw(w) << ShankS()
		      << '\n';
	    term *= -(a + k - 1) * (b + k - 1) / z / k; // 2F0(1,1;;-1/z)
	    BS += term;
	    ABS += term;
	    AKS += term;
	    WBS += term;
	    WKS += term;
	    BTS += term;
	    LTS += term;
	    LDS += term;
	    LVS += term;
	    WTS += term;
	    WDS += term;
	    WPS += term;
	    ShankS += term;
	  }
      }

    // 2F1(1,1;2;-z) = log(1+z)
    for (auto z : {Tp{5}, Tp{1}, Tp{-0.9L}})
      {
	auto a = Tp{1};
	auto b = Tp{1};
	auto c = Tp{2};
	std::cout << "\n\n";
	//std::cout << "  2F1(1,1;2;" << -z << ") = " << std::setw(w) << __gnu_cxx::hyperg(a, b, c, z) << "\n";
	std::cout << "  log(1 + " << z << ") / (" << z << ") = " << std::log1p(z) / z << '\n';
	auto term = Tp{1};
	BS.reset(term);
	ABS.reset(term);
	AKS.reset(term);
	WBS.reset(term);
	WKS.reset(term);
	BTS.reset(term);
	LTS.reset(term);
	LDS.reset(term);
	LVS.reset(term);
	WTS.reset(term);
	WDS.reset(term);
	WPS.reset(term);
	ShankS.reset(term);
	std::cout << std::setw(w) << "k"
		  << std::setw(w) << "Basic"
		  << std::setw(w) << "Aitken-Basic"
		  << std::setw(w) << "Aitken-Kahan"
		  << std::setw(w) << "Winn-Basic"
		  << std::setw(w) << "Winn-Kahan"
		  << std::setw(w) << "BrezinskiT-Basic"
		  << std::setw(w) << "LevinT-Basic"
		  << std::setw(w) << "LevinD-Basic"
		  << std::setw(w) << "LevinV-Basic"
		  << std::setw(w) << "WenigerTau-Basic"
		  << std::setw(w) << "WenigerDelta-Basic"
		  << std::setw(w) << "WenigerPhi-Basic"
		  << std::setw(w) << "Shanks"
		  << '\n';
	std::cout << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << std::setw(w) << "------------------"
		  << '\n';
	for (auto k = 1; k < 100; ++k)
	  {
	    std::cout << std::setw(w) << (k - 1)
		      << std::setw(w) << BS()
		      << std::setw(w) << ABS()
		      << std::setw(w) << AKS()
		      << std::setw(w) << WBS()
		      << std::setw(w) << WKS()
		      << std::setw(w) << BTS()
		      << std::setw(w) << LTS()
		      << std::setw(w) << LDS()
		      << std::setw(w) << LVS()
		      << std::setw(w) << WTS()
		      << std::setw(w) << WDS()
		      << std::setw(w) << WPS()
		      << std::setw(w) << ShankS()
		      << '\n';
	    term *= -((a + k - 1) / (c + k - 1)) * ((b + k - 1) / k) * z; // 2F1(1,1;2;-z)
	    BS += term;
	    ABS += term;
	    AKS += term;
	    WBS += term;
	    WKS += term;
	    BTS += term;
	    LTS += term;
	    LDS += term;
	    LVS += term;
	    WTS += term;
	    WDS += term;
	    WPS += term;
	    ShankS += term;
	  }
      }
  }

int
main()
{
  //std::cout << "\nfloat\n=====\n\n";
  //test<float>();

  std::cout << "\ndouble\n======\n";
  test<double>();

  std::cout << "\nlong double\n===========\n";
  test<long double>();

  //std::cout << "\n__float128\n===========\n";
  //test<__float128>();
}
