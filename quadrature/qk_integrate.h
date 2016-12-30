// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.
//
// Ported from GSL by Jason Dick
// Originally written by Brian Gaugh
//
// Implements Gauss-Kronrod integration
// Based upon gsl-2.3/integration/qk.c

#ifndef QK_INTEGRATE_H
#define QK_INTEGRATE_H 1

#include <vector>
#include <cmath>
#include <cassert>
#include <array>
#include <functional>
#include <tuple>
#include <stdexcept>

#include "integration_error.h"

namespace __gnu_test
{

  template<typename _Tp>
    inline bool
    __test_positivity(_Tp __result, _Tp __resabs)
    {
      return(std::abs(__result) >=
		(1 - 50 * std::numeric_limits<_Tp>::epsilon()) * __resabs);
    }

  enum qk_intrule
  {
    QK_15,
    QK_21,
    QK_31,
    QK_41,
    QK_51,
    QK_61
  };

  namespace __detail
  {
    // Class template for internal implementation of
    // each individual integration rule.
    template<typename _Tp, typename _FuncTp, qk_intrule sz>
      class qk_integrator;
  }

  // Integrates func from a to b using integration rule qkintrule
  // returns a tuple with the results of a single Gauss-Kronrod integration
  // Based upon GSL function gsl_integration_qk()
  // Return values are as follows:
  // 0: result (result of integration using Kronrod scheme)
  // 1: abserr (Estimated error as difference between Gauss and Kronrod)
  // 2: resabs (Integral of absolute value of function)
  // 3: resasc (unknown)
  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp, _Tp, _Tp>
    qk_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		 const qk_intrule __qkintrule)
    {
      //Determine the integration function to use for this routine
      typedef
	std::function<void(const _FuncTp&, _Tp, _Tp,
		 _Tp&, _Tp&, _Tp&, _Tp&)>
	qk_int_func_type;

      qk_int_func_type __qk_int_func;

      switch(__qkintrule)
	{
	case QK_15: __qk_int_func =
	  (qk_int_func_type)(&(__detail::qk_integrator<_Tp, _FuncTp, QK_15>::_S_integrate));
	  break;
	case QK_21: __qk_int_func =
	  (qk_int_func_type)(&(__detail::qk_integrator<_Tp, _FuncTp, QK_21>::_S_integrate));
	  break;
	case QK_31: __qk_int_func =
	  (qk_int_func_type)(&(__detail::qk_integrator<_Tp, _FuncTp, QK_31>::_S_integrate));
	  break;
	case QK_41: __qk_int_func =
	  (qk_int_func_type)(&(__detail::qk_integrator<_Tp, _FuncTp, QK_41>::_S_integrate));
	  break;
	case QK_51: __qk_int_func =
	  (qk_int_func_type)(&(__detail::qk_integrator<_Tp, _FuncTp, QK_51>::_S_integrate));
	  break;
	case QK_61: __qk_int_func =
	  (qk_int_func_type)(&(__detail::qk_integrator<_Tp, _FuncTp, QK_61>::_S_integrate));
	  break;
	default:
	  std::__throw_logic_error("qk_integrate: "
				   "Unrecognized Gauss-Kronrod integration size");
	}

      _Tp __result, __abserr, __resabs, __resasc;
      __qk_int_func(__func, __a, __b, __result, __abserr, __resabs, __resasc);
      return std::make_tuple(__result, __abserr, __resabs, __resasc);
    }

namespace __detail
{

  template<typename _Tp, typename _FuncTp, std::size_t __ksz, std::size_t __gsz>
    void
    qk_integrate(const std::array<_Tp, __ksz> &__xgk,
		 const std::array<_Tp, __gsz> &__wg,
		 const std::array<_Tp, __ksz> &__wgk,
		 const _FuncTp& __func, _Tp __a, _Tp __b,
		 _Tp& __result, _Tp& __abserr,
		 _Tp& __resabs, _Tp& __resasc)
    {
      const auto __center = (__a + __b) / _Tp{2};
      const auto __half_length = (__b - __a) / _Tp{2};
      const auto __abs_half_length = std::abs(__half_length);
      const auto __f_center = __func(__center);
      std::array<_Tp, __ksz> __fv1;
      std::array<_Tp, __ksz> __fv2;

      assert((__ksz / 2) == __gsz);

      auto __result_gauss = _Tp{0};
      auto __result_kronrod = __f_center * __wgk[__ksz - 1];
      auto __result_abs = std::abs(__result_kronrod);

      if (__ksz % 2 == 0)
	__result_gauss = __f_center * __wg[__ksz / 2 - 1];

      for (std::size_t __jj = 0; __jj < (__ksz - 1) / 2; ++__jj)
	{
	  const std::size_t __jtw = __jj * 2 + 1;
	  const auto __abscissa = __half_length * __xgk[__jtw];
	  const auto __fval1 = __func(__center - __abscissa);
	  const auto __fval2 = __func(__center + __abscissa);
	  const auto __fsum = __fval1 + __fval2;
	  __fv1[__jtw] = __fval1;
	  __fv2[__jtw] = __fval2;

	  __result_gauss += __wg[__jj] * __fsum;
	  __result_kronrod += __wgk[__jtw] * __fsum;
	  __result_abs += __wgk[__jtw]
			* (std::abs(__fval1) + std::abs(__fval2));
	}

      for (std::size_t __jj = 0; __jj < __ksz / 2; ++__jj)
	{
	  std::size_t __jtwm1 = __jj * 2;
	  const auto __abscissa = __half_length * __xgk[__jtwm1];
	  const auto __fval1 = __func(__center - __abscissa);
	  const auto __fval2 = __func(__center + __abscissa);
	  __fv1[__jtwm1] = __fval1;
	  __fv2[__jtwm1] = __fval2;

	  __result_kronrod += __wgk[__jtwm1] * (__fval1 + __fval2);
	  __result_abs += __wgk[__jtwm1]
			* (std::abs(__fval1) + std::abs(__fval2));
	}

      auto __mean = __result_kronrod / _Tp{2};
      auto __result_asc = __wgk[__ksz - 1] * std::abs(__f_center - __mean);

      for (std::size_t __jj = 0; __jj < __ksz - 1; ++__jj)
	__result_asc += __wgk[__jj]
		      * (std::abs(__fv1[__jj] - __mean)
		       + std::abs(__fv2[__jj] - __mean));

      auto __err = (__result_kronrod - __result_gauss) * __half_length;

      __result_kronrod *= __half_length;
      __result_abs *= __abs_half_length;
      __result_asc *= __abs_half_length;

      __result = __result_kronrod;
      __resabs = __result_abs;
      __resasc = __result_asc;
      __abserr = __rescale_error(__err, __result_abs, __result_asc);
    }

  template<typename _Tp, typename _FuncTp>
    class qk_integrator<_Tp, _FuncTp, QK_15>
    {

    private:

      qk_integrator() = delete;

    public:

      static void
      _S_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		   _Tp& __result, _Tp& __abserr,
		   _Tp& __resabs, _Tp& __resasc)
      {
	using namespace __detail;

	// Abscissae of the 15-point Kronrod rule
	const std::array<_Tp, 8>
	_S_xgk =
	{
	  _Tp{0.991455371120812639206854697526329L},
	  _Tp{0.949107912342758524526189684047851L},
	  _Tp{0.864864423359769072789712788640926L},
	  _Tp{0.741531185599394439863864773280788L},
	  _Tp{0.586087235467691130294144838258730L},
	  _Tp{0.405845151377397166906606412076961L},
	  _Tp{0.207784955007898467600689403773245L},
	  _Tp{0.000000000000000000000000000000000L}
	};
	// Weights of the 7-point Gauss rule
	const std::array<_Tp, 4>
	_S_wg =
	{
	  _Tp{0.129484966168869693270611432679082L},
	  _Tp{0.279705391489276667901467771423780L},
	  _Tp{0.381830050505118944950369775488975L},
	  _Tp{0.417959183673469387755102040816327L}
	};
	// Weights of the 15-point Kronrod rule
	const std::array<_Tp, 8>
	_S_wgk =
	{
	  _Tp{0.022935322010529224963732008058970L},
	  _Tp{0.063092092629978553290700663189204L},
	  _Tp{0.104790010322250183839876322541518L},
	  _Tp{0.140653259715525918745189590510238L},
	  _Tp{0.169004726639267902826583426598550L},
	  _Tp{0.190350578064785409913256402421014L},
	  _Tp{0.204432940075298892414161999234649L},
	  _Tp{0.209482141084727828012999174891714L}
	};

	qk_integrate(_S_xgk, _S_wg, _S_wgk, __func, __a, __b, __result,
		     __abserr, __resabs, __resasc);
      }
    };

  template<typename _Tp, typename _FuncTp>
    class qk_integrator<_Tp, _FuncTp, QK_21>
    {

    private:

      qk_integrator() = delete;

      // Abscissae of the 21-point Kronrod rule
      static constexpr std::array<_Tp, 11>
      _S_xgk =
      {
	_Tp{0.995657163025808080735527280689003L},
	_Tp{0.973906528517171720077964012084452L},
	_Tp{0.930157491355708226001207180059508L},
	_Tp{0.865063366688984510732096688423493L},
	_Tp{0.780817726586416897063717578345042L},
	_Tp{0.679409568299024406234327365114874L},
	_Tp{0.562757134668604683339000099272694L},
	_Tp{0.433395394129247190799265943165784L},
	_Tp{0.294392862701460198131126603103866L},
	_Tp{0.148874338981631210884826001129720L},
	_Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 10-point Gauss rule
      static constexpr std::array<_Tp, 5>
      _S_wg =
      {
	_Tp{0.066671344308688137593568809893332L},
	_Tp{0.149451349150580593145776339657697L},
	_Tp{0.219086362515982043995534934228163L},
	_Tp{0.269266719309996355091226921569469L},
	_Tp{0.295524224714752870173892994651338L}
      };
      // Weights of the 21-point Kronrod rule
      static constexpr std::array<_Tp, 11>
      _S_wgk =
      {
	_Tp{0.011694638867371874278064396062192L},
	_Tp{0.032558162307964727478818972459390L},
	_Tp{0.054755896574351996031381300244580L},
	_Tp{0.075039674810919952767043140916190L},
	_Tp{0.093125454583697605535065465083366L},
	_Tp{0.109387158802297641899210590325805L},
	_Tp{0.123491976262065851077958109831074L},
	_Tp{0.134709217311473325928054001771707L},
	_Tp{0.142775938577060080797094273138717L},
	_Tp{0.147739104901338491374841515972068L},
	_Tp{0.149445554002916905664936468389821L}
      };

    public:

      static void
      _S_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		   _Tp& __result, _Tp& __abserr,
		   _Tp& __resabs, _Tp& __resasc)
      {
	using namespace __detail;

	qk_integrate(_S_xgk, _S_wg, _S_wgk, __func, __a, __b, __result,
		     __abserr, __resabs, __resasc);
      }
    };

  template<typename _Tp, typename _FuncTp>
    class qk_integrator<_Tp, _FuncTp, QK_31>
    {

    private:

      qk_integrator() = delete;

      // Abscissae of the 31-point Kronrod rule
      static constexpr std::array<_Tp, 16>
      _S_xgk =
      {
	_Tp{0.998002298693397060285172840152271L},
	_Tp{0.987992518020485428489565718586613L},
	_Tp{0.967739075679139134257347978784337L},
	_Tp{0.937273392400705904307758947710209L},
	_Tp{0.897264532344081900882509656454496L},
	_Tp{0.848206583410427216200648320774217L},
	_Tp{0.790418501442465932967649294817947L},
	_Tp{0.724417731360170047416186054613938L},
	_Tp{0.650996741297416970533735895313275L},
	_Tp{0.570972172608538847537226737253911L},
	_Tp{0.485081863640239680693655740232351L},
	_Tp{0.394151347077563369897207370981045L},
	_Tp{0.299180007153168812166780024266389L},
	_Tp{0.201194093997434522300628303394596L},
	_Tp{0.101142066918717499027074231447392L},
	_Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 15-point Gauss rule
      static constexpr std::array<_Tp, 8>
      _S_wg =
      {
	_Tp{0.030753241996117268354628393577204L},
	_Tp{0.070366047488108124709267416450667L},
	_Tp{0.107159220467171935011869546685869L},
	_Tp{0.139570677926154314447804794511028L},
	_Tp{0.166269205816993933553200860481209L},
	_Tp{0.186161000015562211026800561866423L},
	_Tp{0.198431485327111576456118326443839L},
	_Tp{0.202578241925561272880620199967519L}
      };
      // Weights of the 31-point Kronrod rule
      static constexpr std::array<_Tp, 16>
      _S_wgk =
      {
	_Tp{0.005377479872923348987792051430128L},
	_Tp{0.015007947329316122538374763075807L},
	_Tp{0.025460847326715320186874001019653L},
	_Tp{0.035346360791375846222037948478360L},
	_Tp{0.044589751324764876608227299373280L},
	_Tp{0.053481524690928087265343147239430L},
	_Tp{0.062009567800670640285139230960803L},
	_Tp{0.069854121318728258709520077099147L},
	_Tp{0.076849680757720378894432777482659L},
	_Tp{0.083080502823133021038289247286104L},
	_Tp{0.088564443056211770647275443693774L},
	_Tp{0.093126598170825321225486872747346L},
	_Tp{0.096642726983623678505179907627589L},
	_Tp{0.099173598721791959332393173484603L},
	_Tp{0.100769845523875595044946662617570L},
	_Tp{0.101330007014791549017374792767493L}
      };

    public:

      static void
      _S_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		   _Tp& __result, _Tp& __abserr,
		   _Tp& __resabs, _Tp& __resasc)
      {
	using namespace __detail;

	qk_integrate(_S_xgk, _S_wg, _S_wgk, __func, __a, __b, __result,
		     __abserr, __resabs, __resasc);
      }
    };

  template<typename _Tp, typename _FuncTp>
    class qk_integrator<_Tp, _FuncTp, QK_41>
    {

    private:

      qk_integrator() = delete;

      // Abscissae of the 41-point Kronrod rule
      static constexpr std::array<_Tp, 21>
      _S_xgk =
      {
	_Tp{0.998859031588277663838315576545863L},
	_Tp{0.993128599185094924786122388471320L},
	_Tp{0.981507877450250259193342994720217L},
	_Tp{0.963971927277913791267666131197277L},
	_Tp{0.940822633831754753519982722212443L},
	_Tp{0.912234428251325905867752441203298L},
	_Tp{0.878276811252281976077442995113078L},
	_Tp{0.839116971822218823394529061701521L},
	_Tp{0.795041428837551198350638833272788L},
	_Tp{0.746331906460150792614305070355642L},
	_Tp{0.693237656334751384805490711845932L},
	_Tp{0.636053680726515025452836696226286L},
	_Tp{0.575140446819710315342946036586425L},
	_Tp{0.510867001950827098004364050955251L},
	_Tp{0.443593175238725103199992213492640L},
	_Tp{0.373706088715419560672548177024927L},
	_Tp{0.301627868114913004320555356858592L},
	_Tp{0.227785851141645078080496195368575L},
	_Tp{0.152605465240922675505220241022678L},
	_Tp{0.076526521133497333754640409398838L},
	_Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 20-point Gauss rule
      static constexpr std::array<_Tp, 10>
      _S_wg =
      {
	_Tp{0.017614007139152118311861962351853L},
	_Tp{0.040601429800386941331039952274932L},
	_Tp{0.062672048334109063569506535187042L},
	_Tp{0.083276741576704748724758143222046L},
	_Tp{0.101930119817240435036750135480350L},
	_Tp{0.118194531961518417312377377711382L},
	_Tp{0.131688638449176626898494499748163L},
	_Tp{0.142096109318382051329298325067165L},
	_Tp{0.149172986472603746787828737001969L},
	_Tp{0.152753387130725850698084331955098L}
      };
      // Weights of the 41-point Kronrod rule
      static constexpr std::array<_Tp, 21>
      _S_wgk =
      {
	_Tp{0.003073583718520531501218293246031L},
	_Tp{0.008600269855642942198661787950102L},
	_Tp{0.014626169256971252983787960308868L},
	_Tp{0.020388373461266523598010231432755L},
	_Tp{0.025882133604951158834505067096153L},
	_Tp{0.031287306777032798958543119323801L},
	_Tp{0.036600169758200798030557240707211L},
	_Tp{0.041668873327973686263788305936895L},
	_Tp{0.046434821867497674720231880926108L},
	_Tp{0.050944573923728691932707670050345L},
	_Tp{0.055195105348285994744832372419777L},
	_Tp{0.059111400880639572374967220648594L},
	_Tp{0.062653237554781168025870122174255L},
	_Tp{0.065834597133618422111563556969398L},
	_Tp{0.068648672928521619345623411885368L},
	_Tp{0.071054423553444068305790361723210L},
	_Tp{0.073030690332786667495189417658913L},
	_Tp{0.074582875400499188986581418362488L},
	_Tp{0.075704497684556674659542775376617L},
	_Tp{0.076377867672080736705502835038061L},
	_Tp{0.076600711917999656445049901530102L}
      };

    public:

      static void
      _S_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		   _Tp& __result, _Tp& __abserr,
		   _Tp& __resabs, _Tp& __resasc)
      {
	using namespace __detail;

	qk_integrate(_S_xgk, _S_wg, _S_wgk, __func, __a, __b, __result,
		     __abserr, __resabs, __resasc);
      }
    };

  template<typename _Tp, typename _FuncTp>
    class qk_integrator<_Tp, _FuncTp, QK_51>
    {

    private:

      qk_integrator() = delete;

      // Abscissae of the 51-point Kronrod rule
      static constexpr std::array<_Tp, 26>
      _S_xgk =
      {
	_Tp{0.999262104992609834193457486540341L},
	_Tp{0.995556969790498097908784946893902L},
	_Tp{0.988035794534077247637331014577406L},
	_Tp{0.976663921459517511498315386479594L},
	_Tp{0.961614986425842512418130033660167L},
	_Tp{0.942974571228974339414011169658471L},
	_Tp{0.920747115281701561746346084546331L},
	_Tp{0.894991997878275368851042006782805L},
	_Tp{0.865847065293275595448996969588340L},
	_Tp{0.833442628760834001421021108693570L},
	_Tp{0.797873797998500059410410904994307L},
	_Tp{0.759259263037357630577282865204361L},
	_Tp{0.717766406813084388186654079773298L},
	_Tp{0.673566368473468364485120633247622L},
	_Tp{0.626810099010317412788122681624518L},
	_Tp{0.577662930241222967723689841612654L},
	_Tp{0.526325284334719182599623778158010L},
	_Tp{0.473002731445714960522182115009192L},
	_Tp{0.417885382193037748851814394594572L},
	_Tp{0.361172305809387837735821730127641L},
	_Tp{0.303089538931107830167478909980339L},
	_Tp{0.243866883720988432045190362797452L},
	_Tp{0.183718939421048892015969888759528L},
	_Tp{0.122864692610710396387359818808037L},
	_Tp{0.061544483005685078886546392366797L},
	_Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 25-point Gauss rule
      static constexpr std::array<_Tp, 13>
      _S_wg =
      {
	_Tp{0.011393798501026287947902964113235L},
	_Tp{0.026354986615032137261901815295299L},
	_Tp{0.040939156701306312655623487711646L},
	_Tp{0.054904695975835191925936891540473L},
	_Tp{0.068038333812356917207187185656708L},
	_Tp{0.080140700335001018013234959669111L},
	_Tp{0.091028261982963649811497220702892L},
	_Tp{0.100535949067050644202206890392686L},
	_Tp{0.108519624474263653116093957050117L},
	_Tp{0.114858259145711648339325545869556L},
	_Tp{0.119455763535784772228178126512901L},
	_Tp{0.122242442990310041688959518945852L},
	_Tp{0.123176053726715451203902873079050L}
      };
      // Weights of the 51-point Kronrod rule
      static constexpr std::array<_Tp, 26>
      _S_wgk =
      {
	_Tp{0.001987383892330315926507851882843L},
	_Tp{0.005561932135356713758040236901066L},
	_Tp{0.009473973386174151607207710523655L},
	_Tp{0.013236229195571674813656405846976L},
	_Tp{0.016847817709128298231516667536336L},
	_Tp{0.020435371145882835456568292235939L},
	_Tp{0.024009945606953216220092489164881L},
	_Tp{0.027475317587851737802948455517811L},
	_Tp{0.030792300167387488891109020215229L},
	_Tp{0.034002130274329337836748795229551L},
	_Tp{0.037116271483415543560330625367620L},
	_Tp{0.040083825504032382074839284467076L},
	_Tp{0.042872845020170049476895792439495L},
	_Tp{0.045502913049921788909870584752660L},
	_Tp{0.047982537138836713906392255756915L},
	_Tp{0.050277679080715671963325259433440L},
	_Tp{0.052362885806407475864366712137873L},
	_Tp{0.054251129888545490144543370459876L},
	_Tp{0.055950811220412317308240686382747L},
	_Tp{0.057437116361567832853582693939506L},
	_Tp{0.058689680022394207961974175856788L},
	_Tp{0.059720340324174059979099291932562L},
	_Tp{0.060539455376045862945360267517565L},
	_Tp{0.061128509717053048305859030416293L},
	_Tp{0.061471189871425316661544131965264L},
	_Tp{0.061580818067832935078759824240066L}
      };

    public:

      static void
      _S_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		   _Tp& __result, _Tp& __abserr,
		   _Tp& __resabs, _Tp& __resasc)
      {
	using namespace __detail;

	qk_integrate(_S_xgk, _S_wg, _S_wgk, __func, __a, __b, __result,
		     __abserr, __resabs, __resasc);
      }
    };

  template<typename _Tp, typename _FuncTp>
    class qk_integrator<_Tp, _FuncTp, QK_61>
    {

    private:

      qk_integrator() = delete;

      // Abscissae of the 61-point Kronrod rule
      static constexpr std::array<_Tp, 31>
      _S_xgk =
      {
	_Tp{0.999484410050490637571325895705811L},
	_Tp{0.996893484074649540271630050918695L},
	_Tp{0.991630996870404594858628366109486L},
	_Tp{0.983668123279747209970032581605663L},
	_Tp{0.973116322501126268374693868423707L},
	_Tp{0.960021864968307512216871025581798L},
	_Tp{0.944374444748559979415831324037439L},
	_Tp{0.926200047429274325879324277080474L},
	_Tp{0.905573307699907798546522558925958L},
	_Tp{0.882560535792052681543116462530226L},
	_Tp{0.857205233546061098958658510658944L},
	_Tp{0.829565762382768397442898119732502L},
	_Tp{0.799727835821839083013668942322683L},
	_Tp{0.767777432104826194917977340974503L},
	_Tp{0.733790062453226804726171131369528L},
	_Tp{0.697850494793315796932292388026640L},
	_Tp{0.660061064126626961370053668149271L},
	_Tp{0.620526182989242861140477556431189L},
	_Tp{0.579345235826361691756024932172540L},
	_Tp{0.536624148142019899264169793311073L},
	_Tp{0.492480467861778574993693061207709L},
	_Tp{0.447033769538089176780609900322854L},
	_Tp{0.400401254830394392535476211542661L},
	_Tp{0.352704725530878113471037207089374L},
	_Tp{0.304073202273625077372677107199257L},
	_Tp{0.254636926167889846439805129817805L},
	_Tp{0.204525116682309891438957671002025L},
	_Tp{0.153869913608583546963794672743256L},
	_Tp{0.102806937966737030147096751318001L},
	_Tp{0.051471842555317695833025213166723L},
	_Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 30-point Gauss rule
      static constexpr std::array<_Tp, 15>
      _S_wg =
      {
	_Tp{0.007968192496166605615465883474674L},
	_Tp{0.018466468311090959142302131912047L},
	_Tp{0.028784707883323369349719179611292L},
	_Tp{0.038799192569627049596801936446348L},
	_Tp{0.048402672830594052902938140422808L},
	_Tp{0.057493156217619066481721689402056L},
	_Tp{0.065974229882180495128128515115962L},
	_Tp{0.073755974737705206268243850022191L},
	_Tp{0.080755895229420215354694938460530L},
	_Tp{0.086899787201082979802387530715126L},
	_Tp{0.092122522237786128717632707087619L},
	_Tp{0.096368737174644259639468626351810L},
	_Tp{0.099593420586795267062780282103569L},
	_Tp{0.101762389748405504596428952168554L},
	_Tp{0.102852652893558840341285636705415L}
      };
      // Weights of the 61-point Kronrod rule
      static constexpr std::array<_Tp, 31>
      _S_wgk =
      {
	_Tp{0.001389013698677007624551591226760L},
	_Tp{0.003890461127099884051267201844516L},
	_Tp{0.006630703915931292173319826369750L},
	_Tp{0.009273279659517763428441146892024L},
	_Tp{0.011823015253496341742232898853251L},
	_Tp{0.014369729507045804812451432443580L},
	_Tp{0.016920889189053272627572289420322L},
	_Tp{0.019414141193942381173408951050128L},
	_Tp{0.021828035821609192297167485738339L},
	_Tp{0.024191162078080601365686370725232L},
	_Tp{0.026509954882333101610601709335075L},
	_Tp{0.028754048765041292843978785354334L},
	_Tp{0.030907257562387762472884252943092L},
	_Tp{0.032981447057483726031814191016854L},
	_Tp{0.034979338028060024137499670731468L},
	_Tp{0.036882364651821229223911065617136L},
	_Tp{0.038678945624727592950348651532281L},
	_Tp{0.040374538951535959111995279752468L},
	_Tp{0.041969810215164246147147541285970L},
	_Tp{0.043452539701356069316831728117073L},
	_Tp{0.044814800133162663192355551616723L},
	_Tp{0.046059238271006988116271735559374L},
	_Tp{0.047185546569299153945261478181099L},
	_Tp{0.048185861757087129140779492298305L},
	_Tp{0.049055434555029778887528165367238L},
	_Tp{0.049795683427074206357811569379942L},
	_Tp{0.050405921402782346840893085653585L},
	_Tp{0.050881795898749606492297473049805L},
	_Tp{0.051221547849258772170656282604944L},
	_Tp{0.051426128537459025933862879215781L},
	_Tp{0.051494729429451567558340433647099L}
      };

    public:

      static void
      _S_integrate(const _FuncTp& __func, _Tp __a, _Tp __b,
		   _Tp& __result, _Tp& __abserr,
		   _Tp& __resabs, _Tp& __resasc)
      {
	using namespace __detail;

	qk_integrate(_S_xgk, _S_wg, _S_wgk, __func, __a, __b, __result,
		     __abserr, __resabs, __resasc);
      }
    };

} // namespace __detail

} // namespace __gnu_test

#endif // QK_INTEGRATE_H
