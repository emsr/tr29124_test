
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00702_legendre_pnm_vary_02 : public TestCaseReal
    {
    public:
      TestCase_case_00702_legendre_pnm_vary_02() { }
      virtual ~TestCase_case_00702_legendre_pnm_vary_02() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00702_legendre_pnm_vary_02");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        static const e_float delta = ef::euler_gamma() / 100;
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float xk = ef::one_minus() +  (e_float(k) / 5);
          const e_float x  = k <= 5 ? xk + delta : xk - delta;
          const INT32   n  = -((10 - k) * 13);
          const INT32   m  = -((10 - k) * 11);
          data[k] = ef::legendre_p(n, m, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-9.8929251834237557806480930267841486778911263455702574821463607207539185077144789296299985598613280736845816981451481571324527838821681321869603049831295159363095136416535304528544986391627242684507878312828198364728185404084950395724159518261448444314770666363444120789800285991220344441641160165346079430096531986199027691226684681842764478924089504614123623930092407029251604376561582244293222649737e-319"),
           e_float("-8.4145206051996139758581077519068319798447769124596789240402955046669623152975624084794464308757508003478843646323147726121272779225352103878985089652118040505123577731850905899774330091582690077777476923970403777181486669252336956345947768914498715702423921330872085112003165406342108633656367688351797822732564170482480309091362009024768081378700892253726425128280333430835276056257877134893226212351e-210"),
           e_float("-9.284925818058635968135851344358823783539584242198008911654103055398195003352465293226034139108620586930102635683561241802945410940032043261095868611888527587820309673068813270866100144466545752495696357158024175064133024887113636385356747058960567732969579046926350545414426069819946269307246514263195557854298592089225259967124975076528169128361112304517651100045101743363050838363504126383569697941e-174"),
           e_float("8.4134184508782922596132914461705117521219752609452936941238522531303385938136820563342481381751111974035352680062476640000474404282473240937180955045956061990144829606366942309570236393391358197545530817782937691040836735781619342566455699318794116985151098704714915113192098293991964411564140683816400643013124941142693756499204959444021235004103713401589208920220823866511388550185555649505213201326e-147"),
           e_float("1.3199003945914887995296630997854770457983825725306902907742036390553166358589586288302703359676445947032034094833495283984523877423891205712486979509064110117354739577402069457320963052597660239101613242848077289825015589016260268210836666339381721110919511884624007515497527212453417976778205307528256522592718751897218601849855004888310331847828164908346857085057539197792776978950361471977713795045e-121"),
           e_float("6.780173185478943390645199761054432386521350785013634947979444982862578115232119800955513903359629593796441547692956745412385292932756343245534390417241724581441252071625301450044750954015533444287138219784553009678270247001932736648314884020952733760545685346801408186635575435439509511487804227322938990333789628230660336753918006983579074805386679620292940937023267551498423896581017924332368701058e-98"),
           e_float("1.0144650401061410783661716024372569034763085299560360100335813810278135646017859174775097802775135063157236435699079960455625278885249622171112117824705420294461959783523765545991718955153138262878417551744464076696437708408033546960204541987794932228761647519383716871648457392025007294423755185738023015012460075596077538353179218047003411411006581908172396888336650416606817091889129285094916934522e-73"),
           e_float("2.3113245350586958751691579888283854815343408244102181844152168502676369986383519825773951449362029132142666633896305572208718280356902160060827285791769559585832949192267283755950901968882647613044251876607145496929003531527176890410433578570952891721906113253088626440213975400519997353874710855062935256130300648465288339931481679204925248865465994776377129236456658127830345724724740103217894088482e-51"),
           e_float("3.2533544592147699976659070874200924511451001687671521125874122034270486524393387787689487045787404546469385902969105932025602561709079083477623457642135072772766539463386378493926005211164520835457717344082072945777974074232772579459594758933733330773074829454869449030116598160334411796610680172194154122507919627002653188895093878594655113560321193234102933381925228190048742318212290540360741543465e-31"),
           e_float("4.0496038862743418037824823602468766681813396674875783636139676543159809440946914304476440313632570947074475212624458172950544562271490123505965874111147095305662973169460151195990834079143484161490331400737462424273565190873883829803535520432029687970761352772933358118197311374428034552121066686410613581155061425562725003366220518673849507514551037693751271780598477003166458949673022475488080027631e-14"),
           e_float("1."),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00702_legendre_pnm_vary_02(const bool b_write_output)
    {
      return TestCase_case_00702_legendre_pnm_vary_02().execute(b_write_output);
    }
  }
}