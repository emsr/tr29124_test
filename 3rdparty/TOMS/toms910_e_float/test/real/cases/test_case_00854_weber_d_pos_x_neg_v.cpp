
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00854_weber_d_pos_x_neg_v : public TestCaseReal
    {
    public:
      TestCase_case_00854_weber_d_pos_x_neg_v() { }
      virtual ~TestCase_case_00854_weber_d_pos_x_neg_v() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00854_weber_d_pos_x_neg_v");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        const e_float pi_over_ten = ef::pi() / static_cast<INT32>(10);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float v = -ef::euler_gamma() - (k * static_cast<INT32>(10));
          const e_float x = +pi_over_ten       + (k * static_cast<INT32>(20));
          data[k] = ef::weber_d(v, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("1.0303385390113676234556998248540219488151749438105484202095134627399646920409971750331452352825633897850625876845013779493280888775148293341954496558956149357825839302353908435727299886191784099451570732903094820912470812100123886944782606367777138819501123961542534697246358305200732575449708687771364049051475626244850217918368658456918768299221820625410649430348571686468407359593839681047927419549"),
           e_float("1.9950285398788925803536593885271205683465170894626018290846912209631762827091126982020956756198667699301617353781153342234016072579118387841372720534328921818645618211652259872262904499842207651824664198645647252441408669125138984540964054182449016886852417194378229756421393244411859223214021121146887324898386409004208038562166902359489330512174482930827051474308989618833675633679572205463206501587e-59"),
           e_float("2.8077367368684046140814515748537188405188298129547985126416871000178133011798339288235702875171472550096191346178825613766416728833241406096189191854969309706518452249316592979419517488017707631048114189880934867527042818595648482768013500976557787964215255926290577457709979541209580059192749081413141184092085426156194775833125797779925971373193880216873505169978223626032448072692659162453145327332e-210"),
           e_float("3.4176124578317105484353043081965986859733289670250607265517259514436284770628347346250850444773380999440384544404161881212327355556850457606840423588289845865253410683804548968483366945332258861029083133360689298711448565986293418377989455423823873689755104031352583033358579806120250274394359172453577379064594074868177284154516323040343814992449342640444001932668630984370046695155309853995543598515e-450"),
           e_float("2.0559804912098796134634635972500415050380591494640684246601851611951400880043126684403135835676809976702610288507011060226248046551184695686307863746729655018835351825587158252123856134641643693937020862722241295788517960172530740643559787526577307577657754148149390573174899076106076770497587572539021176785367692454484026223341224546298253529560210502918964941507324202344528233259825053755285653723e-778"),
           e_float("1.4186181173471666856212658222935468274231074815020024863830320706400701433058345231252046329965124392133834729934675576034777303462488252418329070580008151621504939116069311246556423466598097167083959013514120395264771881483714640821557629818998963546906105463271800423871907989393576751984200026847247201984382343994286183714799120879103900046179147149212317751253127534522298022432409112863028179922e-1194"),
           e_float("1.850875429623245589218685284827738842964873228104236127439266273077517800708855449480713386447862039271628461867665954948991401853807118259222297134604659903859761866994536136763857340962355755614211255167678902990404092657368135221490502075620514917757019177547307159012532390185346851426437762804398571108931027419665433153370738229319534573645317681829688297267231006453863503956567093432935819226e-1698"),
           e_float("6.3639720893523032021289669085144367646037530701083565580633122639667513712581639197486475538885663693319593352394329948992813111397403622271795367641258872838005455533756579827916304162945551894605145432539811917951171303595362828446664048343290757824805692352630727847349293849842781190649261449886026045552708092220728576363643852012705897732663449611817671506550446648737408253991240472755716986772e-2290"),
           e_float("7.3068902629068803669915036675768790601550116822985459224222031737027547430079190583779712976688853121620346255708696818453463835223194774148371834147006433764151146476505310566979625949003404920880057666656001341362463616250986614693985263491307241004185693879304630396750709415780539502644300981941038658718727132311083704023593511199843637244236615082795735683639675514274531327682626266141108333714e-2969"),
           e_float("3.3454173145711549410028712508248117755293861723168019590532681517911234582674476194260795413210976199678392670670241945442414072920477214468174156681030491772608169279286571848153629873103227282215693698747613746958137381882422674774527875697154855902768941471151636220037639498215805879637183673480768198197226559576708852735395123433164920580523155791522770024303702771416545567656788325323337955167e-3735"),
           e_float("7.0114092335356433459756680878851327312096206279060814234770832902193345236032893750926874790424426758107544281534597896578316801855035102834998101107394402489268899947290546153437151878980179036084569257740244508863411959147487264781120946112799870412033367033523303945274086358041737737026816832328583538363248738990403512983660003211752934464076114838084837260401250428927239114217373743003809140845e-4589"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00854_weber_d_pos_x_neg_v(const bool b_write_output)
    {
      return TestCase_case_00854_weber_d_pos_x_neg_v().execute(b_write_output);
    }
  }
}