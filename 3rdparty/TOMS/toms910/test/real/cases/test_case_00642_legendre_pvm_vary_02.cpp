
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00642_legendre_pvm_vary_02 : public TestCaseReal
    {
    public:
      TestCase_case_00642_legendre_pvm_vary_02() { }
      virtual ~TestCase_case_00642_legendre_pvm_vary_02() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00642_legendre_pvm_vary_02");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        static const e_float delta = ef::euler_gamma() / 100;
        static const e_float sqrt_1_3 = ef::sqrt(ef::third());
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float xk = ef::one_minus() +  (e_float(k) / 5);
          const e_float x  = k <= 5 ? xk + delta : xk - delta;
          const e_float v  = -((10 - k) * 13) - sqrt_1_3;
          const INT32   m  = -((10 - k) * 11);
          data[k] = ef::legendre_p(v, m, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("4.884852646539346890787888813049477838392600369531766545486570027437227965594868915905615561935379934191941545860616790547015818269512291830397327933202989473351747342507068050074026554846727707455033549418802356174655925862686542774601932527109129846178060616558385498507462656354804439623059603876701082287835581375558084316385250471902327786051647966554414590023893762026741717309312923307801167269e-135"),
           e_float("1.1402845648058791321163640990998939380648222599154235417458284826114909749377805666053667166199291819122595352203953830949763328803686538935045942727347284328617899004495474056789667798361313757872193113242568803030179835904280412552866448343028210439320476149722528903065093215241355219633333488145395557627794295520235281785241321319862394580815358623894011231500545868071800619990658049677918801582e-189"),
           e_float("3.2258408778311567093517744119163640774430849036136330867798154495158703245468555721062626637374647801285493082747741122088698033985858606420905904786603236113287283683717097978902691975447337516622353465134473062775016697693515744117079807275930605455037079984065581918727371952950561482340425520730979339108003093916557117546548107151022645759090669482042651940178797432784375057410624348485567837355e-172"),
           e_float("-9.4326051174975253358939527868740397108896380929419587685438542911396381060417564663733673085252844939105506490826006662329479237737659582206292503159957638303061468541593366596737812013083632316281685248971833522616737741114997747633887972391876237730004877709801109137007813907517790035247456984114146025178413200690170236096141910738664495686187509680225680754472180309735418734011382494114389137343e-148"),
           e_float("2.9985052873470709587683551302915040462212503517108537364221664770598261434886943806246676940442204428595222602081414155573341163766230688526883055936404450439979337709356189858357509154934717248891680574905533803767217529467635858508012240759484318892785867688817877575967072440590553773616565963438457963493797470160987589963723367548195386754510973668724085446590508799761077877176550884504500327506e-122"),
           e_float("-1.0865257203858196419851016701202450040849320469391993108451493012170736056745380278505232995673480629771269987431016782972507203817949369236929154746815925234617897887964486659260777401844592508390494685266858234123062648300208550791636145202445882774757113795774350631781060147391503220436686644679575004375107390726761557289719923671769813337557237673203974226777494596856306672520324758421272487568e-97"),
           e_float("5.1236373100634802076276713656428602010378201050491333143027581278550269372820799681372428092248080385964972203264625417361978299848211556797980150167947676131969398707397519217210724900878823978063386726563803911793822463809717347952315801607827287853913897224832185873039934777030007457366530189379684414970484645788129145477943681257817253415248758580742467919869117199904794937352853843804998338653e-74"),
           e_float("7.367772444571398143377098938309686965894770589511767442635971138227777566879749054272574260145349019074251309473919652814054393513539409103843074801843444091972862801888610923278096145855560329346010830189934697874182201314401948469333971978922863520421384316468622037689923917060734189340078210116294157757383739129137056254915311271204193725138414815146655200846397635571647809752176116352093805991e-52"),
           e_float("2.2378489939959600964994197748008627876833944645129068985081497058797393349401167684384565303921716419603039357085966013769508240842636705378170772732909261193712889973961201966012279858130481022838998648899201090777682445615964143593139980850337980508249692447590592829956908241966713657530857675612336324671756431603710828838744139910990058098994909393668664248487468260395615548859057569393798297174e-31"),
           e_float("3.506002840847608054965032550019595204196680097365661189843345229285125404721175179787918653650048921037253739943890832156260029706597108197681296974608827603903980506952464298854700859417696400516265974109775876510626868416529380041001269390628763349682390586779482490429853616221073705352869984019393709193162417112284881559698946409902351021619448891069658525616210536906140799446302968947775443104e-14"),
           e_float("1.0007053945334960366255749083305905652993693614001548027423116259096385859296367372874412211983564861520849560998619397586530545852475248171737023144061034406504056283811182460543648403833072498205264822064177195489618645497996999370920755549793806991088637232221244798018167402046648854302501459201324519316809968628951618976771264310817239814145331322242393294293588148133700889016608020153606105609"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00642_legendre_pvm_vary_02(const bool b_write_output)
    {
      return TestCase_case_00642_legendre_pvm_vary_02().execute(b_write_output);
    }
  }
}