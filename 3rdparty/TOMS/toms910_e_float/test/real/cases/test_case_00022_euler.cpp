
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00022_euler : public TestCaseReal
    {
    public:
      TestCase_case_00022_euler() { }
      virtual ~TestCase_case_00022_euler() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00022_euler");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(101u);
        for(UINT32 k = static_cast<UINT32>(0u); k < static_cast<UINT32>(data.size()); k++)
        {
          data[static_cast<std::size_t>(k)] = ef::euler(10u * k);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 101u> a =
        {{
           e_float("1."),
           e_float("-50521."),
           e_float("3.70371188237525e14"),
           e_float("-4.41543893249023104553682821e26"),
           e_float("1.4851150718114980017877156781405826684425e40"),
           e_float("-6.053285248188621896314383785111649088103498225146815121e54"),
           e_float("1.8108911496579230496545807741652158688733487349236314106008095454231325e70"),
           e_float("-2.85051783223697718732198729556739339504255241778255239879353211106980427546235397447421e86"),
           e_float("1.8622915758412697044482492303043126011920010194518556063577101095681956123546201442832293837005396878225e103"),
           e_float("-4.227240686139909064705589929214593102933845388672369082676644542650248228369590525634078984302153217507945782396923579721e120"),
           e_float("2.903528346661097497054603834764435875077553006646158945080492319146997643370625023889353447129967354174648294748510553528692457632980625125e138"),
           e_float("-5.403078659795293205619115494263476990474882718296420197408516010249814711845963278802717760330206627628886347827381068652979312079075986736179581399343212021e156"),
           e_float("2.4883915747829871631690245540848940823728670907090814055549996853018422439857255460434636907179279971030115914025391078487144429408300462747699810654037377064816073847531472025e175"),
           e_float("-2.63038464627282201918918005755736144830032766236489151905754681577981916252171316354973105486693742840216041344801990544825902436415863797115667159863613239597675722668781930227771509440274344321e194"),
           e_float("5.987386904215954780609340300928990507458814865030864762440959127489603790251637145749604626534662068704764148470130069641596881837375711869240420745498089109837751540966769902924250247702626992182004595330375418925e213"),
           e_float("-2.77857404780457414987248665136951661385802997430099285750035845293762404084709969434016870389430000091054914408044748463535290871776626229195061384515613583189227656542091580147189163312909279754107124154692415552717409356305334976621e233"),
           e_float("2.5071830005737144960191522234762834423570096759890580082823874229833831690163136732846143787945517516134004403918012093742445562612207991873892705782881588541832155633594558800737824587847195977958539176361385963043761572318045891668199089009014838465825e253"),
           e_float("-4.220005513130260808256874149121608869388378226451438744419711766672691776484360898579943327229799071566569220799418193734460802563701559836740994854234132201397018022945474085831250467336792289975504710405831631768298471319957990591150241176345336500257174834567448143108921e273"),
           e_float("1.277331663671980642072877732151869284598756140342756862203699204780475335159504968151690204213663447105620949434717932676191266292193883436288813754387946456540481384238565159064871619288178168239418776112961210567471496320773832611346566378550845636826548410755837970746928546952030526246612725e294"),
           e_float("-6.73012788703472215217655714452463930422605593618476965163010790827494348121050458458808167113687985498382728742497229884244318774989604658830071458046921180500294567182145032339021048987567078430867278430624720984803716761303615031141015161171583843163932861106497453752681124616309213188697054892019876070316741221e314"),
           e_float("5.99547163501087372474968414188065044926728618369425054868091366455415875680150799832459576446626662489154039676247901501453585520247873369633664981039166488273449468831053594931966683906734130221865944184112020908528590877181590543768716367090579274260976307066311777407073312889369750169548196419149504941815733845016897635621925859625e335"),
           e_float("-8.79667520899326528281240524434424310039378504656364299522943462667049171626018849425247625965505612776609078395637261953504294888673423278632880520938768842605783949016944000322188563272074322298186424686479631042754734120997056785729485897741691441995653456106698489480860557570013156424430377672202172196435017546859915449876350667757152200935201473873521e356"),
           e_float("2.075901363395202598099009907137830385539540057964019998173615402655453242515181198090400460199321282414115737226696897084063381172391102049751021657302512954176323297024271251631074962119421852023425211197366278476483991606001923412791507360459499774276435829835587509389397228281624555930682496740279047222874448945155240511472099961732981097618052346907435046143071631202206525e378"),
           e_float("-7.711183769464086426121942550551260111440150236283392297018676581121626964226573714877386417993045640550733395448580631838877494624266367194531545639043083878663568255665641723241848232539283209844775524669712356989838802810271301399636950256005800528470580199678605618653589192901551045149407585044188508348439249373067866591947018800911570010938897809349884507685841574359288848328046749779232505821e399"),
           e_float("4.4208760344171939790055258615701141234600245067934881026055300188084956655939375649533499240396777330172512340673647848245515919995417757243195486674822528085455729734945857825400714230363976612751926667013131057566777975300203484945246312109613746709709521382998865371432062440104638857372097631009039640905904858160867221067713580774575562240050092497541021849510599474200175326586153953799125115116e421"),
           e_float("-3.8417336935123322473718109050035240772730700170603370840008992734751336101683106197337696859074130235723402080505277342076817579539674043006350460172272342909386017389294597747006388405456509545302198874584980742692396098309100698111503209783112427429465759827405445033915439032919589452254966593100950819647284974288209513301120507856396223896914971360825195751044095952257843641108043158362186561768e443"),
           e_float("4.9769476587423519385901224641545544634915692906639077549038923690898726913149484752400778103263331587202410060655718571576080361487469560389973901184830248560988117380446429479096496295602440986485614154666032855832678823730596655282348571009036334329031206127730592155754837773230932132316755149781145751279516450829927459373718060622128162156970344642371683713768849343028293994854983976572624503765e465"),
           e_float("-9.4657577836131743704237836343893513706013774118920995771103023023336490127592379021834908916415895075429581396602194752784011825700181280317615611108086232781335381791908903384108251839450162082722220430357934387282286989030572661711969180253201137553432998706197759930659510818516034633017403825412513126860787425807173586328588011926079758862762246714942131639679414850760607916491630107380892578373e487"),
           e_float("2.6057731045653757617164524096330801925386295356849712408056680605607574730485945379342317551293939788659719153597482399330193714264959391309952684101725494475182363572477963858037924051046970697352926060890746408614980094651193137447837312568376505279068337168624292278142481525116600254961435155728516510540035594116803172539955493683375627392377712989191217982674261547181012013880347760560170447689e510"),
           e_float("-1.024658504519623622872960299762326497174231942858766529385764948104914775013166512585410358364043383813459577627548171101878676477854614769487268816578528758287443848899786446404586022055695352839943781948752940434015796876402129997602737780771041550017605994136556021910006379137393938414054394884770484992169793364241743280116641943889410082769098890616313651468358944266360218077605917126180501497e533"),
           e_float("5.6852563686050582715789855530042087626546497108957923947456306633037454242223164692660463310486845949879941738638619394532581064314311727579362078288887093936753672105014711780188473923202164711381764187934127638590727240305611085738040798981873907165544254238668046181552827686778704604664787163688922140571977298415638208962689185051347676391981219529503486862975224833150672917063747203773810477862e555"),
           e_float("-4.4001903554289422810449168732377484867841409301889907189982979117227378279809290332683053685395109129311613132552914298358857832805166807643187949313536953326317679583162629388552466281759151592497764765859250543878667598397563653691920442332388533225497082409681189284966829200846963347120337700490859678801690758821228149009753633033689470168591771921878817080008669367379982318916575472042039215613e578"),
           e_float("4.6998753867827493538999957758290927215270979927209837954291305512096652114245008616330415180043382813156903894697076754416340691851500533384095663974789147658357954925656412360186195072349921142282838254171511827999687510806739227668604150782509525596700736010859321233924534981204780298026109921171367653694784797714116229524184117564047329259285018727063138506183115594560188688446615095010220465483e601"),
           e_float("-6.8584883413367882826872988492280413263535248100694537566491171922235557726361186230386017761972720055405482374201009082351011307687447048078881745340180339205436041384674382419632703846582870198428484637074405150480393652487356660928214383860187714788725686366246695333127482114963829068425933446688388025107548227576233767243124193504578829974305419131869793195511444777945382022798341743501620274543e624"),
           e_float("1.3545520972660796203898272873543390446840827965402470343388799891166483186158869555816593460621145050033450500267321088187093833136741557338289901742361123628538329177644880662597493756058795155207795137155953418166852178855431335571370595684867837704691366153176263997542963378622014951651273437139926691636284409447619445007330882733648721804904861524839454146148706384419781189323774382970918820458e648"),
           e_float("-3.5886191828791675992325416328693584695168958096194310872897134844256215704844075752672453428730652226132904149358123054381360769249676035611633339251811329844951482960256304630085969573290608201644234366954813247746047857085852576407714463165313427677709244870777953327936504658890539246486269355015216796607696255326634853488728077644056171320514566457023050861442753779063622534158969732192839989199e671"),
           e_float("1.2646865344426782466850882943634747942219482593493367008248951575295026687452869550129049630742161735529809517901321339094583206885434900960374803869800477677913790837210297758651563809486024790570704721087306247704701987584564553622451761791492582028575835818898158470079386259334852686044827289664807241813212841748255252617394378726473476853399144692757990283426857852160780226289609239466032050936e695"),
           e_float("-5.8819838578714861337254486765361141818689224148399302955831461071337143693963415534515673545527607462142631461482781817674644358672007398317752341029461493081645151446107768100943707267312100782619829670530911618938227594237229737439106375152129826259118857844850570113528326629273754893916805570300678261527671050728378023666637631502061542586865546133666455124924627803922135409059569276637175851959e718"),
           e_float("3.5834175043637737164351737597601495847313719731290609217919967211683621094467608218440864316379601083733031616182834568622131929408727869368383220004732952196859616620142353086856831963914376824878744719370412908741315266076682750610260111368819902476050625917091851717151836815552335134401154279920499829154208352742566088202959102719466079133583834864505077458686820081578080489761914672893268837626e742"),
           e_float("-2.8393694847632029317769931319985456543467741590611966585875320007784497397890513502730857912898905313393941618066336189075910281045452849926091201840379454443047208860415679609919489974387946825982281619940211542425012799667657352186507006952091148830885284670171308330943865467984366318129976641123027947237108477519727620578989190351726099606958191542781347296487818563662760846764415516622001560983e766"),
           e_float("2.906521128223345839274838644343293460141781007086153757250387052639712492717724218909276139829054008705786159227281078056342467273714654840123020311632703281011267978419397071630994975368207024797466867142677788112633438613449906486765372025412893331518415756573407426341894396127273961282659185196837209012791004962059724468099888809452127762811155812671844262747789886818518668516417279532060905529e790"),
           e_float("-3.8192170879386787347468341086551432694793103144941323877466353504333365732730837435191628661237411199173058733675533744002478946650801468165655876083184195352454485039783957156770639606170428855289594919892111621400995661030396523677002368642739601094439914051293219109585402144001201683031578078579685135055274615332109830188228246157892184809757433338675726067510204907857107886318117580058346170676e814"),
           e_float("6.4029790761230149924075936237736684048730725436623608271373986472168888169216739372951431793598543416949823837258967183947456136330109869888406427717231463502919467597692262427913374075720537155976828924552725384639153812729974925819386834220860836044189914538579384739161498150525439402019908865383001262692652513716134609709148360830368359872111541252383260213170091667474234286554407517456762018725e838"),
           e_float("-1.3616955056208649494060799793833395526592966533813501489086788340453896582826148367962737159792626077014212151102847115476103661574588113989020771396799202650900954319364383178112422031897993872356014986948657109237496945277101296887111832949026639297825168560150450647236044228708173732789362743718563872708859948168684059261817109158132289862706544856593747869818832864745915695158386727532404385477e863"),
           e_float("3.6531572517423857947105318562288112934532198987563650422093650093410111527644103243668833132741027773950275462888517444044753902667302832214112368119882277460123633590071914626769508440259429103691179511900778866196149939910721599106911211643532019675274660495713701166460077587452744058495311432689030402655573847876473937871191463101837130891393302021087038011423864226925839905832887564844661709208e887"),
           e_float("-1.2298620863167756868471197569688423845089029927846646557905515461111913551620316404822539907496033743619218416528596534764993963223410405574194374968557970057434696631747956630550882641504645348448216911455337438747178593487657285323920341534302615091318440490982599962585672149801527736684507715218482810065993110003742654867620739192868975669967411916891573631800575992874678743257481149617967053693e912"),
           e_float("5.169578200820360663499665315884023783235452093616030232772462390628384219264771191541008407116634738333127243508171116699387520822704526435060944097313826634004721750158267425359552602594693784826748665923378165810136858094740632868323205154041005307970292136641066003576072732658806101689891759665968194549658702268634778417751158672786193782420199696160431402543060659952297318539096483070802796733e936"),
           e_float("-2.7000419728037337047009740797333405049756099152259438974531672850289106328313297954998868856632657052961811223861290612449022093629985208532550961787545770013776211308854084042820505928899045674422874136215654161668260386069408801281338595430460221727390857356587202033666242209445405446343596790146004273814200558961274979342642000772666217191628049200649914792995158670676389508659041137496122514002e961"),
           e_float("1.7442057830155690253675820276132845313258940162135301336372556216878258878275441337004023294815531859032387174843274474182656175962361254656058934263926486403622783733506587059895163031214611500283462467391860221306791988262403665236732808811385835351494961310891490181763075915750633477708444485842831010859697653494796159716772318612725020235673925497656499110288505678924567647553227528820324806098e986"),
           e_float("-1.3874434436959243383355352218858253810974360195779280507594571404076023827740658264733405288812459304050533986029417065084369869353614483627100917823356554284755044867889006256624639521781308024304947855215163825884671254883758118226762755664896942523384354924838185647185793030036099753496902755090184393044112660426964418748444754740091814559944544413437804216528905260218907174904774833001970518749e1011"),
           e_float("1.3532565526862922976467457455320763429600201782554378467845576173227233566229358463036131554687947193504772041827869228682774077856472175498097785446450997143502718438106477343492793342065285778253844352091148138517542246863498766491489435681529508138639781203475993798284348447568031556580565245237959521498108537968995747899916393486616171274840059626122110926261930554679525925135353230806503757267e1036"),
           e_float("-1.6118434469171764007192875818393917359078321286821804351839758727224771919451090700064638457330657765096697678674503618056898863840388535702542289156659906903370715763347144810541327633316972073956849715434803063861822096556135621504128397806162915958803077371187060949327918575859950046441655552773725416781583269694114763645916338194565450081168672664810841709283627741537961415109052003646804513085e1061"),
           e_float("2.3353037022095918685605296719776995156779388172557229547258694311939863978192251489623266747541829456105249611045643896173104176456760915946676049960815755217691821296211681025017499208743374324930798975608302686389340094768494646199303444670894384096808009746691282570487982047278875284608887059619891497361054033930728800085927651977412092188431718369291957693885103256481849233855408405722169274215e1086"),
           e_float("-4.1002177994111345215451650884218470203607573736745100085328295894067395903776814854298319647118565480570165001140437569051307420796953765052326785483034124801510369767536604986656792786351418920854216419718865811298944017098028088168121655334657191864345971315460899038010868255861987559065982015667967823779130658932805075702530141633409803686037875513727570388515007513584598796171389962018798991867e1111"),
           e_float("8.6924146250813954086227299312583439954947952608355223895320145803030911697213925998866911349584983990676704483523586452790202956078087009455253343761402598026422577280272228190314391783753667356327132514191658658140246797567309693730995240769346689920916449904642880939207321625934430865685414413399317677556605663172204939600453673880185400180791873441247802551925974907943504466981292725292373804076e1136"),
           e_float("-2.2173232364264121214148271516676494160059348744662735558546147832292437972284782394891210362208234265289613600836634558104733363909662305559987843095337274548346475445677615585731884363038498356160177298060900157148022686844813318952643773387874190563797497399809044048020720253364818985117189523761323471212090228847143242186896227952646152955248352186996641502421777548160256523959450426801717937426e1162"),
           e_float("6.7828593912647856866560693526640793323157251402235869176297855547213347008759212909688567128027197055847067230486962083076761047848973728075303722148219152643877294628571981927060066695875973512220317610573639891843825895545900825572965516040006346813815515051901627265141393972889583463681579139572556039074445038960488805691576364743458878309427705592856622795460294983392228953627501931607733911948e1187"),
           e_float("-2.48018522090742489825600695394978078232359467376479164623323781148642172848306173997571493106534006160784267134994608432648709967149148920390488287873083410371979403572086607553258138811196627720693237062218299422377804502801051633875766907636779089665129502624002494207709857301672414275736505500636086751315556023040949872677539389061461961502437683928771743525124569062073975935156568168601112713e1213"),
           e_float("1.0806486950160349518989765066030891733585396425176883578077419743671092049256585228350142444575761511120029347609368641469131261446524427679664846470762630386445188739698105743763432230644049940059903569539649400315807987284980511548571382883789004562894602566781122728178472307343772394039393340136414290169478018361238587436856768293360476304388699646913174353632216475238962324942803639768632901545e1239"),
           e_float("-5.5937428535283503906870036089644137034039001838750466747411328578904193790291481165138156673780807947068806906769025094571328182031296155590002876206075846200066450490818818655945392367794985934077479798046764791314078079513031851855705914197478289440240315239164926833654665366220835206079626889283530451820148191102816958588339230306476586705319819275293474474744555541362722110609734594772413812391e1264"),
           e_float("3.4298163540639137017165890976698098555639956990728738593323244854836829275962485041630125282588264368052749525315539439815531589590463691448937755112966754488829332696459141505487454835810660102857051107032399251399266497266799229532555180451286769108260195632270386732469817691082940135816825041231490184244460779957724479514036253868554024746158012057852718397089617366673206940712487865343465012333e1290"),
           e_float("-2.4840750188105643163754219060740751527029826369545700058332172114699149721822122215373513744774203607844960880969950321232474804743502450136419412557195338481755849685066180795065433147981725774593804996367826264272621221016007502892173136756945910284723955755869950564584837765520141303921314261076136289708511023661052453737151442983731122991631765594809292547963993001276764063353914035668171327286e1316"),
           e_float("2.1193330982113064642352827061475705931605041503232761548338520478446604762449138949773901056492332565066765636258454131919731786389824367442729233256733422225470255874017173542335081740831949060606080737621487933227497113652026611812257380005377039452041657585485406695781683161144125985981742507556284988447085964468641368642565935652315595104742089920023487670696951138493045149883849201317058076891e1342"),
           e_float("-2.1243584650429235571239118568321704993858701246731550551073516066591136850421491962222452564897524894089787345318610097300340285209259109510363250983483420959994718230958118203333533944172357637949450230176937559170367982223527022462566766847418745055151877352760639425792285372220776966671519328822803163127379952603842338553276356502259254919182159229163360464754266455680040267654181106962740967045e1368"),
           e_float("2.4954002879901478814703754368044615678867182005091749366443344519883017084854459195066716208769721052670985614340770765982568028537336057747907956422202256145119312480267233694658492610425561041253793997447054795647665091992687186345807075354383187090293073858441715344854905987974084755263036909776797082294685995865553484375071877814857742793936469183888291931041641471616016677505207727662218181282e1394"),
           e_float("-3.4265803578415585218272498509927807655899966191860144880005687320373245188032464967382651037835374713848778900827018000785934352216085648842722903704059166304359713806315973155389104659616010712131386437730415027146823882957867844617747775235657677541801179460268191976062771557307435589015456993017240702865768928364035177088601446292107831637769576601064915575444720215394266928489629990386986865124e1420"),
           e_float("5.4871573892028205407615258862759216291547273047566082009495061137348183709852903647646638671539427066626125858260554341660232226975814731439797121042668850505902717949447143346025508104234728803361122924785120749257102229123795973352195751130102174636933612745514765106635058616709884128170613786757585145618274378779801309833732468867817335040617258901975134179136809931724356399527256812793724913987e1446"),
           e_float("-1.0223245618508421174649787586072326054420678960748161179052062984920027001604954419875009614272467193189606345965171727809972159420550916278785624972064354858101841778494862881673139470650781564086004138845641976900063122025543342705976518226220961600583746076305047044191632407197526939038256773830526732431566189178603360555267147763818043264987762450238873969876823595332436051395303137842890539768e1473"),
           e_float("2.211079189478726541698113620562631965858613071067937771064928160745071301184178564347620626996681344631269518791606987978527018251719943983579638197533778580470238893280690757310157896594967426483021641161892394521498450874434404622226259771418860296803182271447434981815690274365754677982930525183378493133477842242283377975758984715705536210892156049708205191535630811064658717920582762085104170347e1499"),
           e_float("-5.539135560459639292633845218186167746342997127746478413336838108510396329201042469802408792303876852479251402600473120656541527717714670235182177568216502572423679605562408062868379846085248882524072749657762691283246705322172821556804811605043111853897717770949946816368462715525285190959804295235340283702905040959148691258561755458083285669918763542680522840260874299731722400041433773097934880437e1525"),
           e_float("1.6039009001932716708669024528780490236648348018338686395780665491090152417269690535655334908479495492061027625430229157586545147090643434912275650854936317710987677250068157090605600845336234513049605133273128492661760039073080307751125973680031591472409965128876064725600769567486649443496510532405428714193844926288389477229729798827608406728136774251778650657379530953525876695586510687180402486976e1552"),
           e_float("-5.3568927135833829582795875868900318855730191746960501422041850680975626394802825698236265795447262998596710778817609924338266606985446532029084694042846901827860230125969537768081933275937827276715393465048399088032088615988639639331245375385784345566484048261223987726661664837053246081844397628329467307273453738285376827112445772097375750545968362613675396262675840127203208351257886172958594994051e1578"),
           e_float("2.0595652948904640091127460992219149296721431866001447059522723120017662552042844542441593900327728590932796614023799748143237674378147957080803026107039239147404421345279659976521781853422634192186009808602760407811556719953108035457617418696172857174818592352734844875435396512987153436346826164190999919534466557705902465946575642700557259356168111570457676459634929708776114081293762919154466063555e1605"),
           e_float("-9.0973923742505389022264764799792899594839747234322342576911798297395366091909301388569284571512555112578714457704023524597841730565305374917490155597969413911659391548960108501854881929592004501818202016909938940909405120895168061013812343019757936135041509478800148920148879918352027652637977510828957842624358398044396711358084322265458563939097013795667508038852893062120159063648747345962812217139e1631"),
           e_float("4.6079929918894462469545731133143364388465854948931651717882667991838185733908439588656485177048874457160121215940530311312389779671661869912840161874981906489517637465456959192828816901746224353269242673180446414742921080666706745631292140510657384321511572410174120594534903118072221709165079536793065162626101078034172840746585978355000189125214490516606613861133350781726623681067020182724565187032e1658"),
           e_float("-2.6715129464894282293206938524854969823117583331825507558088027901276849638251602799770551930366114691884036309778171134383763998861600735287217642098394204001299007294783780390558758934693015444155599355050674159491287396163867061720486015910979867560636626204885816934487350604452592683303093375436676695943093799050625396357541445652158386070448514798707866048019072087778819795157659533495921701858e1685"),
           e_float("1.7695865673365894593680694547641429556687735428425606284165099727400650529323691302563251569097729360978492788386745927068909566339069808493137248960917205061522945660310652219341707505716772400413383392399270996615877462909347222389359174186657147966970065424721432927350254933186040042936501320517620049609029872312907089216293772529291217545926412055799520642824407098728821527970062335018850213084e1712"),
           e_float("-1.3368861407240290895374458788425699751615006635484963085553016579943631774124757345258458638152040532696818345666585674163116881830619586322717659610548610200342710557684681941968760598515444710939209407705298169230323888629864597968283071307930447088965995331516191122420701417015309848720656677669025402873388319516171187888972054461011622169532414969332501700951235450488872300697124653253444958825e1739"),
           e_float("1.1499630436424284306942320077755088114119461068018681348705402011511630819061816623786786218553051574709172668919222538865844876845305054122733956915700853085583088663011674918486011205506601792828320003153838637484853277722335543029780207147183229428245210446567403304749535613199258452675601323702901928327538820947775671539215894368143643322134417062174384143251894947883601427526089469781641201194e1766"),
           e_float("-1.1243925830679534791400170289383855205865492724974847022463079215216916392043290656120855361100642839815963615239901741477486143700875148493694499260980388870071250293440130538627744444362891199707107005142081498122534479434221241223665325924504105272316190087121650620023232872653606423758961274629546232507909655718784566591674663829195073745323324541134677909201241027334287118716565939786850470588e1793"),
           e_float("1.2476498885711072597461559338954936716273317385815918390300494494368746808924986453077888482778045272450549933902163359182269087407491968062777123591668691253633060466626426994088351395757676592007957929374675015171862873574345978922444606072009185340860196050425542760373539405743614444882815965933001978787967767408166516486242309659695197004107324577371636000822755825187616797754656288275646302463e1820"),
           e_float("-1.5686346296692016075481229340051168789635106611247675285107333483207982006853162710079490029078113447278860859438388142145907618579903394011872726596900155033186786933561970058015925812710976040623791372657790964618195428230006260749438666356425661709566516725394700902166484169473477738684365248236260026576710579760886512148773413020499316633472888325832470757402397453492410787211483623354304256034e1847"),
           e_float("2.2311858339214610243130956198341392597979433203804001347478542305784448142787193094500129053376113361061809705107227398682848731778569123511824530165197282157746349949232352006182946475306498905107371694520935746331543194023739527852160222652711753084182052159116243058319849027188590561122889691735677752993552306644761333822865163163307602817641822500065967228972619318775673003951663460145306336295e1874"),
           e_float("-3.584936441117142342646184324716790156504078189839231855589377459300125799418230682712149452783406643029881929067912898074955087946666794443585790556509871819291049981515521531403877008540153357746109400471211304029468407540419084571909102308497772596165098141475727073326559216904140871219191777918832950895073758749342254764136973060278704978954835942843352745959582721743546039020038185946351682502e1901"),
           e_float("6.4971301788319819573245809144844457037500238679435389737302500899439845583921921978317098682957704889810259068731603372639984854606519100706216934340529542611673489640941827162038041605191007084627631935768733716788895813456246575303023881473692890541687064080333590137197864306637189512037155665304343341699850767250029138457871825373077153874009661181257700678078815560804469527610468527751702933194e1928"),
           e_float("-1.3262758366630254966513949537182869491903189227340688334404289896859303409850885718913316601677071011818613214203160377981106180874929891102366962787163472255906440529816842415707070818133629947936894101292864612502283362926450861901341020048551045465522192408425205956149949813770222797515352952793947034700277508199034644337918674468833307430942702717348502314325036515066984445720614863180597265026e1956"),
           e_float("3.045164428135231270585874593816304541456362692103719181320504314571904995038517005853787991472353785583573498275197640818998368175589175538293450435997816424138446304846674427322865852345215021777686309666376142033075556083491943803151245352408048348042157772189952271659233010725501207973755290962338503595743725991917364982582243073087463406890042073573879324274698720751298906562360130205242047016e1983"),
           e_float("-7.8534192812382401317285809386895324973949362028794377125911709442000717620089721056993515071505912646372670871187141528412719439015105896822692670010760944466628110504132342002766639620994286471112318959807378130711487171750197545907599079367710178293911609902024270243894957187750774592662888967248325542037907620132942150141678023578288313450039167135420805831644823155985371643459258820430482772476e2010"),
           e_float("2.2719467388428790671352158729281669836748433005270784004737191372639898218980903823967844278711499499816872566394448787646832566100547404871788248402729771031705098252037901248359591163592975167291892864326778116223141429695742117158073414296928209555225288933959850048518553733447995916698249122992868454227272946909060449825351144218322602324085414917424704401931373150867307768146565568756255896538e2038"),
           e_float("-7.3631251749538768267282450220916958504534933118641687435413731786706323339404272742573778003366552165177947316724446055924386489454002598632559031110481508161259861925009354946565313839601644592557378138667409098627522530542584135020557502595624342085665942921379480392576945807351950847916164686545536220187548364569119293314241549854363406668571671387290268154972185240673787770601696782796165030545e2065"),
           e_float("2.6699121721860340358113356621447130423688368786902750240814077640694288337534486418116141061987926827688860004641798424510490543716108703427319420331109469773719360002641603819301196455503793971516360133400794906529069437475349655181319030065261065075895094419468389643006206233928237502520101790887313060063743754301842888454984892102442663995353830201644124835468006612567890846621227255903006693919e2093"),
           e_float("-1.0818346873951973608007469379888374183557895372475516763507148365475099980069821863093805536369849063920865445443700274318365227892102562773604839261920946467884401370944513638252918784431939950242787906248622510104873963973215062145039657464607476937229802738569452639240471898791712050686391377072309070462128006501117174416181797023472412796686063618027415124172613524060374718453380809147543189434e2121"),
           e_float("4.8924257570993999081714523774284131312354973319995027619574323604529509486770973673874033790041097527186211182156642644164816490484285028097119572558146978192432733524119998927613656548285502741681594402745347774969390062239501192323537068834557947593254456301382560741409550593639471767342145879524546853069181137886985139176453848182133810421080951765337049182527090635926507029342934109791707833053e2148"),
           e_float("-2.4664299729709884685104286330706073996347689691066910638599856178613214124219718573235280844973133756749155862214946680439198210281523319311066120936811233331209828623520288436181818897252195789248111424710355109735314135681617295148320672856357335439341230542595072397452220871908499386982653679459798417688024726905081603677286499379924149795696720412021734627000016162155166682291298277464427809016e2176"),
           e_float("1.3844824892310532289238590098715603019087906184452440288514782680921510685490279515059358934112457116956988913964707400729816396921516939414293816489897035341745344808683945587736112534224091128234992335132678187764142941276880333036557529512825345198594363610965829042233933453488816924317281607061632237601960350412770547486267133892863828896689652519619109276219888375236016869376798325648685337437e2204"),
           e_float("-8.6433876079683325310348817321381887265729143117222386014229145793340115272957371895070233233200344458306173131186879462635230609273077602781569242460396795654766963195595458736735980914956620160290732851846605697872762757618982646105047460605590278445266275550039572394838708611929688842940518716985436704919904262810885823870905781019115286156473708467406647369699432480839101566129636173573919994322e2231"),
           e_float("5.9947697949800996331266307418079436085383788141739929874716955009062637555622807394749347628975369525105343690907482256852506288635451836244394823961081180237192687476978432618600318820747395573319901582599508819267882343816739954229137151362913284607102168597646230255617189984780615453987066686469716317132642707395424655275955175627662781207717152921305967679491708852818764430251957102121284994337e2259"),
           e_float("-4.6139969523347149452592865374129224592275199567004082791687481934146619315478502536269753737678048193722018029091845910599998989317651233876361544173908739548164107426831094425402334704611169175203389981572138150191546930921171394498063943323212189055522711180914588355092916460752112538572348839104267545912508637927900635943731898027170488788143541383911177305952572713238344605346594606340511622686e2287"),
           e_float("3.9367011253295967975379541454675241771617830882290464046644997960610451705100780497550332927052611063371954848430712528584369662328446272844159491128540709105840791542926615940928621379776285223369952552654625827910300272038387241571169896447435144272720349350569050676608204929401382038500088188217699997244913606048624539940326372013697109178631354881903392251630381864508896961982924640949303214791e2315"),
           e_float("-3.7194738899880770523041913140237029823650892754796213991262257575474327868103416716197569395192367840414308952108146524977947269965955676069618478087310894575754068032696440201232840907284902776918795783095838009693817752042688265366258723120890299703071846625778848545453863447257992462406169423267594087119254818494963740881700366425303247367427185124301132942348049635378836398950066506860693316512e2343"),
           e_float("3.8875618412530706152573358625276255615246331032580348669343051402663220978440828770504159149316049043004153699954965411433169176331319280903231644370432146294956433595873706677392090821991724972792901547407942099273406642444453559681941409552739964098984841270729763245921788073846424409401047153726665557977779902588645947770626410590395989841863017326348632823617405225454796925245391582908750003719e2371"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00022_euler(const bool b_write_output)
    {
      return TestCase_case_00022_euler().execute(b_write_output);
    }
  }
}