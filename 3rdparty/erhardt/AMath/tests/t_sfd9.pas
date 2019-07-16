{Part 9 of regression test for SPECFUN unit  (c) 2013-2018  W.Ehrhardt}

unit t_sfd9;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_hyperg_1F1;
procedure test_hyperg_1F1r;
procedure test_hyperg_2F1;
procedure test_hyperg_2F1r;
procedure test_hyperg_2F0;


implementation

uses
  AMath,SpecFun,t_sfd0;


{---------------------------------------------------------------------------}
procedure test_hyperg_1F1;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
  NE1 = 4;
  NE2 = 20;
  NE3 = 80;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hyperg_1F1');

  {GSL integer a,b}
  y := hyperg_1F1(1, 1, 0.5);
  f := 1.6487212707001281468;
  testrel(1, NE, y, f, cnt,failed);

  y := hyperg_1F1(1, 2, 500.0);
  f := 2.8071844357056748215e+214;
  testrel(2, NE, y, f, cnt,failed);

  y := hyperg_1F1(1, 2, -500.0);
  f := 0.002;
  testrel(3, NE, y, f, cnt,failed);

  y := hyperg_1F1(8, 1, 0.5);
  f := 13.108875178030540372;
  testrel(4, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 1, 1.0);
  f := 131.63017574352619931;
  testrel(5, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 1, 10.0);
  f := 8.514625476546280796e+09;
  testrel(6, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 1, 100.0);
  f := 1.5671363646800353320e+56;
  testrel(7, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 20, 1.0);
  f := 1.6585618002669675465;
  testrel(8, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 20, 10.0);
  f := 265.26686430340188871;
  testrel(9, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 20, 100.0);
  f := 3.640477355063227129e+34;
  testrel(10, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 100, 1.0);
  f := 1.1056660194025527099;
  testrel(11, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 100, 10.0);
  f := 2.8491063634727594206;
  testrel(12, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 100, 100.0);
  f := 8.032171336754168282e+07;
  testrel(13, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 100, 500.0);
  {f := 7.6339612025287314263e+123;}
  f := 0.763396120252873142563e124; {Maple}
  testrel(14, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 1, 1.0);
  f := 6.892842729046469965e+07;
  testrel(15, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 1, 10.0);
  f := 2.4175917112200409098e+28;
  testrel(16, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 1, 100.0);
  f := 1.9303216896309102993e+110;
  testrel(17, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 200, 1.0);
  f := 1.6497469106162459226;
  testrel(18, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 200, 10.0);
  f := 157.93286197349321981;
  testrel(19, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 200, 100.0);
  f := 2.1819577501255075240e+24;
  testrel(20, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 200, 400.0);
  f := 3.728975529926573300e+119;
  testrel(21, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 400, 10.0);
  f := 12.473087623658878813;
  testrel(22, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 400, 100.0);
  f := 9.071230376818550241e+11;
  testrel(23, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 400, 600.0);
  f := 2.4076076354888886030e+112;
  testrel(24, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 1, -1.0);
  f := 0.11394854824644542810;
  testrel(25, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 1, -10.0);
  f := 0.0006715506365396127863;
  testrel(26, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 1, -100.0);
  f := -4.208138537480269868e-32;
  testrel(27, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 50, -1.0);
  {f:= 0.820006196079380;}
  f := 0.8200061960793795459;  {Maple}
  testrel(28, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 100, -10.0);
  f := 0.38378859043466243;
  testrel(29, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 100, -100.0);
  f := 0.0008460143401464189061;
  testrel(30, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 100, -500.0);
  f := 1.1090822141973655929e-08;
  testrel(31, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 100, -10000.0);
  f := 5.173783508088272292e-21;
  testrel(32, NE, y, f, cnt,failed);

  y := hyperg_1F1(50, 1, -100.0);
  f := 4.069661775122048204e-24;
  testrel(33, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 10, -100.0);
  f := -2.7819353611733941962e-37;
  testrel(34, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 1, -90.0);
  f := 7.501705041159802854e-22;
  testrel(35, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 1, -110.0);
  f := -7.007122115422439755e-26;
  testrel(36, NE, y, f, cnt,failed);

  y := hyperg_1F1(200, 50, -1.0);
  f := 0.016087060191732290813;
  testrel(37, NE, y, f, cnt,failed);

  y := hyperg_1F1(200, 50, -300.0);
  f := -4.294975979706421471e-121;
  testrel(38, NE, y, f, cnt,failed);

  y := hyperg_1F1(200, 100, -1.0);
  f := 0.13397521083325179687;
  testrel(39, NE, y, f, cnt,failed);

  y := hyperg_1F1(200, 100, -10.0);
  f := 5.835134393749807387e-10;
  testrel(40, NE, y, f, cnt,failed);

  y := hyperg_1F1(200, 100, -100.0);
  f := 4.888460453078914804e-74;
  testrel(41, NE, y, f, cnt,failed);

  y := hyperg_1F1(200, 100, -500.0);
  f := -1.4478509059582015053e-195;
  testrel(42, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1, 1, 2.0);
  f := -1.0;
  testrel(43, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1, -2, 2.0);
  f := 2.0;
  testrel(44, NE, y, f, cnt,failed);

  y := hyperg_1F1(-2, -3, 2.0);
  f := 3.0;
  testrel(45, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, 1, 1.0);
  f := 0.4189459325396825397;
  testrel(46, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, 1, 10.0);
  f := 27.984126984126984127;
  testrel(47, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, 1, 100.0);
  f := 9.051283795429571429e+12;
  testrel(48, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, 20, 1.0);
  f := 0.20203016320697069566e-2;
  testrel(49, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -20, 1.0);
  f := 1.6379141878548080173;
  testrel(50, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -20, 10.0);
  f := 78.65202404521289970;
  testrel(51, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -20, 100.0);
  f := 4.416169713262624315e+08;
  testrel(52, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -100, 1.0);
  f := 1.1046713999681950919;
  testrel(53, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -100, 10.0);
  f := 2.6035952191039006838;
  testrel(54, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -100, 100.0);
  f := 1151.6852040836932392;
  testrel(55, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, -200, 1.0);
  f := 1.6476859702535324743;
  testrel(56, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, -200, 10.0);
  f := 139.38026829540687270;
  testrel(57, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, -200, 100.0);
  f := 1.1669433576237933752e+19;
  testrel(58, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -20, -1.0);
  f := 0.6025549561148035735;
  testrel(59, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -20, -10.0);
  {GSL 0.00357079636732993491:  err = 0.34e-3 !!!!}
  f := 0.3572011827945780474e-2; {Maple}
  testrel(60, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -20, -100.0);
  {f := 1.64284868563391159e-35; ??? Totally wrong}
  f := 0.4930327226240536973e8; {Maple}
  testrel(61, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -100, -1.0);
  f := 0.9044239725031389133; {Maple}
  testrel(62, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -100, -10.0);
  f := 0.35061515251367213;
  testrel(63, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, -200, -1.0);
  f := 0.6061497939628952629;
  testrel(64, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, -200, -10.0);
  f := 0.6327854390887766225e-2;   {Maple}
  testrel(65, NE, y, f, cnt,failed);

  {No GSL case, GSL case with x=-100 is a problem case after}
  {disallowing Kummer transformation for negative integer b.}
  y := hyperg_1F1(-100, -200, -50.0);
  f := 0.2923923054547819835e-11;  {Maple}
  testrel(66, NE, y, f, cnt,failed);

  {GSL general a,b}
  y := hyperg_1F1(1, 1.5, 1);
  f := 2.0300784692787049755;
  testrel(67, NE, y, f, cnt,failed);

  y := hyperg_1F1(1, 1.5, 10);
  f := 6172.859561078406855;
  testrel(68, NE, y, f, cnt,failed);

  y := hyperg_1F1(1, 1.5, 500);
  f := 5.562895351723513581e+215;
  testrel(69, NE, y, f, cnt,failed);

  y := hyperg_1F1(1.5, 2.5, 1);
  f := 1.8834451238277954398;
  testrel(70, NE, y, f, cnt,failed);

  y := hyperg_1F1(1.5, 2.5, 10);
  f := 3128.7352996840916381;
  testrel(71, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 1.1, 1);
  f := 110.17623733873889579;
  testrel(72, NE1, y, f, cnt,failed);

  y := hyperg_1F1(10, 1.1, 10);
  f := 6.146657975268385438e+09;
  testrel(73, NE1, y, f, cnt,failed);

  y := hyperg_1F1(10, 1.1, 100);
  f := 9.331833897230312331e+55;
  testrel(74, NE1, y, f, cnt,failed);

  y := hyperg_1F1(10, 50.1, 2);
  f := 1.5001295507968071788;
  testrel(75, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 50.1, 10);
  f := 8.713385849265044908;
  testrel(76, NE, y, f, cnt,failed);

  y := hyperg_1F1(10, 50.1, 100);
  f := 5.909423932273380330e+18;
  testrel(77, NE2, y, f, cnt,failed);

  y := hyperg_1F1(10, 50.1, 500);
  f := 9.740060618457198900e+165;
  testrel(78, NE2, y, f, cnt,failed);

  y := hyperg_1F1(100, 1.1, 1);
  f := 5.183531067116809033e+07;
  testrel(79, NE1, y, f, cnt,failed);

  y := hyperg_1F1(100, 1.1, 10);
  f := 1.6032649110096979462e+28;
  testrel(80, NE1, y, f, cnt,failed);

  y := hyperg_1F1(100, 1.1, 100);
  f := 1.1045151213192280064e+110;
  testrel(81, NE1, y, f, cnt,failed);

  y := hyperg_1F1(100, 50.1, 1);
  f := 7.222953133216603757;
  testrel(82, NE, y, f, cnt,failed);

  y := hyperg_1F1(100, 50.1, 10);
  f := 1.0998696410887171538e+08;
  testrel(83, NE1, y, f, cnt,failed);

  y := hyperg_1F1(100, 50.1, 100);
  f := 7.235304862322283251e+63;
  testrel(84, NE2, y, f, cnt,failed);

  y := hyperg_1F1(1, 1.5, -1);
  f := 0.5380795069127684191;
  testrel(85, NE, y, f, cnt,failed);

  y := hyperg_1F1(1, 1.5, -10);
  f := 0.05303758099290164485;
  testrel(86, NE, y, f, cnt,failed);

  y := hyperg_1F1(1, 1.5, -500);
  f := 0.0010010030151059555322;
  testrel(87, NE, y, f, cnt,failed);

  y := hyperg_1F1(1, 1.1, -500);
  f := 0.00020036137599690208265;
  testrel(88, NE2, y, f, cnt,failed);

  y := hyperg_1F1(10, 1.1, -1);
  f := 0.07227645648935938168;
  testrel(89, NE1, y, f, cnt,failed);

  y := hyperg_1F1(10, 1.1, -10);
  f := 0.0003192415409695588126;
  testrel(90, NE1, y, f, cnt,failed);

  {End of regular GSL, start fixes/hacks/bugs}

  y := hyperg_1F1(100, 1.1, -1);
  f := 0.0811637344096042096;
  testrel(91, NE1, y, f, cnt,failed);

  y := hyperg_1F1(100, 1.1, -10);
  f := 0.00025945610092231574387;
  testrel(92, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-10, -10.1, 10.0);
  f := 10959.603204633058116;
  testrel(93, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-10, -10.1, 1000.0);
  f := 2.0942691895502242831e+23;
  testrel(94, NE2, y, f, cnt,failed);

  y := hyperg_1F1(-10, -100.1, 10.0);
  f := 2.6012036337980078062;
  testrel(95, NE, y, f, cnt,failed);

  y := hyperg_1F1(-8.1, -10.1, -10.0);
  f := 0.00018469685276347199258;
  testrel(96, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10, -5.1, 1);
  f := 16.936141866089601635;
  testrel(97, NE2, y, f, cnt,failed);

  y := hyperg_1F1(-10, -5.1, 10);
  f := 771534.0349543820541;
  testrel(98, NE2, y, f, cnt,failed);

  y := hyperg_1F1(-10, -5.1, 100);
  f := 2.2733956505084964469e+17;
  testrel(99, NE2, y, f, cnt,failed);

  y := hyperg_1F1(-100, -50.1, -1);
  f := 0.13854540373629275583;
  testrel(100, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, -50.1, -10);
  f := -9.142260314353376284e+19;
  testrel(101, NE3, y, f, cnt,failed);

  y := hyperg_1F1(-10.5, -8.1, 0.1);
  f := 1.1387201443786421724;
  testrel(102, NE, y, f, cnt,failed);

  y := hyperg_1F1(-10.5, -11.1, 1);
  f := 2.5682766147138452362;
  testrel(103, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100.5, -80.1, 10);
  f := 355145.4517305220603;
  testrel(104, NE2, y, f, cnt,failed);

  y := hyperg_1F1(-100.5, -102.1, 10);
  f := 18678.558725244365016;
  testrel(105, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100.5, -500.1, 100);
  f := 1.2077443075367177662e+8;
  testrel(106, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100, 1.1, 1);
  f := 0.21519810496314438414;
  testrel(107, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, 1.1, 10);
  f := 8.196123715597869948;
  testrel(108, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100, 1.1, 100);
  f := -1.4612966715976530293e+20;
  testrel(109, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100, 20.1, 1);
  f := 0.0021267655527278456412;
  testrel(110, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100, 20.1, 10);
  f := 2.0908665169032186979e-11;
  testrel(111, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100, 20.1, 100);
  f := -0.04159447537001340412;
  testrel(112, NE3, y, f, cnt,failed);

  y := hyperg_1F1(-100, 1.1, -1);
  f := 2.1214770215694685282e+07;
  testrel(113, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, 1.1, -10);
  f := 1.0258848879387572642e+24;
  testrel(114, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100, 1.1, -100);
  f := 1.1811367147091759910e+67;
  testrel(115, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100, 50.1, -1);
  f := 6.965259317271427390;
  testrel(116, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100, 50.1, -10);
  f := 1.0690052487716998389e+07;
  testrel(117, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-100, 50.1, -100);
  f := 6.889644435777096248e+36;
  testrel(118, NE2, y, f, cnt,failed);

  y := hyperg_1F1(-2.05, 1.0, 5.05);
  {f := 3.79393389516785e+00;}
  f := 3.793933895167851280;  {Maple}
  testrel(119, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-26, 2.0, 100.0);
  {f := 1.444786781107436954e+19;}
  f := 1.444786781107436987e19; {Maple}
  testrel(120, NE, y, f, cnt,failed);

  y := hyperg_1F1(1.2, 1.1e-15, 1.5);
  f := 8254503159672429.02;
  testrel(121, NE, y, f, cnt,failed);


  y := hyperg_1F1(-1.5, 1.5, -100.0);
  f := 456.44010011787485545;
  testrel(122, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1.5, 1.5, 100.0);
  f := 1.0893724312430935129254e37;
  testrel(123, NE, y, f, cnt,failed);

  y := hyperg_1F1(1.0, 1000000.5, 800000.5);
  {f := 4.999922505099443804e+00;}
  f := 4.999922505099447552;  {Maple}
  testrel(124, NE, y, f, cnt,failed);

  y := hyperg_1F1(1.5, 1000000.5, 800000.5);
  f := 11.18001288977894650;
  testrel(125, NE, y, f, cnt,failed);

  {-------- end of GSL cases ------------}


  y := hyperg_1F1(-1,1.5,-7);
  f := 5.666666666666666667;
  testrel(126, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1.5,1.5,60);
  f := 0.2263598817446298582e21;
  testrel(127, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1.1,1.5,-701);
  f := 838.8838298222824863;
  testrel(128, NE1, y, f, cnt,failed);

  y := hyperg_1F1(100,150,50);
  f := 0.1637175264548271410e16;
  testrel(129, NE, y, f, cnt,failed);

  y := hyperg_1F1(100,250,50);
  f := 0.1646092714468226994e10;
  testrel(130, NE, y, f, cnt,failed);

  y := hyperg_1F1(2,400,500);
  f := 295912696.5840234021;
  testrel(131, NE, y, f, cnt,failed);

  y := hyperg_1F1(100,400,500);
  f := 0.1140128409643727474e88;
  testrel(132, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1000,-2000,0.5);
  f := 1.284005343911193494;
  testrel(133, NE, y, f, cnt,failed);

  y := hyperg_1F1(-900,-1800,0.5);
  f := 1.284003112382361919;
  testrel(134, NE, y, f, cnt,failed);

  y := hyperg_1F1(-100,6.8,1.2);
  f := 0.0001121335924327979849; {Wolfram}
  testrel(135, NE1, y, f, cnt,failed);

  y := hyperg_1F1(-1000,6.8,1.2);
  f := -1.009676652981919674e-7; {Wolfram}
  testrel(136, NE1, y, f, cnt,failed);

  y := hyperg_1F1(2,1,3);
  f := 80.34214769275067096;
  testrel(137, NE, y, f, cnt,failed);

  y := hyperg_1F1(-3.5,-4.5,3);
  f := 6.695178974395889247;
  testrel(138, NE, y, f, cnt,failed);

  y := hyperg_1F1(3.5,2.5,-3);
  f := -0.9957413673572788596e-2;
  testrel(139, NE, y, f, cnt,failed);

  f := 15957.87215834274658;
  y := hyperg_1F1(-10,-20,-50);
  testrel(140, NE, y, f, cnt,failed);

  f := 0.1853807070382061169e-15;
  y := hyperg_1F1(-10.5,-20.5,-50);
  testrel(141, NE, y, f, cnt,failed);

  f := 1.488987369555254432e25; {wolfram}
  y := hyperg_1F1(-10,-20,5000);
  testrel(142, NE, y, f, cnt,failed);

  f := 16674978349.498007760;
  y := hyperg_1F1(0.5,1e-10,1.25);
  testrel(143, NE, y, f, cnt,failed);

  f := 1667497834983136.1463;
  y := hyperg_1F1(0.5,1e-15,1.25);
  testrel(144, NE, y, f, cnt,failed);

  f := 7/3;
  y := hyperg_1F1(-2,-4,2);
  testrel(145, NE, y, f, cnt,failed);

  y := hyperg_1F1(10,8000.0,11500.0);
  f := 0.1081895784081134874e289;
  testrel(146, NE, y, f, cnt,failed);

  {b = 2a case}
  y := hyperg_1F1(2000,4000,10);
  f := 148.8775590092748498;
  testrel(147, NE, y, f, cnt,failed);

  y := hyperg_1F1(0.125,0.25,5);
  f := 60.52162276513522068;
  testrel(148, NE, y, f, cnt,failed);

  y := hyperg_1F1(-0.125,-0.25,5);
  f := 98.89577758359759834;
  testrel(149, NE, y, f, cnt,failed);

  y := hyperg_1F1(0.125,0.25,-5);
  f := 0.4077914864901255494;
  testrel(150, NE, y, f, cnt,failed);

  y := hyperg_1F1(-0.125,-0.25,-5);
  f := 0.6663545077916252443;
  testrel(151, NE, y, f, cnt,failed);

  y := hyperg_1F1(1.25, 2.5, 2);
  f := 3.124656439439520633;
  testrel(152, NE, y, f, cnt,failed);

  y := hyperg_1F1(1.25, 2.5, -2);
  f := 0.4228762642486532576;
  testrel(153, NE, y, f, cnt,failed);

  f := 0.4074857472142321668e303;
  y := hyperg_1F1(-1500.125,-3000.25, 1600);
  testrel(154, NE, y, f, cnt,failed);

  f := 0.1784570650015316724e-1;
  y := hyperg_1F1(0.5,1,-1000);
  testrel(155, NE, y, f, cnt,failed);

  {large |a|, some problems if a*x < 0}
  y := hyperg_1F1(-1000.5, 1, 2);
  f := 0.1710838395033404881;
  testrel(156, NE, y, f, cnt,failed);

  y := hyperg_1F1(1000.5, 1, 2);
  f := 0.8085466047849697461e38;
  testrel(157, NE, y, f, cnt,failed);

  y := hyperg_1F1(1000.5, 1, -2);
  f := 0.2406013532478147622e-1;
  testrel(158, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1000.5, 1, -2);
  f := 0.1143992271759473033e38;
  testrel(159, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1000.5, -1.5, 2);
  f := -563.1980796259134972;
  testrel(160, NE, y, f, cnt,failed);

  y := hyperg_1F1(1000.5, -1.5, 2);
  f := 0.2612545595507873052e43;
  testrel(161, NE, y, f, cnt,failed);

  y := hyperg_1F1(1000.5, -1.5, -2);
  f := -10.53625857198753570;
  testrel(162, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1000.5, -1.5, -2);
  f := 0.3301330588826002370e42;
  testrel(163, NE, y, f, cnt,failed);

  y := hyperg_1F1(-1000, -1000.1, 200.0);
  f := 7.066514294896245043e+86;
  testrel(164, NE3, y, f, cnt,failed);

  {Related to igammat}
  y := hyperg_1F1(1,-9.5, -234.5);
  f := -0.4709619719604298796e-1;
  testrel(165, NE, y, f, cnt,failed);

  {GSL, 'fixed' in AMath 1.29.01}
  y := hyperg_1F1(-1000, -1000.1, 10.0);
  f := 22004.341698908631636;
  testrel(166, NE1, y, f, cnt,failed);

  {very small a and b}
  y := hyperg_1F1(1e-20, 1e-22, 1.0);
  f := 172.8281828459045235;
  testrel(167, NE, y, f, cnt,failed);

  y := hyperg_1F1(1e-19, -1e-20, -2.5);
  f := 10.17915001376101205;
  testrel(168, NE, y, f, cnt,failed);


(*
  ----------------------------------------
  Problem test cases from GSL:
  ----------------------------------------

  Try to fix in next version(s) !!

  y := hyperg_1F1(100, 1.1, -100);
  f := 1.8315497474210138602e-24;
  writeln(1-y/f);

  y := hyperg_1F1(10, 1.1, -500);
  f := -3.400379216707701408e-23;    {!!!!!!!!!}
  writeln(1-y/f);

  y := hyperg_1F1(50, 1.1, -100);
  f := 4.632883869540640460e-24;     {!!!!!}
  writeln(1-y/f);

  y := hyperg_1F1(50, 1.1, -110.0);     {!!!!!}
  f := 5.642684651305310023e-26;
  writeln(1-y/f);

  y := hyperg_1F1(-100, -50.1, 10);
  f := 1.0551632286359671976e+11;
  writeln(1-y/f);

  y := hyperg_1F1(-26.1, 2.0, 100.0);
  {f := 1.341557199575986995e+19;}
  f := 1.341557199575986910e19;
  writeln(1-y/f);

  {GSL: 'Horrible', AMath too}
  y := hyperg_1F1(100, 10.1, -220);
  f := -4.296130300021696573e-64;
  writeln(1-y/f);

  {??????? Maple/Alpha without Kummer,  GSL with Kummer; both ~ 1e-12}
  y := hyperg_1F1(-10, -100, -100.0);
  f := 4.414825020540336352e-7;  {Maple/Alpha}
  f := 8.19512187960476424e-09;
  writeln(1-y/f);
*)

  {GSL, fixed in AMath 1.29.01}
  {AMath uses h1f1_tricomi, DAMath h1f1_sum}
  y := hyperg_1F1(-1000, -1000.1, 10.0);
  f := 22004.341698908631636;
  testrel(166, NE1, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_hyperg_1F1r;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hyperg_1F1r');

  {functions.wolfram.com}
  y := hyperg_1F1r(-6, -4, -0.5);
  f := 0.203125;
  testrel(1, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-6, -4, 0.5);
  f := -0.171875;
  testrel(2, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-6, -4+ldexpd(1,-40), 0.5);
  f := -0.1718749999519580651;
  testrel(3, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-6, -6, 0.5);
  f := 0;
  testrel(4, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-6, -6-ldexpd(1,-32), 0.5);
  f := -2.763881640260390814e-7;
  testrel(5, NE, y, f, cnt,failed);

  y := hyperg_1F1r(5, -1, 0.5);
  f := 11.12994196346981039;
  testrel(6, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-5, -1, -0.5);
  f := 3.911458333333333333;
  testrel(7, NE, y, f, cnt,failed);

  y := hyperg_1F1r(5, -1, -0.5);
  f := 1.023915364996652648;
  testrel(8, NE, y, f, cnt,failed);

  y := hyperg_1F1r(5, -10, 0.5);
  f := 1.292845512309552162;
  testrel(9, NE, y, f, cnt,failed);

  y := hyperg_1F1r(5, -100, 0.5);
  f := 3.170010090411819724e-24;
  testrel(10, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-5.5, -10, -0.5);
  f := -4.676806175400637187e-8;
  testrel(11, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-5.5, -10, -8);
  f := -46733.02438628188571;
  testrel(12, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-5.5, -10, 4);
  f := 3675.897286648060413;
  testrel(13, NE, y, f, cnt,failed);

  y := hyperg_1F1r(5, -1, -2);
  f := -1.263129310208385124;
  testrel(14, NE, y, f, cnt,failed);

  y := hyperg_1F1r(5, -10, -2);
  f := -178957.2735176149567;
  testrel(15, NE, y, f, cnt,failed);

  y := hyperg_1F1r(5, -100, -2);
  f := -1.515230067091366973e36;
  testrel(16, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-1.5, -2, 2.5);
  f := 2.908235393264730177;
  testrel(17, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-1.5, 3.5, 2.5);
  f := 0.02726382546057175238;
  testrel(18, NE, y, f, cnt,failed);

  y := hyperg_1F1r(-1.5, -3.5, 2.5);
  f := 12.88731786424190386;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_hyperg_2F1;
var
  x,y,b,c,f: double;
  n, cnt, failed: integer;
const
  NE  = 1;
  NE1 = 25;
  NE2 = 100;
  NE3 = 200;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hyperg_2F1');

  {hyp2f1 := (a,b,c,x) -> hypergeom( [a,b],[c], x);}

  {Chu-Vandermonde NIST 15.4.24}
  b := 1.5;
  c := 0.25;
  for n:=1 to 20 do begin
    y := hyperg_2F1(-n,b,c,1);
    f := pochhammer(c-b,n)/pochhammer(c,n);
    testrel(n, NE, y, f, cnt,failed);
  end;

  {-- GSL -----------------------}
  y := hyperg_2F1(1, 1, 1, 0.5);
  f := 2.0;
  testrel(21, NE, y, f, cnt,failed);

  y := hyperg_2F1(8, 8, 1, 0.5);
  f := 12451584.0;
  testrel(22, NE, y, f, cnt,failed);

  y := hyperg_2F1(8, -8, 1, 0.5);
  f := 0.13671875;
  testrel(23, NE, y, f, cnt,failed);

  y := hyperg_2F1(8, -8.1, 1, 0.5);
  f := 0.1414738537889993042;
  testrel(24, NE1, y, f, cnt,failed);

  y := hyperg_2F1(8, -8, 1, -0.5);
  f := 4945.136718750000000;
  testrel(25, NE, y, f, cnt,failed);

  y := hyperg_2F1(8, -8, -5.5, 0.5);
  f := -906.6363636363636364;
  testrel(26, NE, y, f, cnt,failed);

  y := hyperg_2F1(8, -8, -5.5, -0.5);
  f := 24565.363636363636364;
  testrel(27, NE, y, f, cnt,failed);

  y := hyperg_2F1(8, 8, 1, -0.5);
  f := -0.006476312098196747669;
  testrel(28, NE, y, f, cnt,failed);

  y := hyperg_2F1(8, 8, 5, 0.5);
  f := 4205.714285714285714;
  testrel(29, NE, y, f, cnt,failed);

  y := hyperg_2F1(8, 8, 5, -0.5);
  f := 0.0028489656290296436616;
  testrel(30, NE, y, f, cnt,failed);

  y := hyperg_2F1(9, 9, 1, 0.99);
  f := 1.2363536673577259280e+38;
  testrel(31, NE2, y, f, cnt,failed);

  y := hyperg_2F1(9, 9, -1.5, 0.99);
  f := 3.796186436458346579e+46;
  testrel(32, NE2, y, f, cnt,failed);

  y := hyperg_2F1(9, 9, -1.5, -0.99);
  f := 0.14733409946001025146;
  testrel(33, NE, y, f, cnt,failed);

  y := hyperg_2F1(9, 9, -8.5, 0.99);
  f := -1.1301780432998743440e+65;
  testrel(34, NE3, y, f, cnt,failed);

  y := hyperg_2F1(9, 9, -8.5, -0.99);
  f := -0.885646260657534448e1;
  testrel(35, NE2, y, f, cnt,failed);

  y := hyperg_2F1(9, 9, -21.5, 0.99);
  f := 2.0712920991876073253e+95;
  testrel(36, NE3, y, f, cnt,failed);

  y := hyperg_2F1(9, 9, -21.5, -0.99);
  f := -74.30517015382249216;
  testrel(37, NE, y, f, cnt,failed);

  y := hyperg_2F1(9, 9, -100.5, 0.99);
  f := -3.186778061428268980e+262;
  testrele(38, 2e-13, y, f, cnt,failed);    {!!!}

  y := hyperg_2F1(9, 9, -100.5, -0.99);
  f := 2.4454358338375677520;
  testrel(39, NE, y, f, cnt,failed);

  y := hyperg_2F1(25, 25, 1, -0.5);
  f := -2.9995530823639545027e-06;
  testrele(40, NE3, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 0.5, 2.0, 1.0-1.0/64.0);
  f := 3.17175539044729373926;
  testrel(41, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 0.5, 2.0, 1.0-1.0/1024.0);
  f := 4.90784159359675398250;
  testrel(42, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 0.5, 2.0, 1.0-1.0/65536.0);
  f := 7.552266033399683914;
  testrel(43, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 0.5, 2.0, 1.0-1.0/16777216.0);
  f := 11.08235454026043830;
  testrel(44, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 0.5, 2.0, -1.0+1.0/1024.0);
  f := 0.7629109409099549745;
  testrel(45, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 0.5, 2.0, -1.0+1.0/65536.0);
  f := 0.7627621249088454244;
  testrel(46, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 0.5, 2.0, -1.0+1.0/1048576.0);
  f := 0.7627599110890647380;
  testrel(47, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 0.5, 3.0, 1.0);
  f := 1.697652726313550248;
  testrel(48, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, -4.2, 3.0, 1.0);
  f := 0.15583601560025710650;
  testrel(49, NE, y, f, cnt,failed);

  y := hyperg_2F1(-7.4, 0.7, -1.5, 1.0);
  f := -0.3447886695924658500;
  testrel(50, NE, y, f, cnt,failed);

  y := hyperg_2F1(0.1, -2.7, -1.5, 1.0);
  f := 1.059766766063610123;
  testrel(51, NE, y, f, cnt,failed);

  y := hyperg_2F1(0, -2, -4, 0.5);
  f := 1.0;
  testrel(52, NE, y, f, cnt,failed);

  y := hyperg_2F1(-10.34, 2.05, 3.05, 0.1725);
  f := 0.3104735522131300104;
  testrel(53, NE, y, f, cnt,failed);

  y := hyperg_2F1(-9.99999999999, 2.05, 3.05, 0.1725);
  f := 0.3214193463019748754;
  testrel(54, NE, y, f, cnt,failed);

  y := hyperg_2F1(11, -1, 11.0/2.0, 0.125 );
  f := 0.75;
  testrel(55, NE, y, f, cnt,failed);

  {--- Forrey --------------------------}
  x := 0.2;
  y := hyperg_2F1(0.5,1.0,1.5,-x*x);
  f := arctan(x)/x;
  testrel(56, NE, y, f, cnt,failed);

  x := 1;
  y := hyperg_2F1(0.5,1.0,1.5,-x*x);
  f := arctan(x)/x;
  testrel(57, NE, y, f, cnt,failed);

  x := 2;
  y := hyperg_2F1(0.5,1.0,1.5,-x*x);
  f := arctan(x)/x;
  testrel(58, NE, y, f, cnt,failed);

  y := hyperg_2F1(1, 2+1e-10, 3, -2);
  f := 0.4506938556510494928;
  testrel(59, NE, y, f, cnt,failed);

  y := hyperg_2F1(1, 2+1e-15, 3, -2);
  f := 0.4506938556659450053;
  testrel(60, NE, y, f, cnt,failed);

  y := hyperg_2F1(1, 1, 2, 0.5);
  f := -2*ln(0.5);
  testrel(61, NE, y, f, cnt,failed);

  {-------------------------------------}
  y := hyperg_2F1(1.5,-1,0.5,1);
  f := -2;
  testrel(62, NE, y, f, cnt,failed);

  y := hyperg_2F1(2.5,-1,0.5,1);
  f := -4;
  testrel(63, NE, y, f, cnt,failed);

  y := hyperg_2F1(3.5,-2,0.5,1);
  f := 8;
  testrel(64, NE, y, f, cnt,failed);

  y := hyperg_2F1(3.5,-3,0.5,1);
  f := -3.2;
  testrel(65, NE, y, f, cnt,failed);

  y := hyperg_2F1(2,-2,2,3);
  f := 4;
  testrel(66, NE, y, f, cnt,failed);

  y := hyperg_2F1(-500,6.25,-600,-0.5);
  f := 0.1133165935032993343;
  testrel(67, NE1, y, f, cnt,failed);

  y := hyperg_2F1(-5000,6.25,-6000.5,0.5);
  f := 29.01957411332703673;
  testrel(68, NE, y, f, cnt,failed);

  y := hyperg_2F1(-5000,6.25,1e10,0.5);
  f := 0.9999984375014157314;
  testrel(69, NE, y, f, cnt,failed);

  y := hyperg_2F1(-50, 1.75,-60,-2);
  f := 3770.651384581751732;
  testrel(70, NE3, y, f, cnt,failed);

  y := hyperg_2F1(-50, 1.75,-60,-1.25);
  f := 0.2861691872517557611;
  testrel(71, NE1, y, f, cnt,failed);

  x := 0.5*sqrt(2);
  y := hyperg_2F1(-100,-200,-300 + 1e-9,x);
  f := 0.2653635302903697294e-30;
  testrel(72, NE1, y, f, cnt,failed);

  y := hyperg_2F1(-50,7/4,-60, 2);
  f := 0.5273362225866704285e9;
  testrel(73, NE, y, f, cnt,failed);

  y := hyperg_2F1(5.5,-300.125,10.25,0.5);
  f := 0.3402263489853102604e-7;
  testrel(74, NE, y, f, cnt,failed);

  y := hyperg_2F1(90+1/2,1,3/2,9/10);
  f := 0.9860651611217264268e89;
  testrel(75, NE2, y, f, cnt,failed);

  y := hyperg_2F1(-50,6.25,-60,-2);
  f := 0.1499207946299185538e10;
  testrel(76, NE3, y, f, cnt,failed);

  y := hyperg_2F1(-36,6.25,-40,-2);
  f := 0.1512306062973965135e12;
  testrel(77, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,2,3.5,-20);
  f := 0.1014182085412905503;
  testrel(78, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5,6.25,1.5,-2);
  f := 0.1042298608575572767e-2;
  testrel(79, NE, y, f, cnt,failed);

  y := hyperg_2F1(-500,6.25,-600,-0.75);
  f := 0.4804874307780050919e-1;
  testrel(80, NE, y, f, cnt,failed);

  y := hyperg_2F1(-50,6.25,-51,-1);
  f := 33279.23823517919823;
  testrel(81, NE, y, f, cnt,failed);

  y := hyperg_2F1(-50,6.25,-51,-0.999999999999999);
  f := 33279.23823517750261;
  testrel(82, NE1, y, f, cnt,failed);

  y := hyperg_2F1(1,2,1e9,0.5);
  f := 1.000000001000000002;
  testrel(83, NE, y, f, cnt,failed);

  y := hyperg_2F1(2.5,6.25,1e10,0.5);
  f := 1.000000000781250000;
  testrel(84, NE, y, f, cnt,failed);

  y := hyperg_2F1(-50,1.25,-60,-0.75);
  f := 0.5446668263717358280;
  testrel(85, NE, y, f, cnt,failed);

  y := hyperg_2F1(-500.2,1.25,-600.1, -0.75);
  f := 0.5449472884383774076;
  testrel(86, NE, y, f, cnt,failed);

  y := hyperg_2F1(-5,1.25,-6,-5);
  f := -623.9669392903645833;
  testrel(87, NE, y, f, cnt,failed);

  y := hyperg_2F1(2.5,6.25,1e10,-0.9990234375);
  f := 0.9999999984390258809;
  testrel(88, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1,2.5,1.25,0.5);
  f := 0;
  testrel(89, NE, y, f, cnt,failed);

  y := hyperg_2F1(-2,-2.5,-3,1);
  f := -1/24;
  testrel(90, NE, y, f, cnt,failed);

  y := hyperg_2F1(-0.5,-1.0,-1.5,0.9375);
  f := 11/16;
  testrel(91, NE, y, f, cnt,failed);

  y := hyperg_2F1(-0.5,-1.5, 3,-2);
  f := 0.4712023748980449996;
  testrel(92, NE, y, f, cnt,failed);

  y := hyperg_2F1(0,-1, 0, 0.5);
  f := 1; {Maple=0.5, Wolfram = 1}
  testrel(93, NE, y, f, cnt,failed);

  y := hyperg_2F1(-0.5,-1.0,-1.5,1);
  f := 2/3;
  testrel(94, NE, y, f, cnt,failed);

  y := hyperg_2F1(-0.4,-1.1,-1.5,-10);
  f := 3.947566948222521467;
  testrel(95, NE, y, f, cnt,failed);

  y := hyperg_2F1(-0.5,-1.1,-1.5,1);
  f := 0;
  testrel(96, NE, y, f, cnt,failed);

  y := hyperg_2F1(-0.5,-1.1,-1.5,-1);
  f := 1.357579719212638008;
  testrel(97, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.1,0,-1,1);
  f := 1;
  testrel(98, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.2,-1,0.5,1);
  f := -1.4;
  testrel(99, NE, y, f, cnt,failed);

  y := hyperg_2F1(0.1,0.2,0.3,0.5);
  f := 1.046432811217352074;
  testrel(100, NE, y, f, cnt,failed);

  y := hyperg_2F1(-0.1,0.2,0.3,0.5);
  f := 0.9564342109682142140;
  testrel(101, NE, y, f, cnt,failed);

  y := hyperg_2F1(1e-8 ,1e-8 ,1e-8 ,1e-6);
  f := 1.000000000000010000;
  testrel(102, NE, y, f, cnt,failed);

  y := hyperg_2F1(2+1e-9, 3,5, -0.75);
  f := 0.4922388588526510234;
  testrel(103, NE, y, f, cnt,failed);

  y := hyperg_2F1(-2,-3,-5+1e-9 ,0.5);
  f := 0.4749999999137500000;
  testrel(104, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1,-1.5,-2-1e-15,0.5);
  f := 0.6250000000000001875;
  testrel(105, NE, y, f, cnt,failed);

  y := hyperg_2F1(500,-500,500,0.75);
  f := 0.9332636185032188790e-301;
  testrel(106, NE, y, f, cnt,failed);

  y := hyperg_2F1(500,-500, -500, 0.75);
  f := 0.1071508607186267321e302;
  testrel(107, NE, y, f, cnt,failed);

  y := hyperg_2F1(500,500,500,-0.625);
  f := 0.3743840538585102283e-105;
  testrel(108, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1000,-2000,-4000.1,-0.5);
  f := 0.5233580403196932348e95;
  testrel(109, NE1, y, f, cnt,failed);

  y := hyperg_2F1(300,10,5,0.5);
  f := 0.3912238919961547169e99;
  testrel(110, NE, y, f, cnt,failed);

  y := hyperg_2F1(5,-300,10,-0.5);
  f := 0.7931780452678445379e47;
  testrel(111, NE, y, f, cnt,failed);

  y := hyperg_2F1(-5,-300,10,-0.5);
  f := -184448.0300324675325;
  testrel(112, NE, y, f, cnt,failed);

  y := hyperg_2F1(5,-300,10,0.5);
  f := 0.1661006238211309107e-6;
  testrel(113, NE, y, f, cnt,failed);

  y := hyperg_2F1(10,5,-300.5,0.5);
  f := 0.9211827166328477894;
  testrel(114, NE, y, f, cnt,failed);

  y := hyperg_2F1(2.25, 3.75, -0.5, -1);
  f := -0.6312206769497026561;
  testrel(115, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5,3,3,-2);
  f := 0.1924500897298752548;
  testrel(116, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5,3,3,0.75);
  f := 8;
  testrel(117, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1,2,3.5,-1000);
  f := 572.4285714285714286;
  testrel(118, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1,2,3.5,10);
  f := -4.714285714285714286;
  testrel(119, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1,2,1,10);
  f := -19;
  testrel(120, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1,2,3.5,-1e20);
  f := 0.5714285714285714286e20;
  testrel(121, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1.125,2.75,3.5,-1e20);
  f := 0.2422430230688187783e23;
  testrel(122, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1,2.75,3.5,-1e50);
  f := 0.7857142857142857143e50;
  testrel(123, NE, y, f, cnt,failed);

  y := hyperg_2F1(-1,2.75,3.5,1e50);
  f := -0.7857142857142857143e50;
  testrel(124, NE, y, f, cnt,failed);

  y := hyperg_2F1(1.5, 1000, 2.25, 0.5);
  f := 0.7705550948378541170e299;
  testrel(125, NE1, y, f, cnt,failed);

  y := hyperg_2F1(1.125,2.75,3.5,-1e50);
  f := 0.8523483860210403255e-56;
  testrel(126, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,2,3.5,-20);
  f := 0.1014182085412905503;
  testrel(127, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,2,0.5,0.75);
  f := 65.53118474162122841;
  testrel(128, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,2,3.5,0.75);
  f := 1.938669492292365108;
  testrel(129, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,2,3.5,-0.75);
  f := 0.7111642172774965541;
  testrel(130, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,0.5,3.5,-0.75);
  f := 0.9130746471820776254;
  testrel(131, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,0.5,3.5,0.9375);
  f := 1.217759238675065672;
  testrel(132, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,2,3.5,-20);
  f := 0.1014182085412905503;
  testrel(133, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,2,3.5,0.9375);
  f := 2.962465226592878248;
  testrel(134, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,0.5,3.5,-1.75);
  f := 0.8343569941403365761;
  testrel(135, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,-1,3.5,-1.75);
  f := 1.5;
  testrel(136, NE, y, f, cnt,failed);

  y := hyperg_2F1(1,2,3.5,-2.5);
  f := 0.4408735093563970956;
  testrel(137, NE, y, f, cnt,failed);

  y := hyperg_2F1(-50,6.25,-51,-2);
  f := 0.6449742538251207026e20;
  testrel(138, NE, y, f, cnt,failed);

  y := hyperg_2F1(-50,6.25,-51.5,-2);
  f := 0.9886773025921599926e19;
  testrel(139, NE, y, f, cnt,failed);

(*
  {Problem cases}

  y := hyperg_2F1(-500,6.25,-600,-1.5);
  f := 0.627841772681655481946316866475e-2;
  writeln(1-y/f);

  y := hyperg_2F1(-500, 6.25,-600,-1.25);
  f := 0.115277178300449647772660141090e-1;
  writeln(1-y/f);

  y := hyperg_2F1(-500,6.25,-600,-2);
  f := 5.59957326424345415374046211996e28;
  writeln(1-y/f);

  y := hyperg_2F1(-125,6.25,-150,-2);
  f := -0.506691674075524839352943058984e14;
  writeln(1-y/f);

  y := hyperg_2F1(-500,6.25,-600,-5);
  f := 0.271363074642659379575067867518e238;
  writeln(1-y/f);
*)

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_hyperg_2F1r;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
  NE1 = 2;
  NE2 = 100;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hyperg_2F1r');

  {functions.wolfram.com}
  y := hyperg_2F1r(0, 5, 0, 0.5);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := hyperg_2F1r(1, 5, 0, 0.5);
  f := 160.0;
  testrel(2, NE, y, f, cnt,failed);

  y := hyperg_2F1r(1, 5, ldexpd(1,-40), 0.5);
  f := 159.9999999997519665;
  testrel(3, NE, y, f, cnt,failed);

  y := hyperg_2F1r(-6, 5, -2, 1);
  f := 840.0;
  testrel(4, NE, y, f, cnt,failed);

  y := hyperg_2F1r(-6, 5, -5, 0.5);
  f := 2362.50;
  testrel(5, NE, y, f, cnt,failed);

  y := hyperg_2F1r(-6, 5, -6, 0.5);
  f := 0;
  testrel(6, NE, y, f, cnt,failed);

  y := hyperg_2F1r(1.25, -6, -6, 0.5);
  f := 0;
  testrel(7, NE, y, f, cnt,failed);

  y := hyperg_2F1r(-6, 5, -4, 0.5);
  f := -472.5;
  testrel(8, NE, y, f, cnt,failed);

  y := hyperg_2F1r(-6, 5, -4+ldexpd(1,-40), 0.5);
  f := -472.5000000001784741;
  testrel(9, NE, y, f, cnt,failed);

  y := hyperg_2F1r(-6, 5, -4, -2);
  f := 12579840.0;
  testrel(10, NE, y, f, cnt,failed);

  y := hyperg_2F1r(-6, 5, -4, 2);
  f := 6773760.0;
  testrel(11, NE, y, f, cnt,failed);

  y := hyperg_2F1r(-6, 5, -4, Pi);
  f := 1.175999402648652974e8;
  testrel(12, NE1, y, f, cnt,failed);

  {Corrected GSL test case with values from Wolfram Function site.}
  {Note: GSL non-trivial case are totally buggy and incomplete!   }
  {Missing x^(1-c) and a,b negative integers are not allowed!     }
  y := hyperg_2F1r(1, 1, 1, 0.5);
  f := 2.0;
  testrel(13, NE, y, f, cnt,failed);

  y := hyperg_2F1r(8, 8, 1, 0.5);
  f := 12451584.0;
  testrel(14, NE, y, f, cnt,failed);

  y := hyperg_2F1r(8, -8, 1, 0.5);
  f := 0.13671875;
  testrel(15, NE, y, f, cnt,failed);

  y := hyperg_2F1r(8, -8, 1, -0.5);
  f := 4945.13671875;
  testrel(16, NE, y, f, cnt,failed);

  y := hyperg_2F1r(8, -8, -5.5, 0.5);
  f := -83081.19167659493609;
  testrel(17, NE, y, f, cnt,failed);

  y := hyperg_2F1r(8, -8, -5.5, -0.5);
  f := 2.2510895952730178518e+06;
  testrel(18, NE, y, f, cnt,failed);

  y := hyperg_2F1r(8, 8, 5, 0.5);
  f := 175.2380952380952381;
  testrel(19, NE, y, f, cnt,failed);

  y := hyperg_2F1r(9, 9, -1.5, 0.99);
  f := 1.6063266334913066551e+46;
  testrel(20, NE2, y, f, cnt,failed);  {!!! inaccuracy near 1}

  y := hyperg_2F1r(9, 9, -1.5, -0.99);
  f := 0.06234327316254516616;
  testrel(21, NE, y, f, cnt,failed);

  y := hyperg_2F1r(5, 5, -1, 0.5);
  f := 1237440;
  testrel(22, NE, y, f, cnt,failed);

  y := hyperg_2F1r(5, 5, -10, 0.5);
  f := 6.8070553334784e16;
  testrel(23, NE, y, f, cnt,failed);

  y := hyperg_2F1r(5, 5, -100, 0.5);
  f := 1.191184208770032773e176;
  testrel(24, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_hyperg_2F0;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','hyperg_2F0');

  {Computed with MPMath 0.17}
  y := hyperg_2F0(1.5, 2.5, -3.5);
  f := 0.06232362765198486065;
  testrel(1, NE, y, f, cnt,failed);

  y := hyperg_2F0(-2, 1.5, 3.5);
  f := 36.4375;
  testrel(2, NE, y, f, cnt,failed);

  y := hyperg_2F0(1.5, -2, 3.5);
  f := 36.4375;
  testrel(3, NE, y, f, cnt,failed);

  y := hyperg_2F0(1.5, 0, 3.5);
  f := 1;
  testrel(4, NE, y, f, cnt,failed);

  y := hyperg_2F0(1.5, 2.5, -100);
  f := 0.7092776453693892794e-3;
  testrel(5, NE, y, f, cnt,failed);

  y := hyperg_2F0(15, 2.5, -100);
  f := 0.1566288728955602417e-7;
  testrel(6, NE, y, f, cnt,failed);

  y := hyperg_2F0(-15, 2.5, -100);
  f := 6.500742875742012396e+43;
  testrel(7, NE, y, f, cnt,failed);

  y := hyperg_2F0(-15, 2.5, 100);
  f := -6.383615572317773340e+43;
  testrel(8, NE, y, f, cnt,failed);

  y := hyperg_2F0(-2.5, 1.5, 0);
  f := 1;
  testrel(9, NE, y, f, cnt,failed);

  y := hyperg_2F0(1.5, 2.5, -0.01);
  f := 0.9640594689959897311;
  testrel(10, NE, y, f, cnt,failed);

  y := hyperg_2F0(50, 50, -0.03125);
  f := 2.367887757397116595e-17;
  testrel(11, NE, y, f, cnt,failed);

  y := hyperg_2F0(50, 30, -0.03125);
  f := 1.519977201318170588e-11;
  testrel(12, NE, y, f, cnt,failed);

  y := hyperg_2F0(5, 10, -0.03125);
  f := 0.2778878346657118121;
  testrel(13, NE, y, f, cnt,failed);

  {Asyptotic sum}
  y := hyperg_2F0(1.5, 0.5, 0.01);
  f := 1.007644896568125858;
  testrel(14, NE, y, f, cnt,failed);

  y := hyperg_2F0(1.5, 2.5, 0.01);
  f := 1.039232502438289330;
  testrel(15, NE, y, f, cnt,failed);

  y := hyperg_2F0(1.5, 0.5, 0.01953125);
  f := 1.015218046373990184;
  testrel(16, NE, y, f, cnt,failed);

  y := hyperg_2F0(1.5, -2.5, 0.01953125);
  f := 0.9294091153233346690;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


end.

