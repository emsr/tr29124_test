{Part 4 of regression test for SPECFUNX unit  (c) 2018  W.Ehrhardt}

unit t_sfx2c;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_mathematicax;

procedure test_detaix;
procedure test_emlambdax;
procedure test_KleinJx;
procedure test_wplx;
procedure test_wpex;
procedure test_wpe_derx;
procedure test_wpe_invx;
procedure test_wpgx;
procedure test_wpg_derx;
procedure test_wpg_invx;


implementation

uses
  amath, specfunx, t_sfx0;

{---------------------------------------------------------------------------}
procedure test_mathematicax;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 3;
  NE1 = 6;
  NE2 = 10;
  NE3 = 20;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Elliptic integrals (Mathematica)');

  {----------------------------------}
  y := M_EllipticKx(0.0);
  f := Pi_2;
  testrel(1, NE, y, f, cnt,failed);

  y := M_EllipticKx(0.01);
  f := 1.574745561517355953;
  testrel(2, NE, y, f, cnt,failed);

  y := M_EllipticKx(-0.01);
  f := 1.566891273068196358;
  testrel(3, NE, y, f, cnt,failed);

  y := M_EllipticKx(0.5);
  f := 1.854074677301371918;
  testrel(4, NE1, y, f, cnt,failed);

  y := M_EllipticKx(-0.5);
  f := 1.415737208425956199;
  testrel(5, NE, y, f, cnt,failed);

  y := M_EllipticKx(0.99);
  f := 3.695637362989874678;
  testrel(6, NE, y, f, cnt,failed);

  y := M_EllipticKx(-0.99);
  f := 1.312813847340966180;
  testrel(7, NE, y, f, cnt,failed);

  y := M_EllipticKx(3);
  f := 1.001077380456106236;
  testrel(8, NE, y, f, cnt,failed);

  y := M_EllipticKx(-3);
  f := 1.078257823749821618;
  testrel(9, NE, y, f, cnt,failed);

  y := M_EllipticKx(100);
  f := 0.1574745561517355953;
  testrel(10, NE, y, f, cnt,failed);

  y := M_EllipticKx(-100);
  f := 0.3682192486091410329;
  testrel(11, NE, y, f, cnt,failed);

  {----------------------------------}
  y := M_EllipticECx(0);
  f := Pi_2;
  testrel(12, NE, y, f, cnt,failed);

  y := M_EllipticECx(0.01);
  f := 1.566861942021668291;
  testrel(13, NE, y, f, cnt,failed);

  y := M_EllipticECx(-0.01);
  f := 1.574715985016988413;
  testrel(14, NE, y, f, cnt,failed);

  y := M_EllipticECx(0.5);
  f := 1.350643881047675503;
  testrel(15, NE, y, f, cnt,failed);

  y := M_EllipticECx(-0.5);
  f := 1.751771275694817862;
  testrel(16, NE, y, f, cnt,failed);

  y := M_EllipticECx(0.75);
  f := 1.211056027568459525;
  testrel(17, NE, y, f, cnt,failed);

  y := M_EllipticECx(-0.75);
  f := 1.833204967048622144;
  testrel(18, NE, y, f, cnt,failed);

  y := M_EllipticECx(1);
  f := 1;
  testrel(19, NE, y, f, cnt,failed);

  y := M_EllipticECx(-1);
  f := 1.910098894513856009;
  testrel(20, NE, y, f, cnt,failed);

  y := M_EllipticECx(3);
  f := 0.4752239353510171110;
  testrel(21, NE, y, f, cnt,failed);

  y := M_EllipticECx(-3);
  f := 2.422112055136919050;
  testrel(22, NE, y, f, cnt,failed);

  y := M_EllipticECx(100);
  f := 0.07863836119485898078;
  testrel(23, NE, y, f, cnt,failed);

  y := M_EllipticECx(-100);
  f := 10.20926091981457201;
  testrel(24, NE, y, f, cnt,failed);

  {----------------------------------}
  y := M_EllipticPiCx(-1, 0.5);
  f := 1.273127366749682458;
  testrel(25, NE1, y, f, cnt,failed);

  y := M_EllipticPiCx(-1, -0.5);
  f := 1.018679933006156207;
  testrel(26, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(-1, 0.75);
  f := 1.440034318657550564;
  testrel(27, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(-1, -0.75);
  f := 0.9843963918340774593;
  testrel(28, NE1, y, f, cnt,failed);

  y := M_EllipticPiCx(-1, 3);
  f := 0.8606964344717086121;
  testrel(29, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(-1, -3);
  f := 0.8086933678123662133;
  testrel(30, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(-1, 100);
  f := 0.1566920605183800089;
  testrel(31, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(-1, -100);
  f := 0.3071840408785169673;
  testrel(32, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(-1, -1);
{$ifdef BIT16}
  f := 1.910098894513856009*0.5;
{$else}
  f := 0.9550494472569280045;
{$endif}
  testrel(33, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(0.5, 0.5);
  f := 2.701287762095351005;
  testrel(34, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(0.5, -0.5);
  f := 1.967853602214696611;
  testrel(35, NE1, y, f, cnt,failed);

  y := M_EllipticPiCx(0.5, 0.75);
  f := 3.234773471249464853;
  testrel(36, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(0.5, -0.75);
  f := 1.876863274442115336;
  testrel(37, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(0.5, 3);
  f := 1.101695283805732866;
  testrel(38, NE1, y, f, cnt,failed);

  y := M_EllipticPiCx(0.5, -3);
  f := 1.440034318657550564;
  testrel(39, NE1, y, f, cnt,failed);

  y := M_EllipticPiCx(0.5, 100);
  f := 0.1578702221146095725;
  testrel(40, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(0.5, -100);
  f := 0.4457715465792184755;
  testrel(41, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(2, 0.5);
  f := -0.3135446834651840415;
  testrel(42, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(2, -0.5);
  f := 0.1422952035345872497;
  testrel(43, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(2, 0.75);
  f := -0.6824791393936852392;
  testrel(44, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(2, -0.75);
  f := 0.1871128632591755178;
  testrel(45, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(2, 3);
  f := -0.4395535965185157402;
  testrel(46, NE1, y, f, cnt,failed);

  y := M_EllipticPiCx(2, -3);
  f := 0.3412395696968426196;
  testrel(47, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(2, 100);
  f := -0.06532413297793291762;
  testrel(48, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(2, -100);
  f := 0.2416687652198913641;
  testrel(49, NE, y, f, cnt,failed);

  y := M_EllipticPiCx(2, 2);
  f := -0.5990701173677961037;
  testrel(50, NE, y, f, cnt,failed);

  {----------------------------------}
  y := M_EllipticFx(1,1);
  f := 1.226191170883517071;
  testrel(51, NE, y, f, cnt,failed);

  y := M_EllipticFx(1,-3);
  f := 0.7807065662256886254;
  testrel(52, NE, y, f, cnt,failed);

  y := M_EllipticFx(10,-3);
  f := 6.979796589346306317;
  testrel(53, NE, y, f, cnt,failed);

  y := M_EllipticFx(100,-0.5);
  f := 90.08736322283863687;
  testrel(54, NE, y, f, cnt,failed);

  y := M_EllipticFx(10,2);
  f := 8.526932128869660257;
  testrel(55, NE, y, f, cnt,failed);

  y := M_EllipticFx(10,-10);
  f := 5.179035988855946898;
  testrel(56, NE, y, f, cnt,failed);

  y := M_EllipticFx(1000,1/16);
  f := 1016.191702342398549;
  testrel(57, NE1, y, f, cnt,failed);

  y := M_EllipticFx(Pi,2);
  f := 2.622057554292119810;
  testrel(58, NE1, y, f, cnt,failed);

  y := M_EllipticFx(Pi_2,2);
  f := 1.311028777146059905;
  testrel(59, NE1, y, f, cnt,failed);

  y := M_EllipticFx(0.25*Pi,2);
  f := 1.311028777146059905;
  testrel(60, NE, y, f, cnt,failed);

  {----------------------------------}
  y := M_EllipticEx(1,1);
  f := 0.8414709848078965067;
  testrel(61, NE, y, f, cnt,failed);

  y := M_EllipticEx(1,-3);
  f := 1.325663197579998112;
  testrel(62, NE, y, f, cnt,failed);

  y := M_EllipticEx(10,-3);
  f := 15.18760531092766668;
  testrel(63, NE1, y, f, cnt,failed);

  y := M_EllipticEx(100,-0.5);
  f := 111.5708277453651198;
  testrel(64, NE1, y, f, cnt,failed);

  y := M_EllipticEx(10,1);
  f := 6.544021110889369813;
  testrel(65, NE, y, f, cnt,failed);

  y := M_EllipticEx(10,2);
  f := 4.103219034329508954;
  testrel(66, NE, y, f, cnt,failed);

  y := M_EllipticEx(10,-10);
  f := 22.63533049120160840;
  testrel(67, NE, y, f, cnt,failed);

  y := M_EllipticEx(1000,1/16);
  f := 984.1943502924860745;
  testrel(68, NE1, y, f, cnt,failed);

  y := M_EllipticEx(Pi,2);
  f := 1.1981402347355922074;
  testrel(69, NE1, y, f, cnt,failed);

  y := M_EllipticEx(Pi_2,2);
  f := 0.5990701173677961037;
  testrel(70, NE, y, f, cnt,failed);

  y := M_EllipticEx(0.25*Pi,2);
  f := 0.5990701173677961037;
  testrel(71, NE, y, f, cnt,failed);

  {----------------------------------}
  {max(n,m) < 1}
  y := M_EllipticPix(0.5, 1, 0);
  f := 1.178815078927437390;
  testrel(72, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.5, 1, 0);
  f := 1.178815078927437390;
  testrel(73, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.5, pi/4, 0.5);
  f := 0.9190227391656969904;
  testrel(74, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.5, pi/4, 2);
  f := 1.532738349029281103;
  testrel(75, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.25, 10, 0.5);
  f := 13.61321845241576999;
  testrel(76, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.5, 1, 0.5);
  f := 1.288978174244979232;
  testrel(77, NE, y, f, cnt,failed);

  y := M_EllipticPix(-1, 2, 0.5);
  f := 1.577304510600466028;
  testrel(78, NE1, y, f, cnt,failed);

  y := M_EllipticPix(-2, 2, 0.5);
  f := 1.227640864405237886;
  testrel(79, NE1, y, f, cnt,failed);

  y := M_EllipticPix(-0.5, 2, 0.5);
  f := 1.889259502525944879;
  testrel(80, NE1, y, f, cnt,failed);

  y := M_EllipticPix(-0.5, 2, -2);
  f := 1.155284601870055215;
  testrel(81, NE, y, f, cnt,failed);

  y := M_EllipticPix(-1, -1, -1);
  f := -0.7358811532333216979;
  testrel(82, NE, y, f, cnt,failed);

  y := M_EllipticPix(0, -10, 0.5);
  f := -11.71562231566589297;
  testrel(83, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.5, 0, 0.5);
  f := 0;
  testrel(84, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.5, Pi, 0.5);
  f := 5.402575524190702010;
  testrel(85, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.5, 3*Pi_2, 0.5);
  f := 8.103863286286053015;
  testrel(86, NE, y, f, cnt,failed);

  {m < 1 <= n}
  y := M_EllipticPix(1, 10, 0.5);
  f := -4.415211526644026984;
  testrel(87, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, 1, 0);
  f := 0.7617262217813367605;
  testrel(88, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, -1, 0.5);
  f := -0.704583746768798274;
  testrel(89, NE1, y, f, cnt,failed);

  y := M_EllipticPix(2, 1, 0.5);
  f := 0.7045837467687982743;
  testrel(90, NE, y, f, cnt,failed);

  y := M_EllipticPix(3, 1, 0.5);
  f := 0.2715305188238305158;
  testrel(91, NE1, y, f, cnt,failed);

  y := M_EllipticPix(3, -1, -2);
  f := -0.4406586065413746729;
  testrel(92, NE, y, f, cnt,failed);

  y := M_EllipticPix(3, 3*Pi_2, -2);
  f := 0.6961536892272037603;
  testrel(93, NE, y, f, cnt,failed);

  {n < 1 <= m}
  y := M_EllipticPix(0.5, 1, 1);
  f := 1.483099873420077333;
  testrel(94, NE, y, f, cnt,failed);

  y := M_EllipticPix(0.25, 10, 2);
  f := 9.138538826251456719;
  testrel(95, NE, y, f, cnt,failed);

  y := M_EllipticPix(-1, 1/8, 5);
  f := 0.1260234435925222502;
  testrel(96, NE, y, f, cnt,failed);

  {1 <= n <= m}
  y := M_EllipticPix(1, 1, 1);
  f := 2.054332933256248669;
  testrel(97, NE, y, f, cnt,failed);

  y := M_EllipticPix(1, 0.5, 1.5);
  f := 0.5882230881053119906;
  testrel(98, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, -1/8, 3);
  f := -0.1273364223077788701;
  testrel(99, NE, y, f, cnt,failed);

  y := M_EllipticPix(1.5, 1/8, 2);
  f := 0.1266550462017603176;
  testrel(100, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, 1/8, 2);
  f := 0.1269936782273970098;
  testrel(101, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, 3.0*Pi_2, 2);
  f := -1.797210352103388311;
  testrel(102, NE1, y, f, cnt,failed);

  {1 <= m <= n}
  y := M_EllipticPix(3, 1, 1);
  f := 0.11487473777461448788;
  testrel(103, NE2, y, f, cnt,failed);

  y := M_EllipticPix(3, 0.5, 2);
  f := 0.8272068645186950740;
  testrel(104, NE, y, f, cnt,failed);

  y := M_EllipticPix(5, 0.5, 3);
  f := 0.8469758520611422444;
  testrel(105, NE, y, f, cnt,failed);

{$define UseMathematicaMode}
{$ifdef UseMathematicaMode}
  y := M_EllipticPix(2, 2, 0);
  f := -0.4943441954895304141;
  testrel(106, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, 1, 0.5);
  f := 0.7045837467687982743;
  testrel(107, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, 1+8*Pi, 0.5);
  f := 0.7045837467687982743;
  testrel(108, NE3, y, f, cnt,failed);

  y := M_EllipticPix(2, 1-4*Pi, 0.5);
  f := 0.7045837467687982743;
  testrel(109, NE3, y, f, cnt,failed);

  y := M_EllipticPix(5, 8, -2);
  f := -0.1798487636094046750;
  testrel(110, NE, y, f, cnt,failed);

  y := M_EllipticPix(5, 0.5, 2);
  f := 0.8365425120605017677;
  testrel(111, NE, y, f, cnt,failed);

  y := M_EllipticPix(10, 0.5, 2);
  f := 0.1957109568289841125;
  testrel(112, NE1, y, f, cnt,failed);

  y := M_EllipticPix(2, -5, 0.5);
  f := 0.1117793671573657844;
  testrel(113, NE1, y, f, cnt,failed);

  y := M_EllipticPix(3, 2, 0.5);
  f := -0.1327520769500431896;
  testrel(114, NE1, y, f, cnt,failed);
{$else}
  y := M_EllipticPix(2, 2, 0);
  f := -0.4943441954895304141;
  testrel(106, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, 1, 0.5);
  f := 0.7045837467687982743;          {Maple/MPMath}
  testrel(107, NE, y, f, cnt,failed);

  y := M_EllipticPix(2, 1+8*Pi, 0.5);
  f := -4.312131188674146389;          {Maple/MPMath}
  testrel(108, NE1, y, f, cnt,failed);

  y := M_EllipticPix(2, 1-4*Pi, 0.5);
  f := 3.212941214490270606;           {Maple/MPMath}
  testrel(109, NE1, y, f, cnt,failed);

  y := M_EllipticPix(5, 8, -2);
  f := 0.7713463514647902789;          {Maple/MPMath}
  testrel(110, NE, y, f, cnt,failed);

  y := M_EllipticPix(5, 0.5, 2);
  f := 0.8365425120605017677;
  testrel(111, NE, y, f, cnt,failed);

  y := M_EllipticPix(10, 0.5, 2);
  f := 0.1957109568289841125;
  testrel(112, NE1, y, f, cnt,failed);

  y := M_EllipticPix(2, -5, 0.5);
  f := 1.365958101018101950;           {Maple/MPMath}
  testrel(113, NE1, y, f, cnt,failed);

  y := M_EllipticPix(3, 2, 0.5);
  f := -0.5187284939303218479;         {Maple/MPMath}
  testrel(114, NE, y, f, cnt,failed);
{$endif}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_detaix;
var
  x,y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 8;
  NE2 = 100;
  NE3 = 500;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','detaix');

  x := 0.5;
  f := 0.8377557634765980579;
  y := detaix(x);
  testrel(1, NE, y, f, cnt,failed);

  x := 1;
  f := 0.7682254223260566590;
  y := detaix(x);
  testrel(2, NE, y, f, cnt,failed);

  x := 2;
  f := 0.5923827813324158853;
  y := detaix(x);
  testrel(3, NE, y, f, cnt,failed);

  x := 5;
  f := 0.2700908381440137190;
  y := detaix(x);
  testrel(4, NE, y, f, cnt,failed);

  x := 10;
  f := 0.7294906084933912960e-1;
  y := detaix(x);
  testrel(5, NE, y, f, cnt,failed);

  x := 100;
  f := 4.267731135455224686e-12;
  y := detaix(x);
  testrel(6, NE, y, f, cnt,failed);

  x := 1000;
  f := 2.004335936527980120e-114;
  y := detaix(x);
  testrel(7, NE2, y, f, cnt,failed);

  x := 2500;
  f := 0.5687563820662116353e-284;
  y := detaix(x);
  testrel(8, NE2, y, f, cnt,failed);

  x := 1/8;
  f := 0.3483058410597969855;
  y := detaix(x);
  testrel(9, NE, y, f, cnt,failed);

  x := 1/16;
  f := 0.6065847945818627981e-1;
  y := detaix(x);
  testrel(10, NE, y, f, cnt,failed);

  x := 1/64;
  f := 4.230737693555325865e-7;
  y := detaix(x);
  testrel(11, NE1, y, f, cnt,failed);

  x := 1/1024;
  f := 1.197754424502496408e-115;
  y := detaix(x);
  testrel(12, NE2, y, f, cnt,failed);

  x := 1/2048;
  f := 6.340165391144919936e-232;
  y := detaix(x);
  testrel(13, NE2, y, f, cnt,failed);

  x := 1/8192;
  f := 3.486883326152397239e-930;
  y := detaix(x);
  testrel(14, NE3, y, f, cnt,failed);

  x := 40000;
  f := 0.11990037107871639319e-4547;
  y := detaix(x);
  testrel(15, NE3, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_emlambdax;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 6;
  NE1 = 120;
  NE2 = 5200; {extreme value}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','emlambdax');

  {N[ModularLambda[y*I], 30]}
  y := emlambdax(1/15);
  f := 1;
  testrel(1, NE, y, f, cnt,failed);

  y := emlambdax(1/10);
  f := 1.999999999999273248*0.5;
  testrel(2, NE, y, f, cnt,failed);

  y := emlambdax(0.125);
  f := 1.999999999610830185*0.5;
  testrel(3, NE, y, f, cnt,failed);

  y := emlambdax(0.25);
  f := 1.999888408157900115*0.5;
  testrel(4, NE, y, f, cnt,failed);

  y := emlambdax(0.5);
  f := 1.941125496954281171*0.5;
  testrel(5, NE, y, f, cnt,failed);

  y := emlambdax(2/3);
  f := 0.8661058727342564978;
  testrel(6, NE, y, f, cnt,failed);

  y := emlambdax(3/4);
  f := 0.7845002957080673029;
  testrel(7, NE, y, f, cnt,failed);

  y := emlambdax(1);
  f := 0.5;
  testrel(8, NE, y, f, cnt,failed);

  y := emlambdax(2);
  f := 0.2943725152285941438e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := emlambdax(4);
  f := 0.5579592104994237345e-4;
  testrel(10, NE, y, f, cnt,failed);

  y := emlambdax(10);
  f := 3.633761709317889931e-13;
  testrel(11, NE, y, f, cnt,failed);

  y := emlambdax(100);
  f := 5.840964927192880684e-136;
  testrel(12, NE1, y, f, cnt,failed);

  y := emlambdax(200);
  f := 2.132304455043583372e-272;
  testrel(13, NE1, y, f, cnt,failed);

  y := emlambdax(3600);
  f := 2.813495028131470324e-4911;
  testrel(14, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_KleinJx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 8;
  NE1 = 32;
  NE2 = 64;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','KleinJx');

  f := 1.121975409927588378e24;
  y := KleinJx(1/10);
  testrel(1, NE, y, f, cnt,failed);

  y := KleinJx(10);
  testrel(2, NE, y, f, cnt,failed);

  f := 2.527141464766729677e171;
  y := KleinJx(1/64);
  testrel(3, NE2, y, f, cnt,failed);

  y := KleinJx(64);
  testrel(4, NE2, y, f, cnt,failed);

  f := 3.912712369665429648e18;
  y := KleinJx(1/8);
  testrel(5, NE, y, f, cnt,failed);

  y := KleinJx(8);
  testrel(6, NE, y, f, cnt,failed);

  f := 7.306766275995648519e15;
  y := KleinJx(1/7);
  testrel(7, NE1, y, f, cnt,failed);

  y := KleinJx(7);
  testrel(8, NE, y, f, cnt,failed);

  f := 88862.08298421810451;
  y := KleinJx(1/3);
  testrel(9, NE1, y, f, cnt,failed);

  y := KleinJx(3);
  testrel(10, NE, y, f, cnt,failed);

  f := 166.375;
  y := KleinJx(1/2);
  testrel(11, NE, y, f, cnt,failed);

  y := KleinJx(2);
  testrel(12, NE, y, f, cnt,failed);

  f := 2.973899579991082458;
  y := KleinJx(3/4);
  testrel(13, NE, y, f, cnt,failed);

  y := KleinJx(4/3);
  testrel(14, NE, y, f, cnt,failed);

  y := KleinJx(1);
  f := 1;
  testrel(15, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_wplx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 2;
  NE1 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wplx');

  y := wplx(1/1024);
  f := 1.048576000000047684e6;
  testrel(1, NE, y, f, cnt,failed);

  y := wplx(1/16);
  f := 256.0001953125496705;
  testrel(2, NE, y, f, cnt,failed);

  y := wplx(0.5);
  f := 4.012513027096227404;
  testrel(3, NE, y, f, cnt,failed);

  y := wplx(1);
  f := 1.050839791040237018;
  testrel(4, NE1, y, f, cnt,failed);

  y := wplx(1.854075);
  f := 0.5000000000000520672;
  testrel(5, NE, y, f, cnt,failed);

  y := wplx(2);
  f := 0.5107614339783495499;
  testrel(6, NE, y, f, cnt,failed);

  y := wplx(3);
  f := 2.019294402732143640;
  testrel(7, NE1, y, f, cnt,failed);

  y := wplx(4);
  f := 11.744545549928836655;
  testrel(8, NE1, y, f, cnt,failed);

  y := wplx(-1000);
  f := 0.7686394652675544522;
  testrel(9, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_wpex;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 5;
  NE1 = 20;
  NE2 = 1800;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wpex');

  y := wpex(1,2,2);    {e1=e2, mc=0}
  f := 2.181597706907869689;  {Mathematica}
  testrel(1, NE, y, f, cnt,failed);

  y := wpex(4,2,2);    {e1=e2, mc=0}
  f := 2.000000074098977346; {Mathematica}
  testrel(2, NE, y, f, cnt,failed);

  y := wpex(1,-2,-2);  {e1=e2, mc=1}
  f := 12.733139842529623306;
  testrel(3, NE, y, f, cnt,failed);

  y := wpex(4,-2,-2);  {e1=e2, mc=1}
  f := 43.14058060679109643; {Mathematica}
  testrel(4, NE, y, f, cnt,failed);

  y := wpex(1,2,-1);  {e2=e3, mc=1}
  f := 2.079381535373778827; {Mathematica}
  testrel(5, NE, y, f, cnt,failed);

  y := wpex(4,2,-1);  {e2=e3, mc=1}
  f := 7.2997473979957452309; {Mathematica}
  testrel(6, NE, y, f, cnt,failed);

  y := wpex(2,0,0);  {e1=e2=e3=0}
  f := 0.25;
  testrel(7, NE, y, f, cnt,failed);

  y := wpex(10,0,0); {e1=e2=e3=0}
  f := 1/100;
  testrel(8, NE, y, f, cnt,failed);

  y := wpex(2,2,0);
  f := 46.97818219971534662;
  testrel(9, NE, y, f, cnt,failed);

  y := wpex(100,2,0);
  f := 69.41828497727013392;
  testrel(10, NE, y, f, cnt,failed);

  y := wpex(1,3,2);
  f := 3.132015966788728659;
  testrel(11, NE, y, f, cnt,failed);

  y := wpex(-8,3,2);
  f := 3.140119797518689891;
  testrel(12, NE, y, f, cnt,failed);

  y := wpex(4,3,2);
  f := 4.718869292341830266;
  testrel(13, NE, y, f, cnt,failed);

  {near pole}
  y := wpex(0.0078125,3,2);
  f := 16384.00023191762931;
  testrel(14, NE, y, f, cnt,failed);

  y := wpex(3.5,3,2);
  f := 306342.4066537604033;
  testrel(15, NE2, y, f, cnt,failed);

  {Bad range reduction}
  y := wpex(10,3,2);
  f := 4.818826329339614295;
  testrel(16, NE1, y, f, cnt,failed);

  y := wpex(100,3,2);
  f := 11.31473482986812039;
  testrel(17, NE2, y, f, cnt,failed);

  y := wpex(1000,3,2);
  f := 4.985976612159323772;
  testrel(18, NE2, y, f, cnt,failed);

  {negative argument}
  y := wpex(-1,-2,-3);
  f := 45.67407124634771492;
  testrel(19, NE, y, f, cnt,failed);

  y := wpex(-4,-2,-3);
  f := 5.018929622118494273;
  testrel(20, NE, y, f, cnt,failed);

  {a few test for wpex_im}
  y := wpe_imx(1, 0.5, 0);
  f := -1.050839791040237018;
  testrel(21, NE, y, f, cnt,failed);

  y := wpe_imx(1, 1, 0.5);
  f := -1.512365392932099994;
  testrel(22, NE, y, f, cnt,failed);

  y := wpe_imx(1, 0, 2);
  f := -2.043045735913398200;
  testrel(23, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_wpe_derx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 5;
  NE1 = 20;
  NE2 = 60;
  NE3 = 1400;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wpe_derx');

  y := wpe_derx(-8,0,0);
  f := 0.003906250;
  testrel(1, NE, y, f, cnt,failed);

  y := wpe_derx(4,0,0);
  f := -0.03125;
  testrel(2, NE, y, f, cnt,failed);

  y := wpe_derx(1,2,2);    {e1=e2, mc=0}
  f := -0.9030061850406446371;
  testrel(3, NE, y, f, cnt,failed);

  y := wpe_derx(4,2,2);    {e1=e2, mc=0}
  f := -3.630093721594877499e-7;
  testrel(4, NE, y, f, cnt,failed);

  y := wpe_derx(1,-2,-2);  {e1=e2, mc=1}
  f := 87.07841471187595541;
  testrel(5, NE, y, f, cnt,failed);

  y := wpe_derx(4,-2,-2);  {e1=e2, mc=1}
  f := -564.8209126860439458;
  testrel(6, NE, y, f, cnt,failed);

  y := wpe_derx(1,2,-1);  {e2=e3, mc=1}
  f := 1.735214804404404261;
  testrel(7, NE, y, f, cnt,failed);

  y := wpe_derx(4,2,-1);  {e2=e3, mc=1}
  f := -38.21399616480677124;
  testrel(8, NE1, y, f, cnt,failed);

  {Bad range reduction}
  y := wpe_derx(10,3,2);
  f := 14.19023412081285186;
  testrel(9, NE2, y, f, cnt,failed);

  y := wpe_derx(1,3,2);
  f := 2.204797021416015720;
  testrel(10, NE1, y, f, cnt,failed);

  y := wpe_derx(1/16,3,2);
  f := -8191.529157857499106;
  testrel(11, NE1, y, f, cnt,failed);

  y := wpe_derx(1/128,3,2);
  f := -4194303.940633173511;
  testrel(12, NE, y, f, cnt,failed);

  y := wpe_derx(1/2048,3,2);
  f := -17179869183.99628906;
  testrel(13, NE1, y, f, cnt,failed);

  y := wpe_derx(1,2,1);
  f := -0.9456289331929805119e-1;
  testrel(14, NE1, y, f, cnt,failed);

  y := wpe_derx(1,2,0);
  f := 1.192582509608769703;
  testrel(15, NE1, y, f, cnt,failed);

  y := wpe_derx(-1,0,-1);
  f := 1.515722527941721910;
  testrel(16, NE, y, f, cnt,failed);

  y := wpe_derx(10,0.5,0.0);
  f := 1.285113200859947769;
  testrel(17, NE, y, f, cnt,failed);

  y := wpe_derx(1000,0.5,0.0);
  f := 1.023635528833534815;
  testrel(18, NE, y, f, cnt,failed);

  {near pole <> 0}
  y := wpe_derx(3.5,3,2);
  f := -339110061.3727949187;
  testrel(19, NE3, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_wpe_invx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 2;
  NE1 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wpe_invx');

  y := wpe_invx(2, 0, 0);
  f := 0.7071067811865475244;
  testrel(1, NE, y, f, cnt,failed);

  y := wpe_invx(0.5, 0.5, 0);
  f := 1.854074677301371918;
  testrel(2, NE1, y, f, cnt,failed);

  y := wpe_invx(0.75, 0.5, 0);
  f := 1.219046159603532136;
  testrel(3, NE1, y, f, cnt,failed);

  y := wpe_invx(1, 0.5, 0);
  f := 1.028056801052126733;
  testrel(4, NE, y, f, cnt,failed);

  y := wpe_invx(10, 0.5, 0);
  f := 0.3163069054282973787;
  testrel(5, NE, y, f, cnt,failed);

  y := wpe_invx(3,3,2);
  f := 0.8745483141883363653;
  testrel(6, NE, y, f, cnt,failed);

  y := wpe_invx(3+1/1024,3,2);
  f := 0.8635017929130682038;
  testrel(7, NE, y, f, cnt,failed);

  y := wpe_invx(5,3,2);
  f := 0.4823632516973266728;
  testrel(8, NE, y, f, cnt,failed);

  y := wpe_invx(20,3,2);
  f := 0.2246279117357761809;
  testrel(9, NE, y, f, cnt,failed);

  y := wpe_invx(1000,3,2);
  f := 0.3162283661767068010e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := wpe_invx(10000,3,2);
  f := 0.1000000018997858647e-1;
  testrel(11, NE1, y, f, cnt,failed);

  y := wpe_invx(10,2,1);
  f := 0.3183638324100305341;
  testrel(12, NE, y, f, cnt,failed);

  y := wpe_invx(2,2,1);
  f := 1.009452909989211608;
  testrel(13, NE, y, f, cnt,failed);

  y := wpe_invx(2,2,0);
  f := 0.9270373386506859592;
  testrel(14, NE1, y, f, cnt,failed);

  y := wpe_invx(10,2,0);
  f := 0.3175142588489198911;
  testrel(15, NE, y, f, cnt,failed);

  y := wpe_invx(2,2,-1);
  f := 0.9068996821171089253;
  testrel(16, NE, y, f, cnt,failed);

  y := wpe_invx(4,4,-1);
  f := 0.6446666159587342331;
  testrel(17, NE, y, f, cnt,failed);

  y := wpe_invx(5,4,-1);
  f := 0.4843987776422161666;
  testrel(18, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_wpgx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 6;
  NE1 = 12;
  NE2 = 25;
  NE3 = 50;
  NE4 = 400;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wpgx');

  y := wpgx(1,48,-64); {e1=e2}
  f := 2.181597706907869689; {Mathematica}
  testrel(1, NE, y, f, cnt,failed);

  y := wpgx(1,48,64); {e1=e2}
  f := 12.733139842529623306;  {Mathematica}
  testrel(2, NE, y, f, cnt,failed);

  {wpl}
  y := wpgx(5, 1, 0);
  f := 0.6866069655681558007;
  testrel(3, NE, y, f, cnt,failed);

  y := wpgx(50, 1, 0);
  f := 0.5018042255598366777;
  testrel(4, NE, y, f, cnt,failed);

  y := wpgx(500, 1, 0);
  f := 2.794319079242628725;
  testrel(5, NE, y, f, cnt,failed);

  y := wpgx(2, 10, 0);
  f := 137.6154827943146513;
  testrel(6, NE, y, f, cnt,failed);

  y := wpgx(50, 10, 0);
  f := 474.7571842412882356;
  testrel(7, NE4, y, f, cnt,failed);

  y := wpgx(500, 10, 0);
  f := 4.853659971857030600;
  testrel(8, NE4, y, f, cnt,failed);

  y := wpgx(2,0,0);
  f := 0.25;
  testrel(9, NE, y, f, cnt,failed);

  y := wpgx(2,2,0);
  f := 0.9315119326225820305;
  testrel(10, NE, y, f, cnt,failed);

  y := wpgx(2,0,2);
  f := 1.916770245528556302;
  testrel(11, NE, y, f, cnt,failed);

  y := wpgx(1,2,3);
  f := 1.2144337093687324757;
  testrel(12, NE, y, f, cnt,failed);

  y := wpgx(2,4,1);
  f := 4.950267751990133826;
  testrel(13, NE, y, f, cnt,failed);

  y := wpgx(1/128,4,1);
  f := 16384.00001220716430;
  testrel(14, NE, y, f, cnt,failed);

  y := wpgx(1/16,4,1);
  f := 256.0007817957519349;
  testrel(15, NE, y, f, cnt,failed);

  y := wpgx(1/4,4,1);
  f := 16.01264279435187014;
  testrel(16, NE, y, f, cnt,failed);

  y := wpgx(2,1,1);
  f := 1.384027348588571807;
  testrel(17, NE, y, f, cnt,failed);

  y := wpgx(1,2,3);
  f := 1.214433709368732476;
  testrel(18, NE, y, f, cnt,failed);

  y := wpgx(2,2,3);
  f := 6.445549353940210717;
  testrel(19, NE, y, f, cnt,failed);

  y := wpgx(4,2,3);
  f := 1.711910269462141462;
  testrel(20, NE, y, f, cnt,failed);

  y := wpgx(2,3,1); {delta = 0}
  f := 3.183284960632405827;
  testrel(21, NE, y, f, cnt,failed);

  y := wpgx(50,3,1);
  f := 1.000850922654195129;
  testrel(22, NE, y, f, cnt,failed);

  y := wpgx(2,1,2);
  f := 2.658541473896010355;
  testrel(23, NE, y, f, cnt,failed);

  y := wpgx(2,1,-2);
  f := -0.5145833864314375643;
  testrel(24, NE1, y, f, cnt,failed);

  y := wpgx(2,1,-1);
  f := -0.09409585140368366320;
  testrel(25, NE2, y, f, cnt,failed);

  y := wpgx(2,6,-1);  {delta > 0}
  f := 3.963583434450874268;
  testrel(26, NE, y, f, cnt,failed);

  y := wpgx(1/64,6,-1);  {delta > 0}
  f := 4096.000073240059199;
  testrel(27, NE, y, f, cnt,failed);

  y := wpgx(1/16,6,-1);  {delta > 0}
  f := 256.0011713318307081;
  testrel(28, NE, y, f, cnt,failed);

  y := wpgx(10,2,3);
  f := 5.630390556306269352;
  testrel(29, NE2, y, f, cnt,failed);

  y := wpgx(40,2,3);
  f := 2.086001291103351735;
  testrel(30, NE3, y, f, cnt,failed);

  y := wpgx(100,2,3);
  f := 3.158699389454833641;
  testrel(31, NE2, y, f, cnt,failed);

  {equianharmonic case}
  y := wpgx(1/64,0,1);
  f := 4096.000000002128737;
  testrel(32, NE, y, f, cnt,failed);

  y := wpgx(1/32,0,1);
  f := 1024.000000034059797;
  testrel(33, NE, y, f, cnt,failed);

  y := wpgx(1/16,0,1);
  f := 256.0000005449567523;
  testrel(34, NE, y, f, cnt,failed);

  y := wpgx(1/8,0,1);
  f := 64.00000871930812709;
  testrel(35, NE, y, f, cnt,failed);

  y := wpgx(1/4,0,1);
  f := 16.00013950902214234;
  testrel(36, NE, y, f, cnt,failed);

  y := wpgx(1,0,1);
  f := 1.035812586617191909;
  testrel(37, NE, y, f, cnt,failed);

  y := wpgx(2,0,1);
  f := 1.870799354985901704*0.5;
  testrel(38, NE, y, f, cnt,failed);

  y := wpgx(3,0,1);
  f := 278.6309039527418159;
  testrel(39, NE1, y, f, cnt,failed);

  y := wpgx(3.2,0,1);
  f := 50.95348647826213691;
  testrel(40, NE2, y, f, cnt,failed);

  y := wpgx(3.5,0,1);
  f := 5.164471393896693086;
  testrel(41, NE, y, f, cnt,failed);

  y := wpgx(20,0,1);
  f := 0.6446362005489623651;
  testrel(42, NE, y, f, cnt,failed);

  {a few test for wge_im}
  y := wpg_imx(1, 0.5, 0);
  f := -1.025209137571125714;
  testrel(43, NE, y, f, cnt,failed);

  y := wpg_imx(1, 1, 0.5);
  f := -1.032761315967392843;
  testrel(44, NE, y, f, cnt,failed);

  y := wpg_imx(1, 0, 2);
  f := -1.857924845335291353*0.5;
  testrel(45, NE, y, f, cnt,failed);

  y := wpg_imx(1,2,3);
  f := -1.988322550401532842*0.5;
  testrel(46, NE, y, f, cnt,failed);

  y := wpg_imx(2.5,2,3);
  f := 1.022944871169224116;
  testrel(47, NE, y, f, cnt,failed);

  y := wpg_imx(2,4,2);
  f := -0.2147404335837060789;
  testrel(48, NE4, y, f, cnt,failed);

  y := wpg_imx(0.5,4,2);
  f := -4.045729572258790408;
  testrel(49, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_wpg_derx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE  = 4;
  NE1 = 6;
  NE2 = 20;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wpg_derx');

  y := wpg_derx(1,0,1);     {equiharmonian}
  f := -1.856158737883988344;
  testrel(1, NE, y, f, cnt,failed);

  y := wpg_derx(100,0,1);   {equiharmonian}
  f := 2.010812126429562713;
  testrel(2, NE, y, f, cnt,failed);

  y := wpg_derx(-4,0,0);
  f := 0.03125;
  testrel(3, NE, y, f, cnt,failed);

  y := wpg_derx(1, -0.5, 2);
  f := -1.763003174044183231;
  testrel(4, NE, y, f, cnt,failed);

  y := wpg_derx(-1/8, -0.5, 2);
  f := 1024.005691927967364;
  testrel(5, NE, y, f, cnt,failed);

  y := wpg_derx(1/16, -0.5, 2);
  f := -8192.003055244358079;
  testrel(6, NE, y, f, cnt,failed);

  y := wpg_derx(1/32, -0.5, 2);
  f := -65536.00155378065482;
  testrel(7, NE, y, f, cnt,failed);

  y := wpg_derx(1/256, -0.5, 2);
  f := -3.355443200019529547e7;
  testrel(8, NE, y, f, cnt,failed);

  y := wpg_derx(-3,0,1);
  f := -9301.948249601728082;
  testrel(9, NE, y, f, cnt,failed);

  y := wpg_derx(-2,3,1); {delta = 0}
  f := -10.88480183898449443;
  testrel(10, NE1, y, f, cnt,failed);

  y := wpg_derx(1,48,-64); {e1=e2}
  f := -0.9030061850406446371;
  testrel(11, NE, y, f, cnt,failed);

  y := wpg_derx(1,48,64); {e1=e2}
  f := 87.07841471187595541;
  testrel(12, NE, y, f, cnt,failed);

  {lemniscate case}
  y := wpg_derx(5, 1, 0);
  f := -0.7798327505795754523;
  testrel(13, NE, y, f, cnt,failed);

  y := wpg_derx(-50, 1, 0);
  f := 6.023293112221675009e-2;
  testrel(14, NE2, y, f, cnt,failed);

  y := wpg_derx(-500, 1, 0);
  f := -9.191316746291769338;
  testrel(15, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_wpg_invx;
var
  y,f: extended;
  cnt, failed: integer;
const
  NE = 2;
  NE1 = 4;
  NE2 = 40;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','wpg_invx');

  (* u /. FindRoot[WeierstrassP[u, {1, -1}] == 0, {u, 1/10}, WorkingPrecision -> 30]*)
  y := wpg_invx(0.8,0,2);
  f := 1.305375885083806792;
  testrel(1, NE, y, f, cnt,failed);

  y := wpg_invx(2,0,2);
  f := 0.7103461478859869147;
  testrel(2, NE, y, f, cnt,failed);

  y := wpg_invx(5,0,2);
  f := 0.4473415776912397712;
  testrel(3, NE, y, f, cnt,failed);

  y := wpg_invx(0,0,-10);
  f := 1.203596969459945195;
  testrel(4, NE, y, f, cnt,failed);

  y := wpg_invx(0.25,0,0);
  f := 2;
  testrel(5, NE, y, f, cnt,failed);

  y := wpg_invx(-0.5,1,-1);
  f := 2.379880191665208495;
  testrel(6, NE1, y, f, cnt,failed);

  y := wpg_invx(0,1,-1);
  f := 1.907984049121513880;
  testrel(7, NE, y, f, cnt,failed);

  y := wpg_invx(1,1,-1);
  f := 1.0072923743245416045;
  testrel(8, NE1, y, f, cnt,failed);

  y := wpg_invx(1,1,0);
  f := 1.028056801052126733;
  testrel(9, NE, y, f, cnt,failed);

  y := wpg_invx(10,1,0);
  f := 0.3163069054282973787;
  testrel(10, NE, y, f, cnt,failed);

  y := wpg_invx(1,0,1);
  f := 1.019969358038128583;
  testrel(11, NE1, y, f, cnt,failed);

  y := wpg_invx(10,0,1);
  f := 0.3162334135114345463;
  testrel(12, NE, y, f, cnt,failed);

  y := wpg_invx(1.125,3,1); {delta = 0}
  f := 1.053085793963456538;
  testrel(13, NE1, y, f, cnt,failed);

  y := wpg_invx(0.75,2,0);  {delta > 0}
  f := 1.355031378626586435;
  testrel(14, NE, y, f, cnt,failed);

  y := wpg_invx(4,4,1);     {delta > 0}
  f := 0.5033580305553346283;
  testrel(15, NE, y, f, cnt,failed);

  y := wpg_invx(2,4,1);     {delta > 0}
  f := 0.7290483255921761057;
  testrel(16, NE, y, f, cnt,failed);

  y := wpg_invx(-2.5,20,-20);  {extreme}
  f := 2.0801623945988917796;
  testrel(17, NE2, y, f, cnt,failed);

  y := wpg_invx(0,20,-20);     {extreme}
  f := 1.6026399727295061399;
  testrel(18, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



end.
