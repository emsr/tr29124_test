{Part 4 of regression test for SPECFUN unit  (c) 2010  W.Ehrhardt}

unit t_sfd4;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_erf;
procedure test_erfc;
procedure test_erfce;
procedure test_inerfc;
procedure test_erfg;
procedure test_erfi;
procedure test_dawson;
procedure test_dawson2;
procedure test_erfinv;
procedure test_erfcinv;
procedure test_erfce_inv;
procedure test_erfi_inv;
procedure test_erf_z;
procedure test_erf_p;
procedure test_erf_q;
procedure test_fresnel;
procedure test_fresnelfg;
procedure test_erfh2;
procedure test_OwenT;


implementation

uses
  amath, specfun, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_erf;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erf');
  {-Inf .. -7 .. -1 .. 1 .. 7 .. Inf}

  y := erf(NegInf_x);
  f := -1;
  testabs( 1, 0, y, f, cnt,failed);

  y := erf(-100);
  f := -1;
  testrel( 2, NE, y, f, cnt,failed);

  y := erf(-6.0);
  f := -0.99999999999999997848;
  testrel( 3, NE, y, f, cnt,failed);

  y := erf(-2.0);
  f := -0.99532226501895273416;
  testrel( 4, NE, y, f, cnt,failed);

  y := erf(-1.1);
  f := -0.88020506957408169977;
  testrel( 5, NE, y, f, cnt,failed);

  y := erf(-1.0);
  f := -0.84270079294971486934;
  testrel( 6, NE, y, f, cnt,failed);

  y := erf(-0.9);
  f := -0.79690821242283212852;
  testrel( 7, NE, y, f, cnt,failed);

  y := erf(-0.5);
  f := -0.52049987781304653768;
  testrel( 8, NE, y, f, cnt,failed);

  y := erf(-0.125);
  f := -0.14031620480133381739;
  testrel( 9, NE, y, f, cnt,failed);

  y := erf(-0.0078125);
  f := -0.88152828951791887128e-2;
  testrel(10, NE, y, f, cnt,failed);

  y := erf(0.0);
  f := 0;
  testrel(11, NE, y, f, cnt,failed);

  y := erf(0.0078125);
  f := 0.88152828951791887128e-2;
  testrel(12, NE, y, f, cnt,failed);

  y := erf(0.125);
  f := 0.14031620480133381739;
  testrel(13, NE, y, f, cnt,failed);

  y := erf(0.5);
  f := 0.52049987781304653768;
  testrel(14, NE, y, f, cnt,failed);

  y := erf(0.9);
  f := 0.79690821242283212852;
  testrel(15, NE, y, f, cnt,failed);

  y := erf(1.0);
  f := 0.84270079294971486934;
  testrel(16, NE, y, f, cnt,failed);

  y := erf(1.1);
  f := 0.88020506957408169977;
  testrel(17, NE, y, f, cnt,failed);

  y := erf(2.0);
  f := 0.99532226501895273416;
  testrel(18, NE, y, f, cnt,failed);

  y := erf(6.0);
  f := 0.99999999999999997848;
  testrel(19, NE, y, f, cnt,failed);

  y := erf(100);
  f := 1;
  testrel(20, NE, y, f, cnt,failed);

  y := erf(PosInf_x);
  f := 1;
  testabs(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfc;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfc');
  {-Inf .. -7 .. -1 .. 1 .. 7 .. 107 .. Inf}

  y := erfc(NegInf_x);
  f := 2;
  testabs( 1, 0, y, f, cnt,failed);

  y := erfc(-100);
  f := 2;
  testrel( 2, NE, y, f, cnt,failed);

  y := erfc(-6.0);
  f := 1.9999999999999999785;
  testrel( 3, NE, y, f, cnt,failed);

  y := erfc(-2.0);
  f := 1.9953222650189527342;
  testrel( 4, NE, y, f, cnt,failed);

  y := erfc(-1.1);
  f := 1.8802050695740816998;
  testrel( 5, NE, y, f, cnt,failed);

  y := erfc(-1.0);
  f := 1.8427007929497148693;
  testrel( 6, NE, y, f, cnt,failed);

  y := erfc(-0.9);
  f := 1.7969082124228321285;
  testrel( 7, NE, y, f, cnt,failed);

  y := erfc(-0.5);
  f := 1.5204998778130465377;
  testrel( 8, NE, y, f, cnt,failed);

  y := erfc(-0.125);
  f := 1.1403162048013338174;
  testrel( 9, NE, y, f, cnt,failed);

  y := erfc(-0.0078125);
  f := 1.0088152828951791887;
  testrel(10, NE, y, f, cnt,failed);

  y := erfc(0.0);
  f := 1;
  testrel(11, NE, y, f, cnt,failed);

  y := erfc(0.0078125);
  f := 0.99118471710482081129;
  testrel(12, NE, y, f, cnt,failed);

  y := erfc(0.125);
  f := 0.85968379519866618261;
  testrel(13, NE, y, f, cnt,failed);

  y := erfc(0.5);
  f := 0.47950012218695346232;
  testrel(14, NE, y, f, cnt,failed);

  y := erfc(0.9);
  f := 0.20309178757716787148;
  testrel(15, NE, y, f, cnt,failed);

  y := erfc(1.0);
  f := 0.15729920705028513066;
  testrel(16, NE, y, f, cnt,failed);

  y := erfc(1.125);
  f := 0.11161176829829223593;
  testrel(17, NE, y, f, cnt,failed);

  y := erfc(2.0);
  f := 0.46777349810472658379e-2;
  testrel(18, NE, y, f, cnt,failed);

  y := erfc(4.0);
  f := 0.15417257900280018852e-7;
  testrel(19, NE, y, f, cnt,failed);

  y := erfc(6.0);
  f := 0.21519736712498913117e-16;
  testrel(20, NE, y, f, cnt,failed);

  y := erfc(10.0);
  f := 0.20884875837625447570e-44;
  testrel(21, NE, y, f, cnt,failed);

  y := erfc(20.0);
  f := 0.53958656116079009289e-175;
  testrel(22, NE, y, f, cnt,failed);

  y := erfc(26.5);
  f := 0.22109076642637342759e-306;
  testrel(23, NE, y, f, cnt,failed);

  y := erfc(PosInf_x);
  f := 0;
  testabs(24, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfce;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfce');

  x := 1e-6;
  f := 0.999998871621832904;
  y := erfce(x);
  testrel( 1, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.8732218450821508096;
  y := erfce(x);
  testrel( 2, NE, y, f, cnt,failed);

  x := 1/3;  {Fix311}
  f := 0.7122528886000577139;
  y := erfce(x);
  testrel( 3, NE, y, f, cnt,failed);

  x := 1;
  f := 0.4275835761558070044;
  y := erfce(x);
  testrel( 4, NE, y, f, cnt,failed);

  x := 2;
  f := 0.2553956763105057439;
  y := erfce(x);
  testrel( 5, NE, y, f, cnt,failed);

  x := 10.0;
  f := 0.5614099274382258586e-1;
  y := erfce(x);
  testrel( 6, NE, y, f, cnt,failed);

  x := 100.0;
  f := 0.5641613782989432904e-2;
  y := erfce(x);
  testrel( 7, NE, y, f, cnt,failed);

  x := 130.0;
  f := 0.4339791484842787238e-2;
  y := erfce(x);
  testrel( 8, NE, y, f, cnt,failed);

  x := 5000.0;
  f := 0.1128379144527930586e-3;
  y := erfce(x);
  testrel( 9, NE, y, f, cnt,failed);

  x := 70000.0;
  f := 0.8059851192716941733e-5;
  y := erfce(x);
  testrel(10, NE, y, f, cnt,failed);

  x := -1e-15;
  f := 1.000000000000001128;
  y := erfce(x);
  testrel(11, NE, y, f, cnt,failed);

  x := -0.25;
  f := 1.358642370104722115;
  y := erfce(x);
  testrel(12, NE, y, f, cnt,failed);

  x := -1;
  f := 5.008980080762283466;
  y := erfce(x);
  testrel(13, NE, y, f, cnt,failed);

  x := -8/7;
  f := 6.992173861723201774;
  y := erfce(x);
  testrel(14, NE, y, f, cnt,failed);

  x := -2;
  f := 108.9409043899779724;
  y := erfce(x);
  testrel(15, NE, y, f, cnt,failed);

  x := -6.75;
  f := 0.1226231102139051807e21;
  y := erfce(x);
  testrel(16, NE, y, f, cnt,failed);

  x := -26.5;
  f := 0.1924553162418568809e306;
  y := erfce(x);
  testrel(17, NE, y, f, cnt,failed);

  x := -106.5637380121098418;
  f := PosInf_d;
  y := erfce(x);
  testabs(18, 0, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfg;
var
  x,y,f,p: double;
  cnt,failed,i: integer;
const
  NE  = 1;
  NE2 = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfg');

  {Maple: gerf := (p,x) -> int(exp(-t^p), t=0..x);}
  {       gerf := (p,x) -> (GAMMA(1/p)-GAMMA(1/p,x^p))/p;}
  {Pari:  gerf(p,x) = intnum(t=0,x,exp(-t^p))}
  {       gerf(p,x) = incgamc(1/p, x^p)/p}

  p := 0.5;
  y := erfg(p,1/1024);
  f := 0.9564538925403311883e-3;
  testrel( 1, NE, y, f, cnt,failed);

  y := erfg(p,0.125);
  f := 0.9910074638765149504e-1;
  testrel( 2, NE, y, f, cnt,failed);

  y := erfg(p,0.5);
  f := 0.3165581866568181266;
  testrel( 3, NE, y, f, cnt,failed);

  y := erfg(p,1);
  f := 0.5284822353142307136;
  testrel( 4, NE, y, f, cnt,failed);

  y := erfg(p,2);
  f := 0.8261285649781240111;
  testrel( 5, NE, y, f, cnt,failed);

  y := erfg(p,10);
  f := 1.647628069579945710;
  testrel( 6, NE, y, f, cnt,failed);

  y := erfg(p,1000);
  f := 1.999999999998795093;
  testrel( 7, NE, y, f, cnt,failed);


  p := 5;
  y := erfg(p,1/1024);
  f := 2*0.4882812499999999277e-3;
  testrel( 8, NE, y, f, cnt,failed);

  y := erfg(p,0.125);
  f := 0.12499936422241396436;
  testrel( 9, NE, y, f, cnt,failed);

  y := erfg(p,0.5);
  f := 0.4974178699312369020;
  testrel(10, NE, y, f, cnt,failed);

  y := erfg(p,1);
  f := 0.8700746676858926573;
  testrel(11, NE, y, f, cnt,failed);

  y := erfg(p,2);
  f := 0.9181687423997604561;
  testrel(12, NE, y, f, cnt,failed);

  y := erfg(p,10);
  f := 0.9181687423997606106;
  testrel(13, NE, y, f, cnt,failed);


  p := 1.5;
  y := erfg(p,1/1024);
  f := 2*0.4882752895923654593e-3;
  testrel(14, NE, y, f, cnt,failed);

  y := erfg(p,0.125);
  f := 0.12282048474733161436;
  testrel(15, NE, y, f, cnt,failed);

  y := erfg(p,0.5);
  f := 0.4364761380885268077;
  testrel(16, NE, y, f, cnt,failed);

  y := erfg(p,1);
  f := 0.6997923277614944777;
  testrel(17, NE, y, f, cnt,failed);

  y := erfg(p,2);
  f := 0.8772524660847196675;
  testrel(18, NE, y, f, cnt,failed);

  y := erfg(p,10);
  f := 0.9027452929509297575;
  testrel(19, NE, y, f, cnt,failed);

  y := erfg(p,1000);
  f := 0.9027452929509336113;
  testrel(20, NE, y, f, cnt,failed);

  {pari}
  p := 0.01;
  y := erfg(p,1/1024);
  f := 0.0003877209165992133333;
  testrel(21, NE, y, f, cnt,failed);

  y := erfg(p,0.125);
  f := 0.04740070277221238905;
  testrel(22, NE, y, f, cnt,failed);

  y := erfg(p,0.5);
  f := 0.1870537286914434240;
  testrel(23, NE, y, f, cnt,failed);

  y := erfg(p,1);
  f := 0.3715578714528098103;
  testrel(24, NE, y, f, cnt,failed);

  y := erfg(p,2);
  f := 0.7380162212946553285;
  testrel(25, NE, y, f, cnt,failed);

  y := erfg(p,10);
  f := 3.630877526710845930;
  testrel(26, NE, y, f, cnt,failed);

  y := erfg(p,1000);
  f := 346.1598375906238578;
  testrel(27, NE, y, f, cnt,failed);

  y := erfg(p,1e100);
  f := 5.038299480618184746e95;
  testrel(28, 3, y, f, cnt,failed);    {!!??, extended error < 4.4e-19!!}

  {selected values}
  y := erfg(20,0);
  f := 0;
  testrel(29, NE, y, f, cnt,failed);

  y := erfg(0,20);
  f := 7.357588823428846432;
  testrel(30, NE, y, f, cnt,failed);

  y := erfg(1,20);
  f := 2*0.4999999989694231888;
  testrel(31, NE, y, f, cnt,failed);

  y := erfg(0.25,10);
  f := 2.525854435942086756;
  testrel(32, NE, y, f, cnt,failed);

  y := erfg(10,0.125);
  f := 0.1249999999894167889;
  testrel(33, NE, y, f, cnt,failed);

  y := erfg(20,0.25);
  f := 0.2499999999999891727;
  testrel(34, NE, y, f, cnt,failed);

  y := erfg(40,0.5);
  f := 0.4999999999999889086;
  testrel(35, NE, y, f, cnt,failed);

  y := erfg(100,1-1/1024);
  f := 0.991745400651105407;
  testrel(36, NE, y, f, cnt,failed);

  y := erfg(10000,1-1/1024);
  f := 2*0.4995117158972298877;
  testrel(37, NE, y, f, cnt,failed);

  {expint3 tests}
  for i:=0 to 63 do begin
    x := i/16;
    y := erfg(3,x);
    f := expint3(x);
    testrel(100+i, NE2, y, f, cnt,failed);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_inerfc;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','inerfc');

  y := inerfc(3,1e-40);
  f := 0.9403159725795938116e-1;
  testrel(1, NE, y, f, cnt,failed);

  y := inerfc(250,1e-40);
  f := 0.2935791617973459071e-284;
  testrel(2, NE, y, f, cnt,failed);

  y := inerfc(3,0.5);
  f := 0.2775132137764753636e-1;
  testrel(3, NE, y, f, cnt,failed);

  y := inerfc(10,20);
  f := 0.2480339545677166807e-17;
  testrel(4, NE, y, f, cnt,failed);

  y := inerfc(10,2);
  f := 0.4519702487098470796e-8;
  testrel(5, NE, y, f, cnt,failed);

  y := inerfc(10,1);
  f := 0.1305302483913995620e-6;
  testrel(6, NE, y, f, cnt,failed);

  y := inerfc(10,0.5);
  f := 0.9245065129805797129e-6;
  testrel(7, NE, y, f, cnt,failed);

  y := inerfc(10,0.25);
  f := 0.26649953573332452744e-5;
  testrel(8, NE, y, f, cnt,failed);

  y := inerfc(10,0.125);
  f := 0.4622606908901843597e-5;
  testrel(9, NE, y, f, cnt,failed);

  y := inerfc(10,1/1024);
  f := 0.8101666489669277473e-5;
  testrel(10, NE, y, f, cnt,failed);

  y := inerfc(100,8);
  f := 0.3327723273692355744e-132;
  testrel(11, NE, y, f, cnt,failed);

  y := inerfc(100,2);
  f := 0.8417913154612625873e-106;
  testrel(12, NE, y, f, cnt,failed);

  y := inerfc(100,0.5);
  f := 0.2448087104431131554e-97;
  testrel(13, NE, y, f, cnt,failed);

  y := inerfc(100,0.25);
  f := 0.7728184538695565628e-96;
  testrel(14, NE, y, f, cnt,failed);

  y := inerfc(100,1/1024);
  f := 0.2558072521643705004e-94;
  testrel(15, NE, y, f, cnt,failed);

  y := inerfc(100,1e-10);
  f := 0.2593734749448854344e-94;
  testrel(16, NE, y, f, cnt,failed);

  y := inerfc(5,0.125);
  f := 0.6252579526727030871e-2;
  testrel(17, NE, y, f, cnt,failed);

  y := inerfc(20,0.125);
  f := 0.1189350337542980588e-12;
  testrel(18, NE, y, f, cnt,failed);

  y := inerfc(100,0.125);
  f := 0.4442699570076624848e-95;
  testrel(19, NE, y, f, cnt,failed);

  y := inerfc(250,0.125);
  f := 0.1803034248451498754e-285;
  testrel(20, NE, y, f, cnt,failed);

  y := inerfc(5,0.75);
  f := 0.9967894387400569614e-3;
  testrel(21, NE, y, f, cnt,failed);

  y := inerfc(20,0.75);
  f := 0.2815855137500476563e-14;
  testrel(22, NE, y, f, cnt,failed);

  y := inerfc(50,0.75);
  f := 0.4007073297196590972e-43;
  testrel(23, NE, y, f, cnt,failed);

  y := inerfc(125,0.75);
  f := 0.8592345382450687318e-129;
  testrel(24, NE, y, f, cnt,failed);

  y := inerfc(250,0.75);
  f := 0.1984987852267995133e-291;
  testrel(25, NE, y, f, cnt,failed);

  y := inerfc(5,50);     {AE}
  f := 0.1123656973372330445e-11;
  testrel(26, NE, y, f, cnt,failed);

  y := inerfc(20,50);    {AE}
  f := 0.1077656273122095687e-41;
  testrel(27, NE, y, f, cnt,failed);

  y := inerfc(100,50);   {AE}
  f := 0.4110694394268175203e-202;
  testrel(28, NE, y, f, cnt,failed);

  y := inerfc(150,50);   {CF}
  f := 0.1215088638985536402e-302;
  testrel(29, NE, y, f, cnt,failed);

  y := inerfc(10,1000);
  f := 0.5509482087061094431e-36;
  testrel(30, NE, y, f, cnt,failed);

  y := inerfc(90,1000);
  f := 0.4547958756280478593e-300;
  testrel(31, NE, y, f, cnt,failed);

  {x<0}
  y := inerfc(-1,-2);
  f := 0.2066698535409205386e-1;
  testrel(32, NE, y, f, cnt,failed);

  y := inerfc(0,-2);
  f := 1.9953222650189527342;
  testrel(33, NE, y, f, cnt,failed);

  y := inerfc(1,-2);
  f := 4.000978022714951495;
  testrel(34, NE, y, f, cnt,failed);

  y := inerfc(5,-0.5);
  f := 0.4374437542410532658e-1;
  testrel(35, NE, y, f, cnt,failed);

  y := inerfc(5,-2);
  f := 1.325001048378169994;
  testrel(36, NE, y, f, cnt,failed);

  y := inerfc(5,-10);
  f := 1750.625;
  testrel(37, NE, y, f, cnt,failed);

  y := inerfc(5,-300);
  f := 40502250018.75;
  testrel(38, NE, y, f, cnt,failed);

  y := inerfc(20,-0.5);
  f := 0.5711213876890956573e-11;
  testrel(39, NE, y, f, cnt,failed);

  y := inerfc(20,-2);
  f := 0.1557830391267366734e-7;
  testrel(40, NE, y, f, cnt,failed);

  y := inerfc(20,-10);
  f := 196.8746175638803787;
  testrel(41, NE, y, f, cnt,failed);

  y := inerfc(20,-300);
  f := 0.2869385160990788617e32;
  testrel(42, NE, y, f, cnt,failed);

  y := inerfc(1000,-500);
  f := 0.1257148984018811968e133;
  testrel(43, NE, y, f, cnt,failed);

  y := inerfc(100,-1000);
  f := 0.2148330859419028936e143;
  testrel(44, NE, y, f, cnt,failed);

  y := inerfc(5000,-2000);
  f := 0.3183024088162901648e181;
  testrel(45, NE, y, f, cnt,failed);

  y := inerfc(2000,-1000);
  f := 0.1636906036218412122e266;
  testrel(46, NE, y, f, cnt,failed);

  y := inerfc(30000,-11000);
  f := 0.2806686227497745515e-44;
  testrel(47, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_erfi;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfi');

  {erfi := x -> erf(I*x)/I;}
  x := 1e-20;
  y := erfi(x);
  f := 0.1128379167095512574e-19;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1e-5;
  y := erfi(x);
  f := 0.1128379167133125213e-4;
  testrel( 2, NE, y, f, cnt,failed);

  x := -0.125;
  y := erfi(x);
  f := -0.14178547413026538933;
  testrel( 3, NE, y, f, cnt,failed);

  x := 0.5;
  y := erfi(x);
  f := 0.6149520946965109808;
  testrel( 4, NE, y, f, cnt,failed);

  x := 0.75;
  y := erfi(x);
  f := 1.035757284411962968;
  testrel( 5, NE, y, f, cnt,failed);

  x := -1;
  y := erfi(x);
  f := -1.650425758797542876;
  testrel( 6, NE, y, f, cnt,failed);

  x := 2;
  y := erfi(x);
  f := 18.56480241457555260;
  testrel( 7, NE, y, f, cnt,failed);

  x := 4;
  y := erfi(x);
  f := 1296959.730717639232;
  testrel( 8, NE, y, f, cnt,failed);

  x := -10;
  y := erfi(x);
  f := -0.1524307422708669699e43;
  testrel( 9, NE, y, f, cnt,failed);

  x := 26.5;
  y := erfi(x);
  f := 2.0501652832248793153e303;
  testrel(10, NE, y, f, cnt,failed);

  x := 100;
  y := erfi(x);
  f := PosInf_d;
  testabs(11, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_dawson;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','dawson');

  x := 0;
  y := dawson(x);
  f := 0;
  testrel( 1, NE, y, f, cnt,failed);

  x := succd(0); {x=0.5^16445}
  y := dawson(x);
  f := 0.49406564584124654418e-323;
  testrel( 2, NE, y, f, cnt,failed);

  x := 1e-100;
  y := dawson(x);
  f := 1e-100;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1e-6;
  y := dawson(x);
  f := 0.99999999999933333333e-6;
  testrel( 4, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := dawson(x);
  f := 0.97656187911852043720e-3;
  testrel( 5, NE, y, f, cnt,failed);

  x := 0.125;
  y := dawson(x);
  f := 0.12370601848283973394;
  testrel( 6, NE, y, f, cnt,failed);

  x := 0.5;
  y := dawson(x);
  f := 0.42443638350202229593;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.92413887300459176701;
  y := dawson(x);
  f := 0.54104422463518169847;
  testrel( 8, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := dawson(x);
  f := 0.53815344009396028924;
  testrel( 9, NE, y, f, cnt,failed);

  x := 1.0;
  y := dawson(x);
  f := 0.53807950691276841914;
  testrel(10, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := dawson(x);
  f := 0.53800469268824945169;
  testrel(11, NE, y, f, cnt,failed);

  x := 3.9990234375;
  y := dawson(x);
  f := 0.12938197933297510231;
  testrel(12, NE, y, f, cnt,failed);

  x := 4.0;
  y := dawson(x);
  f := 0.12934800123600511559;
  testrel(13, NE, y, f, cnt,failed);

  x := 4.0009765625;
  y := dawson(x);
  f := 0.12931404180823832109;
  testrel(14, NE, y, f, cnt,failed);

  x := 3.0e9;
  y := dawson(x);
  f := 0.16666666666666666667e-9;
  testrel(15, NE, y, f, cnt,failed);

  x := 4.0e9;
  y := dawson(x);
  f := 0.12500000000000000000e-9;
  testrel(16, NE, y, f, cnt,failed);

  x := 1.0e307;
  y := dawson(x);
  f := 5e-308;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_dawson2;
var
  p,x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','dawson2');

  {Pari: gd(p,x)=intnum(t=0,x,exp(t^p))*exp(-x^p)}

  {Standard Dawson}
  p := 2;
  y := dawson2(p,1);
  f := 0.5380795069127684191;
  testrel(1, NE, y, f, cnt,failed);

  y := dawson2(p,3);
  f := 0.1782710306105582873;
  testrel(2, NE, y, f, cnt,failed);

  y := dawson2(p,5);
  f := 0.1021340744242768355;
  testrel(3, NE, y, f, cnt,failed);

  y := dawson2(p,100);
  f := 0.5000250037509378283e-2;
  testrel(4, NE, y, f, cnt,failed);

  y := dawson2(p,1e300);
  f := 0.5e-300;
  testrel(5, NE, y, f, cnt,failed);

  y := dawson2(p,1/16);
  f := 0.6233749361289893996e-1;
  testrel(6, NE, y, f, cnt,failed);

  y := dawson2(p,1e-5);
  f := 0.9999999999333333333e-5;
  testrel(7, NE, y, f, cnt,failed);

  y := dawson2(p,1e-10);
  f := 1e-10;
  testrel(8, NE, y, f, cnt,failed);

  y := dawson2(p,1e-15);
  f := 1e-15;
  testrel(9, NE, y, f, cnt,failed);

  {Special case x=1}
  x := 1;
  y := dawson2(1+1/1024,x);
  f := 0.6319767998075155365;
  testrel(10, NE, y, f, cnt,failed);

  y := dawson2(10,x);
  f := 0.4125034088017899048;
  testrel(11, NE, y, f, cnt,failed);

  y := dawson2(100,x);
  f := 0.3726859445420274430;
  testrel(12, NE, y, f, cnt,failed);

  y := dawson2(1000,x);
  f := 0.3683638488980294743;
  testrel(13, NE, y, f, cnt,failed);

  y := dawson2(1e10,x);
  f := 0.3678794412199252323;
  testrel(14, NE, y, f, cnt,failed);

  y := dawson2(1e300,x);
  f := 0.3678794411714423216;
  testrel(15, NE, y, f, cnt,failed);

  y := dawson2(1/1024,x);
  f := 0.9990253402056163500;
  testrel(16, NE, y, f, cnt,failed);

  y := dawson2(1e-9,x);
  f := 0.999999999000000002;
  testrel(17, NE, y, f, cnt,failed);

  y := dawson2(1e-10,x);
  f := 0.9999999999;
  testrel(18, NE, y, f, cnt,failed);

  y := dawson2(1e-11,x);
  f := 0.99999999999;
  testrel(19, NE, y, f, cnt,failed);

  y := dawson2(1e-16,x);
  f := 0.9999999999999999;
  testrel(20, NE, y, f, cnt,failed);

  {Other cases}
  y := dawson2(1/80,80);
  f := 78.97000764508719285;
  testrel(21, NE, y, f, cnt,failed);

  y := dawson2(1.5,100);
  f := 0.6668891858788577765e-1;
  testrel(22, NE, y, f, cnt,failed);

  y := dawson2(1,4);
  f := 0.981684361111265820;
  testrel(23, NE, y, f, cnt,failed);

  y := dawson2(1,1e-10);
  f := 0.99999999995e-10;
  testrel(24, NE, y, f, cnt,failed);

  y := dawson2(1/2,100);
  f := 18.00009079985952497;
  testrel(25, NE, y, f, cnt,failed);

  y := dawson2(1/2,1e10);
  f := 199998.0;
  testrel(26, NE, y, f, cnt,failed);

  y := dawson2(0.1,10);
  f := 8.964926863611165808;
  testrel(27, NE, y, f, cnt,failed);

  y := dawson2(0.01,100);
  f := 98.9737752892581009;
  testrel(28, NE, y, f, cnt,failed);

  y := dawson2(0.01,1e10);
  f := 9876873773.46623806;
  testrel(29, NE, y, f, cnt,failed);

  y := dawson2(0.01,1e100);
  f := 9.09837270131716558e99;
  testrel(30, NE, y, f, cnt,failed);

  y := dawson2(1/1024,1/2);
  f := 0.49951299954466997501;
  testrel(31, NE, y, f, cnt,failed);

  y := dawson2(1e-5,1e20);
  f := 9.99989995593902475e19;
  testrel(32, NE, y, f, cnt,failed);

  y := dawson2(1e-10,1e10);
  f := 9999999999.0;
  testrel(33, NE, y, f, cnt,failed);

  y := dawson2(25,2);
  f := 0.2384185859227731617e-8;
  testrel(34, NE, y, f, cnt,failed);

  y := dawson2(5,50);
  f := 0.3200000008192000047e-7;
  testrel(35, NE, y, f, cnt,failed);

  y := dawson2(4,100);
  f := 0.2500000018750000328e-6;
  testrel(36, NE, y, f, cnt,failed);

  y := dawson2(40,1.25);
  f := 0.4154375964424212495e-5;
  testrel(37, NE, y, f, cnt,failed);

  y := dawson2(40,1.5);
  f := 0.3391415055475334289e-8;
  testrel(38, NE, y, f, cnt,failed);

  y := dawson2(300,1+1/16);
  f := 4.4722500002707202404e-11;
  testrel(39, NE, y, f, cnt,failed);

  y := dawson2(300,1-1/16);
  f := 0.937499996353225120;
  testrel(40, NE, y, f, cnt,failed);

  y := dawson2(500,3);
  f := 0.1650151773721882811e-240;
  testrel(41, NE, y, f, cnt,failed);

  y := dawson2(200,10);
  f := 5.0e-202;
  testrel(42, NE, y, f, cnt,failed);

  y := dawson2(200,8);
  f := 0.9639679460411536471e-182;
  testrel(43, NE, y, f, cnt,failed);

  y := dawson2(1000,1-1/16);
  f := 0.9375;
  testrel(44, NE, y, f, cnt,failed);

  y := dawson2(1000,1+1/16);
  f := 0.4981845059124721260e-29;
  testrel(45, NE, y, f, cnt,failed);

  y := dawson2(10000,1+1/16);
  f := 0.545684622222586877e-267;
  testrel(46, NE, y, f, cnt,failed);

  y := dawson2(1000,2);
  f := 1.866527237006437758e-304;
  testrel(47, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfinv;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erf_inv');

  x := 0;
  y := erf_inv(x);
  f := 0;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1e-300;
  y := erf_inv(x);
  f := 0.886226925452758013649e-300;
  testrel( 2, NE, y, f, cnt,failed);

  x := -1e-100;
  y := erf_inv(x);
  f := -0.886226925452758013649e-100;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1e-40;
  y := erf_inv(x);
  f := 0.886226925452758013649e-40;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1e-6;
  y := erf_inv(x);
  f := 0.8862269254529900273156e-6;
  testrel( 5, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := erf_inv(x);
  f := 0.86545619796713755345e-3;
  testrel( 6, NE, y, f, cnt,failed);

  x := -0.125;
  y := erf_inv(x);
  f := -0.1112354518409499591471840;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.5;
  y := erf_inv(x);
  f := 0.47693627620446987338;
  testrel( 8, NE, y, f, cnt,failed);

  x := -0.75;
  y := erf_inv(x);
  f := -0.8134198475976185417;
  testrel( 9, NE, y, f, cnt,failed);

  x := 0.9375;
  y := erf_inv(x);
  f := 1.31715033498613074888;
  testrel(10, NE, y, f, cnt,failed);

  x := -0.9990234375;
  y := erf_inv(x);
  f := -2.3314677736219476723;
  testrel(11, NE, y, f, cnt,failed);

  x := 1.0-ldexp(1,-20);
  y := erf_inv(x);
  f := 3.465505025803330717256;
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfcinv;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfc_inv');

  x := succd(0); {x=0.5^1074}
  y := erfc_inv(x);
  f := 27.21329321081294881531;
  testrel( 1, NE, y, f, cnt,failed);

  x := 1e-300;
  y := erfc_inv(x);
  f := 26.209469960516123885998;
  testrel( 2, NE, y, f, cnt,failed);

  x := 1e-200;
  y := erfc_inv(x);
  f := 21.3747830490262576317;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1e-100;
  y := erfc_inv(x);
  f := 15.06557470259264570440;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1e-40;
  y := erfc_inv(x);
  f := 9.448789766720858262503;
  testrel( 5, NE, y, f, cnt,failed);

  x := 1e-10;
  y := erfc_inv(x);
  f := 4.572824967389485278741;
  testrel( 6, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := erfc_inv(x);
  f := 2.33146777362194767232;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.125;
  y := erfc_inv(x);
  f := 1.08478704006928314131;
  testrel( 8, NE, y, f, cnt,failed);

  x := 0.2;
  y := erfc_inv(x);
  f := 0.90619380243682322007;
  testrel( 9, NE, y, f, cnt,failed);

  x := 0.5;
  y := erfc_inv(x);
  f := 0.47693627620446987338;
  testrel(10, NE, y, f, cnt,failed);

  x := 0.75;
  y := erfc_inv(x);
  f := 0.225312055012178104725;
  testrel(11, NE, y, f, cnt,failed);

  x := 0.9375;
  y := erfc_inv(x);
  f := 0.55445948772782020299e-1;
  testrel(12, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := erfc_inv(x);
  f := 0.86545619796713755345e-3;
  testrel(13, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := erfc_inv(x);
  f := -0.86545619796713755345e-3;
  testrel(14, NE, y, f, cnt,failed);

  x := 1.125;
  y := erfc_inv(x);
  f := -0.111235451840949959147;
  testrel(15, NE, y, f, cnt,failed);

  x := 1.25;
  y := erfc_inv(x);
  f := -0.225312055012178104725;
  testrel(16, NE, y, f, cnt,failed);

  x := 1.5;
  y := erfc_inv(x);
  f := -0.47693627620446987338;
  testrel(17, NE, y, f, cnt,failed);

  x := 1.75;
  y := erfc_inv(x);
  f := -0.81341984759761854169;
  testrel(18, NE, y, f, cnt,failed);

  x := 1.875;
  y := erfc_inv(x);
  f := -1.08478704006928314131;
  testrel(19, NE, y, f, cnt,failed);

  x := 1.9375;
  y := erfc_inv(x);
  f := -1.31715033498613074888;
  testrel(20, NE, y, f, cnt,failed);

  x := 1.9990234375;
  y := erfc_inv(x);
  f := -2.3314677736219476723;
  testrel(21, NE, y, f, cnt,failed);

  x := 2.0-ldexp(1,-20);
  y := erfc_inv(x);
  f := -3.4655050258033307173;
  testrel(22, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfce_inv;
var
  x,y,f,h,u,w: double;
  n,cnt, failed: integer;
const
  NE  = 1;
  NE1 = 16;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfce_inv');

  y := erfce_inv(0.01);
  f := 56.41009747672935547;
  testrel(1, NE, y, f, cnt,failed);

  y := erfce_inv(0.1);
  f := 5.554585892541128945;
  testrel(2, NE, y, f, cnt,failed);

  y := erfce_inv(0.5);
  f := 0.7690797710613142052;
  testrel(3, NE, y, f, cnt,failed);

  y := erfce_inv(0.875);
  f := 0.1230493217177805782;
  testrel(4, NE, y, f, cnt,failed);

  y := erfce_inv(1);
  f := 0;
  testabs(5, 0, y, f, cnt,failed);

  y := erfce_inv(1.75);
  f := -0.4289449529676623929;
  testrel(6, NE, y, f, cnt,failed);

  y := erfce_inv(2);
  f := -0.5151980774824833669;
  testrel(7, NE, y, f, cnt,failed);

  y := erfce_inv(3);
  f := -0.7494365560122995770;
  testrel(8, NE, y, f, cnt,failed);

  y := erfce_inv(4);
  f := -0.8952923401402964116;
  testrel(9, NE, y, f, cnt,failed);

  y := erfce_inv(100);
  f := -1.978533994229951307;
  testrel(10, NE, y, f, cnt,failed);

  y := erfce_inv(10000);
  f := -2.918426210428600120;
  testrel(11, NE, y, f, cnt,failed);

  y := erfce_inv(1e100);
  f := -15.15141452534530266;
  testrel(12, NE, y, f, cnt,failed);

  y := erfce_inv(1e300);
  f := -26.26941911648702140;
  testrel(13, NE, y, f, cnt,failed);

  y := erfce_inv(1e-5);
  f := 56418.95834591335944;
  testrel(14, NE, y, f, cnt,failed);

  y := erfce_inv(1e-9);
  f := 564189583.5477562861;
  testrel(15, NE, y, f, cnt,failed);

  y := erfce_inv(1e-300);
  f := 0.5641895835477562869e300;
  testrel(16, NE, y, f, cnt,failed);

  inc(cnt);
  u := 0;
  for n:=-25 to 25 do begin
    x := ldexp(1,n);
    y := erfce_inv(x);
    w := erfce(y);
    h := abs(1-w/x);
    if h>u then u := h;
  end;
  u := u/eps_d;
  if u > NE1 then begin
    inc(failed);
    writeln(' Fail: Max error erfce(erfce_inv(x)) = ', u:1:2, ' eps_d.');
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_erfi_inv;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfi_inv');

  {Maple: erfi_inv := x->fsolve(erfi(y)=x, y=sqrt(ln(x)));}
  x := 0;
  f := 0;
  y := erfi_inv(x);
  testrel(1, NE, y, f, cnt,failed);

  x := -1e-100;
  f := -0.8862269254527580136e-100;
  y := erfi_inv(x);
  testrel(2, NE, y, f, cnt,failed);

  x := sqrt_epsh;
  f := 0.2063407854765555853e-9;
  y := erfi_inv(x);
  testrel(3, NE, y, f, cnt,failed);

  x := 1e-6;
  f := 0.8862269254525260000e-6;
  y := erfi_inv(x);
  testrel(4, NE, y, f, cnt,failed);

  x := -0.5;
  f := -0.4175276255915004598;
  y := erfi_inv(x);
  testrel(5, NE, y, f, cnt,failed);

  x := 1.0;
  f := 0.7316971534684923963;
  y := erfi_inv(x);
  testrel(6, NE, y, f, cnt,failed);

  x := 1+1/1024;
  f := 0.7322036437239068925;
  y := erfi_inv(x);
  testrel(7, NE, y, f, cnt,failed);

  x := 1.25;
  f := 0.8499314171548616670;
  y := erfi_inv(x);
  testrel(8, NE, y, f, cnt,failed);

  x := -10;
  f := -1.8003080434321563762;
  y := erfi_inv(x);
  testrel(9, NE, y, f, cnt,failed);

  x := 100;
  f := 2.439855131304673795;
  y := erfi_inv(x);
  testrel(10, NE, y, f, cnt,failed);

  x := 1e10;
  f := 5.019014073055847358;
  y := erfi_inv(x);
  testrel(11, NE, y, f, cnt,failed);

  x := 1e300;
  f := 26.35562280355559192;
  y := erfi_inv(x);
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_erf_z;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','erf_z');
  {Maple: statevalf[pdf,normald](x)}

  y := erf_z(NegInf_x);
  f := 0;
  testabs( 1, 0, y, f, cnt,failed);

  y := erf_z(-35);
  f := 0.3940396277136024331e-266;
  testrel( 2, NE, y, f, cnt,failed);

  y := erf_z(-6.0);
  f := 0.60758828498232854869e-8;
  testrel( 3, NE, y, f, cnt,failed);

  y := erf_z(-2.0);
  f := 0.53990966513188051949e-1;
  testrel( 4, NE, y, f, cnt,failed);

  y := erf_z(-1.1);
  f := 0.21785217703255053138;
  testrel( 5, NE, y, f, cnt,failed);

  y := erf_z(-1.0);
  f := 0.24197072451914334980;
  testrel( 6, NE, y, f, cnt,failed);

  y := erf_z(-0.9);
  f := 0.26608524989875482182;
  testrel( 7, NE, y, f, cnt,failed);

  y := erf_z(-0.5);
  f := 0.35206532676429947777;
  testrel( 8, NE, y, f, cnt,failed);

  y := erf_z(-0.125);
  f := 0.39583768694474948414;
  testrel( 9, NE, y, f, cnt,failed);

  y := erf_z(-0.0078125);
  f := 0.39893010583499324766;
  testrel(10, NE, y, f, cnt,failed);

  y := erf_z(0.0);
  f := 0.39894228040143267794;
  testrel(11, NE, y, f, cnt,failed);

  y := erf_z(0.0078125);
  f := 0.39893010583499324766;
  testrel(12, NE, y, f, cnt,failed);

  y := erf_z(0.125);
  f := 0.39583768694474948414;
  testrel(13, NE, y, f, cnt,failed);

  y := erf_z(0.5);
  f := 0.35206532676429947777;
  testrel(14, NE, y, f, cnt,failed);

  y := erf_z(0.9);
  f := 0.26608524989875482182;
  testrel(15, NE, y, f, cnt,failed);

  y := erf_z(1.0);
  f := 0.24197072451914334980;
  testrel(16, NE, y, f, cnt,failed);

  y := erf_z(1.1);
  f := 0.21785217703255053138;
  testrel(17, NE, y, f, cnt,failed);

  y := erf_z(2.0);
  f := 0.53990966513188051949e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := erf_z(4.0);
  f := 0.5*0.26766045152977070355e-3; {0.133830225764885351774074488441e-3}
  testrel(19, NE, y, f, cnt,failed);

  y := erf_z(6.0);
  f := 0.60758828498232854869e-8;
  testrel(20, NE, y, f, cnt,failed);

  y := erf_z(10.0);
  f := 0.76945986267064193463e-22;
  testrel(21, NE, y, f, cnt,failed);

  y := erf_z(20.0);
  f := 0.5520948362159763190e-87;
  testrel(22, NE, y, f, cnt,failed);

  y := erf_z(30);
  f := 0.1473646134878547519e-195;
  testrel(23, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erf_p;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','erf_p');
  {Maple: statevalf[cdf,normald](x)}

  y := erf_p(NegInf_x);
  f := 0;
  testabs( 1, 0, y, f, cnt,failed);

  y := erf_p(-35);
  f := 0.1124910706472406243978835e-267;
  testrel( 2, NE, y, f, cnt,failed);

  y := erf_p(-30);
  f := 0.4906713927148187059531905e-197;
  testrel( 3, NE, y, f, cnt,failed);

  y := erf_p(-20);
  f := 0.27536241186062336951e-88;
  testrel( 4, NE, y, f, cnt,failed);

  y := erf_p(-10);
  f := 0.76198530241605260660e-23;
  testrel( 5, NE, y, f, cnt,failed);

  y := erf_p(-6.0);
  f := 0.98658764503769814070e-9;
  testrel( 6, NE, y, f, cnt,failed);

  y := erf_p(-2.0);
  f := 0.227501319481792072003e-1;
  testrel( 7, NE, y, f, cnt,failed);

  y := erf_p(-1.1);
  f := 0.13566606094638267518;
  testrel( 8, NE, y, f, cnt,failed);

  y := erf_p(-1.0);
  f := 0.15865525393145705142;
  testrel( 9, NE, y, f, cnt,failed);

  y := erf_p(-0.9);
  f := 0.18406012534675948856;
  testrel(10, NE, y, f, cnt,failed);

  y := erf_p(-0.5);
  f := 0.30853753872598689637;
  testrel(11, NE, y, f, cnt,failed);

  y := erf_p(-0.125);
  f := 0.45026177516988710702;
  testrel(12, NE, y, f, cnt,failed);

  y := erf_p(-0.0078125);
  f := 0.49688329513915741955;
  testrel(13, NE, y, f, cnt,failed);

  y := erf_p(0.0);
  f := 0.5;
  testrel(14, NE, y, f, cnt,failed);

  y := erf_p(0.0078125);
  f := 0.50311670486084258045;
  testrel(15, NE, y, f, cnt,failed);

  y := erf_p(0.125);
  f := 0.54973822483011289298;
  testrel(16, NE, y, f, cnt,failed);

  y := erf_p(0.5);
  f := 0.69146246127401310364;
  testrel(17, NE, y, f, cnt,failed);

  y := erf_p(0.9);
  f := 0.81593987465324051145;
  testrel(18, NE, y, f, cnt,failed);

  y := erf_p(1.0);
  f := 0.84134474606854294859;
  testrel(19, NE, y, f, cnt,failed);

  y := erf_p(1.1);
  f := 0.86433393905361732483;
  testrel(20, NE, y, f, cnt,failed);

  y := erf_p(2.0);
  f := 0.97724986805182079280;
  testrel(21, NE, y, f, cnt,failed);

  y := erf_p(4.0);
  f := 0.99996832875816688008;
  testrel(22, NE, y, f, cnt,failed);

  y := erf_p(6.0);
  f := 0.99999999901341235496;
  testrel(23, NE, y, f, cnt,failed);

  y := erf_p(10.0);
  f := 1;
  testrel(24, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erf_q;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;

  writeln('Function: ','erf_q');
  {Maple: statevalf[cdf,normald](-x)}
  y := erf_q(0);
  f := 0.5;
  testrel(1, NE, y, f, cnt,failed);

  y := erf_q(0.5);
  f := 0.3085375387259868964;
  testrel(2, NE, y, f, cnt,failed);

  y := erf_q(1);
  f := 0.1586552539314570514;
  testrel(3, NE, y, f, cnt,failed);

  y := erf_q(2);
  f := 0.2275013194817920720e-1;
  testrel(4, NE, y, f, cnt,failed);

  y := erf_q(5);
  f := 0.2866515718791939117e-6;
  testrel(5, NE, y, f, cnt,failed);

  y := erf_q(10);
  f := 0.7619853024160526066e-23;
  testrel(6, NE, y, f, cnt,failed);

  y := erf_q(20);
  f := 0.2753624118606233695e-88;
  testrel(7, NE, y, f, cnt,failed);

  y := erf_q(35);
  f := 0.1124910706472406244e-267;
  testrel(8, NE, y, f, cnt,failed);

  y := erf_q(-0.25);
  f := 0.5987063256829237242;
  testrel(9, NE, y, f, cnt,failed);

  y := erf_q(-1.5);
  f := 0.9331927987311419340;
  testrel(10, NE, y, f, cnt,failed);

  y := erf_q(-4);
  f := 0.9999683287581668801;
  testrel(11, NE, y, f, cnt,failed);

  y := erf_q(-8);
  f := 2*0.4999999999999996890; {0.999999999999999377903942572822;}
  testrel(12, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_fresnel;
var
  x,yc,ys,fs,fc: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Fresnel(CS)');

  x  := 0;
  fs := 0;
  fc := 0;
  Fresnel(x,ys,yc);
  testrel( 1, NE, ys, fs, cnt,failed);
  testrel( 2, NE, yc, fc, cnt,failed);

  x  := 1e-300;
  fs := 0.5235987755982988731e-900;
  fc := 1e-300;;
  Fresnel(x,ys,yc);
  testrel( 3, NE, ys, fs, cnt,failed);
  testrel( 4, NE, yc, fc, cnt,failed);

  x  := 1e-10;
  fs := 0.5235987755982988731e-30;
  fc := 1e-10;
  Fresnel(x,ys,yc);
  testrel( 5, NE, ys, fs, cnt,failed);
  testrel( 6, NE, yc, fc, cnt,failed);

  x  := 0.00048828125;
  fs := 0.6095491996946437605e-10;
  fc := 0.4882812499999931516e-3;
  Fresnel(x,ys,yc);
  testrel( 7, NE, ys, fs, cnt,failed);
  testrel( 8, NE, yc, fc, cnt,failed);

  x  := -0.125;
  fs := -0.1022609856621743054e-2;
  fc := -0.1249924702994110995;
  Fresnel(x,ys,yc);
  testrel( 9, NE, ys, fs, cnt,failed);
  testrel(10, NE, yc, fc, cnt,failed);

  x  := 0.5;
  fs := 0.6473243285999927761e-1;
  fc := 0.4923442258714463929;
  Fresnel(x,ys,yc);
  testrel(11, NE, ys, fs, cnt,failed);
  testrel(12, NE, yc, fc, cnt,failed);

  x  := 1.0;
  fs := 0.4382591473903547661;
  fc := 0.7798934003768228295;
  Fresnel(x,ys,yc);
  testrel(13, NE, ys, fs, cnt,failed);
  testrel(14, NE, yc, fc, cnt,failed);

  x  := -1.5;
  fs := -0.6975049600820930131;
  fc := -0.4452611760398215351;
  Fresnel(x,ys,yc);
  testrel(15, NE, ys, fs, cnt,failed);
  testrel(16, NE, yc, fc, cnt,failed);

  x  := 2.0;
  fs := 0.3434156783636982422;
  fc := 0.4882534060753407545;
  Fresnel(x,ys,yc);
  testrel(17, NE, ys, fs, cnt,failed);
  testrel(18, NE, yc, fc, cnt,failed);

  x  := 10.0;
  fs := 0.4681699785848822404;
  fc := 0.4998986942055157236;
  Fresnel(x,ys,yc);
  testrel(19, NE, ys, fs, cnt,failed);
  testrel(20, NE, yc, fc, cnt,failed);

  x  := -100.0;
  fs := -0.4968169011478375533;
  fc := -0.4999998986788178976;
  Fresnel(x,ys,yc);
  testrel(21, NE, ys, fs, cnt,failed);
  testrel(22, NE, yc, fc, cnt,failed);

  x  := 1000.0;
  fs := 0.4996816901138163061;
  fc := 0.4999999998986788164;
  Fresnel(x,ys,yc);
  testrel(23, NE, ys, fs, cnt,failed);
  testrel(24, NE, yc, fc, cnt,failed);

  x := 5000.10009765625;
  fs := 0.4999986583294711585;
  fc := 0.5000636465631322710;
  Fresnel(x,ys,yc);
  testrel(25, NE, ys, fs, cnt,failed);
  testrel(26, NE, yc, fc, cnt,failed);

  x := 6000;
  fs := 0.4999469483523027016;
  fc := 0.4999999999995309204;
  Fresnel(x,ys,yc);
  testrel(27, NE, ys, fs, cnt,failed);
  testrel(28, NE, yc, fc, cnt,failed);

  x := 1e5;
  fs := 0.4999968169011381621;
  fc := 0.4999999999999998987;
  Fresnel(x,ys,yc);
  testrel(29, NE, ys, fs, cnt,failed);
  testrel(20, NE, yc, fc, cnt,failed);

  x := 2e6;
  fs := 0.4999998408450569081;
  fc := 0.5;
  Fresnel(x,ys,yc);
  testrel(31, NE, ys, fs, cnt,failed);
  testrel(32, NE, yc, fc, cnt,failed);

  x  := 1e7;
  fs := 0.4999999681690113816;
  fc := 0.5;
  Fresnel(x,ys,yc);
  testrel(33, NE, ys, fs, cnt,failed);
  testrel(34, NE, yc, fc, cnt,failed);

  x  := 2e19;
  fs := 0.5;
  fc := 0.5;
  Fresnel(x,ys,yc);
  testrel(35, NE, ys, fs, cnt,failed);
  testrel(36, NE, yc, fc, cnt,failed);

  {one test for separate functions}
  x  := 0.75;
  fs := 0.2088771112333835702;
  fc := 0.6935259907871358975;
  ys := FresnelS(x);
  yc := FresnelC(x);
  testrel(37, NE, ys, fs, cnt,failed);
  testrel(38, NE, yc, fc, cnt,failed);

  x  := 30;
  fs := 0.4893896744421937968;
  fc := 0.4999962473706098869;
  Fresnel(x,ys,yc);
  testrel(39, NE, ys, fs, cnt,failed);
  testrel(40, NE, yc, fc, cnt,failed);

  x  := 50;
  fs := 0.4936338025859387415;
  fc := 0.4999991894307279680;
  Fresnel(x,ys,yc);
  testrel(41, NE, ys, fs, cnt,failed);
  testrel(42, NE, yc, fc, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_fresnelfg;
var
  x,f,g,fm,gm: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','Fresnel(FG)');

  x  := 0;
  fm := 0.5;
  gm := 0.5;
  FresnelFG(x,f,g);
  testrel( 1, NE, f, fm, cnt,failed);
  testrel( 2, NE, g, gm, cnt,failed);

  FresnelFG(1e-10, f,g);
  fm := 0.5;
  gm := 0.4999999999;
  testrel( 3, NE, f, fm, cnt,failed);
  testrel( 4, NE, g, gm, cnt,failed);

  FresnelFG(1e-8, f,g);
  fm := 0.4999999999999999215;
  gm := 0.4999999900000000785;
  testrel( 5, NE, f, fm, cnt,failed);
  testrel( 6, NE, g, gm, cnt,failed);

  FresnelFG(0.125, f,g);
  fm := 0.4896239619982552894;
  gm := 0.3871401026031354497;
  testrel( 7, NE, f, fm, cnt,failed);
  testrel( 8, NE, g, gm, cnt,failed);

  FresnelFG(1.0, f,g);
  fm := 0.2798934003768228295;
  gm := 0.06174085260964523392;
  testrel( 9, NE, f, fm, cnt,failed);
  testrel(10, NE, g, gm, cnt,failed);

  FresnelFG(1.5, f,g);
  fm := 0.2034184312260139559;
  gm := 0.2500979694279809422e-1;
  testrel(11, NE, f, fm, cnt,failed);
  testrel(12, NE, g, gm, cnt,failed);

  FresnelFG(1.5625, f,g);
  fm := 0.1962747584000039877;
  gm := 0.2257740933297765834e-1;
  testrel(13, NE, f, fm, cnt,failed);
  testrel(14, NE, g, gm, cnt,failed);

  FresnelFG(2, f,g);
  fm := 0.1565843216363017578;
  gm := 0.1174659392465924550e-1;
  testrel(15, NE, f, fm, cnt,failed);
  testrel(16, NE, g, gm, cnt,failed);

  FresnelFG(10, f,g);
  fm := 0.03183002141511775956;
  gm := 1.0130579448427638585e-4;
  testrel(17, NE, f, fm, cnt,failed);
  testrel(18, NE, g, gm, cnt,failed);

  FresnelFG(30, f,g);
  fm := 0.1061032555780620321e-1;
  gm := 0.3752629390113089832e-5;
  testrel(19, NE, f, fm, cnt,failed);
  testrel(20, NE, g, gm, cnt,failed);

  FresnelFG(50, f,g);
  fm := 0.6366197414061258547e-2;
  gm := 0.8105692720320441898e-6;
  testrel(21, NE, f, fm, cnt,failed);
  testrel(22, NE, g, gm, cnt,failed);

  x := 100;
  fm := 0.3183098852162446729e-2;
  gm := 0.1013211821024405315e-6;
  f  := FresnelF(x);
  g  := FresnelG(x);
  testrel(23, NE, f, fm, cnt,failed);
  testrel(24, NE, g, gm, cnt,failed);

  FresnelFG(10000, f,g);
  fm := 0.3183098861837906619e-4;
  gm := 0.1013211836423377560e-12;
  testrel(25, NE, f, fm, cnt,failed);
  testrel(26, NE, g, gm, cnt,failed);

  FresnelFG(65000, f,g);
  fm := 0.4897075172058318024e-5;
  gm := 0.3689437729352308473e-15;
  testrel(27, NE, f, fm, cnt,failed);
  testrel(28, NE, g, gm, cnt,failed);

  FresnelFG(-1, f,g);
  fm := -1.279893400376822829;
  gm :=  0.9382591473903547661;
  testrel(29, NE, f, fm, cnt,failed);
  testrel(30, NE, g, gm, cnt,failed);

  FresnelFG(-5, f,g);
  fm := -1.063631188704012231;
  gm :=  0.9991913819171168868;
  testrel(31, NE, f, fm, cnt,failed);
  testrel(32, NE, g, gm, cnt,failed);

  FresnelFG(-10, f,g);
  fm := 0.9681699785848822404;
  gm := 0.9998986942055157236;
  testrel(33, NE, f, fm, cnt,failed);
  testrel(34, NE, g, gm, cnt,failed);

  FresnelFG(-100, f,g);
  fm := 0.9968169011478375533;
  gm := 0.9999998986788178976;
  testrel(35, NE, f, fm, cnt,failed);
  testrel(36, NE, g, gm, cnt,failed);

  FresnelFG(-1000, f,g);
  fm := 0.9996816901138163061;
  gm := 0.9999999998986788164;
  testrel(37, NE, f, fm, cnt,failed);
  testrel(38, NE, g, gm, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_erfh2;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','erfh/2');

  y := erfh(1e-9, 1e-10);          {case 1}
  f := 0.2256758334191025146e-9;
  testrel( 1, NE, y, f, cnt,failed);

  y := erfh(1e-4, 1e-6);           {case 2} {double: case 1}
  f := 0.2256758311622689666e-5;
  testrel( 2, NE, y, f, cnt,failed);

  y := erfh(1e-3, 1e-6);           {maincase 1} {double: case 2}
  f := 0.2256756077433067085e-5;
  testrel( 3, NE, y, f, cnt,failed);

  y := erfh(1e-3, 5e-10);          {case 3}
  f := 0.1128378038716909668e-8;
  testrel( 4, NE, y, f, cnt,failed);

  y := erfh(0.5,2.25);             {case 4}
  f := 1.986571049297062807;
  testrel( 5, NE, y, f, cnt,failed);

  y := erfh(8, 9);                 {case 5}
  f := 1.842700792949714869;
  testrel( 6, NE, y, f, cnt,failed);

  y := erfh(0.25, 8.75);           {case 6}
  f := 2.0;
  testrel( 7, NE, y, f, cnt,failed);

  y := erfh(110, 1);               {case 7}
  f := 0;
  testabs( 8, 1, y, f, cnt,failed);

  y := erfh(4, 3);                 {case 8}
  f := 0.1572992070502851307;
  testrel( 9, NE, y, f, cnt,failed);

  y := erfh(0.2, 0.1);             {case 9}
  f := 0.2161638434408425354;
  testrel(10, NE, y, f, cnt,failed);

  y := erfh(5, 2);                 {case 10}
  f := 0.2209049699858544133e-4;
  testrel(11, NE, y, f, cnt,failed);

  y := erfh(2.0, 0.5);             {maincase 2: ax > 0.4}
  f := 0.3348790150724431399e-1;
  testrel(12, NE, y, f, cnt,failed);

  y := erfh(6, 1.5);               {maincase 2: ax > 0.4}
  f := 0.1966160441542887199e-9;
  testrel(13, NE, y, f, cnt,failed);

  y := erfh(10, 1);                {maincase 2: ax > 0.4, x-h >= 9} {double: case 8}
  f := 0.4137031746513810224e-36;
  testrel(14, NE, y, f, cnt,failed);

  y := erfh(9.5, 0.5);             {maincase 2: ax > 0.4, x-h >= 9}
  f := 0.4137031725628934400e-36;
  testrel(15, NE, y, f, cnt,failed);

  y := erfh(1/4, 1/16384);         {maincase 1: ax <= 0.4}
  f := 0.1293962558867630992e-3;
  testrel(16, NE, y, f, cnt,failed);

  y := erfh(0.2, 0.01);            {maincase 1: ax <= 0.4}
  f := 0.2168203082435382999e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := erfh(0.25, -0.1);           {case 9}
  f := -0.2113860821349468255;
  testrel(18, NE, y, f, cnt,failed);

  y := erfh(1.5, -1);              {case 10}
  f := -0.4790931701695085034;
  testrel(19, NE, y, f, cnt,failed);

  y := erfh(1, -1.5);              {case 4}
  f := -1.520092925795601579;
  testrel(20, NE, y, f, cnt,failed);

  y := erfh(6.5,0.5);              {maincase 2: ax > 0.4}
  f := 0.2151969487424283532e-16;
  testrel(21, NE, y, f, cnt,failed);

  y := erf2(4,5);
  f := 0.1541572044048559082e-7;
  testrel(22, NE, y, f, cnt,failed);

  y := erf2(6,5);
  f := -0.1537438274691322351e-11;
  testrel(23, NE, y, f, cnt,failed);

  y := erf2(-0.1, 0.2);
  f := 0.3351655052287633463;
  testrel(24, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_OwenT;
var
  a,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','OwenT');

  {Test values from [72]}
  y := OwenT(0.0625, 0.25);
  f := 3.8911930234701366897e-2;
  testrel(1, NE, y, f, cnt,failed);

  y := OwenT(6.5, 0.4375);
  f := 2.000577304850831541e-11;
  testrel(2, NE, y, f, cnt,failed);

  y := OwenT(7,0.96875);
  f := 6.399062719389868531e-13;
  testrel(3, NE, y, f, cnt,failed);

  y := OwenT(4.78125, 0.0625);
  f := 1.0632974804687463806e-7;
  testrel(4, NE, y, f, cnt,failed);

  y := OwenT(2, 0.5);
  f := 8.625077985521507131e-3;
  testrel(5, NE, y, f, cnt,failed);

  y := OwenT(1, 0.9999975);
  f := 6.674180897822859277e-2;
  testrel(6, NE, y, f, cnt,failed);

  {--------------------------------------}
  a := 1/8;
  y := OwenT(0,a);
  f := 0.01979171208028277101;
  testrel(7, NE, y, f, cnt,failed);

  y := OwenT(1/8,a);
  f := 0.01963689602674331987;
  testrel(8, NE, y, f, cnt,failed);

  y := OwenT(1/2,a);
  f := 0.01745480704924987218;
  testrel(9, NE, y, f, cnt,failed);

  y := OwenT(3/4,a);
  f := 0.01491780303257272623;
  testrel(10, NE, y, f, cnt,failed);

  y := OwenT(1,a);
  f := 0.01197322093905603733;
  testrel(11, NE, y, f, cnt,failed);

  y := OwenT(2,a);
  f := 0.002650989165102699101;
  testrel(12, NE, y, f, cnt,failed);

  y := OwenT(10,a);
  f := 3.021698334728800304e-24;
  testrel(13, NE, y, f, cnt,failed);

  y := OwenT(30,a);
  f := 2.452930207593037796e-198;
  testrel(14, NE, y, f, cnt,failed);

  {--------------------------------------}
  a := 1/4;
  y := OwenT(0,a);
  f := 0.03898956518868466273;
  testrel(15, NE, y, f, cnt,failed);

  y := OwenT(1/8,a);
  f := 0.03867995195352368921;
  testrel(16, NE, y, f, cnt,failed);

  y := OwenT(1/2,a);
  f := 0.03432021712709420963;
  testrel(17, NE, y, f, cnt,failed);

  y := OwenT(3/4,a);
  f := 0.02926208986224979475;
  testrel(18, NE, y, f, cnt,failed);

  y := OwenT(1,a);
  f := 0.02340824574252598197;
  testrel(19, NE, y, f, cnt,failed);

  y := OwenT(2,a);
  f := 0.005068174277162414139;
  testrel(20, NE, y, f, cnt,failed);

  y := OwenT(10,a);
  f := 3.765695914803479760e-24;
  testrel(21, NE, y, f, cnt,failed);

  y := OwenT(30,a);
  f := 2.453356963573946311e-198;
  testrel(22, NE, y, f, cnt,failed);

  {--------------------------------------}
  a := 3/4;
  y := OwenT(0,a);
  f := 0.1024163823495667258;
  testrel(23, NE, y, f, cnt,failed);

  y := OwenT(1/8,a);
  f := 0.1014881460301621857;
  testrel(24, NE, y, f, cnt,failed);

  y := OwenT(1/2,a);
  f := 0.08854921402151783342;
  testrel(25, NE, y, f, cnt,failed);

  y := OwenT(3/4,a);
  f := 0.07386590149693262928;
  testrel(26, NE, y, f, cnt,failed);

  y := OwenT(1,a);
  f := 0.05736472678167914879;
  testrel(27, NE, y, f, cnt,failed);

  y := OwenT(2,a);
  f := 0.01042929792412484401;
  testrel(29, NE, y, f, cnt,failed);

  y := OwenT(10,a);
  f := 3.809926512080107827e-24;
  testrel(29, NE, y, f, cnt,failed);

  y := OwenT(30,a);
  f := 2.453356963574093530e-198;
  testrel(30, NE, y, f, cnt,failed);

  {--------------------------------------}
  a := 1;
  y := OwenT(0,a);
  f := 0.125;
  testrel(31, NE, y, f, cnt,failed);

  y := OwenT(1/8,a);
  f := 0.1237630544953745706;
  testrel(32, NE, y, f, cnt,failed);

  y := OwenT(1/2,a);
  f := 0.1066710629614485163;
  testrel(33, NE, y, f, cnt,failed);

  y := OwenT(3/4,a);
  f := 0.08763369776575950553;
  testrel(34, NE, y, f, cnt,failed);

  y := OwenT(1,a);
  f := 0.06674188216570096662;
  testrel(35, NE, y, f, cnt,failed);

  y := OwenT(2,a);
  f := 0.01111628172225982148;
  testrel(36, NE, y, f, cnt,failed);

  y := OwenT(10,a);
  f := 3.809926512080263033e-24;
  testrel(37, NE, y, f, cnt,failed);

  y := OwenT(30,a);
  f := 2.453356963574093530e-198;
  testrel(38, NE, y, f, cnt,failed);

  {--------------------------------------}
  a := 2;
  y := OwenT(0,a);
  f := 0.1762081911747833629;
  testrel(39, NE, y, f, cnt,failed);

  y := OwenT(1/8,a);
  f := 0.1737438887583106110;
  testrel(40, NE, y, f, cnt,failed);

  y := OwenT(1/2,a);
  f := 0.1415806036539783935;
  testrel(41, NE, y, f, cnt,failed);

  y := OwenT(3/4,a);
  f := 0.1095704813486006094;
  testrel(42, NE, y, f, cnt,failed);

  y := OwenT(1,a);
  f := 0.07846818699308409634;
  testrel(43, NE, y, f, cnt,failed);

  y := OwenT(2,a);
  f := 0.01137490879318756594;
  testrel(44, NE, y, f, cnt,failed);

  y := OwenT(10,a);
  f := 3.809926512080263033e-24;
  testrel(45, NE, y, f, cnt,failed);

  y := OwenT(30,a);
  f := 2.453356963574093530e-198;
  testrel(46, NE, y, f, cnt,failed);

  {--------------------------------------}
  a := 10;
  y := OwenT(0,a);
  f := 0.2341372412847232152;
  testrel(47, NE, y, f, cnt,failed);

  y := OwenT(1/8,a);
  f := 0.2231420381435428845;
  testrel(48, NE, y, f, cnt,failed);

  y := OwenT(1/2,a);
  f := 0.1542687674982540145;
  testrel(49, NE, y, f, cnt,failed);

  y := OwenT(3/4,a);
  f := 0.1133136761884339769;
  testrel(50, NE, y, f, cnt,failed);

  y := OwenT(1,a);
  f := 0.0793276269657285257;
  testrel(51, NE, y, f, cnt,failed);

  y := OwenT(2,a);
  f := 0.0113750659740896036;
  testrel(52, NE, y, f, cnt,failed);

  y := OwenT(10,a);
  f := 3.809926512080263033e-24;
  testrel(53, NE, y, f, cnt,failed);

  y := OwenT(30,a);
  f := 2.453356963574093530e-198;
  testrel(54, NE, y, f, cnt,failed);

  {--------------------------------------}
  {test symmetry and special values}
  y := OwenT(2, 3);
  f := 0.01137506597154504061;
  testrel(55, NE, y, f, cnt,failed);

  y := OwenT(-2, 3);
  f := 0.01137506597154504061;
  testrel(56, NE, y, f, cnt,failed);

  y := OwenT(2, -3);
  f := -0.01137506597154504061;
  testrel(57, NE, y, f, cnt,failed);

  y := OwenT(-2, -3);
  f := -0.01137506597154504061;
  testrel(58, NE, y, f, cnt,failed);

  y := OwenT(0,-3);
  f := -0.1987918088252166371;
  testrel(59, NE, y, f, cnt,failed);

  y := OwenT(0,3);
  f := 0.19879180882521663709;
  testrel(60, NE, y, f, cnt,failed);

  y := OwenT(0,0);
  f := 0;
  testabs(61, NE, y, f, cnt,failed);

  y := OwenT(0,100);
  f := 0.2484085036175458724254;
  testrel(62, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;

end.
