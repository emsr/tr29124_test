{Part 8a of regression test for SPECFUN unit  (c) 2012  W.Ehrhardt}

unit t_sfd8a;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface


procedure test_primezeta;
procedure test_harmonic;
procedure test_harmonic2;
procedure test_zetah;
procedure test_cl2;
procedure test_dilog;
procedure test_trilog;
procedure test_ti2;
procedure test_lobachevsky_c;
procedure test_lobachevsky_s;
procedure test_euler;


implementation

uses
  amath, specfun, t_sfd0;

{---------------------------------------------------------------------------}
procedure test_zetah;
var
  y,f,a: double;
  cnt, failed: integer;
const
  NE = 1;
  NE2 = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','zetah');

  y := zetah(1.125,0.125);
  f := 0.1876387417613955387e2;
  testrel(1, NE, y, f, cnt,failed);

  y := zetah(1.5,0.125);
  f := 0.2501730935051015065e2;
  testrel(2, NE, y, f, cnt,failed);

  y := zetah(10,0.125);
  f := 0.1073741824308490817e10;
  testrel(3, NE, y, f, cnt,failed);

  y := zetah(64,0.125);
  f := 0.6277101735386680764e58;
  testrel(4, NE, y, f, cnt,failed);

  y := zetah(1.125,0.25);
  f := 0.1298083692721034518e2;
  testrel(5, NE, y, f, cnt,failed);

  y := zetah(10,0.25);
  f := 0.1048576107683114754e7;
  testrel(6, NE, y, f, cnt,failed);

  y := zetah(64,0.25);
  f := 0.3402823669209384635e39;
  testrel(7, NE, y, f, cnt,failed);

  y := zetah(10000,0.25);
  f := PosInf_x;
  testabs(8, 1, y, f, cnt,failed);

  y := zetah(1.125,0.5);
  f := 0.1014048375787649331e2;
  testrel(9, NE, y, f, cnt,failed);

  y := zetah(10,0.5);
  f := 0.1024017450355757901e4;
  testrel(10, NE, y, f, cnt,failed);

  y := zetah(60,0.5);
  f := 0.1152921504606846976e19;
  testrel(11, NE, y, f, cnt,failed);

  y := zetah(21165,0.5);
  f := PosInf_x;
  testabs(12, 1, y, f, cnt,failed);

  y := zetah(1.125, 1);
  f := 0.8586241294510575300e1;
  testrel(13, NE, y, f, cnt,failed);

  y := zetah(3, 1);
  f := 0.1202056903159594285e1;
  testrel(14, NE, y, f, cnt,failed);

  y := zetah(30, 1);
  f := 1.000000000931327432;
  testrel(15, NE, y, f, cnt,failed);

  y := zetah(1.125,1.25);
  f := 8.22400846719946091484921086761;
  testrel(16, NE, y, f, cnt,failed);

  y := zetah(10,1.25);
  f := 0.1076831147539267451;
  testrel(17, NE, y, f, cnt,failed);

  y := zetah(64,1.25);
  f := 0.6277101735386681052e-6;
  testrel(18, NE, y, f, cnt,failed);

  y := zetah(1.125,3.5);
  f := 6.969036597548218952;
  testrel(19, NE, y, f, cnt,failed);

  y := zetah(10,3.5);
  f := 0.3968242068686308214e-5;
  testrel(20, NE, y, f, cnt,failed);

  y := zetah(500,3.5);
  f := 0.9246509599237502930e-272;
  testrel(21, NE, y, f, cnt,failed);

  y := zetah(10000,3.5);
  f := 0;
  testrel(22, NE, y, f, cnt,failed);

  y := zetah(9,9.25);
  f := 0.3500985287686790903e-8;
  testrel(23, NE, y, f, cnt,failed);

  y := zetah(10,10.25);
  f := 0.1342531440842145817e-9;
  testrel(24, NE, y, f, cnt,failed);

  y := zetah(2,1/1024);
  f := 0.1048577642589392153e7;
  testrel(25, NE, y, f, cnt,failed);

  y := zetah(70,1-1/1024);
  f := 1.0707858041599152997;
  testrel(26, NE, y, f, cnt,failed);

  y := zetah(60,1-1/1024);
  f := 1.060374745217192150;
  testrel(27, NE, y, f, cnt,failed);

  y := zetah(60,0.875);
  f := 0.3016593692709887877e4;
  testrel(28, NE, y, f, cnt,failed);

  {extension to s < 1}
  y := zetah(0.3,0.5);
  f := -0.2090838188536853394;
  testrel(29, NE, y, f, cnt,failed);

  y := zetah(0.25,0.75);
  f := -0.5092007213592424391;
  testrel(30, NE, y, f, cnt,failed);

  y := zetah(0.125,0.25);
  f := 0.2941047482442775943;
  testrel(31, NE, y, f, cnt,failed);

  y := zetah(1-1/128,2);
  f := -128.4233535044541389;
  testrel(33, NE, y, f, cnt,failed);

  y := zetah(1-1/128,1e300);
  f := -28246.19608428274944;
  testrel(34, NE, y, f, cnt,failed);

  y := zetah(1-1/128,1e-300);
  f := 0.4531583637600817883e298;
  testrel(35, NE, y, f, cnt,failed);

  y := zetah(1e-6,0.125);
  f := 0.3750011004809786947;
  testrel(36, NE, y, f, cnt,failed);

  y := zetah(1e-6,1/128);
  f := 0.4921914286430100681;
  testrel(37, NE, y, f, cnt,failed);

  y := zetah(1e-6,2);
  f := -1.500000918939536384;
  testrel(38, NE, y, f, cnt,failed);

  y := zetah(0,1.25);
  f := -0.75;
  testrel(39, NE, y, f, cnt,failed);

  y := zetah(10.25, 1e-10);
  f := 0.3162277660168379332e103;
  testrel(40, NE2, y, f, cnt,failed);   {!!!}

  {Negative s}
  y := zetah(-10.25, 1e-10);
  f := 0.5259072290844016389e-2;
  testrel(41, NE, y, f, cnt,failed);

  y := zetah(-10.125,3);
  f := -0.1117677409454477505e4;
  testrel(42, NE, y, f, cnt,failed);

  y := zetah(-3,2.5);
  f := -3.507291666666666667;
  testrel(43, NE, y, f, cnt,failed);

  y := zetah(-7.875, 1234.5);
  f := -0.3070809320067946656e27;
  testrel(44, NE, y, f, cnt,failed);

  y := zetah(-8.0, 1234.5);
  f := -0.7372131481430697100e27;
  testrel(45, NE, y, f, cnt,failed);

  y := zetah(-8.125, 1234.5);
  f := -0.1770178566902331096e28;
  testrel(46, NE, y, f, cnt,failed);

  y := zetah(-10.125, 234.5);
  f := -0.2047336215461645541e26;
  testrel(47, NE, y, f, cnt,failed);

  y := zetah(-100.125, 234.5);
  f := -0.3790996246227456684e238;
  testrel(48, NE, y, f, cnt,failed);

  y := zetah(-500, 2);
  f := -1;
  testrel(49, NE, y, f, cnt,failed);

  y := zetah(-100, 2.5);
  f := -0.4065611775352152374e18;
  testrel(50, NE, y, f, cnt,failed);

  y := zetah(-100.25, 2.75);
  f :=-0.8348610836307083910e78;
  testrel(51, NE, y, f, cnt,failed);

  y := zetah(-200, 30);
  f := -0.3019871986270282730e293;
  testrel(52, NE, y, f, cnt,failed);

  y := zetah(-250.75, 16.25);
  f := -0.5070678782430602848e297;
  testrel(53, NE, y, f, cnt,failed);

  {a near 1, s < 0}
  a := 1+ldexp(1,-32);
  y := zetah(-10,a);
  f := -0.1763868512529315370e-10;
  testrel(54, NE, y, f, cnt,failed);

  a := 1-ldexp(1,-32);

  y := zetah(-2.5,a);
  f := 0.8516928792684670420e-2;
  testrel(55, NE, y, f, cnt,failed);

  y := zetah(-2,a);
  f := 0.3880510724853988384e-10;
  testrel(56, NE, y, f, cnt,failed);

  y := zetah(-9,a);
  f := -0.7575757575757575749e-2;
  testrel(57, NE, y, f, cnt,failed);

  y := zetah(-99,a);
  f := 0.2838224957069370693e77;
  testrel(58, NE, y, f, cnt,failed);

  y := zetah(-259,a);
  f := 0.8760156344622921506e307;
  testrel(59, NE, y, f, cnt,failed);

  a := 1-ldexp(1,-16);
  y := zetah(-9,a);
  f := -0.7575757540832979237e-2;
  testrel(60, NE, y, f, cnt,failed);

  y := zetah(-0.25,1234.5);
  f := -5851.0516641390874229;
  testrel(61, NE, y, f, cnt,failed);

  y := zetah(-1/1024, 20);
  f := -19.537565883074590864;
  testrel(62, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_primezeta;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','primezeta');

  {Values computed with mpmath 0.17, http://code.google.com/p/mpmath}
  {function primezeta(); parameters mp.prec=200, mpf.precision = 100}

  y := primezeta(1.125);
  f := 1.913164591315230872;
  testrel(1, NE, y, f, cnt,failed);

  y := primezeta(1.5);
  f := 0.8495626836215664464;
  testrel(2, NE, y, f, cnt,failed);

  y := primezeta(2);
  f := 0.4522474200410654985;
  testrel(3, NE, y, f, cnt,failed);

  y := primezeta(3);
  f := 0.1747626392994435364;
  testrel(4, NE, y, f, cnt,failed);

  y := primezeta(10);
  f := 0.9936035744369802179e-3;
  testrel(5, NE, y, f, cnt,failed);

  y := primezeta(36);
  f := 0.1455192189083022590e-10;
  testrel(6, NE, y, f, cnt,failed);

  y := primezeta(40);
  f := 0.9094947840255617476e-12;
  testrel(7, NE, y, f, cnt,failed);

  y := primezeta(50);
  f := 0.8881784210930808014e-15;
  testrel(8, NE, y, f, cnt,failed);

  y := primezeta(60);
  f := 0.8673617380119933721e-18;
  testrel(9, NE, y, f, cnt,failed);

  y := primezeta(120);
  f := 0.7523163845262640051e-36;
  testrel(10, NE, y, f, cnt,failed);

  y := primezeta(120);
  f := 0.7523163845262640051e-36;
  testrel(11, NE, y, f, cnt,failed);

  y := primezeta(1000);
  f := 0.9332636185032188794e-301;
  testrel(12, NE, y, f, cnt,failed);

  y := primezeta(1+1/1024);
  f := 6.617053486643312845;
  testrel(13, NE, y, f, cnt,failed);

  y := primezeta(1+ldexpd(1,-20));
  f := 13.54722642999334505;
  testrel(14, NE, y, f, cnt,failed);

  y := primezeta(1+ldexpd(1,-50));
  f := 34.34164057594337658;
  testrel(15, NE, y, f, cnt,failed);

  y := primezeta(1+ldexpd(1,-52));
  f := 35.72793493706326631;
  testrel(16, NE, y, f, cnt,failed);

  y := primezeta(27);
  f := 0.7450711734323300780e-8;
  testrel(17, NE, y, f, cnt,failed);

  {Arguments < 1}
  y := primezeta(0.203125);
  f := -1.310594771640394252;
  testrel(18, NE, y, f, cnt,failed);

  y := primezeta(0.328125);
  f := -1.968247421867130864;
  testrel(19, NE, y, f, cnt,failed);

  y := primezeta(0.34375);
  f := -1.757028476951072126;
  testrel(20, NE, y, f, cnt,failed);

  y := primezeta(0.375);
  f := -1.353421131778164875;
  testrel(21, NE, y, f, cnt,failed);

  y := primezeta(0.46875);
  f := -1.496570620236068512;
  testrel(22, NE, y, f, cnt,failed);

  y := primezeta(0.53125);
  f := -1.255414114971863083;
  testrel(23, NE, y, f, cnt,failed);

  y := primezeta(0.625);
  f := -0.2291137470699003535;
  testrel(24, NE, y, f, cnt,failed);

  y := primezeta(0.65673828125);
  f := 3.748098974927893859e-5;
  {only absolute error for P(x) near the zero}
  {for this x relative error is about 3.4e-16}
  testabs(25, 1, y, f, cnt,failed);

  y := primezeta(0.75);
  f := 0.6149705292500156942;
  testrel(26, NE, y, f, cnt,failed);

  y := primezeta(0.984375);
  f := 3.822024460290017266;
  testrel(27, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cl2;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NR = 1;
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','cl2');

  {The first 20 test are from MISCUN[22]}
  x := 1.0/512.0;  y := cl2(x);  f := +0.14137352886760576684e-1;
  testrel( 1, NR, y, f, cnt,failed);

  x :=  1.0/32.0;  y := cl2(x);  f := +0.13955467081981281934e0;
  testrel( 2, NR, y, f, cnt,failed);

  x :=  -1.0/8.0;  y := cl2(x);  f := -0.38495732156574238507e0;
  testrel( 3, NR, y, f, cnt,failed);

  x :=   1.0/2.0;  y := cl2(x);  f := +0.84831187770367927099e0;
  testrel( 4, NR, y, f, cnt,failed);

  x :=   1.0/1.0;  y := cl2(x);  f := +0.10139591323607685043e1;
  testrel( 5, NR, y, f, cnt,failed);

  x :=  -3.0/2.0;  y := cl2(x);  f := -0.93921859275409211003e0;
  testrel( 6, NR, y, f, cnt,failed);

  x :=   2.0/1.0;  y := cl2(x);  f := +0.72714605086327924743e0;
  testrel( 7, NR, y, f, cnt,failed);

  x :=   5.0/2.0;  y := cl2(x);  f := +0.43359820323553277936e0;
  testrel( 8, NR, y, f, cnt,failed);

  x :=  -3.0/1.0;  y := cl2(x);  f := -0.98026209391301421161e-1;
  testrel( 9, NR, y, f, cnt,failed);

  x :=   4.0/1.0;  y := cl2(x);  f := -0.56814394442986978080e0;
  testrel(10, NR, y, f, cnt,failed);

  x :=  17.0/4.0;  y := cl2(x);  f := -0.70969701784448921625e0;
  testrel(11, NR, y, f, cnt,failed);

  x :=  -5.0/1.0;  y := cl2(x);  f := +0.99282013254695671871e0;
  testrel(12, NR, y, f, cnt,failed);

  x :=  11.0/2.0;  y := cl2(x);  f := -0.98127747477447367875e0;
  testrel(13, NR, y, f, cnt,failed);

  x :=   6.0/1.0;  y := cl2(x);  f := -0.64078266570172320959e0;
  testrel(14, NR, y, f, cnt,failed);

  x :=   8.0/1.0;  y := cl2(x);  f := +0.86027963733231192456e0;
  testrel(15, NR, y, f, cnt,failed);

  x := -10.0/1.0;  y := cl2(x);  f := +0.39071647608680211043e0;
  testrel(16, NR, y, f, cnt,failed);

  x :=  15.0/1.0;  y := cl2(x);  f := +0.47574793926539191502e0;
  testrel(17, NR, y, f, cnt,failed);

  x :=  20.0/1.0;  y := cl2(x);  f := +0.10105014481412878253e1;
  testrel(18, NR, y, f, cnt,failed);

  x := -30.0/1.0;  y := cl2(x);  f := +0.96332089044363075154e0;
  testrel(19, NR, y, f, cnt,failed);

  x :=  50.0/1.0;  y := cl2(x);  f := -0.61782699481929311757e0;
  testrel(20, NR, y, f, cnt,failed);

  {cl2   := x -> -I*(polylog(2,exp(I*x))-polylog(2,exp(-I*x)))/2;}
  {cl2(x) = -I*(dilog(exp(I*x))-dilog(exp(-I*x)))/2}
  x := 1e-9;
  y := cl2(x);
  f := 0.2172326583694641116e-7;
  testrel(21, NE, y, f, cnt,failed);

  x := 199.0/64.0;
  y := cl2(x);
  f := 0.2233018233094811502e-1;
  testrel(22, NE, y, f, cnt,failed);

  x := 3198.0/1024.0;
  y := cl2(x);
  f := 0.1285468835749490331e-1;
  testrel(23, NE, y, f, cnt,failed);

  x := 3.125;
  y := cl2(x);
  f := 0.1150096070973137289e-1;
  testrel(24, NE, y, f, cnt,failed);

  x := 403.0/128.0;
  y := cl2(x);
  f := -0.4744472628326165504e-2;
  testrel(25, NE, y, f, cnt,failed);

  x := 3.140625;
  y := cl2(x);
  f := 0.6707263197711507115e-3;
  testrel(26, NE, y, f, cnt,failed);

  x := -3.140625;
  y := cl2(x);
  f := -0.6707263197711507115e-3;
  testrel(27, NE, y, f, cnt,failed);

  x := 3.14453125;
  y := cl2(x);
  f := -0.2036878759212949702e-2;
  testrel(28, NE, y, f, cnt,failed);

  x := 3.1416015625;
  y := cl2(x);
  f :=-0.6175185991649017015e-5;
  testrel(29, NE, y, f, cnt,failed);

  x := 6.28125;
  y := cl2(x);
  f := -0.1402611801429578260e-1;
  testrel(30, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := cl2(x);
  f :=  0.7745577948090751733e-2;
  testrel(32, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := cl2(x);
  f :=  -0.7745577948090751733e-2;
  testrel(32, NE, y, f, cnt,failed);

  x := 1e-8;
  y := cl2(x);
  f := 0.1942068074395236547e-6;
  testrel(33, NE, y, f, cnt,failed);

  x := 3.1416015625;
  y := cl2(x);
  f := -0.6175185991649017015e-5;
  testrel(34, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_dilog;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
  NA = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','dilog');

  x := 0.0;
  y := dilog(x);
  f := 0;
  testabs( 1, 1, y, f, cnt,failed);

  x := 1e-15;
  y := dilog(x);
  f := 0.10000000000000002500e-14;
  testrel( 2, NE, y, f, cnt,failed);

  x := -1e-15;
  y := dilog(x);
  f := -0.99999999999999975e-15;
  testrel( 3, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := dilog(x);
  f := 0.976801022116266601e-3;
  testrel( 4, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := dilog(x);
  f := -0.97632418484437659009e-3;
  testrel( 5, NE, y, f, cnt,failed);

  x := 0.125;
  y := dilog(x);
  f := 0.12913986010995340567;
  testrel( 6, NE, y, f, cnt,failed);

  x := -0.125;
  y := dilog(x);
  f := -0.12129662872272647871;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.25;
  y := dilog(x);
  f := 0.26765263908273260692;
  testrel( 8, NE, y, f, cnt,failed);

  x := -0.25;
  y := dilog(x);
  f := -0.23590029768626345382;
  testrel( 9, NE, y, f, cnt,failed);

  x := 0.5;
  y := dilog(x);
  f := 0.58224052646501250590;
  testrel(10, NE, y, f, cnt,failed);

  x := -0.5;
  y := dilog(x);
  f := -0.4484142069236462025;
  testrel(11, NE, y, f, cnt,failed);

  x := 0.75;
  y := dilog(x);
  f := 0.9784693929303061040;
  testrel(12, NE, y, f, cnt,failed);

  x := -0.75;
  y := dilog(x);
  f := -0.642761268839978879;
  testrel(13, NE, y, f, cnt,failed);

  x := 0.875;
  y := dilog(x);
  f := 1.2381234817964802347;
  testrel(14, NE, y, f, cnt,failed);

  x := 0.96875;
  y := dilog(x);
  f := 1.503403870941998162;
  testrel(15, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := dilog(x);
  f := 1.6371849430542471844;
  testrel(16, NE, y, f, cnt,failed);

  x := 1.0;
  y := dilog(x);
  f := 1.6449340668482264365;
  testrel(17, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := dilog(x);
  f := 1.6526761034351458339;
  testrel(18, NE, y, f, cnt,failed);

  x := 1.125;
  y := dilog(x);
  f := 2.0111536328199940118;
  testrel(19, NE, y, f, cnt,failed);

  x := 1.96875;
  y := dilog(x);
  f := 2.4671517720141088262;
  testrel(20, NE, y, f, cnt,failed);

  x := 1.9990234375;
  y := dilog(x);
  f := 2.4674008616984453258;
  testrel(21, NE, y, f, cnt,failed);

  x := 2.0;
  y := dilog(x);
  f := 2.4674011002723396547;
  testrel(22, NE, y, f, cnt,failed);

  x := 2.0009765625;
  y := dilog(x);
  f := 2.4674008620088863024;
  testrel(23, NE, y, f, cnt,failed);

  x := 2.125;
  y := dilog(x);
  f := 2.4637968182885816094;
  testrel(24, NE, y, f, cnt,failed);

  x := 2.5;
  y := dilog(x);
  f := 2.4207908065659338439;
  testrel(25, NE, y, f, cnt,failed);

  x := -2.5;
  y := dilog(x);
  f := -1.6988958419950141730;
  testrel(26, NE, y, f, cnt,failed);

  x := 5.0;
  y := dilog(x);
  f := 1.7837191612666306277;
  testrel(27, NE, y, f, cnt,failed);

  x := -5.0;
  y := dilog(x);
  f := -2.7492791260608082900;
  testrel(28, NE, y, f, cnt,failed);

  x := 10.0;
  y := dilog(x);
  f := 0.53630128735786273655;
  testrel(29, NE, y, f, cnt,failed);

  x := -10.0;
  y := dilog(x);
  f := -4.1982778868581038579;
  testrel(30, NE, y, f, cnt,failed);

  x := 12.59517036984501613;
  y := dilog(x);
  f := -0.264629547172758555217e-18;
  testabs(31, NA, y, f, cnt,failed);

  x := 12.595170369845;
  y := dilog(x);
  f := 0.3138080782974712711461e-14;
  testabs(32, NA, y, f, cnt,failed);

  x := 12.59517;
  y := dilog(x);
  f := 0.71959170440116884065945e-7;
  testabs(33, NA, y, f, cnt,failed);

  x := 12.59518;
  y := dilog(x);
  f := -0.1873697847668815973913e-5;
  testabs(34, NA, y, f, cnt,failed);

  x := 12.59517;
  y := dilog(x);
  f := 0.71959170440116884065945e-7;
  testabs(35, NA, y, f, cnt,failed);

  x := 100.0;
  y := dilog(x);
  f := -7.3239531990004822427;
  testrel(36, NE, y, f, cnt,failed);

  x := -100.0;
  y := dilog(x);
  f := -12.238755177314938922;
  testrel(37, NE, y, f, cnt,failed);

  x := 1e10;
  y := dilog(x);
  f := -261.80503739032344766;
  testrel(38, NE, y, f, cnt,failed);

  x := -1e10;
  y := dilog(x);
  f := -266.73983959066812696;
  testrel(39, NE, y, f, cnt,failed);

  x := MaxDouble;
  y := dilog(x);
  f := -251892.45989301229208;
  testrel(40, NE, y, f, cnt,failed);

  x := -MaxDouble;
  y := dilog(x);
  f := -251897.39469521283676;
  testrel(41, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_trilog;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','trilog');

  y := trilog(1);
  f := 1.202056903159594285;
  testrel(1, NE, y, f, cnt,failed);

  y := trilog(1.5);
  f := 2.060877507320280871;
  testrel(2, NE, y, f, cnt,failed);

  y := trilog(2);
  f := 2.762071906228924136;
  testrel(3, NE, y, f, cnt,failed);

  y := trilog(5);
  f := 4.805344102965590473;
  testrel(4, NE, y, f, cnt,failed);

  y := trilog(10);
  f := 5.641811414751341251;
  testrel(5, NE, y, f, cnt,failed);

  {absolute error near positive zero x=85.1716733428841653527728949292}
  x := 85.171875;
  y := trilog(x);
  f := -0.1562520906191944201e-4;
  testabs(5, NE, y, f, cnt,failed);

  y := trilog(128);
  f := -3.067549375920255031;
  testrel(7, NE, y, f, cnt,failed);

  y := trilog(4294967296.0);
  f := -1745.787022447498125;
  testrel(8, NE, y, f, cnt,failed);

  y := trilog(1e10);
  f := -1958.926579067720252;
  testrel(9, NE, y, f, cnt,failed);

  y := trilog(1e19);
  f := -13811.93163941220694;
  testrel(10, NE, y, f, cnt,failed);

  y := trilog(MaxDouble);
  f := -0.5959474778573073133e8;
  testrel(11, NE, y, f, cnt,failed);

  y := trilog(1+1/1024);
  f := 1.203666516900211776;
  testrel(12, NE, y, f, cnt,failed);

  y := trilog(1+1/1024);
  f := 1.203666516900211776;
  testrel(13, NE, y, f, cnt,failed);

  x := 1+ldexp(1,-20);
  y := trilog(x);
  f := 1.202058471897204243;
  testrel(14, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-20);
  y := trilog(x);
  f := 1.202055334434460785;
  testrel(15, NE, y, f, cnt,failed);

  x := -123;
  y := trilog(x);
  f := -26.49656819677117693;
  testrel(16, NE, y, f, cnt,failed);

  y := trilog(-4294967296.0);
  f := -1855.244437869677651;
  testrel(17, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_ti2;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','ti2');

  {Calculated with Pari/GP and the definition}
  {ti2(x) = 1/2*I*(dilog(-I*x)-dilog(I*x))}

  {Maple ti2 := x -> 1/2*I*(polylog(2,-I*x)-polylog(2,I*x))}
  {works for |x|<1000, for |x|>=1000 Release 4.00b hangs. This}
  {is the same for ti2 := x ->int(arctan(t)/t, t=0..x);}

  x := 0.0;
  y := ti2(x);
  f := 0;
  testrel( 1, 1, y, f, cnt,failed);

  x := 4e-10;
  y := ti2(x);
  f := 3.9999999999999999999e-10;
  testrel( 2, NE, y, f, cnt,failed);

  x := 6e-10;
  y := ti2(x);
  f := 5.9999999999999999998e-10;
  testrel( 3, NE, y, f, cnt,failed);

  x := -3e-5;
  y := ti2(x);
  f := -0.2999999999700000000e-4;
  testrel( 4, NE, y, f, cnt,failed);

  x := 0.5;
  y := ti2(x);
  f := 0.48722235829452235711;
  testrel( 5, NE, y, f, cnt,failed);

  x := -1.0;
  y := ti2(x);
  f := -0.9159655941772190151;
  testrel( 6, NE, y, f, cnt,failed);

  x := 2.0;
  y := ti2(x);
  f := 1.576015403446323422;
  testrel( 7, NE, y, f, cnt,failed);

  x := -10.0;
  y := ti2(x);
  f := -3.716781493068068590;
  testrel( 8, NE, y, f, cnt,failed);

  x := 100.0;
  y := ti2(x);
  f := 7.243784301308353497;
  testrel( 9, NE, y, f, cnt,failed);

  x := 1000;
  y := ti2(x);
  f := 10.85167661851208615;
  testrel(10, NE, y, f, cnt,failed);

  x := 1e5;
  y := ti2(x);
  f := 18.08447103103866192;
  testrel(11, NE, y, f, cnt,failed);

  x := 3.0e9;
  y := ti2(x);
  f := 34.27772600381452632;
  testrel(12, NE, y, f, cnt,failed);

  x := 3.04e9;
  y := ti2(x);
  f := 34.29853155733663828;
  testrel(13, NE, y, f, cnt,failed);

  x := 3.1e9;
  y := ti2(x);
  f := 34.32923213705038316;
  testrel(14, NE, y, f, cnt,failed);

  x := 1e20;
  y := ti2(x);
  f := 72.33784412415464813;
  testrel(15, NE, y, f, cnt,failed);

  x := 1e308;
  y := ti2(x);
  f := 1114.002799511981581;
  testrel(16, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_harmonic;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
  NA  = 1;
  NA2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','harmonic');

  y := harmonic(0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := harmonic(1);
  f := 1;
  testrel(2, NE, y, f, cnt,failed);

  y := harmonic(2);
  f := 1.5;
  testrel(3, NE, y, f, cnt,failed);

  y := harmonic(12);
  f := 3.103210678210678211;
  testrel(4, NE, y, f, cnt,failed);

  y := harmonic(13);
  f := 3.180133755133755134;
  testrel(5, NE, y, f, cnt,failed);

  y := harmonic(0.125);
  f := 0.1887230016056779928;
  testrel(6, NE, y, f, cnt,failed);

  y := harmonic(0.25);
  f := 0.3497621315252674525;
  testrel(7, NE, y, f, cnt,failed);

  y := harmonic(1/3);
  f := 0.4451818848807265376;
  testrel(8, NE, y, f, cnt,failed);

  y := harmonic(0.4);
  f := 0.5158311203164167149;
  testrel(9, NE, y, f, cnt,failed);

  y := harmonic(0.5);
  f := 0.6137056388801093812;
  testrel(10, NE, y, f, cnt,failed);

  y := harmonic(-0.125);
  f := -0.2268014066461625217;
  testrel(11, NE, y, f, cnt,failed);

  y := harmonic(-0.25);
  f := -0.5086452148849393090;
  testrel(12, NE, y, f, cnt,failed);

  y := harmonic(-1/3);
  f := -0.7410187508850556118;
  testrel(13, NE, y, f, cnt,failed);

  y := harmonic(-2/7);
  f := -0.6029657754379658942;
  testrel(14, NE, y, f, cnt,failed);

  y := harmonic(-0.5);
  f := -1.386294361119890619;
  testrel(15, NE, y, f, cnt,failed);

  y := harmonic(-0.75);
  f := -3.650237868474732547;
  testrel(16, NE, y, f, cnt,failed);

  y := harmonic(100.0);
  f := 5.187377517639620261;
  testrel(17, NE, y, f, cnt,failed);

  y := harmonic(1000.0);
  f := 7.485470860550344913;
  testrel(18, NE, y, f, cnt,failed);

  y := harmonic(1234567890.0);
  f := 21.51120252446859658;
  testrel(19, NE, y, f, cnt,failed);

  y := harmonic(1e20);
  f := 46.62891752478244654;
  testrel(20, NE, y, f, cnt,failed);

  {cot(Pi*x)=0}
  y := harmonic(-1234.5);
  f := 7.695231896729530704;
  testrel(21, NE, y, f, cnt,failed);

  {Extreme: large negative x with cot(Pi*x)<>0}
  y := harmonic(-12345.675);
  f := 8.073065788498956869;
  testrele(22, 1.5E-12, y, f, cnt,failed);  {!!!!!!!!}

  {absolute near zero}
  y := harmonic(-2.625);
  f := 0.3860123114100866946e-1;
  testabs(23, NA, y, f, cnt,failed);

  y := harmonic(-4.6875);
  f := -0.8747240023647986249e-1;
  testabs(24, NA, y, f, cnt,failed);

  y := harmonic(-4.68118577225);
  f := 0.2464763120598810321e-10;
  testabs(25, NA2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_harmonic2;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','harmonic2');

  {integer x}
  y := harmonic2(4, -10);
  f := 1108650.0;
  testrel(1, NE, y, f, cnt,failed);

  y := harmonic2(5, -11);
  f := 53201625.0;
  testrel(2, NE, y, f, cnt,failed);

  y := harmonic2(5, -11.25);
  f := 79182293.86892715208;
  testrel(3, NE, y, f, cnt,failed);

  y := harmonic2(5, 11.25);
  f := 1.000415065530641588;
  testrel(4, NE, y, f, cnt,failed);

  y := harmonic2(6, -2.5);
  f := 198.3286416953014298;
  testrel(5, NE, y, f, cnt,failed);

  y := harmonic2(20, 30);
  f := 1.000000000931327432;
  testrel(6, NE, y, f, cnt,failed);

  y := harmonic2(10, 8.5);
  f := 1.002859248026791897;
  testrel(7, NE, y, f, cnt,failed);

  {x<1}
  y := harmonic2(0.5, 10);
  {f:= 0.9835442193699167854;}
  f := 1.967088438739833571*0.5;
  testrel(8, NE, y, f, cnt,failed);

  y := harmonic2(0.5, 1);
  f := 0.6137056388801093812;
  testrel(9, NE, y, f, cnt,failed);

  y := harmonic2(0.5, 0.9375);
  f := 0.6070681530436531038;
  testrel(10, NE, y, f, cnt,failed);

  y := harmonic2(0.96875, 0.9375);
  f := 0.9790742810492755181;
  testrel(11, NE, y, f, cnt,failed);

  y := harmonic2(0.53125, 0.9375);
  f := 0.6360771119883223991;
  testrel(12, NE, y, f, cnt,failed);

  y := harmonic2(1e-9, 2.5);
  f := 0.2816834663678296263e-8;
  testrel(13, NE, y, f, cnt,failed);

  y := harmonic2(0.125, 2.5);
  f := 0.29122756828655908281;
  testrel(14, NE, y, f, cnt,failed);

  y := harmonic2(0.125, -0.5);
  {f:= 0.9622063502817188448e-1;}
  f := 0.1924412700563437690*0.5;
  testrel(15, NE, y, f, cnt,failed);

  y := harmonic2(0.125, -2);
  f := 0.29296875e-1;
  testrel(16, NE, y, f, cnt,failed);

  y := harmonic2(0.125, -1);
  f := 0.703125e-1;
  testrel(17, NE, y, f, cnt,failed);

  y := harmonic2(0.5, 2.5);
  f := 0.7512308721744858897;
  testrel(18, NE, y, f, cnt,failed);

  y := harmonic2(0.125, 1e-19);
  f := 0.125;
  testrel(19, NE, y, f, cnt,failed);

  y := harmonic2(0.125, 1e-15);
  f := 0.1250000000000000600;
  testrel(20, NE, y, f, cnt,failed);

  y := harmonic2(0.125, 1e-9);
  f := 0.1250000000600231841;
  testrel(21, NE, y, f, cnt,failed);

  y := harmonic2(0.125, -1 + ldexpd(1,-30));
  f := 0.7031250004521736822e-1;
  testrel(22, NE, y, f, cnt,failed);

  y := harmonic2(-0.125, 2.5);
  f := -0.4399995003575235926;
  testrel(23, NE, y, f, cnt,failed);

  y := harmonic2(-0.5, 0.75);
  f := -1.095041682395986023;
  testrel(24, NE, y, f, cnt,failed);

  y := harmonic2(-0.75, 0.75);
  f := -2.503693422758070408;
  testrel(25, NE, y, f, cnt,failed);

  y := harmonic2(-0.5, 1e-9);
  f := -0.5000000005723649432;
  testrel(26, NE, y, f, cnt,failed);

  y := harmonic2(-0.75, 1e-9);
  f := -0.7500000012880225257;
  testrel(27, NE, y, f, cnt,failed);

  y := harmonic2(-0.75,1);
  f := -3.650237868474732547;
  testrel(28, NE, y, f, cnt,failed);

  y := harmonic2(-0.75, -3.125);
  f := 0.6556786052139600106e-2;
  testrel(29, NE, y, f, cnt,failed);

  y := harmonic2(0.75, -3.125);
  f := 0.4163030422854981335;
  testrel(30, NE, y, f, cnt,failed);

  y := harmonic2(0.96875, -3.125);
  f := 0.9060665748486035285;
  testrel(31, NE, y, f, cnt,failed);

  y := harmonic2(-0.5, -3.125);
  f := 0.1425234865536558298e-1;
  testrel(32, NE, y, f, cnt,failed);

  y := harmonic2(0.5, -3.125);
  f := 0.1288778540559494870;
  testrel(33, NE, y, f, cnt,failed);

  y := harmonic2(-0.5, 3.5);
  f := -10.49407079612483783;
  testrel(34, NE, y, f, cnt,failed);

  y := harmonic2(0.5, 3.5);
  f := 0.8196377028599225592;
  testrel(35, NE, y, f, cnt,failed);

  y := harmonic2(-0.625, 3.5);
  f := -30.24352029038063375;
  testrel(36, NE, y, f, cnt,failed);

  y := harmonic2(0.625, 3.5);
  f := 0.8874329085593194128;
  testrel(37, NE, y, f, cnt,failed);

  y := harmonic2(-0.96875, 3.5);
  f := -185363.6921486060383;
  testrel(38, NE, y, f, cnt,failed);

  y := harmonic2(0.96875, 3.5);
  f := 0.9938171087571505753;
  testrel(39, NE, y, f, cnt,failed);

  {x >= 1}
  y := harmonic2(2.5, 1.75);
  f := 1.380842090808948150;
  testrel(40, NE, y, f, cnt,failed);

  y := harmonic2(2.5, -0.5);
  f := 3.244215792104376753;
  testrel(41, NE, y, f, cnt,failed);

  y := harmonic2(2.5, 1e-9);
  f := 2.499999998799026398;
  testrel(42, NE, y, f, cnt,failed);

  y := harmonic2(1, -1 + 1e-9);
  f := 1.0;
  testrel(43, NE, y, f, cnt,failed);

  y := harmonic2(2, -1 + 1e-9);
  f := 2.999999998613705639;
  testrel(44, NE, y, f, cnt,failed);

  y := harmonic2(1000,-2);
  f := 333833500.0;
  testrel(45, NE, y, f, cnt,failed);

  y := harmonic2(1000,-3);
  f := 250500250000.0;
  testrel(46, NE, y, f, cnt,failed);

  y := harmonic2(1000,-3.125);
  f := 0.5760647784978448461e12;
  testrel(47, NE, y, f, cnt,failed);

  y := harmonic2(2.5, 3.5);
  f := 1.102039769777552195;
  testrel(48, NE, y, f, cnt,failed);

  y := harmonic2(2.5, -0.75);
  f := 3.750362305074930384;
  testrel(49, NE, y, f, cnt,failed);

  y := harmonic2(1.25, -0.75);
  f := 1.360786713750249277;
  testrel(50, NE, y, f, cnt,failed);

  y := harmonic2(1+1/1024, -0.75);
  f := 1.001328406180784655;
  testrel(51, NE, y, f, cnt,failed);

  y := harmonic2(2.5, -3.5);
  f := 28.93568599767752340;
  testrel(51, NE, y, f, cnt,failed);

  y := harmonic2(1.25, -3.5);
  f := 2.191909416239121815;
  testrel(53, NE, y, f, cnt,failed);

  y := harmonic2(1+1/1024, -3.5);
  f := 1.003393139272674024;
  testrel(54, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_lobachevsky_c;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lobachevsky_c');

  x := ldexp(1,63); {1/eps_x}
  f := 0.6393154322601327830e19;
  y := lobachevsky_c(x);
  testrel(1, NE, y, f, cnt,failed);

  x := ldexp(1,52); {1/eps_d}
  f := 0.31216573840826801117e16;
  y := lobachevsky_c(x);
  testrel(2, NE, y, f, cnt,failed);

  x := -0.125;
  f := -0.3260309790120720032e-3;
  y := lobachevsky_c(x);
  testrel(3, NE, y, f, cnt,failed);

  x := 1e-5;
  f := 0.1666666666683333333e-15;
  y := lobachevsky_c(x);
  testrel(4, NE, y, f, cnt,failed);

  x := 1e-10;
  f := 0.1666666666666666667e-30;
  y := lobachevsky_c(x);
  testrel(5, NE, y, f, cnt,failed);

  x := Pi_4;
  f := 0.8641372548729102510e-1;
  y := lobachevsky_c(x);
  testrel(6, NE, y, f, cnt,failed);

  x := Pi_2;
  {computed with x = asd(Pi/2) = 1.570796326794896557998981734272...}
  f := 1.088793045151798718;
  y := lobachevsky_c(x);
  testrel(7, NE, y, f, cnt,failed);

  x := -5;
  f := -3.913719502878449559;
  y := lobachevsky_c(x);
  testrel(8, NE, y, f, cnt,failed);

  x := -25;
  f := -17.42029821132434980;
  y := lobachevsky_c(x);
  testrel(9, NE, y, f, cnt,failed);

  x := 1e10;
  f := 0.6931471805929815760e10;
  y := lobachevsky_c(x);
  testrel(10, NE, y, f, cnt,failed);

  x := 1e20;
  f := 0.6931471805599453094e20;
  y := lobachevsky_c(x);
  testrel(11, NE, y, f, cnt,failed);

  {Miscfun test cases}
  x := 0.001953125;
  f := 0.1241763906516139386e-8;
  y := lobachevsky_c(x);
  testrel(12, NE, y, f, cnt,failed);

  x := 0.0078125;
  f := 0.7947334477000108823e-07;
  y := lobachevsky_c(x);
  testrel(13, NE, y, f, cnt,failed);

  x := 0.03125;
  f := 0.5086759818620883420e-05;
  y := lobachevsky_c(x);
  testrel(14, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.3260309790120720032e-03;
  y := lobachevsky_c(x);
  testrel(16, NE, y, f, cnt,failed);

  x := 0.5;
  f := 0.2138053681540821442e-01;
  y := lobachevsky_c(x);
  testrel(16, NE, y, f, cnt,failed);

  x := 1.0;
  f := 0.1875381690208382405;
  y := lobachevsky_c(x);
  testrel(17, NE, y, f, cnt,failed);

  x := 1.5;
  f := 0.83051199971883645115;
  y := lobachevsky_c(x);
  testrel(18, NE, y, f, cnt,failed);

  x := 2.0;
  f := 0.1885436242667903490e+01;
  y := lobachevsky_c(x);
  testrel(19, NE, y, f, cnt,failed);

  x := 2.5;
  f := 0.2131598898651641105e+01;
  y := lobachevsky_c(x);
  testrel(20, NE, y, f, cnt,failed);

  x := 3.0;
  f := 0.2177112018561342722e+01;
  y := lobachevsky_c(x);
  testrel(21, NE, y, f, cnt,failed);

  x := 4.0;
  f := 0.2292102792189665085e+01;
  y := lobachevsky_c(x);
  testrel(22, NE, y, f, cnt,failed);

  x := 5.0;
  f := 0.3913719502878449559e+01;
  y := lobachevsky_c(x);
  testrel(23, NE, y, f, cnt,failed);

  x := 6.0;
  f := 0.4351356398383642790e+01;
  y := lobachevsky_c(x);
  testrel(24, NE, y, f, cnt,failed);

  x := 7.0;
  f := 0.4420064496847818590e+01;
  y := lobachevsky_c(x);
  testrel(25, NE, y, f, cnt,failed);

  x := 10.0;
  f := 0.6565601313362382916e+01;
  y := lobachevsky_c(x);
  testrel(26, NE, y, f, cnt,failed);

  x := 15.0;
  f := 0.1082550466150459948e+02;
  y := lobachevsky_c(x);
  testrel(27, NE, y, f, cnt,failed);

  x := 20.0;
  {f:= 0.13365512855474227325e+02;}
  f := 0.2673102571094845465e+02*0.5;
  y := lobachevsky_c(x);
  testrel(28, NE, y, f, cnt,failed);

  x := 30.0;
  f := 0.2113100268563995993e+02;
  y := lobachevsky_c(x);
  testrel(29, NE, y, f, cnt,failed);

  x := 50.0;
  f := 0.3483823658944911739e+02;
  y := lobachevsky_c(x);
  testrel(30, NE, y, f, cnt,failed);

  x := 100.0;
  f := 0.6965706243783739428e+02;
  y := lobachevsky_c(x);
  testrel(31, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lobachevsky_s;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','lobachevsky_s');

  x := 0;
  f := 0;
  y := lobachevsky_s(x);
  testabs(1, 0, y, f, cnt,failed);

  x := Pi;
  {f computed with actual float value for Pi}
  f := -0.4524526404044344276e-14;
  {extended f := 0.2244699373958346916e-17;}
  y := lobachevsky_s(x);
  testabs(2, NE, y, f, cnt,failed);

  x := Pi/3;
  f := 0.3383138688032178750;
  y := lobachevsky_s(x);
  testrel(3, NE, y, f, cnt,failed);

  x := Pi/4;
  f := 0.4579827970886095075;
  y := lobachevsky_s(x);
  testrel(4, NE, y, f, cnt,failed);

  x := Pi/6;
  f := 0.5074708032048268125; {Maximum value}
  y := lobachevsky_s(x);
  testrel(5, NE, y, f, cnt,failed);

  x := -0.125;
  {f:= -0.2983953360169009094}
  f := -0.596790672033801819*0.5;
  y := lobachevsky_s(x);
  testrel(6, NE, y, f, cnt,failed);

  x := 1;
  f := 0.3635730254316396237;
  y := lobachevsky_s(x);
  testrel(7, NE, y, f, cnt,failed);

  x := -2;
  f := 0.2840719722149348904;
  y := lobachevsky_s(x);
  testrel(8, NE, y, f, cnt,failed);

  x := 3.125;
  f := -0.7310164582093342783e-1;
  y := lobachevsky_s(x);
  testrel(9, NE, y, f, cnt,failed);

  x := -1e10;
  f := 0.5072900264888338291;
  y := lobachevsky_s(x);
  testrel(10, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_euler;
var
  y,f: double;
  n,k,cnt,failed: integer;
const
  NE = 1;
const
  nmax = 180;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','euler');
  for k:=0 to nmax div 2 do begin
    n := 2*k;
    y := euler(n);
    f := zetah(-n,0.25);
    f := ldexpd(f,n+n+2);
    testrel(k, NE, y, f, cnt,failed);
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;

end.
