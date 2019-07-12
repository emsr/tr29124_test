{Part 2b of regression test for SPECFUND unit  (c) 2013  W.Ehrhardt}

unit t_sfd2b;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface


procedure test_jacobi_arcsn;
procedure test_jacobi_arccn;
procedure test_jacobi_arcdn;

procedure test_jacobi_arcsc;
procedure test_jacobi_arccs;
procedure test_jacobi_arcnc;
procedure test_jacobi_arcns;
procedure test_jacobi_arcnd;
procedure test_jacobi_arccd;
procedure test_jacobi_arcdc;
procedure test_jacobi_arcsd;
procedure test_jacobi_arcds;

implementation

uses
  damath, specfund, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_jacobi_arcsn;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE = 4;
  NE2 = 12;    {D64, otherwise NE2=8}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcsn');

  k := 0.75;
  x := 1e-8;
  f := 1.000000000000000026e-8;
  y := jacobi_arcsn(x,k);
  testrel(1, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.1255125343116088719;
  y := jacobi_arcsn(x,k);
  testrel(2, NE, y, f, cnt,failed);

  x := 0.25;
  f := 0.2541977440434425896;
  y := jacobi_arcsn(x,k);
  testrel(3, NE, y, f, cnt,failed);

  x := 0.5;
  f := 0.5372335944267905037;
  y := jacobi_arcsn(x,k);
  testrel(4, NE, y, f, cnt,failed);

  x := 0.75;
  f := 0.9067855002649757614;
  y := jacobi_arcsn(x,k);
  testrel(5, NE, y, f, cnt,failed);

  x := 0.875;
  f := 1.182037116036673001;
  y := jacobi_arcsn(x,k);
  testrel(6, NE, y, f, cnt,failed);

  x := 0.9375;
  f := 1.386909274425256848;
  y := jacobi_arcsn(x,k);
  testrel(7, NE, y, f, cnt,failed);

  x := 1.0;
  f := 1.910989780751829197;
  y := jacobi_arcsn(x,k);
  testrel(8, NE, y, f, cnt,failed);

  k := 1.5;
  x := 1e-8;
  f := 1.000000000000000054e-8;
  y := jacobi_arcsn(x,k);
  testrel(9, NE, y, f, cnt,failed);

  x := 0.125;
  f := 0.1260756176604139572;
  y := jacobi_arcsn(x,k);
  testrel(10, NE, y, f, cnt,failed);

  x := 0.25;
  f := 0.2590679342828028516;
  y := jacobi_arcsn(x,k);
  testrel(11, NE, y, f, cnt,failed);

  x := 0.5;
  f := 0.5950912525404356888;
  y := jacobi_arcsn(x,k);
  testrel(12, NE, y, f, cnt,failed);

  x := 0.625;
  f := 0.8935505613066573745;
  y := jacobi_arcsn(x,k);
  testrel(13, NE, y, f, cnt,failed);

  x := 0.666015625;
  f := 1.166923593833933571;
  y := jacobi_arcsn(x,k);
{$ifdef CPUARM}
  testrel(14, 20, y, f, cnt,failed);
{$else}
  testrel(14, NE, y, f, cnt,failed);
{$endif}

  k := 4.0;
  x := 1/4;
  f := 0.3990605555329458775;
  y := jacobi_arcsn(x,k);
  testrel(15, NE, y, f, cnt,failed);

  x := 0.75;
  y := jacobi_arcsn(x,0);
  f := arcsin(x);
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arcsn(x,1);
  f := arctanh(x);
  testrel(17, NE, y, f, cnt,failed);

  {x/k near 0 or 1}
  x := 1-ldexpd(1,-20);
  k := 1-ldexpd(1,-14);
  f := 5.767236895237212015;
  y := jacobi_arcsn(x,k);
{$ifdef CPUARM}
  testrel(18, 20, y, f, cnt,failed);
{$else}
  testrel(18, NE2, y, f, cnt,failed);
{$endif}

  x := 1-ldexpd(1,-20);
  k := 1+ldexpd(1,-20);
  f := 7.970212451787251854;
  y := jacobi_arcsn(x,k);
  testrel(19, NE, y, f, cnt,failed);

  x := 1-ldexpd(1,-30);
  k := ldexpd(1,-20);
  f := 1.570753168422375252;
  y := jacobi_arcsn(x,k);
  testrel(20, NE, y, f, cnt,failed);

  x := ldexpd(1,-30);
  k := ldexpd(0.875,30);
  f := 0.1134016488823230046e-8;
  y := jacobi_arcsn(x,k);
  testrel(21, NE, y, f, cnt,failed);

  {negative x}
  y := jacobi_arcsn(-0.75,0.5);
  f := -0.8716588079216072607;
  testrel(22, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-0.75,1.25);
  f := -1.138143885639268067;
  testrel(23, NE, y, f, cnt,failed);


  k := 0.5;
  y := jacobi_arcsn(0, k);
  f := 0.0;
  testrel(24, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-1/1024, k);
  f := -9.765626940256182572e-4;
  testrel(25, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.0625, k);
  f := 0.6255095074767487577e-1;
  testrel(26, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-0.125, k);
  f := -0.1254097402642616810;
  testrel(27, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.25, k);
  f := 0.2533486577214804679;
  testrel(28, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-0.5, k);
  f := -0.5294286270519058177;
  testrel(29, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.75, k);
  f := 0.8716588079216072607;
  testrel(30, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-0.875, k);
  f := -1.109808675860371512;
  testrel(31, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.9921875, k);
  f := 1.541443506748182545;
  testrel(32, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-1, k);
  f := -1.685750354812596043;
  testrel(33, NE, y, f, cnt,failed);

  k := 2;
  y := jacobi_arcsn(0, k);
  f := 0.0;
  testrel(34, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-1/1024, k);
  f := -9.765632761034555790e-4;
  testrel(35, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.0625, k);
  f := 0.6270487013213084050e-1;
  testrel(36, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-0.125, k);
  f := -0.1266743288607402339;
  testrel(37, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.25, k);
  f := 0.2647143135259529090;
  testrel(38, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-0.375, k);
  f := -0.4358294039608036304;
  testrel(39, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.4375, k);
  f := 0.5549043379301857559;
  testrel(40, NE, y, f, cnt,failed);

  y := jacobi_arcsn(-0.4921875, k);
  f := -0.7408556987000292771;
  testrel(41, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.5, k);
  f := 0.8428751774062980214;
  testrel(42, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.25, 4);
  f := 0.3990605555329458775;
  testrel(43, NE, y, f, cnt,failed);

  y := jacobi_arcsn(0.19921875, 5);
  f := 0.2993264734350559155;
  testrel(44, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arccn;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 16;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arccn');

  k := 0.75;
  x := 0;
  f := 1.910989780751829197;
  y := jacobi_arccn(x,k);
  testrel(1, NE, y, f, cnt,failed);

  x := 1e-8;
  f := 1.910989765633250276;
  y := jacobi_arccn(x,k);
  testrel(2, NE, y, f, cnt,failed);

  x := 0.125;
  f := 1.722141969722491233;
  y := jacobi_arccn(x,k);
  testrel(3, NE, y, f, cnt,failed);

  x := 0.25;
  f := 1.533954617380815167;
  y := jacobi_arccn(x,k);
  testrel(4, NE, y, f, cnt,failed);

  x := 0.5;
  f := 1.157958094010720378;
  y := jacobi_arccn(x,k);
  testrel(5, NE, y, f, cnt,failed);

  x := 0.75;
  f := 0.7589125065434269319;
  y := jacobi_arccn(x,k);
  testrel(6, NE, y, f, cnt,failed);

  x := 0.875;
  f := 0.5176098467406406717;
  y := jacobi_arccn(x,k);
  testrel(7, NE, y, f, cnt,failed);

  x := 0.9375;
  f := 0.3596574867018126771;
  y := jacobi_arccn(x,k);
  testrel(8, NE, y, f, cnt,failed);

  x := 1.0;
  f := 0;
  y := jacobi_arccn(x,k);
  testrel(9, NE, y, f, cnt,failed);

  k := 1.5;
  x := 0.75;
  f := 1.0945810939666114305;
  y := jacobi_arccn(x,k);
  testrel(10, NE, y, f, cnt,failed);

  x := 0.875;
  f := 0.5680522359771010260;
  y := jacobi_arccn(x,k);
  testrel(11, NE, y, f, cnt,failed);

  x := 0.9375;
  f := 0.3742727651159121512;
  y := jacobi_arccn(x,k);
  testrel(12, NE, y, f, cnt,failed);

  x := 255/256;
  f := 0.8867800933123157563e-1;
  y := jacobi_arccn(x,k);
  testrel(13, NE, y, f, cnt,failed);

  {x or k near 1 or 0}
  x := 1-ldexpd(1,-20);
  f := 0.1381069029580945514e-2;
  y := jacobi_arccn(x,1.5);
  testrel(14, NE, y, f, cnt,failed);

  x := 1;
  f := 0;
  y := jacobi_arccn(x,k);
  testrel(15, NE, y, f, cnt,failed);

  x := 0.75;
  y := jacobi_arccn(x,0);
  f := arccos(x);
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arccn(x,1);
  f := arcsech(x);
  testrel(17, NE, y, f, cnt,failed);

  k := 0.875;
  x := ldexpd(1,-30);
  f := 2.185488467354492049;
  y := jacobi_arccn(x,k);
  testrel(18, NE, y, f, cnt,failed);

  k := 1-ldexpd(1,-14);
  x := 0;
  f := 5.891915583987215421;
  y := jacobi_arccn(x,k);
  testrel(19, NE, y, f, cnt,failed);

  x := ldexpd(1,-15);
  f := 5.889153409486920880;
  y := jacobi_arccn(x,k);
  testrel(20, NE2, y, f, cnt,failed);

  x := -ldexpd(1,-15);
  f := 5.894677758487509961;
  y := jacobi_arccn(x,k);
  testrel(21, NE2, y, f, cnt,failed);

  {negative x}
  k := 0.75;
  x := -1.0;
  f := 3.821979561503658393;
  y := jacobi_arccn(x,k);
  testrel(22, NE, y, f, cnt,failed);

  k := 0.75;
  x := -0.25;
  f := 2.288024944122843227;
  y := jacobi_arccn(x,k);
  testrel(23, NE, y, f, cnt,failed);

  k := 0.75;
  x := ldexpd(-1,-30);
  f := 1.910989782159856581;
  y := jacobi_arccn(x,k);
  testrel(24, NE, y, f, cnt,failed);

  k := ldexpd(1,-30);
  x := ldexpd(-1,-30);
  f := 1.570796327726219194;
  y := jacobi_arccn(x,k);
  testrel(25, NE, y, f, cnt,failed);

  k := ldexpd(-1,-30);
  x := -0.5;
  f := 2.094395102393195493;
  y := jacobi_arccn(x,k);
  testrel(26, NE, y, f, cnt,failed);

  k := 0.5;
  y := jacobi_arccn(1, k);
  f := 0.0;
  testrel(27, NE, y, f, cnt,failed);

  y := jacobi_arccn(1-1/1024, k);
  f := 0.0442013679399762574185108812419;
  testrel(28, NE, y, f, cnt,failed);

  y := jacobi_arccn(0.75, k);
  f := 0.737672368302886292896633710301;
  testrel(29, NE, y, f, cnt,failed);

  y := jacobi_arccn(0.5, k);
  f := 1.089550670051885409;
  testrel(30, NE, y, f, cnt,failed);

  y := jacobi_arccn(0.25, k);
  f := 1.394992625369530698;
  testrel(31, NE, y, f, cnt,failed);

  y := jacobi_arccn(0.03125, k);
  f := 1.649662045285716059;
  testrel(32, NE, y, f, cnt,failed);

  y := jacobi_arccn(0, k);
  f := 1.685750354812596043;
  testrel(33, NE, y, f, cnt,failed);

  y := jacobi_arccn(-1/1024, k);
  f := 1.686877992176595841;
  testrel(34, NE, y, f, cnt,failed);

  y := jacobi_arccn(-0.25, k);
  f := 1.976508084255661388;
  testrel(35, NE, y, f, cnt,failed);

  y := jacobi_arccn(-0.75, k);
  f := 2.633828341322305793;
  testrel(36, NE, y, f, cnt,failed);

  y := jacobi_arccn(-1+1/1024, k);
  f := 3.327299341685215828;
  testrel(37, NE, y, f, cnt,failed);

  y := jacobi_arccn(-1, k);
  f := 3.371500709625192086;
  testrel(38, NE, y, f, cnt,failed);

  k := 2;
  {x ~ sqrt(3)/2..1}
  y := jacobi_arccn(0.8660888671875, k);
  f := 0.8307680659816238384;
  testrel(39, NE2, y, f, cnt,failed);

  y := jacobi_arccn(0.875, k);
  f := 0.6974963126847653491;
  testrel(40, NE, y, f, cnt,failed);

  y := jacobi_arccn(0.90625, k);
  f := 0.5227981471475147268;
  testrel(41, NE, y, f, cnt,failed);

  y := jacobi_arccn(0.96875, k);
  f := 0.2623795549460648467;
  testrel(42, NE, y, f, cnt,failed);

  y := jacobi_arccn(0.9921875, k);
  f := 0.1264199600557561270;
  testrel(43, NE, y, f, cnt,failed);

  y := jacobi_arccn(1, k);
  f := 0.0;
  testrel(44, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arcdn;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 32;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcdn');

  {-----------------------------------------}
  k := 1.5;
  x := ldexpd(1,-16);
  f := 1.206431349115219658;
  y := jacobi_arcdn(x,k);
  testrel(1, NE, y, f, cnt,failed);

  x := 0.125;
  f := 1.094581093966611431;
  y := jacobi_arcdn(x,k);
  testrel(2, NE, y, f, cnt,failed);

  x := 0.25;
  f := 0.982298793673304697;
  y := jacobi_arcdn(x,k);
  testrel(3, NE, y, f, cnt,failed);

  x := 0.5;
  f := 0.7530028596237113767;
  y := jacobi_arcdn(x,k);
  testrel(4, NE, y, f, cnt,failed);

  x := 0.75;
  f := 0.5003371719928295923;
  y := jacobi_arcdn(x,k);
  testrel(5, NE, y, f, cnt,failed);

  x := 0.875;
  f := 0.3432704279622044052;
  y := jacobi_arcdn(x,k);
  testrel(6, NE, y, f, cnt,failed);

  x := 0.9375;
  f := 0.2391638269349492710;
  y := jacobi_arcdn(x,k);
  testrel(7, NE, y, f, cnt,failed);

  x := 1-ldexpd(1,-20);
  f := 0.9207121579245082433e-3;
  y := jacobi_arcdn(x,k);
  testrel(8, NE, y, f, cnt,failed);

  x := 1;
  f := 0;
  y := jacobi_arcdn(x,k);
  testrel(9, NE, y, f, cnt,failed);

  {-----------------------------------------}
  k := 0.75;
  x := 0.75;
  f := 1.201282778685676692;
  y := jacobi_arcdn(x,k);
  testrel(10, NE, y, f, cnt,failed);

  x := 0.875;
  f := 0.7347510784697389252;
  y := jacobi_arcdn(x,k);
  testrel(11, NE, y, f, cnt,failed);

  x := 0.9375;
  f := 0.4931327726691525676;
  y := jacobi_arcdn(x,k);
  testrel(12, NE, y, f, cnt,failed);

  x := 255/256;
  f := 0.1181638862247553136;
  y := jacobi_arcdn(x,k);
  testrel(13, NE, y, f, cnt,failed);

  {x near 1 or 0, k near 1}
  x := 1-ldexpd(1,-20);
  f := 0.1841425096347588164e-2;
  y := jacobi_arcdn(x,k);
  testrel(14, NE, y, f, cnt,failed);

  x := 1;
  f := 0;
  y := jacobi_arcdn(x,k);
  testrel(15, NE, y, f, cnt,failed);

  x := 0.625;
  k := 1;
  y := jacobi_arcdn(x,k);
  f := 1.046967915003188411;
  testrel(16, NE, y, f, cnt,failed);

  k := 1+ldexpd(1,-16);
  x := 0;
  f := 6.5848517916584914796;
  y := jacobi_arcdn(x,k);
  testrel(17, NE, y, f, cnt,failed);

  k := 1+ldexpd(1,-16);
  x := ldexpd(1,-20);
  f := 6.584679158826389041;
  y := jacobi_arcdn(x,k);
  testrel(18, NE, y, f, cnt,failed);

  k := 1-ldexpd(1,-12);
  x := ldexpd(3,-6);
  f := 3.814124582731960262;
  y := jacobi_arcdn(x,k);
  testrel(19, NE2, y, f, cnt,failed);

  k := 0.5;
  {x = kc .. 1 = sqrt(3)/2 .. }
  y := jacobi_arcdn(0.8662109375, k);
  f := 1.644343100326904530;
  testrel(20, NE+1, y, f, cnt,failed);    {+1 for ARM}

  y := jacobi_arcdn(0.875,0.5);
  f := 1.394992625369530698;
  testrel(21, NE, y, f, cnt,failed);

  y := jacobi_arcdn(0.90625, k);
  f := 1.045596294295029454;
  testrel(22, NE, y, f, cnt,failed);

  y := jacobi_arcdn(0.96875, k);
  f := 0.5247591098921296932;
  testrel(23, NE, y, f, cnt,failed);

  y := jacobi_arcdn(1, k);
  f := 0.0;
  testrel(24, NE, y, f, cnt,failed);

  k := 2;
  y := jacobi_arcdn(-0.9990234375, k);
  f := 1.663649670842607914;
  testrel(25, NE, y, f, cnt,failed);

  y := jacobi_arcdn(-0.875, k);
  f := 1.430444414968407738;
  testrel(26, NE, y, f, cnt,failed);

  y := jacobi_arcdn(-0.25, k);
  f := 2*0.4941270210639153469;
  testrel(27, NE, y, f, cnt,failed);

  y := jacobi_arcdn(0, k);
  f := 0.8428751774062980214;
  testrel(28, NE, y, f, cnt,failed);

  y := jacobi_arcdn(0.5, k);
  f := 0.5447753350259427046;
  testrel(29, NE, y, f, cnt,failed);

  y := jacobi_arcdn(0.875, k);
  f := 0.2553059398441883044;
  testrel(30, NE, y, f, cnt,failed);

  y := jacobi_arcdn(0.96875, k);
  f := 0.1256540933223839525;
  testrel(31, NE, y, f, cnt,failed);

  y := jacobi_arcdn(1, k);
  f := 0;
  testrel(32, NE, y, f, cnt,failed);

  y := jacobi_arcdn(-0.875, 5);
  f := 0.5335108266509599867;
  testrel(33, NE, y, f, cnt,failed);

  y := jacobi_arcdn(-0.75, 3);
  f := 0.8351988606817798565;
  testrel(34, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arcsc;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
{$ifdef FPC271or3}   {64-bit}
  NE2 = 40;
{$else}
  NE2 = 32;
{$endif}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcsc');

  {--------}
  k := 0.5;
  y := jacobi_arcsc(0,k);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arcsc(-0.0625,k);
  f := -0.6242893947277133946e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arcsc(0.125,k);
  f := 0.1244350127068800583;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arcsc(-0.5,k);
  f := -0.4677190513505932145;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arcsc(1,k);
  f := 0.8043661012320655565;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arcsc(-2,k);
  f := -1.156321727760690225;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_arcsc(8,k);
  f := 1.542280274951022840;
  testrel(7, NE, y, f, cnt,failed);

  y := jacobi_arcsc(-1e6,k);
  f := -1.685749200112057664;
  testrel(8, NE, y, f, cnt,failed);

  {--------}
  k := 2;
  y := jacobi_arcsc(0,k);
  f := 0.0;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_arcsc(-0.0625,k);
  f := -0.6258195616920016351e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_arcsc(0.125,k);
  f := 0.1256698453399509173;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arcsc(-0.25,k);
  f := -0.2558646551628484008;
  testrel(12, NE, y, f, cnt,failed);

  y := jacobi_arcsc(0.5,k);
  f := 0.5781608638803451126;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arcsc(-0.55,k);
  f := -0.6877487722016756098;
  testrel(14, NE, y, f, cnt,failed);

  x := 0.57733154296875; {~1.0/sqrt(3)};
  y := jacobi_arcsc(x,k);
  f := 0.8388480744361355186;
  testrel(15, NE2, y, f, cnt,failed);

  {--------}
  k := 1-1/1024;
  x := 100;
  y := jacobi_arcsc(x,k);
  f := 4.282976939350126374;
  testrel(16, NE, y, f, cnt,failed);

  k := 1+1/1024;
  x := 1/1024;
  y := jacobi_arcsc(x,k);
  f := 0.97656234508295035466e-3;
  testrel(17, NE, y, f, cnt,failed);

  k := 1+1/1024;
  x := 16;
  y := jacobi_arcsc(x,k);
  f := 3.623592182529732367;
{$ifdef CPUARM}
  testrel(18, 40, y, f, cnt,failed);
{$else}
  testrel(18, NE2, y, f, cnt,failed);
{$endif}

  k := 1/1024;
  x := 100;
  y := jacobi_arcsc(x,k);
  f := 1.560797029847403445;
  testrel(19, NE, y, f, cnt,failed);

  k := 1/1024;
  x := 0.125;
  y := jacobi_arcsc(x,k);
  f := 0.1243549948514774209;
  testrel(20, NE, y, f, cnt,failed);

  y := jacobi_arcsc(0.0999755859375,10);
  f := 0.1472175842452071406;
  testrel(21, NE, y, f, cnt,failed);

  y := jacobi_arcsc(1/128,100);
  f := 0.8966378520558796833e-2;
  testrel(22, NE, y, f, cnt,failed);

  y := jacobi_arcsc(1e10,1);
  f := 23.71899811050040215;
  testrel(23, NE, y, f, cnt,failed);

  y := jacobi_arcsc(42,0);
  f := 1.546991300609826676;
  testrel(24, NE, y, f, cnt,failed);

  x := 1/3;
  y := jacobi_arcsc(x,3);
  f := 0.4256118745355927045;
  testrel(26, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arccs;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
{$ifdef FPC271or3}   {64-bit}
  NE2 = 40;
{$else}
  NE2 = 16;
{$endif}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arccs');

  {--------}
  k := 0.5;
  y := jacobi_arccs(0, k);
  f := 1.685750354812596043;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arccs(1/1024, k);
  f := 1.684622717986295370;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arccs(0.5, k);
  f := 1.156321727760690225;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arccs(1, k);
  f := 0.8043661012320655565;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arccs(2, k);
  f := 0.4677190513505932145;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arccs(8, k);
  f := 0.1244350127068800583;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_arccs(1024, k);
  f := 2*0.4882811141821932303e-3;
  testrel(7, NE, y, f, cnt,failed);

  y := jacobi_arccs(1e6, k);
  f := 9.999999999997083333e-7;
  testrel(8, NE, y, f, cnt,failed);

  y := jacobi_arccs(1e10, k);
  f := 1e-10;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_arccs(-5, k);
  f := -0.1977149385392947046;
  testrel(10, NE, y, f, cnt,failed);

  {--------}
  k := 2;
  x := 1.732421875; {~ sqrt(3) = boundary}
  y := jacobi_arccs(x, k);
  f := 0.8325261366952091590;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arccs(1.75, k);
  f := 0.7711401374755114200;
{$ifdef CPUARM}
  testrel(12, 8, y, f, cnt,failed);
{$else}
  testrel(12, NE, y, f, cnt,failed);
{$endif}

  y := jacobi_arccs(5, k);
  f := 0.2028726523469851893;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arccs(128, k);
  f := 0.7812658963183475901e-2;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_arccs(1e6, k);
  f := 2*0.5000000000001666667e-6;
  testrel(15, NE, y, f, cnt,failed);

  y := jacobi_arccs(-2, k);
  f := -0.5781608638803451126;
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arccs(-10, k);
  f := -0.1003394354919195643;
  testrel(17, NE, y, f, cnt,failed);

  {--------}
  k := 1-1/1024;
  x := 1/1024;
  y := jacobi_arccs(x,k);
  f := 4.485312916914429815;
{$ifdef CPUARM}
  testrel(18, 8, y, f, cnt,failed);
{$else}
  testrel(18, NE, y, f, cnt,failed);
{$endif}

  y := jacobi_arccs(1024,k);
  f := 2*0.4882811722383103635e-3;
  testrel(19, NE, y, f, cnt,failed);

  k := 1+1/1024;
  x := 0.0625;
  y := jacobi_arccs(x,k);
  f := 3.623592182529732367;
{$ifdef CPUARM}
  testrel(20, 40, y, f, cnt,failed);
{$else}
  testrel(20, NE2, y, f, cnt,failed);
{$endif}

  y := jacobi_arccs(1024,k);
  f := 2*0.4882811725414751773e-3;
  testrel(21, NE, y, f, cnt,failed);

  y := jacobi_arccs(0.125,0);
  f := 1.446441332248135184;
  testrel(22, NE, y, f, cnt,failed);

  y := jacobi_arccs(0.125,1/1024);
  f := 1.446441647762956175;
  testrel(23, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arcnc;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
{$ifdef FPC271or3}   {64-bit}
  NE2 = 48;
{$else}
  NE2 = 16;
{$endif}
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcnc');

  k := 0.5;
  y := jacobi_arcnc(1, k);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arcnc(1.015625, k);
  f := 0.1758624151346494556;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arcnc(1.125, k);
  f := 0.4802796703198898856;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arcnc(2, k);
  f := 1.089550670051885409;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arcnc(8, k);
  f := 1.541159831597026902;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arcnc(16, k);
  f := 1.613550174333460224;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_arcnc(1e6, k);
  f := 1.685749200112057663;
  testrel(7, NE, y, f, cnt,failed);

  k := 2;
  y := jacobi_arcnc(1, k);
  f := 0.0;
  testrel(8, NE, y, f, cnt,failed);

  y := jacobi_arcnc(1+1/1024, k);
  f := 0.04423385703235096325;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_arcnc(1.015625, k);
  f := 0.1794404858826064887;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_arcnc(1.0625, k);
  f := 0.3791057562248457994;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arcnc(1.125, k);
  f := 0.6070171936345210231;
  testrel(12, NE, y, f, cnt,failed);

  x := 1.1546630859375; {~sqrt(4/3)};
  y := jacobi_arcnc(x, k);
  f := 0.8348206448834532702;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arcnc(2, 1);
  f := 1.316957896924816709;
  testrel(14, NE, y, f, cnt,failed);

  k := 1-1/1024;
  y := jacobi_arcnc(2, k);
  f := 1.315912552507747042;
  testrel(15, NE, y, f, cnt,failed);

{$ifdef CPUARM}
  y := jacobi_arcnc(1e6, k);
  f := 4.507390964955750622;
  testrel(16, 6, y, f, cnt,failed);

  k := 1+1/1024;
  y := jacobi_arcnc(2, k);
  f := 1.318009394793904406;
  testrel(17, NE, y, f, cnt,failed);

  y := jacobi_arcnc(16, k);
  f := 3.620832060474116051;
  testrel(18, 48, y, f, cnt,failed);
{$else}
  y := jacobi_arcnc(1e6, k);
  f := 4.507390964955750622;
  testrel(16, NE, y, f, cnt,failed);

  k := 1+1/1024;
  y := jacobi_arcnc(2, k);
  f := 1.318009394793904406;
  testrel(17, NE, y, f, cnt,failed);

  y := jacobi_arcnc(16, k);
  f := 3.620832060474116051;
  testrel(18, NE2, y, f, cnt,failed);
{$endif}


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_jacobi_arcns;
var
  k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcns');

  k := 0.5;
  y := jacobi_arcns(1, k);
  f := 1.685750354812596043;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arcns(-1.0625, k);
  f := -1.290246106568424218;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arcns(1.125, k);
  f := 1.142653490597676858;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arcns(-2, k);
  f := -0.5294286270519058177;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arcns(4, k);
  f := 0.2533486577214804679;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arcns(-8, k);
  f := -0.1254097402642616810;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_arcns(1e6, k);
  f := 1.000000000000208333e-6;
  testrel(7, NE, y, f, cnt,failed);

  k := 2;
  y := jacobi_arcns(2, k);
  f := 0.8428751774062980214;
  testrel(8, NE, y, f, cnt,failed);

  y := jacobi_arcns(-2.5, k);
  f := -0.4788095499820080662;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_arcns(3, k);
  f := 0.3725431893564766129;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_arcns(-4, k);
  f := -0.2647143135259529089;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arcns(8, k);
  f := 0.1266743288607402339;
  testrel(12, NE, y, f, cnt,failed);

  y := jacobi_arcns(-64, k);
  f := -0.1562818028887092862e-1;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arcns(1e6, k);
  f := 1.000000000000833333e-6;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_arcns(3+1/1024, 3);
  f := 0.5301091595130078372;
  testrel(15, NE, y, f, cnt,failed);

  k := 1-1/1024;
  y := jacobi_arcns(1, k);
  f := 4.507413597899042267;
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arcns(64, k);
  f := 0.1562627051053165931e-1;
  testrel(17, NE, y, f, cnt,failed);

  k := 1+1/1024;
  y := jacobi_arcns(1+1/256, k);
  f := 3.188030956268532245;
  testrel(18, NE, y, f, cnt,failed);

  y := jacobi_arcns(64, k);
  f := 0.1562627299478631160e-1;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arcnd;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcnd');

  k := 0.5;
  y := jacobi_arcnd(1, k);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1+1/1024, k);
  f := 0.0884677140647019265;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1.015625, k);
  f := 0.3588809717652129774;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1.0625, k);
  f := 0.7582115124496915987;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1.09375, k);
  f := 0.9762944264382908490;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1.125, k);
  f := 1.214034387269042046;
  testrel(6, NE, y, f, cnt,failed);

  x := 1.154296875; {~sqrt(4/3)};
  y := jacobi_arcnd(x, k);
  f := 1.632840659787347001;
  testrel(7, NE, y, f, cnt,failed);

  {-----}
  k := 2;

  y := jacobi_arcnd(1, k);
  f := 0.0;
  testrel(8, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1+1/1024, k);
  f := 0.02208989700966151846;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1.125, k);
  f := 0.2401398351599449428;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1.5, k);
  f := 0.4320593425516389592;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arcnd(4, k);
  f := 0.6974963126847653491;
  testrel(12, NE, y, f, cnt,failed);

  y := jacobi_arcnd(16, k);
  f := 0.8067750871667301122;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1024, k);
  f := 0.8423113587242981223;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1e6, k);
  f := 0.8428746000560288317;
  testrel(15, NE, y, f, cnt,failed);

  {-----}
{$ifdef CPUARM}
  k := 1+1/1024;
  y := jacobi_arcnd(16, k);
  f := 3.356848125290103039;
  testrel(16, 6, y, f, cnt,failed);

  y := jacobi_arcnd(1e6, k);
  f := 4.503479395715522424;
  testrel(17, 8, y, f, cnt,failed);
{$else}
  k := 1+1/1024;
  y := jacobi_arcnd(16, k);
  f := 3.356848125290103039;
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arcnd(1e6, k);
  f := 4.503479395715522424;
  testrel(17, NE, y, f, cnt,failed);
{$endif}

  k := 1+ldexpd(1,-30);
  y := jacobi_arcnd(2, k);
  f := 1.316957894698462387;
  testrel(18, NE, y, f, cnt,failed);

  k := 1024;
  y := jacobi_arcnd(2, k);
  f := 0.1022654001591530943e-2;
  testrel(19, NE, y, f, cnt,failed);

  {-----}
  x := 2;
  y := jacobi_arcnd(x,2);
  f := 0.5447753350259427046;
  testrel(20, NE, y, f, cnt,failed);

  k := 10;
  y := jacobi_arcnd(x,k);
  f := 0.1048738631962168580;
  testrel(21, NE, y, f, cnt,failed);

  k := 10000;
  y := jacobi_arcnd(x,k);
  f := 0.1047197552732059875e-3;
  testrel(22, NE, y, f, cnt,failed);

  k := 0.875;
  y := jacobi_arcnd(x,k);
  f := 1.892594693704786101;
  testrel(23, NE, y, f, cnt,failed);

  x := 5;
  k := 0.984375;
  y := jacobi_arcnd(x,k);
  f := 2.615385071029827950;
  testrel(24, NE2, y, f, cnt,failed);     {near boundary}

  y := jacobi_arcnd(10,1);
  f := 2.993222846126380898;
  testrel(25, NE, y, f, cnt,failed);

  k := 1-1/1024;
  y := jacobi_arcnd(1.0078125, k);
  f := 0.1250415108313409798;
  testrel(26, NE, y, f, cnt,failed);

  y := jacobi_arcnd(2, k);
  f := 1.319298803421139938;
  testrel(27, NE, y, f, cnt,failed);

  y := jacobi_arcnd(16, k);
  f := 3.624571434598265302;
  testrel(28, NE+1, y, f, cnt,failed);    {+1 for ARM}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arccd;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arccd');

  k := 0.5;
  y := jacobi_arccd(1, k);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arccd(0.9990234375, k);
  f := 0.5102965534724259892e-1;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arccd(0.75, k);
  f := 0.8140915468909887821;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arccd(0.5,k);
  f := 1.156321727760690225;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arccd(0.25,k);
  f := 1.432401697091115575;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arccd(0.125,k);
  f := 1.560340614548334362;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_arccd(0.0078125, k);
  f := 1.677937755468838312;
  testrel(7, NE, y, f, cnt,failed);

  y := jacobi_arccd(0, k);
  f := 1.685750354812596043;
  testrel(8, NE, y, f, cnt,failed);

  y := jacobi_arccd(-0.0078125, k);
  f := 1.693562954156353773;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_arccd(-0.75, k);
  f := 2.557409162734203304;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_arccd(-0.9990234375, k);
  f := 3.320471054277949487;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arccd(-1, k);
  f := 3.371500709625192086;
  testrel(12, NE, y, f, cnt,failed);

  {-----------}
  k := 2;
  y := jacobi_arccd(1, k);
  f := 0.0;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arccd(1.0009765625, k);
  f := 0.2550237905074735039e-1;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_arccd(1.5, k);
  f := 0.4703319880498214086;
  testrel(15, NE, y, f, cnt,failed);

  y := jacobi_arccd(2, k);
  f := 0.5781608638803451126;
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arccd(4, k);
  f := 0.7162008485455577875;
  testrel(17, NE, y, f, cnt,failed);

  y := jacobi_arccd(128, k);
  f := 0.8389688777344191562;
  testrel(18, NE, y, f, cnt,failed);

  f := 0.8428746774062980213;
  y := jacobi_arccd(1e6, k);
  testrel(19, NE, y, f, cnt,failed);

  {-----------}
  x := 16;
  k := 3;
  f := 0.5182804821705849459;
  y := jacobi_arccd(x,k);
  testrel(20, NE, y, f, cnt,failed);

  x := 100;
  k := 3;
  f := 0.5357955168104662044;
  y := jacobi_arccd(x,k);
  testrel(21, NE, y, f, cnt,failed);

  {0 <= x <= 1, k <= 1}
  k := 0.9375;
  x := 0.5;
  f := 1.942585238624454519;
  y := jacobi_arccd(x,k);
  testrel(22, NE, y, f, cnt,failed);

  x := 1/128;
  k := 3/4;
  f := 1.903177156571793799;
  y := jacobi_arccd(x,k);
  testrel(23, NE, y, f, cnt,failed);

  x := 0.875;
  k := 0;
  f := 0.5053605102841573070;
  y := jacobi_arccd(x,k);
  testrel(24, NE, y, f, cnt,failed);

  x := 0.875;
  k := ldexpd(1,-20);
  f := 0.5053605102843685299;
  y := jacobi_arccd(x,k);
  testrel(25, NE, y, f, cnt,failed);

  x := 0.875;
  k := 3/4;
  f := 0.7289526647151561956;
  y := jacobi_arccd(x,k);
  testrel(26, NE, y, f, cnt,failed);

  x := 0.875;
  k := 1-1/1024;
  f := 3.154546630787595406;
  y := jacobi_arccd(x,k);
  testrel(27, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_jacobi_arcdc;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcdc');

  k := 0.5;
  y := jacobi_arcdc(1.0,0.5);
  f := 0.0;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arcdc(1+1/1024,k);
  f := 0.05100475810149470078;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arcdc(1.25,k);
  f := 0.7281312548485799106;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arcdc(2,k);
  f := 1.156321727760690225;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arcdc(10,k);
  f := 1.585541094196814353;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arcdc(1e6,k);
  f := 1.685749354812596043;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_arcdc(-1-1/1024,k);
  f := 3.320495951523697385;
  testrel(7, NE, y, f, cnt,failed);

  y := jacobi_arcdc(-1.5,k);
  f := 2.430836733525549269;
  testrel(8, NE, y, f, cnt,failed);

  {-----------}
  k := 2;
  y := jacobi_arcdc(0,k);
  f := 0.8428751774062980214;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_arcdc(1/128,k);
  f := 0.8389688777344191562;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_arcdc(0.0625,k);
  f := 0.8115997020324605835;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arcdc(0.5,k);
  f := 0.5781608638803451126;
  testrel(12, NE, y, f, cnt,failed);

  y := jacobi_arcdc(1-1/128,k);
  f := 0.07215342403220674891;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arcdc(-0.75,k);
  f := 1.278704581367101652;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_arcdc(-0.5,k);
  f := 1.107589490932250930;
  testrel(15, NE, y, f, cnt,failed);


  {-----------}
  k := 1/4;
  y := jacobi_arcdc(10,k);
  f := 1.496064349943746668;
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arcdc(2,k);
  f := 1.071217826209477658;
  testrel(17, NE, y, f, cnt,failed);

  y := jacobi_arcdc(2,0);
  f := 1.047197551196597746;
  testrel(18, NE, y, f, cnt,failed);

  y := jacobi_arcdc(0.5,3);
  f := 0.3637463872488045010;
  testrel(19, NE, y, f, cnt,failed);

  y := jacobi_arcdc(-0.5,3);
  f := 0.7145114365010171167;
  testrel(20, NE, y, f, cnt,failed);

  y := jacobi_arcdc(0.0625,3);
  f := 0.5182804821705849459;
  testrel(21, NE, y, f, cnt,failed);

  k := 1-1/1024;
  x := 1+1/2048;
  f := 0.6585890808021593225;
  y := jacobi_arcdc(x,k);
  testrel(22, NE, y, f, cnt,failed);

  f := 4.507412597899042266;
  y := jacobi_arcdc(1e6,k);
  testrel(23, NE, y, f, cnt,failed);

  k := 1+1/1024;

  x := 1-1/8192;
  f := 0.3464959454674026608;
  y := jacobi_arcdc(x,k);
  testrel(24, NE, y, f, cnt,failed);

  x := 0;
  f := 4.503502017610268965;
  y := jacobi_arcdc(x,k);
  testrel(25, NE, y, f, cnt,failed);

  x := -1;
  f := 9.007004035220537931;
  y := jacobi_arcdc(x,k);
  testrel(26, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arcsd;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 8;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcsd');

  k := 0.5;
  y := jacobi_arcsd(0,k);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arcsd(-0.0625,k);
  f := -0.06252038087319812021;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arcsd(0.125,k);
  f := 0.1251639123384003270;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arcsd(-0.25,k);
  f := -0.2513396906799018347;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arcsd(0.5,k);
  f := 0.5117293103256968015;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arcsd(-0.75,k);
  f := -0.7970472253705041676;
  testrel(6, NE, y, f, cnt,failed);

  y := jacobi_arcsd(1,k);
  f := 1.156321727760690225;
  testrel(7, NE, y, f, cnt,failed);

  y := jacobi_arcsd(-1.125,k);
  f := -1.457964146174858544;
  testrel(8, NE, y, f, cnt,failed);

  x := 1.154541015625; {~sqrt(4/3)};
  y := jacobi_arcsd(x,k);
  f := 1.669127652279864100;
  testrel(9, NE, y, f, cnt,failed);

  {-----------}
  k := 2;

  y := jacobi_arcsd(0,k);
  f := 0.0;
  testrel(10, NE, y, f, cnt,failed);

  x := 1e-6;
  f := 0.9999999999988333333e-6;
  y := jacobi_arcsd(x,k);
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arcsd(0.125,k);
  f := 0.1227940011442581578;
  testrel(12, NE, y, f, cnt,failed);

  y := jacobi_arcsd(-0.5,k);
  f := -0.4021830506160327783;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arcsd(1,k);
  f := 0.5781608638803451126;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_arcsd(-2,k);
  f := -0.7018985502494669139;
  testrel(15, NE, y, f, cnt,failed);

  y := jacobi_arcsd(4,k);
  f := 0.7711401374755114200;
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arcsd(-8,k);
  f := -0.8068454503662184059;
  testrel(17, NE, y, f, cnt,failed);

  y := jacobi_arcsd(1e6,k);
  f := 0.8428748887311634267;
  testrel(18, NE, y, f, cnt,failed);

  {-----------}
  k := 1.125;
  x := 1.5;
  f := 1.074318100359308022;
  y := jacobi_arcsd(x,k);
  testrel(19, NE, y, f, cnt,failed);

  k := 4/3;
  x := 3;
  f := 1.158876953547957870;
  y := jacobi_arcsd(x,k);
  testrel(20, NE, y, f, cnt,failed);

  x := 1;
  k := 0.0;
  y := jacobi_arcsd(x,k);
  f := Pi_2;
  testrel(21, NE, y, f, cnt,failed);

  x := 2;
  k := 0.875;
  y := jacobi_arcsd(x,k);
  f := 1.930730485989356762;
  testrel(22, NE, y, f, cnt,failed);

  x := 5;
  k := 0.984375;
  y := jacobi_arcsd(x,k);
  f := 2.625095962638367141;
  testrel(23, NE2, y, f, cnt,failed);  {near boundary}

  x := 55;
  k := 1;
  y := jacobi_arcsd(x,k);
  f := 4.700563000177194743;
  testrel(24, NE, y, f, cnt,failed);

  k := 1+1/1024;
  x := -1/128;
  y := jacobi_arcsd(x,k);
  f := -0.7812420218741962711e-2;
  testrel(25, NE, y, f, cnt,failed);

  x := 1e6;
  y := jacobi_arcsd(x,k);
  f := 4.503479417785663647;
  testrel(26, NE, y, f, cnt,failed);

  k := 1-1/1024;
  x := 0.015625;
  y := jacobi_arcsd(x,k);
  f := 0.1562436676891280678e-1;
  testrel(27, NE, y, f, cnt,failed);

  x := 16;
  y := jacobi_arcsd(x,k);
  f := 3.625954867718658024;
  testrel(28, NE+4, y, f, cnt,failed);    {+4 for ARM}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_jacobi_arcds;
var
  x,k,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 64;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','jacobi_arcds');

  k := 0.5;
  x := 0.8671875; {~sqrt(3)/2;}
  y := jacobi_arcds(x, k);
  f := 1.63396862632269164172604217374;
  testrel(1, NE, y, f, cnt,failed);

  y := jacobi_arcds(-0.875, k);
  f := -1.542280274951022840;
  testrel(2, NE, y, f, cnt,failed);

  y := jacobi_arcds(1, k);
  f := 1.156321727760690225;
  testrel(3, NE, y, f, cnt,failed);

  y := jacobi_arcds(-2, k);
  f := -0.5117293103256968015;
  testrel(4, NE, y, f, cnt,failed);

  y := jacobi_arcds(8, k);
  f := 0.1251639123384003270;
  testrel(5, NE, y, f, cnt,failed);

  y := jacobi_arcds(1e6, k);
  f := 1.000000000000083333e-6;
  testrel(6, NE, y, f, cnt,failed);

  {----------}
  k := 2;
  y := jacobi_arcds(0, k);
  f := 0.8428751774062980214;
  testrel(7, NE, y, f, cnt,failed);

  y := jacobi_arcds(1e-6, k);
  f := 0.8428748887311634267;
  testrel(8, NE, y, f, cnt,failed);

  y := jacobi_arcds(1/1024, k);
  f := 0.8425932681213084324;
  testrel(9, NE, y, f, cnt,failed);

  y := jacobi_arcds(0.5, k);
  f := 0.7018985502494669139;
  testrel(10, NE, y, f, cnt,failed);

  y := jacobi_arcds(1, k);
  f := 0.5781608638803451126;
  testrel(11, NE, y, f, cnt,failed);

  y := jacobi_arcds(-2, k);
  f := -0.4021830506160327783;
  testrel(12, NE, y, f, cnt,failed);

  y := jacobi_arcds(4, k);
  f := 0.2338595256752966072;
  testrel(13, NE, y, f, cnt,failed);

  y := jacobi_arcds(-8, k);
  f := -0.1227940011442581578;
  testrel(14, NE, y, f, cnt,failed);

  y := jacobi_arcds(1e6, k);
  f := 9.999999999988333333e-7;
  testrel(15, NE, y, f, cnt,failed);

  {----------}
  k := 1-1/1024;

  y := jacobi_arcds(1e6, k);
  f := 0.9999999999998339841e-6;
  testrel(16, NE, y, f, cnt,failed);

  y := jacobi_arcds(8, k);
  f := 0.1246780061216979402;
  testrel(17, NE, y, f, cnt,failed);

  y := jacobi_arcds(0.5, k);
  f := 1.445655745153180149;
  testrel(18, NE, y, f, cnt,failed);

  y := jacobi_arcds(0.0625, k);
  f := 3.62595486771865802422708362769;
  testrel(19, NE+4, y, f, cnt,failed);    {+4 for ARM}

  {----------}
  k := 1/128; {i.e. x > 0.9999694819562}
  x := 0.999969482421875;
  y := jacobi_arcds(x, k);
  f := 1.570789777558152897;
  testrel(20, NE2, y, f, cnt,failed);   {!!!!}

  y := jacobi_arcds(1, k);
  f := 1.563007716587740338;
  testrel(21, NE, y, f, cnt,failed);

  y := jacobi_arcds(8, k);
  f := 0.1253277910543881852;
  testrel(22, NE, y, f, cnt,failed);

  {----------}
  x := 1;
  k := 0;
  y := jacobi_arcds(x,k);
  f := 1.570796326794896619;
  testrel(23, NE, y, f, cnt,failed);

  x := 0;
  k := 3;
  y := jacobi_arcds(x,k);
  f := 0.5391289118749108089;
  testrel(24, NE, y, f, cnt,failed);

  x := 2;
  k := 0.75;
  y := jacobi_arcds(x,k);
  f := 0.4981805085698164591;
  testrel(25, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


end.
