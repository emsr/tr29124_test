{Part 8 of regression test for SPECFUN unit  (c) 2012  W.Ehrhardt}

unit t_sfd8;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_zeta;
procedure test_zeta1p;
procedure test_zetaint;

{$ifndef VER50}
procedure test_zetam1;
procedure test_etam1;
procedure test_eta;
procedure test_LerchPhi;
procedure test_DirichletBeta;
procedure test_DirichletLambda;
procedure test_LegendreChi;
procedure test_polylog;
procedure test_polylogr;
procedure test_fermi_dirac;
procedure test_fermi_dirac_half;
{$endif}


implementation

uses
  amath, specfun, t_sfd0;



{---------------------------------------------------------------------------}
procedure test_zeta;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','zeta');

  x := 0.0;
  y := zeta(x);
  f := -0.5;
  testrel( 1, NE, y, f, cnt,failed);

  x := 0.5;
  y := zeta(x);
  f := -1.460354508809586813;
  testrel( 2, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := zeta(x);
  f := -1023.422855448942979;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := zeta(x);
  f := 1024.577286769504594;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1.5;
  y := zeta(x);
  f := 2.6123753486854883433;
  testrel( 5, NE, y, f, cnt,failed);

  x := 1.9990234375;
  y := zeta(x);
  f := 1.645850590810321192;
  testrel( 6, NE, y, f, cnt,failed);

  x := 2.0009765625;
  y := zeta(x);
  f := 1.644019440013418371;
  testrel( 7, NE, y, f, cnt,failed);

  x := 3.0;
  y := zeta(x);
  f := 1.202056903159594285;
  testrel( 8, NE, y, f, cnt,failed);

  x := 3.9990234375;
  y := zeta(x);
  f := 1.082390560902667758;
  testrel( 9, NE, y, f, cnt,failed);

  x := 4.0009765625;
  y := zeta(x);
  f := 1.082255968563913698;
  testrel(10, NE, y, f, cnt,failed);

  x := 6;
  y := zeta(x);
  f := power(Pi,6)/945.0;
  testrel(11, NE, y, f, cnt,failed);

  x := 6.9990234375;
  y := zeta(x);
  f := 1.008355171623342114;
  testrel(12, NE, y, f, cnt,failed);

  x := 7.0009765625;
  y := zeta(x);
  f := 1.008343387409452773;
  testrel(13, NE, y, f, cnt,failed);

  x := 14;
  y := zeta(x);
  f := 2.0*power(Pi,14)/18243225.0;
  testrel(14, NE, y, f, cnt,failed);

  x := 14.9990234375;
  y := zeta(x);
  f := 1.000030608976823097;
  testrel(15, NE, y, f, cnt,failed);

  x := 15.0009765625;
  y := zeta(x);
  f := 1.0000305675098559807;
  testrel(16, NE, y, f, cnt,failed);

  x := 30.0;
  y := zeta(x);
  f := 1.0000000009313274324;
  testrel(17, NE, y, f, cnt,failed);

  x := 41.9990234375;
  y := zeta(x);
  f := 1.00000000000022752765;
  testrel(18, NE, y, f, cnt,failed);

  x := 42.0009765625;
  y := zeta(x);
  f := 1.00000000000022721983;
  testrel(19, NE, y, f, cnt,failed);

  x := 60.0;
  y := zeta(x);
  f := 1.00000000000000000087;
  testrel(20, NE, y, f, cnt,failed);

  x := 64.0;
  y := zeta(x);
  f := 1.0;
  testrel(21, NE, y, f, cnt,failed);

  x := -1.0;
  y := zeta(x);
  f := -1.0/12.0;
  testrel(22, NE, y, f, cnt,failed);

  x := -1.5;
  y := zeta(x);
  f := -0.254852018898330359495e-1;
  testrel(23, NE, y, f, cnt,failed);

  x := -2.0;
  y := zeta(x);
  f := 0.0;
  testrel(24, NE, y, f, cnt,failed);

  x := -2.0009765625;
  y := zeta(x);
  f := 0.297034757704571229761e-4;
  testrel(25, NE, y, f, cnt,failed);

  x := -3.0;
  y := zeta(x);
  f := 1.0/120.0;
  testrel(26, NE, y, f, cnt,failed);

  x := -55.5;
  y := zeta(x);
  f := 0.107190862583058341659e30;
  testrel(27, NE, y, f, cnt,failed);

  y := zeta(1.5e-10);
  f := -0.5000000001378407800;
  testrel(28, NE, y, f, cnt,failed);

  y := zeta(-1e-10);
  f := -0.4999999999081061467;
  testrel(29, NE, y, f, cnt,failed);

  {Test values for s < -MaxGAMD}
  y := zeta(-172.5);
  f := -0.1301138852431751653e175;
  testrel(30, NE, y, f, cnt,failed);

  y := zeta(-234.5);
  f := 0.5247027227116272750e268;
  testrel(31, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_zeta1p;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','zeta1p');

  x := 1.0;
  y := zeta1p(x);
  f := 1.64493406684822643647;
  testrel( 1, NE, y, f, cnt,failed);

  x := -0.9990234375;
  y := zeta1p(x);
  f := -0.500898358549607584099;
  testrel( 2, NE, y, f, cnt,failed);

  x := -1e-5;
  y := zeta1p(x);
  f := -99999.4227850632574065;
  testrel( 3, NE, y, f, cnt,failed);

  x := +1e-5;
  y := zeta1p(x);
  f := 100000.577216393059503;
  testrel( 4, NE, y, f, cnt,failed);

  x := -1e-10;
  y := zeta1p(x);
  f := -9999999999.42278433511;
  testrel( 5, NE, y, f, cnt,failed);

  x := +1e-10;
  y := zeta1p(x);
  f := 10000000000.5772156649;
  testrel( 6, NE, y, f, cnt,failed);

  x := -1e-18;
  y := zeta1p(x);
  f := -999999999999999999.423;
  testrel( 7, NE, y, f, cnt,failed);

  x := +1e-18;
  y := zeta1p(x);
  f := 1000000000000000000.58;
  testrel( 8, NE, y, f, cnt,failed);

  x := +1e-25;
  y := zeta1p(x);
  f := 1/x;
  testrel( 9, NE, y, f, cnt,failed);

  x := -1e-25;
  y := zeta1p(x);
  f := 1/x;
  testrel(10, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_zetaint;
var
  x,y,f: double;
  i: integer;
begin
  writeln('Function: ','zetaint');
  {Because zeta now uses zetaint test with zeta(succ(i)) and a suitable tolerance}
  for i:=-260 to 100 do if i<>1 then begin
    if i<0 then x := abs(i*eps_d*10)
    else x := 10*eps_d;
    y := zetaint(i);
    if (i<0) and (i and 1 = 0) then f := 0.0
    else f := zeta(succd(i));
    if IsInf(y) or IsInf(f) then begin
      if f<>y then begin
       writeln('*** failed first for ',i);
       inc(total_failed);
       exit;
     end;
    end
    else if abs(y-f) > abs(f)*x then begin
      writeln('*** failed first for ',i);
      inc(total_failed);
      exit;
    end;
  end;
  writeln(' Tests OK');
end;


{$ifndef VER50}
{---------------------------------------------------------------------------}
procedure test_zetam1;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','zetam1');

  x := 0.0;
  y := zetam1(x);
  f := -1.5;
  testrel( 1, NE, y, f, cnt,failed);

  x := 0.5;
  y := zetam1(x);
  f := -2.4603545088095868129;
  testrel( 2, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := zetam1(x);
  f := -1024.4228554489429787;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := zetam1(x);
  f := 1023.57728676950459406;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1.5;
  y := zetam1(x);
  f := 1.61237534868548834335;
  testrel( 5, NE, y, f, cnt,failed);

  x := 1.9990234375;
  y := zetam1(x);
  f := 0.64585059081032119165;
  testrel( 6, NE, y, f, cnt,failed);

  x := 2.0009765625;
  y := zetam1(x);
  f := 0.64401944001341837056;
  testrel( 7, NE, y, f, cnt,failed);

  x := 3.0;
  y := zetam1(x);
  f := 0.20205690315959428540;
  testrel( 8, NE, y, f, cnt,failed);

  x := 6.0;
  y := zetam1(x);
  f := 0.17343061984449139715e-1;
  testrel( 9, NE, y, f, cnt,failed);

  x := 20;
  y := zetam1(x);
  f := 0.95396203387279611315e-6;
  testrel(10, NE, y, f, cnt,failed);

  x := 40.0;
  y := zetam1(x);
  f := 0.90949478402638892825e-12;
  testrel(11, NE, y, f, cnt,failed);

  x := 64.0;
  y := zetam1(x);
  f := 0.54210108624566454109e-19;
  testrel(12, NE, y, f, cnt,failed);

  x := 100.0;
  y := zetam1(x);
  f := 0.7888609052210118074e-30;
  testrel(13, NE, y, f, cnt,failed);

  x := 1000.0;
  y := zetam1(x);
  f := 9.33263618503218879e-302;
  testrel(14, NE, y, f, cnt,failed);

  x := 12500.0;
  y := zetam1(x);
  f := 0;
  testrel(15, NE, y, f, cnt,failed);

  x := -1.0;
  y := zetam1(x);
  f := -1.083333333333333333;
  testrel(16, NE, y, f, cnt,failed);

  x := -1.5;
  y := zetam1(x);
  f := -1.025485201889833036;
  testrel(17, NE, y, f, cnt,failed);

  x := -2.0;
  y := zetam1(x);
  f := -1.0;
  testrel(18, NE, y, f, cnt,failed);

  x := -2.0009765625;
  y := zetam1(x);
  f := -0.999970296524229543;
  testrel(19, NE, y, f, cnt,failed);

  x := -19.25;
  y := zetam1(x);
  f := 31.49461534567477545;
  testrel(20, NE, y, f, cnt,failed);

  x := -55.5;
  y := zetam1(x);
  f := 0.10719086258305834166e30;
  testrel(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_etam1;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','etam1');

  {etam1 :=  x -> Zeta(x)*(1-2^(1-x)) - 1.0;}

  x := 0.0;
  y := etam1(x);
  f := -0.5;
  testrel( 1, NE, y, f, cnt,failed);

  x := 0.5;
  y := etam1(x);
  f := -0.3951013565783696298;
  testrel( 2, NE, y, f, cnt,failed);

  x := 0.9990234375;
  y := etam1(x);
  f := -0.3070089725899074737;
  testrel( 3, NE, y, f, cnt,failed);

  x := 1.0009765625;
  y := etam1(x);
  f := -0.3066967286343630638;
  testrel( 4, NE, y, f, cnt,failed);

  x := 1.5;
  y := etam1(x);
  f := -0.2348529753745920546;
  testrel( 5, NE, y, f, cnt,failed);

  x := 1.9990234375;
  y := etam1(x);
  f := -0.1776319325704589369;
  testrel( 6, NE, y, f, cnt,failed);

  x := 2.0009765625;
  y := etam1(x);
  f := -0.1774340486232085190;
  testrel( 7, NE, y, f, cnt,failed);

  x := 3.0;
  y := etam1(x);
  f := -0.9845732263030428595e-1;
  testrel( 8, NE, y, f, cnt,failed);

  x := 6.0;
  y := etam1(x);
  f := -0.1444890870256489590e-1;
  testrel( 9, NE, y, f, cnt,failed);

  x := 20;
  y := etam1(x);
  f := -0.9533884184778849492e-6;
  testrel(10, NE, y, f, cnt,failed);

  x := 40.0;
  y := etam1(x);
  f := -0.9094946195211219090e-12;
  testrel(11, NE, y, f, cnt,failed);

  x := 64.0;
  y := etam1(x);
  f := -0.5421010862398398930e-19;
  testrel(12, NE, y, f, cnt,failed);

  x := 100.0;
  y := etam1(x);
  f := -0.7888609052210118035e-30;
  testrel(13, NE, y, f, cnt,failed);

  x := 1000.0;
  y := etam1(x);
  f := -9.33263618503218879e-302;
  testrel(14, NE, y, f, cnt,failed);

  x := 12500.0;
  y := etam1(x);
  f := 0;
  testrel(15, NE, y, f, cnt,failed);

  x := -1.0;
  y := etam1(x);
  f := -0.75;
  testrel(16, NE, y, f, cnt,failed);

  x := -1.5;
  y := etam1(x);
  f := -0.8813191292801597880;
  testrel(17, NE, y, f, cnt,failed);

  x := -2.0;
  y := etam1(x);
  f := -1.0;
  testrel(18, NE, y, f, cnt,failed);

  x := -2.0009765625;
  y := etam1(x);
  f := -1.0002080852354742793;
  testrel(19, NE, y, f, cnt,failed);

  x := -19.25;
  y := etam1(x);
  f := -40519910.275413219422;
  testrel(20, NE, y, f, cnt,failed);

  x := -55.5;
  y := etam1(x);
  f := -0.10923266281825727645e47;
  testrel(21, NE, y, f, cnt,failed);

  x := -0.9e-9;
  y := etam1(x);
  f := -0.5000000002032122174;
  testrel(22, NE, y, f, cnt,failed);

  x := 1e-10;
  y := etam1(x);
  f := -0.4999999999774208647;
  testrel(23, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_LerchPhi;
var
  y,f,z,s,a: double;
  cnt, failed: integer;
const
  NE  = 1;
  NE2 = 3;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','LerchPhi');

  z := -1;
  y := LerchPhi(z, 2, 3);
  f := 0.7246703342411321824e-1;
  testrel( 1, NE, y, f, cnt,failed);

  y := LerchPhi(z, 100, 3);
  f := 1.940325217482010536e-48;
  testrel( 2, NE, y, f, cnt,failed);

  y := LerchPhi(z, 0.125, 5);
  f := 0.4139549024672754777;
  testrel( 3, NE, y, f, cnt,failed);

  y := LerchPhi(z, 16, 0.5);
  f := 65535.99847798871839;
  testrel( 4, NE, y, f, cnt,failed);

  y := LerchPhi(z, 2, 0.5);
  f := 3.663862376708876060;
  testrel( 5, NE, y, f, cnt,failed);

  y := LerchPhi(0.5, 16, 1e-19);
  f := 1e304;
  testrel( 6, NE2, y, f, cnt,failed);

  y := LerchPhi(0.5, 2, 3);
  f := 0.1579242117201000472;
  testrel( 7, NE, y, f, cnt,failed);

  z := -1+ldexp(1,-52);
  y := LerchPhi(z, 1+1/128, 1e10);
  f := 4.176812734999603534e-11;
  testrel( 8, NE, y, f, cnt,failed);

  z := 1-ldexp(1,-52);
  y := LerchPhi(z, 1+1/128, 1e10);
  f := 9.89847701124157791;
  testrel( 9, NE, y, f, cnt,failed);

  y := LerchPhi(0.75, 1+1/128, 200);
  f := 0.1890864747918034799e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := LerchPhi(-0.75, 1+1/128, 200);
  f := 0.2747206889737371425e-2;
  testrel(11, NE, y, f, cnt,failed);

  y := LerchPhi(-0.75, 1+1/128, 30);
  f := 0.1881378370869183099e-1;
  testrel(12, NE, y, f, cnt,failed);

  z := 1-1/16;
  y := LerchPhi(z, 2, 10000);
  f := 1.595222182725875097e-7;
  testrel(13, NE, y, f, cnt,failed);

  y := LerchPhi(z, 2, 1000);
  f := 0.1554103481247648889e-4;
  testrel(14, NE, y, f, cnt,failed);

  y := LerchPhi(0.3,0.1,100);
  f := 0.9009852125614140412;
  testrel(14, NE, y, f, cnt,failed);

  y := LerchPhi(-0.3,0.1,100);
  f := 0.4854634757432184874;
  testrel(15, NE, y, f, cnt,failed);

  y := LerchPhi(-0.3,100,100);
  f := 0.9002485412225464125e-200;
  testrel(16, NE, y, f, cnt,failed);

  y := LerchPhi(0.75,3,1e-10);
  f := 1e30;
  testrel(18, NE, y, f, cnt,failed);

  y := LerchPhi(0.75,3,1e-4);
  f := 0.1000000000000844188e13;
  testrel(19, NE2, y, f, cnt,failed);

  z := -0.96875;
  s := 16;
  y := LerchPhi(z,s,0.25);
  f := 4.294967295972734281e9;
  testrel(20, NE, y, f, cnt,failed);

  z := -0.96875;
  y := LerchPhi(z,s,1/64);
  f := 7.922816251426433759e28;
  testrel(21, NE, y, f, cnt,failed);

  z := -0.96875;
  y := LerchPhi(z,s,1);
  f := 0.9999852396432582781;
  testrel(22, NE, y, f, cnt,failed);

  z := -0.96875;
  s := 4;
  y := LerchPhi(z,s,0.25);
  f := 255.6336018166029893;
  testrel(23, NE, y, f, cnt,failed);

  y := LerchPhi(z,s,0.015625);
  f := 0.1677721513782503849e8;
  testrel(24, NE, y, f, cnt,failed);

  y := LerchPhi(z,s,0.75);
  f := 3.070223141356004934;
  testrel(25, NE, y, f, cnt,failed);

  z := 0.96875;
  s := 16;
  y := LerchPhi(z, s, 1/16);
  f := 1.844674407370955162e19;
  testrel(26, NE, y, f, cnt,failed);

  y := LerchPhi(z, s, 1/8);
  f := 2.814749767106561472e14;
  testrel(27, NE, y, f, cnt,failed);

  y := LerchPhi(z, s, 2);
  f := 0.1528151848585480626e-4;
  testrel(28, NE, y, f, cnt,failed);

  y := LerchPhi(z, s, 20);
  f := 2.813806862758912383e-21;
  testrel(29, NE, y, f, cnt,failed);

  y := LerchPhi(z, s, 1);
  f := 1.000014803971033172;
  testrel(30, NE, y, f, cnt,failed);

  y := LerchPhi(z, s, 0.5);
  f := 65536.00147526752452;
  testrel(31, NE, y, f, cnt,failed);

  y := LerchPhi(z, 2, 1e100);
  f := 32e-200;
  testrel(32, NE, y, f, cnt,failed);

  y := LerchPhi(z, 2, 1000);
  f := 0.3018300795020214889e-4;
  testrel(33, NE, y, f, cnt,failed);

  y := LerchPhi(z, 8, 1/1024);
  f := 1.208925819614629175e24;
  testrel(34, NE, y, f, cnt,failed);

  y := LerchPhi(z, 32, 1/128);
  f := 2.695994666715063979e67;
  testrel(35, NE, y, f, cnt,failed);

  y := LerchPhi(z, 32, 1/64);
  f := 6.277101735386680764e57;
  testrel(36, NE, y, f, cnt,failed);

  y := LerchPhi(z, 32, 1);
  f := 1.000000000225555193;
  testrel(37, NE, y, f, cnt,failed);

  y := LerchPhi(z, 32, 2);
  f := 2.328311664999511981e-10;
  testrel(38, NE, y, f, cnt,failed);

  {Test cases from Aksenov/Jentschura}
  {Run 5}
  z := 0.99999;
  s := 2;
  a := 1;
  y := LerchPhi(z, s, a);
  f := 1.644825385246778980;
  testrel(39, NE2, y, f, cnt,failed);

  {Run 6}
  z := -0.99999;
  s := 2;
  a := 1;
  y := LerchPhi(z, s, a);
  f := 8.224683266259164962E-1;
  testrel(40, NE, y, f, cnt,failed);

  {Run 9}
  z := 1e-5;
  s := 2;
  a := 1;
  y := LerchPhi(z, s, a);
  f := 1.000002500011111174;
  testrel(41, NE, y, f, cnt,failed);

  {Run 10}
  z := -6.3E-5;
  s := 2;
  a := 1;
  y := LerchPhi(z, s, a);
  f := 0.999984250440984373;
  testrel(42, NE, y, f, cnt,failed);

  {Run 11}
  z := 3.4709929976435479E-6;
  s := 1;
  a := 1.5172413793103448;
  y := LerchPhi(z, s, a);
  f := 0.6590922879819636657;
  testrel(43, NE, y, f, cnt,failed);

  {More z=-1 cases}
  y := LerchPhi(-1, 16, 10);
  f := 0.82489438389165980294e-16;
  testrel(44, NE, y, f, cnt,failed);

  y := LerchPhi(-1, 2, 1e-10);
  f := 1e20;
  testrel(45, NE, y, f, cnt,failed);

  y := LerchPhi(-1, 1+1/128, 1e-10);
  f := 0.1197085030426290549e11;
  testrel(46, NE, y, f, cnt,failed);

  f := 0.5000000004e-160;
  y := LerchPhi(-1, 16, 1e10);
  testrel(47, NE, y, f, cnt,failed);

  y := LerchPhi(-1, 16, 1);
  f := 0.9999847642149061064;
  testrel(48, NE, y, f, cnt,failed);

  y := LerchPhi(-1, 16, 1e-10);
  f := 1e160;
  testrel(49, NE2, y, f, cnt,failed);

  y := LerchPhi(-1,0.125,1e24);
  f := 0.5e-3;
  testrel(50, NE, y, f, cnt,failed);

  {Tests for s<0}
  y := LerchPhi(-0.9921875, -0.9375, 5);
  f := 2.057733305230137650;
  testrel(51, NE, y, f, cnt,failed);

  y := LerchPhi(-0.9921875, -0.5, 5);
  f := 1.066663573920271744;
  testrel(52, NE, y, f, cnt,failed);

  y := LerchPhi(-0.9921875, -0.125, 5);
  f := 0.6062199942770585427;
  testrel(53, NE, y, f, cnt,failed);

  y := LerchPhi(-0.25, -0.9375, 7.5);
  f := 5.158129174261488297;
  testrel(54, NE, y, f, cnt,failed);

  y := LerchPhi(-0.25, -0.5, 7.5);
  f := 2.162257387021207087;
  testrel(55, NE, y, f, cnt,failed);

  y := LerchPhi(-0.25, -0.125, 7.5);
  f := 1.025823888398309961;
  testrel(56, NE, y, f, cnt,failed);

  y := LerchPhi(0.375, -0.9375, 2.5);
  f := 4.609912789035665679;
  testrel(57, NE, y, f, cnt,failed);

  y := LerchPhi(0.375, -0.5, 2.5);
  f := 2.788814534753693403;
  testrel(58, NE, y, f, cnt,failed);

  y := LerchPhi(0.375, -0.125, 2.5);
  f := 1.835284527035582744;
  testrel(59, NE, y, f, cnt,failed);

  y := LerchPhi(0.75, -0.9375, 12.5);
  f := 52.17092092836499335;
  testrel(60, NE, y, f, cnt,failed);

  y := LerchPhi(0.75, -0.5, 12.5);
  f := 15.66354877711018779;
  testrel(61, NE, y, f, cnt,failed);

  y := LerchPhi(0.75, -0.125, 12.5);
  f := 5.621603260975934263;
  testrel(62, NE, y, f, cnt,failed);

  y := LerchPhi(-0.875, -0.125, 1000);
  f := 1.264658869110458885;
  testrel(63, NE, y, f, cnt,failed);

  y := LerchPhi(0.875, -0.125, 1000);
  f := 18.98748182387054361;
  testrel(64, NE, y, f, cnt,failed);

  y := LerchPhi(-0.375, -0.125, 2000);
  f := 1.880696209365743222;
  testrel(65, NE, y, f, cnt,failed);

  y := LerchPhi(0.375, -0.125, 2000);
  f := 4.137757266383553387;
  testrel(66, NE, y, f, cnt,failed);

  y := LerchPhi(0.375, -0.125, 1e6);
  f := 8.997461877854526012;
  testrel(67, NE, y, f, cnt,failed);

  y := LerchPhi(-0.75, -0.125, 1e6);
  f := 3.213378828942416160;
  testrel(68, NE, y, f, cnt,failed);

  y := LerchPhi(-0.75, -1, 5);
  f := 2.612244897959183673;
  testrel(69, NE, y, f, cnt,failed);

  y := LerchPhi(0.375, -1, 5);
  f := 8.96;
  testrel(70, NE, y, f, cnt,failed);

  {z=1, no special case}
  y := LerchPhi(1,0.75,1.5);
  f := -4.028036535056665957;
  testrel(71, NE, y, f, cnt,failed);

  y := LerchPhi(1,-0.75,1.5);
  f := -0.5404252521791210354;
  testrel(72, NE, y, f, cnt,failed);

  {z=0}
  y := LerchPhi(0,0.75,1.5);
  f := 0.7377879464668810616;
  testrel(73, NE, y, f, cnt,failed);

  y := LerchPhi(0,-0.75, 0.5);
  f := 0.5946035575013605334;
  testrel(74, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_DirichletBeta;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','DirichletBeta');

  y := DirichletBeta(40.5);
  f := 1.0;
  testrel( 1, NE, y, f, cnt,failed);

  y := DirichletBeta(32);
  f := 0.9999999999999994603;
  testrel( 2, NE, y, f, cnt,failed);

  y := DirichletBeta(27.6);
  f := 0.9999999999999321651;
  testrel( 3, NE, y, f, cnt,failed);

  y := DirichletBeta(22.8);
  f := 0.9999999999867678140;
  testrel( 4, NE, y, f, cnt,failed);

  y := DirichletBeta(22.0);
  f := 0.999999999968134064;
  testrel( 5, NE, y, f, cnt,failed);

  y := DirichletBeta(10.0);
  f := 0.9999831640261968774;
  testrel( 6, NE, y, f, cnt,failed);

  y := DirichletBeta(5.0);
  f := 0.996157828077088064;
  testrel( 7, NE, y, f, cnt,failed);

  y := DirichletBeta(2.0); {Catalan}
  f := 0.9159655941772190151;
  testrel( 8, NE, y, f, cnt,failed);

  y := DirichletBeta(1.0);
  f := Pi/4;
  testrel( 9, NE, y, f, cnt,failed);

  y := DirichletBeta(0.25);
  f := 0.5907230564424947319;
  testrel(10, NE, y, f, cnt,failed);

  y := DirichletBeta(1/128);
  f := 0.5030522278808478693;
  testrel(11, NE, y, f, cnt,failed);

  y := DirichletBeta(ldexp(1,-30));
  f := 0.5000000003647006979;
  testrel(12, NE, y, f, cnt,failed);

  y := DirichletBeta(1e-10);
  f := 0.5000000000391594393;
  testrel(13, NE, y, f, cnt,failed);

  y := DirichletBeta(ldexp(1,-42));
  f := 0.5000000000000890383;
  testrel(14, NE, y, f, cnt,failed);

  y := DirichletBeta(0);
  f := 0.5;
  testrel(15, NE, y, f, cnt,failed);

  y := DirichletBeta(-1/1024);
  f := 0.4996174725733061643;
  testrel(16, NE, y, f, cnt,failed);

  y := DirichletBeta(-0.75);
  f := 0.1425165634878215728;
  testrel(17, NE, y, f, cnt,failed);

  y := DirichletBeta(-3);
  f := 0.0;
  testrel(18, NE, y, f, cnt,failed);

  y := DirichletBeta(-3.25);
  f := 0.4612305190414612658;
  testrel(19, NE, y, f, cnt,failed);

  y := DirichletBeta(-4);
  f := 5/2;
  testrel(20, NE, y, f, cnt,failed);

  y := DirichletBeta(-9.875);
  f := -19551.74572283995210;
  testrel(21, NE, y, f, cnt,failed);

  y := DirichletBeta(-10.125);
  f := -31440.68339647519141;
  testrel(22, NE, y, f, cnt,failed);

  y := DirichletBeta(-13);
  f := 0;
  testrel(23, NE, y, f, cnt,failed);

  y := DirichletBeta(-13.03125);
  f := -586962.5557661498894;
  testrel(24, NE, y, f, cnt,failed);

  y := DirichletBeta(-100 + 1/16);
  f := 1.114115509708245184e138;
  testrel(25, NE, y, f, cnt,failed);

  y := DirichletBeta(-100 - 1/16);
  f := 1.873640276815095868e138;
  testrel(26, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_DirichletLambda;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','DirichletLambda');

  y := DirichletLambda(40);
  f := 1.0;
  testrel( 1, NE, y, f, cnt,failed);

  y := DirichletLambda(20);
  f := 1.000000000286807697;
  testrel( 2, NE, y, f, cnt,failed);

  y := DirichletLambda(10);
  f := 1.000017041363044825;
  testrel( 3, NE, y, f, cnt,failed);

  y := DirichletLambda(2);
  f := 1.233700550136169827;
  testrel( 4, NE, y, f, cnt,failed);

  y := DirichletLambda(1.125);
  f := 4.649432303012021381;
  testrel( 5, NE, y, f, cnt,failed);

  y := DirichletLambda(0.96875);
  f := -15.36847262070367157;
  testrel( 6, NE, y, f, cnt,failed);

  y := DirichletLambda(0.5);
  f := -0.4277279326939782213;
  testrel( 7, NE, y, f, cnt,failed);

  y := DirichletLambda(1e-20);
  f := -3.465735902799726547e-21;
  testrel( 8, NE, y, f, cnt,failed);

  y := DirichletLambda(0);
  f := 0;
  testabs( 9, 0, y, f, cnt,failed);

  y := DirichletLambda(-1e-10);
  f := 3.465735902282880147e-11;
  testrel(10, NE, y, f, cnt,failed);

  y := DirichletLambda(-5.5);
  f := 0.1182249311977602014;
  testrel(11, NE, y, f, cnt,failed);

  y := DirichletLambda(-9);
  f := 511/132;
  testrel(12, NE, y, f, cnt,failed);

  y := DirichletLambda(-9.5);
  f := 4.824496622405952553;
  testrel(13, NE, y, f, cnt,failed);

  y := DirichletLambda(-215);
  f := -0.1914521579845423869e303;
  testrel(14, NE, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_LegendreChi;
var
  s,x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','LegendreChi');

  {lchi := (s,x) -> (polylog(s,x)-polylog(s,-x))/2;}

  y := LegendreChi(27.6,1);
  f := 1.000000000000067835;
  testrel( 1, NE, y, f, cnt,failed);

  y := LegendreChi(27.6,0.25);
  f := 0.2500000000000010599;
  testrel( 2, NE, y, f, cnt,failed);

  y := LegendreChi(22.8,1);
  f := 1.000000000013232417;
  testrel( 3, NE, y, f, cnt,failed);

  y := LegendreChi(22.8,0.25);
  f := 0.2500000000002067548;
  testrel( 4, NE, y, f, cnt,failed);

  y := LegendreChi(41,1);
  f := 1;
  testrel( 5, NE, y, f, cnt,failed);

  y := LegendreChi(41,0.9);
  f := 0.9;
  testrel( 6, NE, y, f, cnt,failed);

  y := LegendreChi(41,0.1);
  f := 0.1;
  testrel( 7, NE, y, f, cnt,failed);

  x := 1-ldexp(1,-40);
  y := LegendreChi(0,x);
  f := 549755813887.7499999;
  testrel( 8, NE, y, f, cnt,failed);

  y := LegendreChi(1,x);
  f := 14.20951720147865147;
  testrel( 9, NE, y, f, cnt,failed);

  y := LegendreChi(1.5,x);
  f := 1.688759496312122536;
  testrel(10, NE, y, f, cnt,failed);

  y := LegendreChi(2,x);
  f := 1.233700550122791599;
  testrel(11, NE, y, f, cnt,failed);

  y := LegendreChi(1.9375,x);
  f := 1.261856101682040955;
  testrel(12, NE, y, f, cnt,failed);

  x := 0.9375;
  y := LegendreChi(3,x);
  f := 0.9773160665376878901;
  testrel(13, NE, y, f, cnt,failed);

  y := LegendreChi(2,x);
  f := 1.090626526643523340;
  testrel(14, NE, y, f, cnt,failed);

  y := LegendreChi(1,x);
  f := 1.716993602242573123;
  testrel(15, NE, y, f, cnt,failed);

  y := LegendreChi(0.125,x);
  f := 5.935102407465141301;
  testrel(16, NE, y, f, cnt,failed);

  y := LegendreChi(0,x);
  f := 7.74193548387096774193548387097;
  testrel(17, NE, y, f, cnt,failed);

  s := 1+ldexp(1,-30);
  y := LegendreChi(s,1);
  f := 536870912.6351814228;
  testrel(18, NE, y, f, cnt,failed);

  s := 1+ldexp(1,-30);
  y := LegendreChi(s,0.99609375);
  f := 3.118184789086591631;
  testrel(19, NE, y, f, cnt,failed);

  y := LegendreChi(2,0.5);
  f := 0.5153273666943293542;
  testrel(20, NE, y, f, cnt,failed);

  y := LegendreChi(Pi,0.5);
  f := 0.5041812993238665247;
  testrel(21, NE, y, f, cnt,failed);

  y := LegendreChi(Pi,1);
  f := 1.042956220681483755;
  testrel(22, NE, y, f, cnt,failed);

  y := LegendreChi(0,3/4);
  f := 12/7;
  testrel(23, NE, y, f, cnt,failed);

  y := LegendreChi(2,1-ldexp(1,-50));
  f := 1.233700550136153684;
  testrel(24, NE, y, f, cnt,failed);

  y := LegendreChi(2,0.1);
  f := 0.1001115131643563575;
  testrel(25, NE, y, f, cnt,failed);

  y := LegendreChi(2,0.99);
  f := 1.202075664776857538;
  testrel(26, NE, y, f, cnt,failed);

  y := LegendreChi(2,1);
  f := sqr(Pi)/8;
  testrel(27, NE, y, f, cnt,failed);

  y := LegendreChi(17,1-ldexp(1,-50));
  f := 1.000000007744838568;
  testrel(28, NE, y, f, cnt,failed);

  y := LegendreChi(17,0.125);
  f := 0.125000000015124111;
  testrel(29, NE, y, f, cnt,failed);

  y := LegendreChi(17,0.99);
  f := 0.990000007514784503;
  testrel(30, NE, y, f, cnt,failed);

  y := LegendreChi(17,1);
  f := 1.000000007744839456;
  testrel(31, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_eta;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','eta');

  {eta :=  x -> Zeta(x)*(1-2^(1-x));}

  x := 0.0;
  y := eta(x);
  f := 0.5;
  testrel( 1, NE, y, f, cnt,failed);

  x := 0.5;
  y := eta(x);
  f := 0.6048986434216303702;
  testrel( 2, NE, y, f, cnt,failed);

  x := 1.0;
  y := eta(x);
  f := ln(2);
  testrel( 3, NE, y, f, cnt,failed);

  x := 1.5;
  y := eta(x);
  f := 0.7651470246254079454;
  testrel( 4, NE, y, f, cnt,failed);

  x := 2.0;
  y := eta(x);
  f := 0.8224670334241132182;
  testrel( 5, NE, y, f, cnt,failed);

  x := 5.0;
  y := eta(x);
  f := 0.9721197704469093060;
  testrel( 6, NE, y, f, cnt,failed);

  x := 10.0;
  y := eta(x);
  f := 0.9990395075982715656;
  testrel( 7, NE, y, f, cnt,failed);

  x := 20;
  y := eta(x);
  f := 0.9999990466115815221;
  testrel( 8, NE, y, f, cnt,failed);

  x := 40.0;
  y := eta(x);
  f := 0.9999999999990905054;
  testrel( 9, NE, y, f, cnt,failed);

  x := 64.0;
  y := eta(x);
  f := 1.0;
  testrel(10, NE, y, f, cnt,failed);

  x := -0.5;
  y := eta(x);
  f := 0.3801048126096840168;
  testrel(11, NE, y, f, cnt,failed);

  x := -1.0;
  y := eta(x);
  f := 0.25;
  testrel(12, NE, y, f, cnt,failed);

  x := -1.5;
  y := eta(x);
  f := 0.1186808707198402120;
  testrel(13, NE, y, f, cnt,failed);

  x := -2.0;
  y := eta(x);
  f := 0;
  testrel(14, NE, y, f, cnt,failed);

  x := -7.5;
  y := eta(x);
  f := -1.180249705900827553;
  testrel(15, NE, y, f, cnt,failed);

  x := -19.25;
  y := eta(x);
  f := -40519909.27541321942;
  testrel(16, NE, y, f, cnt,failed);

  x := -55.5;
  y := eta(x);
  f := -0.10923266281825727645e47;
  testrel(17, NE, y, f, cnt,failed);

  x := -0.9e-9;
  y := eta(x);
  f := 0.4999999997967877826;
  testrel(18, NE, y, f, cnt,failed);

  x := 1e-10;
  y := eta(x);
  f := 0.5000000000225791353;
  testrel(19, NE, y, f, cnt,failed);

  x := 1+1e-6;
  y := eta(x);
  f := 0.6931473404288163656;
  testrel(20, NE, y, f, cnt,failed);

  x := 1-1e-6;
  y := eta(x);
  f := 0.6931470206910088807;
  testrel(21, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;




{---------------------------------------------------------------------------}
procedure test_polylog;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','polylog');

  x := -4;
  y := polylog(-10,x);
  f := -9.179867136;
  testrel( 1, NE, y, f, cnt,failed);

  x := -2;
  y := polylog(-10,x);
  f := 12.98673982624599909;
  testrel( 2, NE, y, f, cnt,failed);

  x := -0.9375;
  y := polylog(-10,x);
  f := -5.513617961192165296;
  testrel( 3, NE, y, f, cnt,failed);

  x := -0.375;
  y := polylog(-10,x);
  f := 2.7526396512211967113;
  testrel( 4, NE, y, f, cnt,failed);

  x := -0.125;
  y := polylog(-10,x);
  f := -0.4948902328188429910;
  testrel( 5, NE, y, f, cnt,failed);

  x := 0.125;
  y := polylog(-10,x);
  f := 1154.3834518967005141;
  testrel( 6, NE, y, f, cnt,failed);

  x := 0.375;
  y := polylog(-10,x);
  f := 4489888.134893568;
  testrel( 7, NE, y, f, cnt,failed);

  x := 0.9375;
  y := polylog(-10,x);
  f := 44849120154417799440.0;
  testrel( 8, NE, y, f, cnt,failed);

  x := 4;
  y := polylog(-10,x);
  f := -99851.129401005944216;
  testrel( 9, NE, y, f, cnt,failed);

  x := -4;
  y := polylog(-5,x);
  f := 0.11648;
  testrel(10, NE, y, f, cnt,failed);

  x := -2;
  y := polylog(-5,x);
  f := -0.57613168724279835391e-1;
  testrel(11, NE, y, f, cnt,failed);

  x := -0.9375;
  y := polylog(-5,x);
  f := -0.24779282014042711334;
  testrel(12, NE, y, f, cnt,failed);

  x := -0.375;
  y := polylog(-5,x);
  f := 0.45505630345215321403e-1;
  testrel(13, NE, y, f, cnt,failed);

  x := -0.125;
  y := polylog(-5,x);
  f := 0.78262685792025831654e-1;
  testrel(14, NE, y, f, cnt,failed);

  x := 0.125;
  y := polylog(-5,x);
  f := 1.4851634948023357615;
  testrel(15, NE, y, f, cnt,failed);

  x := 0.375;
  y := polylog(-5,x);
  f := 134.77632;
  testrel(16, NE, y, f, cnt,failed);

  x := 0.9375;
  y := polylog(-5,x);
  f := 1660608240.0;
  testrel(17, NE, y, f, cnt,failed);

  x := 4;
  y := polylog(-5,x);
  f := 16.905349794238683128;
  testrel(18, NE, y, f, cnt,failed);


  x := -4;
  y := polylog(3,x);
  f := -2.9671076939431944597;
  testrel(19, NE, y, f, cnt,failed);

  x := -2;
  y := polylog(3,x);
  f := -1.6682833639665712120;
  testrel(20, NE, y, f, cnt,failed);

  x := -1.25;
  y := polylog(3,x);
  f := -1.1032795677851719720;
  testrel(21, NE, y, f, cnt,failed);

  x := -0.9375;
  y := polylog(3,x);
  f := -0.84988320611888332528;
  testrel(22, NE, y, f, cnt,failed);

  x := -0.625;
  y := polylog(3,x);
  f := -0.58339381700091951204;
  testrel(23, NE, y, f, cnt,failed);

  x := -0.375;
  y := polylog(3,x);
  f := -0.35911489585411991165;
  testrel(24, NE, y, f, cnt,failed);

  x := -0.125;
  y := polylog(3,x);
  f := -0.12311562602883607226;
  testrel(25, NE, y, f, cnt,failed);

  x := 0.125;
  y := polylog(3,x);
  f := 0.12702954097934859904;
  testrel(26, NE, y, f, cnt,failed);

  x := 0.375;
  y := polylog(3,x);
  f := 0.3949165234332217059;
  testrel(27, NE, y, f, cnt,failed);

  x := 0.9375;
  y := polylog(3,x);
  f := 1.104748926956492455;
  testrel(28, NE, y, f, cnt,failed);

  x := -4;
  y := polylog(10,x);
  f := -3.9852813502030905105;
  testrel(29, NE, y, f, cnt,failed);
  x := -2;
  y := polylog(10,x);
  f := -1.9962164930503395001;
  testrel(30, NE, y, f, cnt,failed);

  x := -1.25;
  y := polylog(10,x);
  f := -1.2485051313747192922;
  testrel(31, NE, y, f, cnt,failed);

  x := -0.9375;
  y := polylog(10,x);
  f := -0.93665497525551672216;
  testrel(32, NE, y, f, cnt,failed);

  x := -0.625;
  y := polylog(10,x);
  f := -0.62462252819071471797;
  testrel(33, NE, y, f, cnt,failed);

  x := -0.375;
  y := polylog(10,x);
  f := -0.37486354581717511487;
  testrel(34, NE, y, f, cnt,failed);

  x := -0.125;
  y := polylog(10,x);
  f := -0.12498477405751377867;
  testrel(35, NE, y, f, cnt,failed);

  x := 0.125;
  y := polylog(10,x);
  f := 0.12501529210142635343;
  testrel(36, NE, y, f, cnt,failed);

  x := 0.375;
  y := polylog(10,x);
  f := 0.37513824183158653697;
  testrel(37, NE, y, f, cnt,failed);

  x := 0.9375;
  y := polylog(10,x);
  f := 0.93837308609808966430;
  testrel(38, NE, y, f, cnt,failed);

  x := -4;
  y := polylog(6,x);
  f := -3.8050745875177576463;
  testrel(39, NE, y, f, cnt,failed);

  x := -2;
  y := polylog(6,x);
  f := -1.9458305048460724536;
  testrel(40, NE, y, f, cnt,failed);

  x := -1.25;
  y := polylog(6,x);
  f := -1.2278089280721467032;
  testrel(41, NE, y, f, cnt,failed);

  x := -0.9375;
  y := polylog(6,x);
  f := -0.9247444157530764757;
  testrel(42, NE, y, f, cnt,failed);

  x := -0.625;
  y := polylog(6,x);
  f := -0.619199203891737450347;
  testrel(43, NE, y, f, cnt,failed);

  x := -0.375;
  y := polylog(6,x);
  f := -0.3728706669691755157;
  testrel(44, NE, y, f, cnt,failed);

  x := -0.125;
  y := polylog(6,x);
  f := -0.12475848082937029109;
  testrel(45, NE, y, f, cnt,failed);

  x := 0.125;
  y := polylog(6,x);
  f := 0.12524688145264086435;
  testrel(46, NE, y, f, cnt,failed);

  x := 0.375;
  y := polylog(6,x);
  f := 0.37727497647999170526;
  testrel(47, NE, y, f, cnt,failed);

  x := 0.9375;
  y := polylog(6,x);
  f := 0.95262262228338317537;
  testrel(48, NE, y, f, cnt,failed);


  x := -2;
  y := polylog(-9,x);
  f := 3.4540466392318244170;
  testrel(49, NE, y, f, cnt,failed);

  y := polylog(-8,x);
  f := -2.0256058527663465935;
  testrel(50, NE, y, f, cnt,failed);

  y := polylog(-7,x);
  f := -0.14540466392318244170;
  testrel(51, NE, y, f, cnt,failed);

  y := polylog(-6,x);
  f := 0.40329218106995884774;
  testrel(52, NE, y, f, cnt,failed);

  y := polylog(-5,x);
  f := -0.57613168724279835391e-1;
  testrel(53, NE, y, f, cnt,failed);

  y := polylog(-4,x);
  f := -0.12345679012345679012;
  testrel(54, NE, y, f, cnt,failed);

  y := polylog(-3,x);
  f := 2/27;
  testrel(55, NE, y, f, cnt,failed);

  y := polylog(-2,x);
  f := 2/27;
  testrel(56, NE, y, f, cnt,failed);

  y := polylog(-1,x);
  f := -2/9;
  testrel(57, NE, y, f, cnt,failed);

  y := polylog( 0,x);
  f := -2/3;
  testrel(58, NE, y, f, cnt,failed);

  y := polylog(1,x);
  f := -1.0986122886681096914;
  testrel(59, NE, y, f, cnt,failed);

  y := polylog(2,x);
  f := -1.4367463668836809464;
  testrel(60, NE, y, f, cnt,failed);

  x := -0.125;
  y := polylog(-9,x);
  f := -1.1418293161051686144;
  testrel(61, NE, y, f, cnt,failed);

  y := polylog(-8,x);
  f := -0.44952366987487850701;
  testrel(62, NE, y, f, cnt,failed);

  y := polylog(-7,x);
  f := -0.86590567490610957336e-2;
  testrel(63, NE, y, f, cnt,failed);

  y := polylog(-6,x);
  f := 0.10864046996750344817;
  testrel(64, NE, y, f, cnt,failed);

  y := polylog(-5,x);
  f := 0.78262685792025831654e-1;
  testrel(65, NE, y, f, cnt,failed);

  y := polylog(-4,x);
  f := 0.14225473759081440837e-1;
  testrel(66, NE, y, f, cnt,failed);

  y := polylog(-3,x);
  f := -0.40237768632830361225e-1;
  testrel(67, NE, y, f, cnt,failed);

  y := polylog(-2,x);
  f := -0.76817558299039780521e-1;
  testrel(68, NE, y, f, cnt,failed);

  y := polylog(-1,x);
  f := -0.987654320987654321e-1;
  testrel(69, NE, y, f, cnt,failed);

  y := polylog( 0,x);
  f := -1/9;
  testrel(70, NE, y, f, cnt,failed);

  y := polylog(1,x);
  f := -0.11778303565638345454;
  testrel(71, NE, y, f, cnt,failed);

  y := polylog(2,x);
  f := -0.12129662872272647871;
  testrel(72, NE, y, f, cnt,failed);

  x := 0.75;
  y := polylog(-9,x);
  f := 93461994732.0;
  testrel(73, NE, y, f, cnt,failed);

  y := polylog(-8,x);
  f := 2987482260.0;
  testrel(74, NE, y, f, cnt,failed);

  y := polylog(-7,x);
  f := 107430636.0;
  testrel(75, NE, y, f, cnt,failed);

  y := polylog(-6,x);
  f := 4415124.0;
  testrel(76, NE, y, f, cnt,failed);

  y := polylog(-5,x);
  f := 211692.0;
  testrel(77, NE, y, f, cnt,failed);

  y := polylog(-4,x);
  f := 12180.0;
  testrel(78, NE, y, f, cnt,failed);

  y := polylog(-3,x);
  f := 876.0;
  testrel(79, NE, y, f, cnt,failed);

  y := polylog(-2,x);
  f := 84.0;
  testrel(80, NE, y, f, cnt,failed);

  y := polylog(-1,x);
  f := 12.0;
  testrel(81, NE, y, f, cnt,failed);

  y := polylog( 0,x);
  f := 3.0;
  testrel(82, NE, y, f, cnt,failed);

  y := polylog(1,x);
  f := 1.3862943611198906188;
  testrel(83, NE, y, f, cnt,failed);

  y := polylog(2,x);
  f := 0.97846939293030610374;
  testrel(84, NE, y, f, cnt,failed);

  {larger negative n, x < -4}
  y := polylog(-16, -20);
  f := -225.6984685459220749;
  testrel(85, NE, y, f, cnt,failed);

  y := polylog(-16, -12);
  f := 2208.471352806776057;
  testrel(86, NE, y, f, cnt,failed);

  y := polylog(-16, -10);
  f := 3714.753659565581804;
  testrel(87, NE, y, f, cnt,failed);

  y := polylog(-16, -5);
  f := -20022.19871067531485;
  testrel(88, NE, y, f, cnt,failed);

  y := polylog(42, -0.75);
  f := -0.7499999999998721023;
  testrel(89, NE, y, f, cnt,failed);

  y := polylog(64, 0.75);
  f := 0.75;
  testrel(90, NE, y, f, cnt,failed);

  y := polylog(16,1.0);
  f := 0.50000764112970432594*2; {=1.000015282259408651872}
  testrel(91, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_polylogr;
var
  s,y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','polylogr');

  y := polylogr(5,0.75);
  f := 0.7697354105997573810;
  testrel( 1, NE, y, f, cnt,failed);

  y := polylogr(42.5,-0.75);
  f := -0.7499999999999095627;
  testrel( 2, NE, y, f, cnt,failed);

  y := polylogr(42.5,0.75);
  f := 0.7500000000000904373;
  testrel( 3, NE, y, f, cnt,failed);

  y := polylogr(42.5,1);
  f := 1.000000000000160777;
  testrel( 4, NE, y, f, cnt,failed);

  y := polylogr(42.5,-1);
  f := -0.9999999999998392225;
  testrel( 5, NE, y, f, cnt,failed);

  y := polylogr(40.5,-1);
  f := -0.9999999999993568902;
  testrel( 6, NE, y, f, cnt,failed);

  y := polylogr(40.5, 1);
  f := 1.000000000000643110;
  testrel( 7, NE, y, f, cnt,failed);

  y := polylogr(40.5,-0.75);
  f := -0.7499999999996382507;
  testrel( 8, NE, y, f, cnt,failed);

  y := polylogr(40.5, 0.75);
  f := 0.7500000000003617493;
  testrel( 9, NE, y, f, cnt,failed);

  y := polylogr(40.25,-1);
  f := -0.9999999999992352092;
  testrel(10, NE, y, f, cnt,failed);

  y := polylogr(40.25, 1);
  f := 1.000000000000764791;
  testrel(11, NE, y, f, cnt,failed);

  y := polylogr(40.25,-0.75);
  f := -0.7499999999995698052;
  testrel(12, NE, y, f, cnt,failed);

  y := polylogr(40.25, 0.75);
  f := 0.7500000000004301949;
  testrel(13, NE, y, f, cnt,failed);

  y := polylogr(40,-0.75);
  f := -0.7499999999994884093;
  testrel(14, NE, y, f, cnt,failed);

  y := polylogr(40, 0.75);
  f := 0.75000000000051159080;
  testrel(15, NE, y, f, cnt,failed);

  y := polylogr(0,0.75);
  f := 3;
  testrel(16, NE, y, f, cnt,failed);

  y := polylogr(2.5,0.5);
  f := 0.5549972787175122932;
  testrel(17, NE, y, f, cnt,failed);

  y := polylogr(2.5,-0.5);
  f := -0.4622977821900634382;
  testrel(18, NE, y, f, cnt,failed);

  y := polylogr(2.5,0.50390625);
  f := 0.5598843647061405067;
  testrel(19, NE, y, f, cnt,failed);

  y := polylogr(2.5,-0.50390625);
  f := -0.4656545654919774045;
  testrel(20, NE, y, f, cnt,failed);

  y := polylogr(2.5,0.9921875);
  f := 1.322594574342688176;
  testrel(21, NE, y, f, cnt,failed);

  y := polylogr(2.5,-0.9921875);
  f := -0.8612172798686255019;
  testrel(22, NE, y, f, cnt,failed);

  y := polylogr(2.5,1);
  f := 1.341487257250917180;
  testrel(23, NE, y, f, cnt,failed);

  y := polylogr(2.5,-1);
  f := -0.8671998890121841382;
  testrel(24, NE, y, f, cnt,failed);

  y := polylogr(0.5,-1);
  f := -0.6048986434216303702;
  testrel(25, NE, y, f, cnt,failed);

  y := polylogr(0.5,0.9990234375);
  f := 55.244520579175209022;
  testrel(26, NE, y, f, cnt,failed);

  y := polylogr(0.5,1-ldexp(1,-40));
  f := 0.1858551108812170978e7;
  testrel(27, NE, y, f, cnt,failed);

  y := polylogr(0.5,1.0);
 {f := PosInf_x;
  testabs(28, 0, y, f, cnt,failed);}
  f := -1.460354508809586813;
  testrel(28, NE, y, f, cnt,failed);    {change: polylog(s,1)=zeta(s)}

  s := 1+ldexp(1,-40);
  y := polylogr(s,1);
  f := 0.1099511627776577216e13;
  testrel(29, NE, y, f, cnt,failed);

  y := polylogr(s,-1);
  f := -0.6931471805600907093;
  testrel(30, NE, y, f, cnt,failed);

  y := polylogr(s,0.75);
  f := 1.386294361119239802;
  testrel(31, NE, y, f, cnt,failed);

  y := polylogr(s,-0.75);
  f := -0.5596157879355182278;
  testrel(32, NE, y, f, cnt,failed);

  y := polylogr(s,0.25);
  f := 0.2876820724517544193;
  testrel(33, NE, y, f, cnt,failed);

  y := polylogr(s,-0.25);
  f := -0.2231435513142252513;
  testrel(34, NE, y, f, cnt,failed);

  y := polylogr(-0.9921875, -0.9375);
  f := -0.2517042218174219566;
  testrel(35, NE, y, f, cnt,failed);

  y := polylogr(-0.9921875, -0.5);
  f := -0.2232513428133884107;
  testrel(36, NE, y, f, cnt,failed);

  y := polylogr(-0.9921875, -0.375);
  f := -0.1990593704344978727;
  testrel(37, NE, y, f, cnt,failed);

  y := polylogr(-0.375, -0.9375);
  f := -0.4014043247268414126;
  testrel(38, NE, y, f, cnt,failed);

  y := polylogr(-0.375, -0.5);
  f := -0.2965680466460651464;
  testrel(39, NE, y, f, cnt,failed);

  y := polylogr(-0.375, -0.125);
  f := -0.1073243322124273833;
  testrel(40, NE, y, f, cnt,failed);

  y := polylogr(-0.75, -1);
  f := -0.3158761453565543129;
  testrel(41, NE, y, f, cnt,failed);

  y := polylogr(0.875, -1);
  f := -0.6726499691730939986;
  testrel(42, NE, y, f, cnt,failed);

  y := polylogr(-0.75, 1);
  f := -0.1336427744365845624; {mpmath 0.17}
  testrel(43, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;

{---------------------------------------------------------------------------}
procedure test_fermi_dirac;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','fermi_dirac');

  {Maple: fd := (n,x) -> -polylog(n+1,-exp(x));}

  y := fermi_dirac(0,1);
  f := 1.313261687518222834;
  testrel(1, NE, y, f, cnt,failed);

  y := fermi_dirac(0,-1);
  f := 0.3132616875182228340;
  testrel(2, NE, y, f, cnt,failed);

  y := fermi_dirac(1,1);
  f := 1.806286070444774257;
  testrel(3, NE, y, f, cnt,failed);

  y := fermi_dirac(-1,1);
  f := 0.7310585786300048793;
  testrel(4, NE, y, f, cnt,failed);

  y := fermi_dirac(1,-1);
  f := 0.3386479964034521798;
  testrel(5, NE, y, f, cnt,failed);

  y := fermi_dirac(-1,-1);
  f := 0.2689414213699951207;
  testrel(6, NE, y, f, cnt,failed);

  y := fermi_dirac(-2,1);
  f := 0.1966119332414818525;
  testrel(7, NE, y, f, cnt,failed);

  y := fermi_dirac(-3,1);
  f := -0.9085774767294840944e-1;
  testrel(8, NE, y, f, cnt,failed);

  y := fermi_dirac(-2,-1);
  f := 0.1966119332414818525;
  testrel(9, NE, y, f, cnt,failed);

  y := fermi_dirac(-3,-1);
  f := 0.9085774767294840944e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := fermi_dirac(-4,1e-5);
  f := -0.1249999999875000000;
  testrel(11, NE, y, f, cnt,failed);

  y := fermi_dirac(-4,-1e-5);
  f := -0.1249999999875000000;
  testrel(12, NE, y, f, cnt,failed);

  y := fermi_dirac(4,1e-5);
  f := 0.9721292408202815493;
  testrel(13, NE, y, f, cnt,failed);

  y := fermi_dirac(4,-1e-5);
  f := 0.9721103001636913303;
  testrel(14, NE, y, f, cnt,failed);

  y := fermi_dirac(8,10);
  f := 7946.319776573429079;
  testrel(15, NE, y, f, cnt,failed);

  y := fermi_dirac(8,-10);
  f := 0.4539992573679893686e-4;
  testrel(16, NE, y, f, cnt,failed);

  y := fermi_dirac(9,10);
  f := 10388.99016731291248;
  testrel(17, NE, y, f, cnt,failed);

  y := fermi_dirac(9,-10);
  f := 0.4539992774964110184e-4;
  testrel(18, NE, y, f, cnt,failed);

  y := fermi_dirac(-3,10);
  f := -0.4539168599011319565e-4;
  testrel(19, NE, y, f, cnt,failed);

  y := fermi_dirac(-3,-10);
  f := 0.4539168599011319565e-4;
  testrel(20, NE, y, f, cnt,failed);

  y := fermi_dirac(8,10000);
  f := 0.2755735186158236597e31;
  testrel(21, NE, y, f, cnt,failed);

  y := fermi_dirac(8,15000);
  f := 0.1059396483983069418e33;
  testrel(22, NE, y, f, cnt,failed);

  y := fermi_dirac(5,20000);
  f := 0.8888889985511638002e23;
  testrel(23, NE, y, f, cnt,failed);

  y := fermi_dirac(1,12000);
  f := 0.7200000164493406684823e8;
  testrel(24, NE, y, f, cnt,failed);

  y := fermi_dirac(0,12000);
  f := 12000;
  testrel(25, NE, y, f, cnt,failed);

  y := fermi_dirac(-1,12000);
  f := 1;
  testrel(26, NE, y, f, cnt,failed);

  y := fermi_dirac(-2,12000);
  f := 0;
  testabs(27, 0, y, f, cnt,failed);

  y := fermi_dirac(100,200);
  f := 0.4216103500655924161e73;
  testrel(28, NE, y, f, cnt,failed);

  y := fermi_dirac(200,100);
  f := 0.2688117141816135437e44;
  testrel(29, NE, y, f, cnt,failed);

  y := fermi_dirac(-100,200);
  f := 0.1383896526736737531e-86;
  testrel(30, NE, y, f, cnt,failed);

  y := fermi_dirac(-200,100);
  f := -0.1111917984507457867e-26;
  testrel(31, NE, y, f, cnt,failed);

  y := fermi_dirac(500,10);
  f := 22026.46579480671652;
  testrel(32, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_fermi_dirac_half;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','fermi_dirac_m05/p05/p15/p25');

  y := fermi_dirac_m05(-20);
  f := 0.2061153619434517730e-8;
  testrel(1, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(-10);
  f := 0.4539847236080549532e-4;
  testrel(2, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(-2);
  f := 0.1236656218012099427;
  testrel(3, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(-0.5);
  f := 0.4312314419263970869;
  testrel(4, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(0);
  f := 0.6048986434216303702;
  testrel(5, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(1);
  f := 1.027057125474350689;
  testrel(6, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(3.5);
  f := 2.026193699100460058;
  testrel(7, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(4);
  f := 2.185869690623811946;
  testrel(8, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(4.5);
  f := 2.334245766142244249;
  testrel(9, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(10);
  f := 3.552779239536617160;
  testrel(10, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(13.5);
  f := 4.136326785587186739;
  testrel(11, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(14);
  f := 4.212933653105178558;
  testrel(12, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(14.5);
  f := 4.288146102878802782;
  testrel(13, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(20);
  f := 5.041018507535328603;
  testrel(14, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(39.5);
  f := 7.089878704263408613;
  testrel(15, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(40);
  f := 7.134657233550764731;
  testrel(16, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(40.5);
  f := 7.179155892384824069;
  testrel(17, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(100);
  f := 11.28332744292768060;
  testrel(18, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(1000);
  f := 35.68246764915936962; {alpha}
  testrel(19, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(10000);
  f := 112.8379162455239043; {alpha}
  testrel(20, NE, y, f, cnt,failed);

  y := fermi_dirac_m05(1000000);
  f := 1128.379167095048547; {alpha}
  testrel(21, NE, y, f, cnt,failed);

  {---------}
  y := fermi_dirac_p05(-20);
  f := 0.2061153620936537778e-8;
  testrel(22, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(-10);
  f := 0.4539920105264132755e-4;
  testrel(23, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(-2);
  f := 0.1292985133200755911;
  testrel(24, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(-0.5);
  f := 0.5075371035546378361;
  testrel(25, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(0);
  f := 0.7651470246254079454;
  testrel(26, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(1);
  f := 1.575640776151300231;
  testrel(27, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(3.5);
  f := 5.458044388745459911;
  testrel(28, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(4);
  f := 6.511567592754791415;
  testrel(29, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(4.5);
  f := 7.642031275793455330;
  testrel(30, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(9);
  f := 20.62401821156795236;
  testrel(31, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(11.5);
  f := 29.61226879150170305;
  testrel(32, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(12);
  f := 31.54020328704424259;
  testrel(33, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(12.5);
  f := 33.50923983932172004;
  testrel(34, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(20);
  f := 67.49151222165892049;
  testrel(35, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(35.5);
  f := 159.2691126191554513;
  testrel(36, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(36);
  f := 162.6413796464579914;
  testrel(37, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(36.5);
  f := 166.0371774151569984;
  testrel(38, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(100);
  f := 752.3455915521961188;
  testrel(39, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(10000);
  f := 752252.7873442217908;  {alpha}
  testrel(40, NE, y, f, cnt,failed);

  y := fermi_dirac_p05(1000000);
  f := 7.522527780646031039e8; {alpha}
  testrel(41, NE, y, f, cnt,failed);

  {---------}
  y := fermi_dirac_p15(-20);
  f := 0.2061153621687547803e-8;
  testrel(42, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(-10);
  f := 0.4539956540456176333e-4;
  testrel(43, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(-2);
  f := 0.1322467822517723668;
  testrel(44, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(-0.5);
  f := 0.5526495259473540854;
  testrel(45, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(0);
  f := 0.8671998890121841382;
  testrel(46, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(1);
  f := 2.002258148778464457;
  testrel(47, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(3.5);
  f := 10.27141115061732259;
  testrel(48, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(4);
  f := 13.26048817729106850;
  testrel(49, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(4.5);
  f := 16.79579730928991532;
  testrel(50, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(9);
  f := 78.66626344736753771;
  testrel(51, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(11.5);
  f := 141.2287997372907052;
  testrel(52, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(12);
  f := 156.5151864279572804;
  testrel(53, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(12.5);
  f := 172.7758530334651074;
  testrel(54, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(20);
  f := 546.5630100657601959;
  testrel(55, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(35.5);
  f := 2270.464576509283294;
  testrel(56, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(36);
  f := 2350.941215700957481;
  testrel(57, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(36.5);
  f := 2433.109877926863218;
  testrel(58, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(100);
  f := 30108.67168135486936;
  testrel(59, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(10000);
  f := 3.009011297865632890e9; {alpha}
  testrel(60, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(20000);
  f := 1.702153755962129317e10; {alpha}
  testrel(61, NE, y, f, cnt,failed);

  y := fermi_dirac_p15(1000000);
  f := 3.00901111227325e14;  {GSL/double}
  testrele(62, 5e-15, y, f, cnt,failed);

  {---------}
  y := fermi_dirac_p25(-20);
  f := 0.2061153622063052815e-8;
  testrel(63, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(-10);
  f := 0.4539974758252285430e-4;
  testrel(64, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(-2);
  f := 0.1337669290459732776;
  testrel(65, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(-0.5);
  f := 0.5779521605410086652;
  testrel(66, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(0);
  f := 0.9275535777739480351;
  testrel(67, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(1);
  f := 2.294832963121518129;
  testrel(68, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(3.5);
  f := 15.60768391014749828;
  testrel(69, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(4);
  f := 21.46870822801142616;
  testrel(70, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(4.5);
  f := 28.959226461407335057;
  testrel(71, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(9);
  f := 221.7903608693186891;
  testrel(72, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(11.5);
  f := 491.9765499165918583;
  testrel(73, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(12);
  f := 566.3723808356910136;
  testrel(74, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(12.5);
  f := 648.654118806843732657;
  testrel(75, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(20);
  f := 3186.735096831265885;
  testrel(76, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(35.5);
  f := 23178.76346773165267;
  testrel(77, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(36);
  f := 24334.04466016358338;
  testrel(78, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(36.5);
  f := 25529.98668772732444;
  testrel(79, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(100);
  f := 860954.9737352777100;
  testrel(80, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(1000);
  f := 0.2718704450106142809e10;
  testrel(81, NE, y, f, cnt,failed);

  y := fermi_dirac_p25(5000);
  f := 0.7598904953966710865e12;
  testrel(82, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;

{$endif}


end.
