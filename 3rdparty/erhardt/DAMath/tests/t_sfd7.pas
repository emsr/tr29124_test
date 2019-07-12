{Part 7 of regression test for SPECFUND unit  (c) 2013+  W.Ehrhardt}

unit t_sfd7;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}

interface

procedure test_agm;
procedure test_bernpoly;
procedure test_lambertw;
procedure test_lambertw1;
procedure test_RiemannR;
procedure test_RiemannR_inv;
procedure test_cosint;
procedure test_sinint;
procedure test_fibfun;
procedure test_fibpoly;
procedure test_lucpoly;
procedure test_catalan;
procedure test_bring;
procedure test_omega;
procedure test_MarcumQ;
procedure test_expreln;


implementation

uses
  DAMath, SpecFunD, t_sfd0;


{---------------------------------------------------------------------------}
procedure test_agm;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','agm');

{$ifdef CPUARM}
  y := agm(MinDouble, Sqrt_MaxDbl);
  f := 1.97815801598132823e151;
  testrel( 1, 4, y, f, cnt,failed);

  y := agm(1, Sqrt_MaxDbl);
  f := 5.91138270923603144e151;
  testrel( 2, 5, y, f, cnt,failed);
{$else}
  y := agm(MinDouble, Sqrt_MaxDbl);
  f := 1.97815801598132823e151;
  testrel( 1, NE, y, f, cnt,failed);

  y := agm(1, Sqrt_MaxDbl);
  f := 5.91138270923603144e151;
  testrel( 2, 2, y, f, cnt,failed);
{$endif}

  y := agm(MinDouble, 1);
  f := 0.0022130664755015594768;
  testrel( 3, NE, y, f, cnt,failed);

  y := agm(0, 1);
  f := 0;
  testabs( 4, 0, y, f, cnt,failed);

  y := agm(1, 0);
  f := 0;
  testabs( 5, 0, y, f, cnt,failed);

  y := agm(1, 1e-100);
  f := 0.0067810557455754508824;
  testrel( 5, NE+1, y, f, cnt,failed);       {+1 for ARM}

  y := agm(1, 1e-10);
  f := 0.064344870476013322929;
  testrel( 7, NE, y, f, cnt,failed);

  y := agm(1, 1e-5);
  f := 0.12177452186538904490;
  testrel( 8, NE, y, f, cnt,failed);

  y := agm(1, 0.125);
  f := 0.45196952219967034359;
  testrel( 9, NE, y, f, cnt,failed);

  y := agm(1, 0.5);
  f := 0.72839551552345343459;
  testrel(10, NE, y, f, cnt,failed);

  y := agm(1, 1);
  f := 1;
  testrel(11, 0, y, f, cnt,failed);

  y := agm(1, 1000);
  f := 189.38830240995087556;
  testrel(12, NE, y, f, cnt,failed);

  y := agm(1, 1e20);
  f := 3311261967046375735.6;
  testrel(13, NE, y, f, cnt,failed);

  y := agm(1, 1e100);
  f := 6.7810557455754508824e97;
  testrel(14, NE+1, y, f, cnt,failed);   {+1 for ARM}

  y := agm(1, 1e150);
  f := 4.52974001126015639814e147;
  testrel(15, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lambertw;
var
  x,y,f: double;
  cnt, failed: integer;
  i: integer;
const
  NE = 2;
  NA = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','LambertW');

  x := -exp(-1);
  y := LambertW(x);
  f := -1;
  testrel( 1, NE, y, f, cnt,failed);

  x := -6027/16384.0;
  y := LambertW(x);
  f := -0.9894660903569092419;
  testrel( 2, 16, y, f, cnt,failed);

  x := -3013/8192.0;
  y := LambertW(x);
  f := -0.9790854111351907084;
  testrel( 3, 8, y, f, cnt,failed);

  x := -753/2048;
  y := LambertW(x);
  f := -0.9670887700916448631;
  testrel( 4, 8, y, f, cnt,failed);

  x := -0.25;
  y := LambertW(x);
  f := -0.3574029561813889031;
  testrel( 5, NE, y, f, cnt,failed);

  x := -0.125;
  y := LambertW(x);
  f := -0.1444213531375097292;
  testrel( 6, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := LambertW(x);
  f := -0.977517573730222695e-3;
  testrel( 7, NE, y, f, cnt,failed);

  x := -1e-10;
  y := LambertW(x);
  f := -0.1000000000100000000e-9;
  testrel( 8, NE, y, f, cnt,failed);

  x := 0.0;
  y := LambertW(x);
  f := 0;
  testrel( 9, 1, y, f, cnt,failed);

  x := 1e-10;
  y := LambertW(x);
  f := 0.9999999999000000000e-10;
  testrel(10, NE, y, f, cnt,failed);

  x := 0.0009765625;
  y := LambertW(x);
  f := 0.9756102202467530500e-3;
  testrel(11, NE, y, f, cnt,failed);

  x := 0.125;
  y := LambertW(x);
  f := 0.1117801089327885068;
  testrel(12, NE, y, f, cnt,failed);

  x := 0.25;
  y := LambertW(x);
  f := 0.2038883547022401644;
  testrel(13, NE, y, f, cnt,failed);

  x := 1;
  y := LambertW(x);
  f := +0.5671432904097838730;
  testrel(14, NE, y, f, cnt,failed);

  x := 3.125;
  y := LambertW(x);
  f := 1.070918030310010008;
  testrel(15, NE, y, f, cnt,failed);

  x := 10;
  y := LambertW(x);
  f := 1.745528002740699383;
  testrel(16, NE, y, f, cnt,failed);

  x := 100;
  y := LambertW(x);
  f := 3.385630140290050185;
  testrel(17, NE, y, f, cnt,failed);

  x := 1e4;
  y := LambertW(x);
  f := 7.231846038093372706;
  testrel(18, NE, y, f, cnt,failed);

  x := 1e10;
  y := LambertW(x);
  f := 20.02868541330495078;
  testrel(19, NE, y, f, cnt,failed);

  x := MaxDouble;  {2^1024}
  y := LambertW(x);
  f := 703.2270331047701870;
  testrel(20, NE, y, f, cnt,failed);

  x := 1e-300;
  y := LambertW(x);
  f := x;
  testrel(21, NE, y, f, cnt,failed);

  x := MaxDouble;
  i := 1000;
  while x>1e-200 do begin
    y := LambertW(x);
    f := ln(x/y);
    if y>=0.5 then testrel(i, NE, y, f, cnt,failed)
    else testabs(i, NA, y, f, cnt,failed);
    x := 0.125*x;
    inc(i);
  end;

  x := -exp(-1);
  i := 10000;
  while x<-MinDouble do begin
    y := LambertW(x);
    f := ln(x/y);
    if y>=0.5 then testrel(i, NE, y, f, cnt,failed)
    else testabs(i, NA, y, f, cnt,failed);
    if x < -1e-7 then x := 0.95*x
    else x := 0.01*x;
    inc(i);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_lambertw1;
var
  x,y,f: double;
  cnt, failed: integer;
  i: integer;
const
  NE = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','LambertW1');

  x := -exp(-1);
  y := LambertW1(x);
  f := -1;
  testrel( 1, NE, y, f, cnt,failed);

  x := -6027/16384.0;
  y := LambertW1(x);
  f := -1.010608408692175856;
  testrel( 2, 16, y, f, cnt,failed);

  x := -3013/8192.0;
  y := LambertW1(x);
  f := -1.021210331605726089434;
  testrel( 3, 8, y, f, cnt,failed);

  x := -753/2048;
  y := LambertW1(x);
  f := -1.033649565301978476;
  testrel( 4, 8, y, f, cnt,failed);

  x := -0.25;
  y := LambertW1(x);
  f := -2.153292364110349649;
  testrel( 5, NE, y, f, cnt,failed);

  x := -0.125;
  y := LambertW1(x);
  f := -3.261685684576488777;
  testrel( 6, NE, y, f, cnt,failed);

  x := -0.0009765625;
  y := LambertW1(x);
  f := -9.144639686625083192;
  testrel( 7, NE, y, f, cnt,failed);

  x := -1e-10;
  y := LambertW1(x);
  f := -26.29523881924692569;
  testrel( 8, NE, y, f, cnt,failed);

  x := -1e-100;
  y := LambertW1(x);
  f := -235.7211588756853137;
  testrel( 9, NE, y, f, cnt,failed);

  x := -1e-299;
  y := LambertW1(x);
  f := -695.0168789367294442;
  testrel(10, NE, y, f, cnt,failed);

  x := -exp(-1);
  i := 1000;
  while x<-MinDouble do begin
    y := LambertW1(x);
    f := ln(x/y);
    testrel(i, 2, y, f, cnt,failed);
    if x < -1e-7 then x := 0.95*x
    else
      x := 0.01*x;
    inc(i);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_RiemannR;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 2;
  NE1 = 6;
  NE2 = 256;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','RiemannR');

  y := RiemannR(ldexpd(1,1023));
  f := 1.269399252238138694e305;
  testrel( 1, NE, y, f, cnt,failed);

  y := RiemannR(1e100);
  f := 4.361971987140703159e97;
  testrel( 2, NE, y, f, cnt,failed);

  y := RiemannR(1e30);
  f := 0.1469239889772043272e29;
  testrel( 3, NE, y, f, cnt,failed);

  y := RiemannR(1e24);
  f := 0.18435599767347541878e23;
  testrel( 4, NE, y, f, cnt,failed);

  y := RiemannR(1e20);
  f := 0.2220819602556027015e19;
  testrel( 5, NE, y, f, cnt,failed);

  y := RiemannR(1e19);
  f := 0.2340576673002289402e18;
  testrel( 6, NE, y, f, cnt,failed);

  y := RiemannR(1e18);
  f := 0.2473995428423949440e17;
  testrel( 7, NE+1, y, f, cnt,failed);  {+1 for ARM}

  y := RiemannR(1e16);
  f := 0.2792383413609771872e15;
  testrel( 8, NE, y, f, cnt,failed);

  y := RiemannR(1e6);
  f := 78527.39942912770486;
  testrel( 9, NE, y, f, cnt,failed);

  y := RiemannR(1000);
  f := 168.3594462811673481;
  testrel(10, NE1, y, f, cnt,failed);

  y := RiemannR(100);
  f := 25.66163326692418259;
  testrel(11, NE, y, f, cnt,failed);

  y := RiemannR(10);
  f := 4.564583141005090240;
  testrel(12, NE1, y, f, cnt,failed);

  y := RiemannR(8);
  f := 3.901186044934149947;
  testrel(13, NE, y, f, cnt,failed);

  y := RiemannR(4);
  f := 2.426657752706807136;
  testrel(14, NE, y, f, cnt,failed);

  y := RiemannR(2);
  f := 1.541009016187131883;
  testrel(15, NE, y, f, cnt,failed);

  y := RiemannR(1);
  f := 1;
  testrel(16, NE, y, f, cnt,failed);

  y := RiemannR(0.5);
  f := 0.6635262381124574212;
  testrel(17, NE, y, f, cnt,failed);

  y := RiemannR(0.1);
  f := 0.2790883528560020161;
  testrel(18, NE, y, f, cnt,failed);

  y := RiemannR(0.125);
  f := 0.3124612249259569001;
  testrel(19, NE, y, f, cnt,failed);

  y := RiemannR(0.0625);
  f := 0.2216077332920197402;
  testrel(20, NE1, y, f, cnt,failed);

  y := RiemannR(1e15);
  f := 2.984457049588692738e13;
  testrel(21, NE1, y, f, cnt,failed);

  {New small argmuments with slightly larger errors}
  y := RiemannR(1/128);
  f := 0.089281671873970447697;
  testrel(22, NE2, y, f, cnt,failed);

  y := RiemannR(1/256);
  f := 0.06839846233073223110;
  testrel(23, NE2, y, f, cnt,failed);

  y := RiemannR(1/512);
  f := 0.05325261152444241778;
  testrel(24, NE2, y, f, cnt,failed);

  f := 0.04208120148226344838;
  y := RiemannR(1/1024);
  testrel(25, 1024, y, f, cnt,failed);


  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_RiemannR_inv;
var
  y,f,a,x: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE1 = 8;
  NE2 = 20;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','RiemannR_inv');

  (*
  InvRi[x_?NumberQ]:=y/.FindRoot[RiemannR[y]==x,{y,x*Log[x]},WorkingPrecision->30,PrecisionGoal->20]
  *)

  y := RiemannR_inv(9/8);
  f := 1.212078769477618119;
  testrel(1, NE2, y, f, cnt,failed);

  y := RiemannR_inv(5/4);
  f := 1.436305732692753600;
  testrel(2, NE1, y, f, cnt,failed);

  y := RiemannR_inv(3/2);
  f := 1.917301426243511654;
  testrel(3, NE1, y, f, cnt,failed);

  y := RiemannR_inv(2);
  f := 2.989364977431284783;
  testrel(4, NE, y, f, cnt,failed);

  y := RiemannR_inv(10);
  f := 29.33556598589806959;
  testrel(5, NE+1, y, f, cnt,failed);   {+1 ARM}

  y := RiemannR_inv(100);
  f := 536.4792680045458659;
  testrel(6, NE, y, f, cnt,failed);

  y := RiemannR_inv(1000);
  f := 7922.569910466747789;
  testrel(7, NE, y, f, cnt,failed);

  y := RiemannR_inv(1e6);
  f := 1.548403976528648830e7;
  testrel(8, NE, y, f, cnt,failed);

  y := RiemannR_inv(1e10);
  f := 2.520977157769366229e11;
  testrel(9, NE1, y, f, cnt,failed);   {SSE}

  y := RiemannR_inv(1e16);
  f := 3.949069137982249754e17;
  testrel(10, NE, y, f, cnt,failed);

  y := RiemannR_inv(1e19);
  f := 4.656754651117253795e20;
  testrel(11, NE, y, f, cnt,failed);

  y := RiemannR_inv(1e36);
  f := 8.633948839782661051e37;
  testrel(12, NE, y, f, cnt,failed);

  y := RiemannR_inv(1e40);
  f := 9.565345409988923395e41;
  testrel(13, NE, y, f, cnt,failed);

  y := RiemannR_inv(1e100);
  f := 2.347125735865764178e102;
  testrel(14, NE, y, f, cnt,failed);

  y := RiemannR_inv(1e200);
  f := 4.656583139411941691e202;
  testrel(15, NE, y, f, cnt,failed);

  y := RiemannR_inv(1e300);
  f := 6.963198968074983689e302;
  testrel(16, NE, y, f, cnt,failed);

  x := 9/8;
  a := 0;
  repeat
    y := RiemannR_inv(x);
    f := RiemannR(y);
    y := abs(1-f/x);
    if y > a then a := y;
    x := x+x/16;
  until x>2e40;
  inc(cnt);
  a := a/eps_d;
  if a > 2*NE2 then begin
    inc(failed);
    writeln(' Test RiemannR(RiemannR_inv) failed: max. rel. error = ', a:1:2, ' eps_d');
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_bernpoly;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE1 = 6;
  NE2 = 12;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bernpoly');

  y := bernpoly(1,-1);
  f := -3/2;
  testrel(1, NE, y, f, cnt,failed);

  y := bernpoly(7,-1);
  f := -7;
  testrel(2, NE, y, f, cnt,failed);

  y := bernpoly(8,-1);
  f := 239/30;
  testrel(3, NE, y, f, cnt,failed);

  y := bernpoly(1,1);
  f := 0.5;
  testrel(4, NE, y, f, cnt,failed);

  y := bernpoly(7,1);
  f := 0;
  testrel(5, NE, y, f, cnt,failed);

  y := bernpoly(8,1);
  f := -0.3333333333333333333e-1;
  testrel(6, NE, y, f, cnt,failed);

  {some 'normal' case}
  y := bernpoly(10,-0.75);
  f := 0.7507730252815015388;
  testrel(7, NE, y, f, cnt,failed);

  y := bernpoly(12,1/3);
  f := 0.1265560621400686431;
  testrel(8, NE, y, f, cnt,failed);

  y := bernpoly(11,-0.125);
  f := -0.9375500085297971964e-1;
  testrel(9, NE, y, f, cnt,failed);

  y := bernpoly(20, 0.375);
  f := 374.14698287953266697;
  testrel(10, NE, y, f, cnt,failed);

  y := bernpoly(15,-1.25);
  f := -343.8455539522692561;
  testrel(11, NE, y, f, cnt,failed);

  y := bernpoly(12,2.5);
  f := 1038.229552462511447;
  testrel(12, NE, y, f, cnt,failed);


  x := 10000+1/3;
  y := bernpoly(15,x);
  f := 0.9997499416835207086e60;
  testrel(13, NE1, y, f, cnt,failed);

  y := bernpoly(10,x);
  f := 0.9998333083377781148e40;
  testrel(14, NE, y, f, cnt,failed);

  x := sqrt(1e21);
  y := bernpoly(15,x);
  f := 0.3162277659418379332e158;
  testrel(15, NE1, y, f, cnt,failed);   {NE1 for ARM}

  x := sqrt(1e21);
  y := bernpoly(10,x);
  f := 0.9999999998418861170e105;
  testrel(16, NE, y, f, cnt,failed);

  y := bernpoly(256,0.5);
  f := 0.7950212504588525285e303;
  testrel(17, NE, y, f, cnt,failed);

  y := bernpoly(100,0.5);
  f := 2.838224957069370696e78;
  testrel(18, NE, y, f, cnt,failed);

  x := -987654321.0/65536.0;  {=15070.4089508056641}
  y := bernpoly(70,x);
  f := 0.2949584500818898822e293;
  testrel(19, NE, y, f, cnt,failed);

  x := 1e-5;
  y := bernpoly(7,x);
  f := 0.1666666665500000000e-5;
  testrel(20, NE, y, f, cnt,failed);

  x := 1e-5;
  y := bernpoly(6,x);
  f := 0.2380952375952380955e-1;
  testrel(21, NE, y, f, cnt,failed);

  x := 0.5e-5;
  y := bernpoly(6,x);
  f := 0.2380952379702380953e-1;
  testrel(22, NE, y, f, cnt,failed);

  x := 1e-11;
  y := bernpoly(2,x);
  f := 0.1666666666566666667;
  testrel(23, NE, y, f, cnt,failed);

  y := bernpoly(3,x);
  f := 0.4999999999850000000e-11;
  testrel(24, NE, y, f, cnt,failed);

  y := bernpoly(4,x);
  f := -0.3333333333333333333e-1;
  testrel(25, NE, y, f, cnt,failed);

  y := bernpoly(5,x);
  f := -0.1666666666666666667e-11;
  testrel(26, NE, y, f, cnt,failed);

  x := 1e-12;
  y := bernpoly(15,x);
  f := 0.175e-10;
  testrel(27, NE, y, f, cnt,failed);

  x := 1e-10;
  y := bernpoly(15,x);
  f := 0.175e-8;
  testrel(28, NE, y, f, cnt,failed);

  x := 1e-10;
  y := bernpoly(20,x);
  f := -529.1242424242424241;
  testrel(29, NE, y, f, cnt,failed);

  x := 2e-10;
  y := bernpoly(51,x);
  f := 0.7650884080998503652e17;
  testrel(30, NE, y, f, cnt,failed);

  x := 1e-5;
  y := bernpoly(51,x);
  f := 0.3825442037982211854e22;
  testrel(31, NE, y, f, cnt,failed);

  x := 1e-5;
  y := bernpoly(101,x);
  f := -0.2866607204753912463e76;
  testrel(32, NE, y, f, cnt,failed);

  x := 3e5;
  y := bernpoly(16,x);
  f := 0.4304557309700593800e88;
  testrel(33, NE, y, f, cnt,failed);

  x := -100;
  y := bernpoly(2,x);
  f := 10100.16666666666667;
  testrel(34, NE, y, f, cnt,failed);

  x := -Pi;
  y := bernpoly(15,x);
  f := -0.1374730009236393778e9;
  testrel(35, NE, y, f, cnt,failed);

  x := sqrt(10);
  y := bernpoly(20,x);
  f := 46168783.47767854148;
  testrel(36, NE+1, y, f, cnt,failed);   {+1 for ARM}

  y := bernpoly(15,x);
  f := 732699.8879814299995;
  testrel(37, NE, y, f, cnt,failed);

  x := 1/4;
  y := bernpoly(68,x);
  f := 0.8896458292761226510e22;
  testrel(38, NE2, y, f, cnt,failed);

  x := 2/3;
  y := bernpoly(70,x);
  f := -0.1606254105135901626e45;
  testrel(39, NE2, y, f, cnt,failed);

  x := -1.75;
  y := bernpoly(75,x);
  f := 0.6794407537645821705e50;
  testrel(40, NE2, y, f, cnt,failed);

  x := 4.75;
  y := bernpoly(120,x);
  f := 0.2450175593271593322e71;
  testrel(41, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_cosint;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 6;
  NE3 = 10;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','cosint');

  {special case}
  y := cosint(10000,0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  f := -1234.5;
  y := cosint(0,f);
  testrel(2, NE, y, f, cnt,failed);

  y := cosint(1,-1234.5);
  f := -0.1453956505229364208;
  testrel(3, NE, y, f, cnt,failed);

  {------------------------------}
  y := cosint(10,4);
  f := 1.158627877632916986;
  testrel(4, NE+1, y, f, cnt,failed);    {+1 for ARM}

  y := cosint(10,5);
  f := 1.159689565748002596;
  testrel(5, NE2, y, f, cnt,failed);

  y := cosint(10,Pi);
  f := 0.7731263170943631798;
  testrel(6, NE, y, f, cnt,failed);

  y := cosint(10,Pi_2);
  f := 0.3865631585471815899;
  testrel(7, NE3, y, f, cnt,failed);

  y := cosint(10,19*Pi_2);
  f := 7.344700012396450208;
  testrel(8, NE2, y, f, cnt,failed);

  y := cosint(10,-19*Pi_2);
  f := -7.344700012396450208;
  testrel(9, NE2, y, f, cnt,failed);

  y := cosint(10,100);
  f := 24.38341351832059336;
  testrel(10, NE2, y, f, cnt,failed);

  y := cosint(5,1);
  f := 0.5286328129112155881;
  testrel(11, NE, y, f, cnt,failed);

  y := cosint(17,Pi_2);
  f := 32768/109395;
  testrel(12, NE, y, f, cnt,failed);

  y := cosint(20,1.25);
  f := 0.2767696820752703606;
  testrel(13, NE, y, f, cnt,failed);

  y := cosint(20,10);
  f := 1.935371243364068614;
  testrel(14, NE, y, f, cnt,failed);

  y := cosint(99,TwoPi);
  f := 0;
  testabs(15, 2, y, f, cnt,failed);

  y := cosint(99,4);
  f := -0.1256451290185490101;
  testrel(16, NE2, y, f, cnt,failed);

  y := cosint(9,30);
  f := -0.4063492055789358166;
  testrel(17, NE, y, f, cnt,failed);

  y := cosint(200,10000);
  f := 563.5558003428517485;
  testrel(18, NE, y, f, cnt,failed);

  y := cosint(600,100000);
  f := 3255.962631924894514;
  testrel(19, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_sinint;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
  NE2 = 8;
  NE3 = 16;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','sinint');

  {special case}
  y := sinint(10000,0);
  f := 0;
  testrel(1, NE, y, f, cnt,failed);

  f := -1234.5;
  y := sinint(0,f);
  testrel(2, NE, y, f, cnt,failed);

  y := sinint(1,-1234.5);
  f := 1.989373592132422007;
  testrel(3, NE, y, f, cnt,failed);

  {------------------------------}
  y := sinint(6,0.5);
  f := 0.9185547761246106100e-3;
  testrel(4, NE, y, f, cnt,failed);

  y := sinint(6,1.5);
  f := 0.4204309461889264874;
  testrel(5, NE, y, f, cnt,failed);

  y := sinint(6,2.5);
  f := 0.9771099714670822183;
  testrel(6, NE, y, f, cnt,failed);

  y := sinint(6,3.0);
  f := 0.9817475437693720085;
  testrel(7, NE, y, f, cnt,failed);

  y := sinint(6,5.0);
  f := 1.737945254534703918;
  testrel(18, NE, y, f, cnt,failed);

  y := sinint(6,8.0);
  f := 2.597326791887639688;
  testrel(9, NE, y, f, cnt,failed);

  y := sinint(6,11.0);
  f := 3.440542590614796164;
  testrel(10, NE, y, f, cnt,failed);

  y := sinint(6,-4.0);
  f := -1.009340246947754459;
  testrel(11, NE, y, f, cnt,failed);

  y := sinint(6,-20.0);
  f := -6.025751555775555279;
  testrel(12, NE, y, f, cnt,failed);

  y := sinint(5,0.5);
  f := 0.2226985853239443664e-2;
  testrel(13, NE, y, f, cnt,failed);

  y := sinint(5,1.5);
  f := 0.4628317450440416392;
  testrel(14, NE, y, f, cnt,failed);

  y := sinint(5,2.5);
  f := 1.057683460168249835;
  testrel(15, NE, y, f, cnt,failed);

  y := sinint(5,3.5);
  f := 1.066340666683688580;
  testrel(16, NE, y, f, cnt,failed);

  y := sinint(5,4.5);
  f := 0.7379679184068533471;
  testrel(17, NE, y, f, cnt,failed);

  y := sinint(5,5.5);
  f := 0.2618445135450898748e-1;
  testrel(18, NE, y, f, cnt,failed);

  y := sinint(5,6.5);
  f := 0.16811939182727524784e-4;
  testrel(19, NE, y, f, cnt,failed);

  y := sinint(10,Pi_2);
  f := 0.3865631585471815899;
  testrel(20, NE, y, f, cnt,failed);

  y := sinint(10,Pi);
  f := 0.7731263170943631798;
  testrel(21, NE+1, y, f, cnt,failed);    {+1 for ARM}

  y := sinint(10,-13*Pi_2);
  f := -5.025321061113360669;
  testrel(22, NE2, y, f, cnt,failed);

  y := sinint(9,Pi_2);
  f := 0.4063492063492063492;
  testrel(23, NE, y, f, cnt,failed);

  y := sinint(9,Pi);
  f := 0.8126984126984126984;
  testrel(24, NE, y, f, cnt,failed);

  y := sinint(9,-TwoPi);
  f := 0;
  testrel(25, NE, y, f, cnt,failed);

  y := sinint(9,-13*Pi_2);
  f := 0.4063492063492063492;
  testrel(26, NE3, y, f, cnt,failed);

  y := sinint(99,0.5);
  f := 0.1341426041012494184e-33;
  testrel(27, NE2, y, f, cnt,failed);

  y := sinint(99,2);
  f := 0.2512885482248272477;
  testrel(28, NE, y, f, cnt,failed);

  y := sinint(99,Pi_2);
  f := 0.1256451290185490101;
  testrel(29, NE, y, f, cnt,failed);

  y := sinint(100,0.5);
  f := 0.6367642770299571293e-34;
  testrel(30, NE, y, f, cnt,failed);

  y := sinint(100,2);
  f := 0.2500354235665634526;
  testrel(31, NE, y, f, cnt,failed);

  y := sinint(100,Pi_2);
  f := 0.1250184817401874538;
  testrel(32, NE2, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_fibfun;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 6;
  NE2 = 12;
  NE3 = 360;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','fibfun');


  y := fibfun(0.4, 1/Pi);
  f := 0.3828881779428055359;
  testrel(1, NE, y, f, cnt,failed);

  y := fibfun(0.4, -1/Pi);
  f := 0.3008780080154496036;
  testrel(2, NE, y, f, cnt,failed);

  y := fibfun(-0.4, 1/Pi);
  f := 0.3008780080154496036;
  testrel(3, NE, y, f, cnt,failed);

  y := fibfun(0.4, 0);
  f := 0.3454915028125262879;
  testrel(4, NE, y, f, cnt,failed);

  y := fibfun(20,1);
  f := 6765.0;
  testrel(5, NE, y, f, cnt,failed);

  y := fibfun(20,0.5);
  f := 68.43055915832519531;
  testrel(6, NE, y, f, cnt,failed);

  y := fibfun(20,Pi);
  f := 0.1387302985715636070e11;
  testrel(7, NE, y, f, cnt,failed);

  y := fibfun(20.5, Pi);
  f := 0.2570404097775231245e11;
  testrel(8, NE, y, f, cnt,failed);

  y := fibfun(20.5, -Pi);
  f := 0.2805010692595192992e-11;
  testrel(9, NE2, y, f, cnt,failed);

  y := fibfun(-Pi, 20.5);
  f := 583.3871056411511285;
  testrel(10, NE, y, f, cnt,failed);

  y := fibfun(-Pi, -20.5);
  f := 646.2795671312729308;
  testrel(11, NE, y, f, cnt,failed);

  y := fibfun(100,1);
  f := 3.542248481792619151e20;
  testrel(12, NE, y, f, cnt,failed);

  y := fibfun(-Pi,0);
  f := 0.9513426809665355331;
  testrel(13, NE, y, f, cnt,failed);

  y := fibfun(1000.5,0.125);
  f := 0.6877821969022931296e27;
  testrel(14, NE3, y, f, cnt,failed);

  y := fibfun(1000.5,0);
  f := 0.5;
  testrel(15, NE, y, f, cnt,failed);

  y := fibfun(1000.5,-0.125);
  f := 0.3620728167662394508e-27;
  testrel(16, NE3, y, f, cnt,failed);

  y := fibfun(1400,1);
  f := 0.1710847690234022724e293;
  testrel(17, NE3, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;



{---------------------------------------------------------------------------}
procedure test_fibpoly;
var
  x,y,f,r: double;
  cnt, failed, k,n: integer;
const
  NE  = 4;
{$ifdef FPC}
  NE1 = 540;   {32-Bit/SSE2 only}
{$else}
  NE1 = 420;
{$endif}

{$ifdef BIT16}
const
  tol=5;
{$else}
const
  tol=2;
{$endif}

{Values directly from Maple, note that the 16-bit compilers are inaccurate}
{for some literals of length > 19, therefore the tolerance is increased !}
const
  fp15: array[-10..10] of double = (
          -450359962737049.599999999999998,
          14073748835532.8000000000000113,
          -439804651110.399999999999636200,
          13743895347.2000000000116415321,
          -429496729.599999999627470970153,
          13421772.8000000119209289550781,
{$ifdef BIT16}  {Fix311}
          -419430.0-0.3999996185302734375,
{$else}
          -419430.3999996185302734375,
{$endif}
          13107.20001220703125,
          -409.599609375,
          12.8125,
          0,
          12.8125,
          409.599609375,
          13107.20001220703125,
{$ifdef BIT16}
          419430.0 + 0.3999996185302734375,
{$else}
          419430.3999996185302734375,
{$endif}
          13421772.8000000119209289550781,
          429496729.599999999627470970153,
          13743895347.2000000000116415321,
          439804651110.399999999999636200,
          14073748835532.8000000000000113,
          450359962737049.599999999999998);

  fm25: array[-10..10] of double = (
          17491671035336584642638.2191390,
          92899200757966079853.9563808855,
          493392625783670135.494300316212,
          2620434634437155.42014947329880,
          13917268549466.5442822966724634,
          73915357907.6294983029365539548,
          392568420.6778049468994140625,
{$ifdef BIT16}
          2084951.0 + 0.88653564453125,
{$else}
          2084951.88653564453125,
{$endif}
          11073.291015625,
          58.8125,
          0,
          58.8125,
          -11073.291015625,
{$ifdef BIT16}
          2084951.0 + 0.88653564453125,
{$else}
          2084951.88653564453125,
{$endif}
          -392568420.6778049468994140625,
          73915357907.6294983029365539548,
          -13917268549466.5442822966724634,
          2620434634437155.42014947329880,
          -493392625783670135.494300316212,
          92899200757966079853.9563808855,
          -17491671035336584642638.2191390);

begin
  {pari: fib(n,x) = ([x,1;1,0]^(n-1))[1,1]}
  cnt := 0;
  failed := 0;
  writeln('Function: ','fibpoly');

  y := fibpoly(10,1.5);
  f := 409.599609375;
  testrel(1, NE, y, f, cnt,failed);

  y := fibpoly(15,-Pi);
  f := 0.2909849191767565688e8; {Fib(15, asd(Pi))}
  testrel(2, NE, y, f, cnt,failed);

  y := fibpoly(123,1.5);
  f := 0.4253529586511730793e37;
  testrel(3, NE, y, f, cnt,failed);

  y := fibpoly(123,-1.5);
  f := 0.4253529586511730793e37;
  testrel(4, NE, y, f, cnt,failed);

  y := fibpoly(-123,1.5);
  f := 0.4253529586511730793e37;
  testrel(5, NE, y, f, cnt,failed);

  y := fibpoly(-123,-1.5);
  f := 0.4253529586511730793e37;
  testrel(6, NE, y, f, cnt,failed);

  y := fibpoly(234,1.5);
  f := 0.1104279415486490206e71;
  testrel(7, NE, y, f, cnt,failed);

  y := fibpoly(234,-1.5);
  f := -0.1104279415486490206e71;
  testrel(8, NE, y, f, cnt,failed);

  y := fibpoly(-234,1.5);
  f := -0.1104279415486490206e71;
  testrel(9, NE, y, f, cnt,failed);

  y := fibpoly(-234,-1.5);
  f := 0.1104279415486490206e71;
  testrel(10, NE, y, f, cnt,failed);

  {Max n for F(n) = fibpoly(n,1) is 23601, 1476 for double}
  y := fibpoly(1476, 1);
  f := 0.1306989223763399318e309;
  testrel(11, NE, y, f, cnt,failed);

  y := fibpoly(32000, 1/32);
  f := 6.875799035044984665e216;
  testrel(12, NE1, y, f, cnt,failed);

  x := 1.5;
  for k:=-10 to 10 do begin
    inc(cnt);
    n := 5*k;
    y := fibpoly(n,x);
    f := fp15[k];
    if f=0.0 then r := y-f
    else r := 1.0-y/f;
    if abs(r) > tol*eps_d then begin
      inc(failed);
      writeln('Test failed: n,x,err= ',n:4, x:8:2, r:30);
    end;
  end;

  x := -2.5;
  for k:=-10 to 10 do begin
    inc(cnt);
    n := 5*k;
    y := fibpoly(n,x);
    f := fm25[k];
    if f=0.0 then r := y-f
    else r := 1.0-y/f;
    if abs(r) > tol*eps_d then begin
      inc(failed);
      writeln('Test failed: n,x,err= ',n:4, x:8:2, r:30);
    end;
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_lucpoly;
var
  y,f: double;
  cnt, failed: integer;
const
  NE  = 4;
begin
  {Pari:  luc(n,x) = trace([x,1;1,0]^n)    }
  {Maple: luc := (n,x) -> fibonacci(n+1,x) + fibonacci(n-1,x); }
  cnt := 0;
  failed := 0;
  writeln('Function: ','lucpoly');

  y := lucpoly(0,0);
  f := 2;
  testrel(1, NE, y, f, cnt,failed);

  y := lucpoly(1,0);
  f := 0;
  testabs(2, 0, y, f, cnt,failed);

  y := lucpoly(2,0);
  f := 2;
  testrel(3, NE, y, f, cnt,failed);

  y := lucpoly(3,0);
  f := 0;
  testabs(4, 0, y, f, cnt,failed);

  y := lucpoly(-1,0);
  f := 0;
  testabs(5, 0, y, f, cnt,failed);

  y := lucpoly(-2,0);
  f := 2;
  testrel(6, NE, y, f, cnt,failed);

  y := lucpoly(-3,0);
  f := 0;
  testabs(7, 0, y, f, cnt,failed);

  y := lucpoly(1,1.5);
  f := 1.5;
  testrel(8, NE, y, f, cnt,failed);

  y := lucpoly(2,1.5);
  f := 4.25;
  testrel(9, NE, y, f, cnt,failed);

  y := lucpoly(9,1.5);
  f := 511.998046875;
  testrel(10, NE, y, f, cnt,failed);

  y := lucpoly(9,-1.5);
  f := -511.998046875;
  testrel(11, NE, y, f, cnt,failed);

  y := lucpoly(-9,1.5);
  f := -511.998046875;
  testrel(12, NE, y, f, cnt,failed);

  y := lucpoly(-9,-1.5);
  f := 511.998046875;
  testrel(13, NE, y, f, cnt,failed);

  y := lucpoly(10,1.5);
  f := 1024.0009765625;
  testrel(14, NE, y, f, cnt,failed);

  y := lucpoly(10,-1.5);
  f := 1024.0009765625;
  testrel(15, NE, y, f, cnt,failed);

  y := lucpoly(-10,1.5);
  f := 1024.0009765625;
  testrel(16, NE, y, f, cnt,failed);

  y := lucpoly(-10,-1.5);
  f := 1024.0009765625;
  testrel(17, NE, y, f, cnt,failed);

  y := lucpoly(125,1.5);
  f := 0.4253529586511730793e38;
  testrel(18, NE, y, f, cnt,failed);

  y := lucpoly(125,-1.5);
  f := -0.4253529586511730793e38;
  testrel(19, NE, y, f, cnt,failed);

  y := lucpoly(-125,1.5);
  f := -0.4253529586511730793e38;
  testrel(20, NE, y, f, cnt,failed);

  y := lucpoly(-125,-1.5);
  f := 0.4253529586511730793e38;
  testrel(21, NE, y, f, cnt,failed);

  y := lucpoly(234,1.5);
  f := 0.2760698538716225515e71;
  testrel(22, NE, y, f, cnt,failed);

  y := lucpoly(234,-1.5);
  f := 0.2760698538716225515e71;
  testrel(23, NE, y, f, cnt,failed);

  y := lucpoly(-234,1.5);
  f := 0.2760698538716225515e71;
  testrel(24, NE, y, f, cnt,failed);

  y := lucpoly(-234,-1.5);
  f := 0.2760698538716225515e71;
  testrel(25, NE, y, f, cnt,failed);

  y := lucpoly(-123,4.5);
  f := -6.407383962038300961e82; {Wolfram Alpha}
  testrel(26, NE, y, f, cnt,failed);

  {Max n for L(n) = lucpoly(n,1) is 23599, 1474 for double}
  y := lucpoly(1474,1);
  f := 1.116302065883468286e308;
  testrel(27, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_catalan;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 2;
  NE1 = 4;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','catalan');

  y := catalan(0);
  f := 1;
  testrel(1, NE, y, f, cnt,failed);

  y := catalan(1e-9);
  f := 1.999999998000000005*0.5;
  testrel(2, NE1, y, f, cnt,failed);  {'large' error for arg 1e-9}

  y := catalan(0.125);
  f := 0.4542281453407636314 * 2;
  testrel(3, NE, y, f, cnt,failed);

  y := catalan(0.5);
  f := 0.8488263631567751241;
  testrel(4, NE, y, f, cnt,failed);

  y := catalan(2.5);
  f := 3.104279270973349025;
  testrel(5, NE+1, y, f, cnt,failed);   {+1 for ARM}

  y := catalan(21);
  f := 24466267020.0;
  testrel(6, NE, y, f, cnt,failed);

  y := catalan(35);
{$ifdef BIT16}
  f := 0.38953568686341265775e18 * 8;
{$else}
  f := 3116285494907301262.0;
{$endif}
  testrel(7, NE+1, y, f, cnt,failed);   {+1 for ARM}

  y := catalan(100);
  f := 0.8965199470901314967e57;
  testrel(8, NE, y, f, cnt,failed);

  y := catalan(500);
  f := 0.5394974869170390609e297;
  testrel(9, NE, y, f, cnt,failed);

  y := catalan(-1);
  f := -0.5;
  testrel(10, NE, y, f, cnt,failed);

  y := catalan(-1.25);
  f := -0.3934468663386987420;
  testrel(11, NE1, y, f, cnt,failed);  {D18/64 on AMD}

  y := catalan(-12.375);
  f := -0.1218624678667657878e-8;
  testrel(12, NE1, y, f, cnt,failed);

  y := catalan(-45.625);
  f := 0.1538947487306520572e-29;
  testrel(13, NE+1, y, f, cnt,failed);  {+1 for ARM}

  y := catalan(-123.875);
  f := 0.4497289298048582343e-78;
  testrel(14, NE, y, f, cnt,failed);

  y := catalan(-456.75);
  f := 0.5916664229531300294e-279;
  testrel(15, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_bring;
var
  y,f: double;
  cnt, failed: integer;
const
  NE = 1;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','bring');

  y := bring(0);
  f := 0;
  testabs(1, 0, y, f, cnt,failed);

  y := bring(-1e-10);
  f := 0.1e-9;
  testrel(2, NE, y, f, cnt,failed);

  y := bring(0.01);
  f := -0.9999999900000005e-2;
  testrel(3, NE, y, f, cnt,failed);

  y := bring(-0.125);
  f := 0.1249695196112396472;
  testrel(4, NE, y, f, cnt,failed);

  y := bring(0.5);
  f := -0.4756527435396047855;
  testrel(5, NE, y, f, cnt,failed);

  y := bring(-1);
  f := 0.7548776662466927600;
  testrel(6, NE+1, y, f, cnt,failed);    {+1 for ARM}

  y := bring(2);
  f := -1;
  testrel(7, NE, y, f, cnt,failed);

  y := bring(5);
  f := -1.299152792190219369;
  testrel(8, NE, y, f, cnt,failed);

  y := bring(-10);
  f := 1.533012798646982508;
  testrel(9, NE, y, f, cnt,failed);

  y := bring(1e10);
  f := -99.99999979999999960;
  testrel(10, NE+1, y, f, cnt,failed);   {+1 for ARM}

  y := bring(-1e25);
  f := 1e5;
  testrel(11, NE, y, f, cnt,failed);

  y := bring(ldexpd(1,1000));
  f := -0.1606938044258990276e61;
  testrel(12, NE+1, y, f, cnt,failed);   {+1 for ARM}

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_omega;
var
  x,y,f: double;
  cnt, failed: integer;
const
  NE  = 1;
  NE2 = 2;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','omega');

  x := -4.25;
  y := omega(x);
  f := 0.1406501160968660630e-1;
  testrel(1, NE2, y, f, cnt,failed);

  x := -2;
  y := omega(x);
  f := 0.1200282389876412295;
  testrel(2, NE2, y, f, cnt,failed);

  x := -1.5;
  y := omega(x);
  f := 0.1853749184489398015;
  testrel(3, NE, y, f, cnt,failed);

  x := -1.125;
  y := omega(x);
  f := 0.2522670460430938558;
  testrel(4, NE, y, f, cnt,failed);

  x := -1;
  y := omega(x);
  f := 0.2784645427610737951;
  testrel(5, NE, y, f, cnt,failed);

  x := -0.5;
  y := omega(x);
  f := 0.4046738485459385206;
  testrel(6, NE, y, f, cnt,failed);

  x := -0.25;
  y := omega(x);
  f := 0.4812884663975864440;
  testrel(7, NE, y, f, cnt,failed);

  x := 0;
  y := omega(x);
  f := 0.567143290409783873;
  testrel(8, NE, y, f, cnt,failed);

  x := 0.1;
  y := omega(x);
  f := 0.6040681900258921248;
  testrel(9, NE, y, f, cnt,failed);

  x := 0.25;
  y := omega(x);
  f := 0.6621950814645123396;
  testrel(10, NE, y, f, cnt,failed);

  x := 0.5;
  y := omega(x);
  f := 0.766248608161750259;
  testrel(11, NE, y, f, cnt,failed);

  x := 1;
  y := omega(x);
  f := 1.0;
  testrel(12, NE, y, f, cnt,failed);

  x := 1.125;
  y := omega(x);
  f := 1.063466316835943815;
  testrel(13, NE, y, f, cnt,failed);

  x := 4;
  y := omega(x);
  f := 2.926271062443500913;
  testrel(14, NE, y, f, cnt,failed);

  x := 4.5;
  y := omega(x);
  f := 3.304664918169325271;
  testrel(15, NE, y, f, cnt,failed);

  x := 5;
  y := omega(x);
  f := 3.693441358960649804;
  testrel(16, NE, y, f, cnt,failed);

  x := -45.0;
  y := omega(x);
  f := 0.2862518580549393644e-19;
  testrel(17, NE, y, f, cnt,failed);

  x := -46.0;
  y := omega(x);
  f := 0.1053061735755381238e-19;
  testrel(18, NE, y, f, cnt,failed);

  x := 1e7;
  y := omega(x);
  f := 9999983.881905960852;
  testrel(19, NE, y, f, cnt,failed);

  x := 2e7;
  y := omega(x);
  f := 19999983.18875800904;
  testrel(20, NE, y, f, cnt,failed);

  x := 2.9e-3;
  y := omega(x);
  f := 0.5681934091514999128;
  testrel(21, NE, y, f, cnt,failed);

  x := -2.9e-3;
  y := omega(x);
  f := 0.5660944109285194683;
  testrel(22, NE, y, f, cnt,failed);

  x := 3e-3;
  y := omega(x);
  f := 0.5682296422435454525;
  testrel(23, NE, y, f, cnt,failed);

  x := -3e-3;
  y := omega(x);
  f := 0.5660582647762504257;
  testrel(24, NE+1, y, f, cnt,failed);   {+1 for ARM}

  f := exp(1.0);
  x := 1.0 + f;
  y := omega(x);
  testrel(25, NE, y, f, cnt,failed);

  {double}
  x := 8e5;
  y := omega(x);
  f := 799986.4076499839318;
  testrel(26, NE, y, f, cnt,failed);

  x := 11e5;
  y := omega(x);
{$ifdef BIT16}
  f := 1000000.0 + 99986.08919190850054;
{$else}
  f := 1099986.08919190850054;
{$endif}
  testrel(27, NE, y, f, cnt,failed);

  x := -37.0;
  y := omega(x);
  f := 0.8533047625744065066e-16;
  testrel(28, NE, y, f, cnt,failed);

  x := -38.0;
  y := omega(x);
  f := 0.3139132792048029530e-16;
  testrel(29, NE, y, f, cnt,failed);

  x := 0.9e-2;
  y := omega(x);
  f := 0.5704063236320519836;
  testrel(30, NE, y, f, cnt,failed);

  x := -0.9e-2;
  y := omega(x);
  f := 0.5638921929704859383;
  testrel(31, NE, y, f, cnt,failed);

  x := 1.1e-2;
  y := omega(x);
  f := 0.5711330624359316365;
  testrel(32, NE, y, f, cnt,failed);

  x := -1.1e-2;
  y := omega(x);
  f := 0.5631713483645817001;
  testrel(33, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{
In[83]:= m:=1
In[84]:= N[MarcumQ[m, 0,1/2],22]
Out[84]= 0.8824969025845954028649
In[85]:= N[MarcumQ[m, 0,3/2],22]
Out[85]= 0.3246524673583497297971
In[86]:= N[MarcumQ[m, 0,5/2],22]
Out[86]= 0.04393693362340741732675
In[87]:= N[MarcumQ[m, 0,10],22]
Out[87]= 1.928749847963917783017E-22

In[88]:= N[MarcumQ[m, 1,1/2],22]
Out[88]= 0.9265273979566479682693
In[89]:= N[MarcumQ[m, 1,3/2],22]
Out[89]= 0.4880399991353009476633
In[90]:= N[MarcumQ[m, 1,5/2],22]
Out[90]= 0.1208541231487328978838
In[91]:= N[MarcumQ[m, 1,10],22]
Out[91]= 3.635319497851740426873E-19

In[92]:= N[MarcumQ[m, 2,1/2],22]
Out[92]= 0.9820693672916649480461
In[93]:= N[MarcumQ[m, 2,3/2],22]
Out[93]= 0.7907677793967700955673
In[94]:= N[MarcumQ[m, 2,5/2],22]
Out[94]= 0.3941039245275404662816
In[95]:= N[MarcumQ[m, 2,10],22]
Out[95]= 1.408336693991234558317E-15

In[96]:= N[MarcumQ[m, 4,1/2],22]
Out[96]= 0.99993782390866670768100
In[97]:= N[MarcumQ[m, 4,3/2],22]
Out[97]= 0.9965615816010135876409
In[98]:= N[MarcumQ[m, 4,5/2],22]
Out[98]= 0.9515004104807076314221
In[99]:= N[MarcumQ[m, 4,10],22]
Out[99]= 1.577105283785697995448E-9


In[117]:= m:=2
In[118]:= N[MarcumQ[m, 0,1/2],22]
Out[118]= 0.9928090154076698282230
In[119]:= N[MarcumQ[m, 0,3/2],22]
Out[119]= 0.6898864931364931758188
In[120]:= N[MarcumQ[m, 0,5/2],22]
Out[120]= 0.1812398511965555964728
In[121]:= N[MarcumQ[m, 0,10],22]
Out[121]= 9.836624224615980693388E-21

In[122]:= N[MarcumQ[m, 1,1/2],22]
Out[122]= 0.9955478351118699124944
In[123]:= N[MarcumQ[m, 1,3/2],22]
Out[123]= 0.7779923705497920894830
In[124]:= N[MarcumQ[m, 1,5/2],22]
Out[124]= 0.2885246636204469823841
In[125]:= N[MarcumQ[m, 1,10],22]
Out[125]= 3.488176713111085295551E-18

In[126]:= N[MarcumQ[m, 2,1/2],22]
Out[126]= 0.9989440246248604872508
In[127]:= N[MarcumQ[m, 2,3/2],22]
Out[127]= 0.9210420030198507456016
In[128]:= N[MarcumQ[m, 2,5/2],22]
Out[128]= 0.5749851966379568028418
In[129]:= N[MarcumQ[m, 2,10],22]
Out[129]= 6.949303115516158763715E-15

In[130]:= N[MarcumQ[m, 4,1/2],22]
Out[130]= 0.9999966863659479068774
In[131]:= N[MarcumQ[m, 4,3/2],22]
Out[131]= 0.9990668346791311089799
In[132]:= N[MarcumQ[m, 4,5/2],22]
Out[132]= 0.9761055534244281993673
In[133]:= N[MarcumQ[m, 4,10],22]
Out[133]= 3.956112991806315661648E-9


In[134]:= m:=10
In[135]:= N[MarcumQ[m, 0,1/2],22]
Out[135]= 0.999999999999999770909
In[136]:= N[MarcumQ[m, 0,3/2],22]
Out[136]= 0.9999996767165933071709
In[137]:= N[MarcumQ[m, 0,5/2],22]
Out[137]= 0.9985150433209065237141
In[138]:= N[MarcumQ[m, 0,10],22]
Out[138]= 1.259608459166090750622E-12

In[139]:= N[MarcumQ[m, 1,1/2],22]
Out[139]= 0.999999999999999860258
In[140]:= N[MarcumQ[m, 1,3/2],22]
Out[140]= 0.9999997937492857668413
In[141]:= N[MarcumQ[m, 1,5/2],22]
Out[141]= 0.9989669143045678057673
In[142]:= N[MarcumQ[m, 1,10],22]
Out[142]= 7.820126791065489343802E-12

In[143]:= N[MarcumQ[m, 2,1/2],22]
Out[143]= 0.9999999999999999682845
In[144]:= N[MarcumQ[m, 2,3/2],22]
Out[144]= 0.9999999465070362442234
In[145]:= N[MarcumQ[m, 2,5/2],22]
Out[145]= 0.9996550641297873841872
In[146]:= N[MarcumQ[m, 2,10],22]
Out[146]= 4.116606905023441216598E-10

In[147]:= N[MarcumQ[m, 4,1/2],22]
Out[147]= 0.9999999999999999999159
In[148]:= N[MarcumQ[m, 4,3/2],22]
Out[148]= 0.9999999997620614508745
In[149]:= N[MarcumQ[m, 4,5/2],22]
Out[149]= 0.9999961437065028322035
In[150]:= N[MarcumQ[m, 4,10],22]
Out[150]= 2.568840381567226608978E-6

In[1]:= m:=-2
In[2]:= N[MarcumQ[m, 0,1/2],22]
Out[2]= 0
In[3]:= N[MarcumQ[m, 0,3/2],22]
Out[3]= 0
In[4]:= N[MarcumQ[m, 0,5/2],22]
Out[4]= 0
In[5]:= N[MarcumQ[m, 0,10],22]
Out[5]= 0

In[6]:= N[MarcumQ[m, 1,1/2],22]
Out[6]= 0.01289149715690371942840
In[7]:= N[MarcumQ[m, 1,3/2],22]
Out[7]= 0.005347220771553924032236
In[8]:= N[MarcumQ[m, 1,5/2],22]
Out[8]= 0.0009129258033471960295218
In[9]:= N[MarcumQ[m, 1,10],22]
Out[9]= 2.210161784811189747775E-22

In[10]:= N[MarcumQ[m, 4,1/2],22]
Out[10]= 0.9824423557649863037565
In[11]:= N[MarcumQ[m, 4,3/2],22]
Out[11]= 0.9351891880933253524691
In[12]:= N[MarcumQ[m, 4,5/2],22]
Out[12]= 0.7609227640223771066394
In[13]:= N[MarcumQ[m, 4,10],22]
Out[13]= 8.615364706315978402765E-11
}


{---------------------------------------------------------------------------}
procedure test_MarcumQ;
var
  a,y,f: double;
  m,cnt, failed: integer;
const
  NE  = 3;
  NE1 = 6;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','MarcumQ');

  {------------------------------------------------}
  m := 1;
  a := 0;
  y := MarcumQ(m,a,1/2);
  f := 0.8824969025845954029;
  testrel(1, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.3246524673583497298;
  testrel(2, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.04393693362340741733;
  testrel(3, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 1.928749847963917783E-22;
  testrel(4, NE, y, f, cnt,failed);

  a := 1;
  y := MarcumQ(m,a,1/2);
  f := 0.9265273979566479683;
  testrel(5, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.4880399991353009477;
  testrel(6, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.1208541231487328979;
  testrel(7, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 3.635319497851740427E-19;
  testrel(8, NE, y, f, cnt,failed);

  a := 2;
  y := MarcumQ(m,a,1/2);
  f := 0.9820693672916649480;
  testrel(9, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.7907677793967700956;
  testrel(10, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.3941039245275404663;
  testrel(11, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 1.408336693991234558E-15;
  testrel(12, NE, y, f, cnt,failed);

  a := 4;
  y := MarcumQ(m,a,1/2);
  f := 0.9999378239086667077;
  testrel(13, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.9965615816010135876;
  testrel(14, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.9515004104807076314;
  testrel(15, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 1.577105283785697995E-9;
  testrel(16, NE, y, f, cnt,failed);

  {------------------------------------------------}
  m := 2;
  a := 0;
  y := MarcumQ(m,a,1/2);
  f := 0.9928090154076698282;
  testrel(17, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.6898864931364931758;
  testrel(18, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.1812398511965555965;
  testrel(19, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 9.836624224615980693E-21;
  testrel(20, NE, y, f, cnt,failed);

  a := 1;
  y := MarcumQ(m,a,1/2);
  f := 0.9955478351118699125;
  testrel(21, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.7779923705497920895;
  testrel(22, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.2885246636204469824;
  testrel(23, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 3.488176713111085296E-18;
  testrel(24, NE, y, f, cnt,failed);

  a := 2;
  y := MarcumQ(m,a,1/2);
  f := 0.9989440246248604873;
  testrel(25, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.9210420030198507456;
  testrel(26, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.5749851966379568028;
  testrel(27, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 6.949303115516158764E-15;
  testrel(28, NE, y, f, cnt,failed);

  a := 4;
  y := MarcumQ(m,a,1/2);
  f := 0.9999966863659479069;
  testrel(29, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.9990668346791311090;
  testrel(30, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.9761055534244281994;
  testrel(31, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 3.956112991806315662E-9;
  testrel(32, NE, y, f, cnt,failed);

  {------------------------------------------------}
  m := 10;
  a := 0;
  y := MarcumQ(m,a,1/2);
  f := 0.9999999999999997709;
  testrel(33, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.9999996767165933072;
  testrel(34, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.9985150433209065237141;
  testrel(35, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 1.259608459166090751E-12;
  testrel(36, NE, y, f, cnt,failed);

  a := 4;
  y := MarcumQ(m,a,1/2);
  f := 1;
  testrel(37, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.9999999997620614509;
  testrel(38, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.9999961437065028322;
  testrel(39, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 2.568840381567226609E-6;
  testrel(40, NE1, y, f, cnt,failed);

  {-------- m negative ------------------------------}

  m := -2;

  y := MarcumQ(m,0,1/2);
  f := 0;
  testrel(41, NE, y, f, cnt,failed);

  {Use MarcumQ[-2, 1, 1e-100], b=0 gives infinity}
  {0.014387677966970686643825755639331763156883029833321197015...}
  y := MarcumQ(m,1,0);
  f := 0.01438767796697068664;
  testrel(42, NE, y, f, cnt,failed);

  a := 1;
  y := MarcumQ(m,a,1/2);
  f := 0.01289149715690371943;
  testrel(43, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.005347220771553924032;
  testrel(44, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.0009129258033471960295;
  testrel(45, NE1, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 2.210161784811189748E-22;
  testrel(46, NE, y, f, cnt,failed);

  a := 4;
  y := MarcumQ(m,a,1/2);
  f := 0.9824423557649863038;
  testrel(47, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,3/2);
  f := 0.9351891880933253525;
  testrel(48, NE1, y, f, cnt,failed);

  y := MarcumQ(m,a,5/2);
  f := 0.7609227640223771066;
  testrel(49, NE, y, f, cnt,failed);

  y := MarcumQ(m,a,10);
  f := 8.615364706315978403E-11;
  testrel(50, NE, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_expreln;
var
  y,f: double;
  cnt, failed: integer;
const
  NE   = 2;
  NE1  = 8;
  NE2  = 16;
  NE3  = 256;
begin
  cnt := 0;
  failed := 0;
  writeln('Function: ','expreln');

  y := expreln(0,1e-6);
  f := 1.0000010000005000002;
  testrel(1, NE, y, f, cnt,failed);

  y := expreln(1,-1e-6);
  f := 0.9999995000001666666;
  testrel(2, NE, y, f, cnt,failed);

  y := expreln(2,-0.5);
  f := 0.8522452777010673888;
  testrel(3, NE, y, f, cnt,failed);

  y := expreln(2,1e-5);
  f := 1.000003333341666683;
  testrel(4, NE, y, f, cnt,failed);

  y := expreln(6,1e-9);
  f := 1.00000000014285714329;
  testrel(5, NE, y, f, cnt,failed);

  y := expreln(50,-500);
  f := 0.9105982792090725795e-1;
  testrel(6, NE+1, y, f, cnt,failed);    {+1 for ARM}

  y := expreln(50,-1000);
  f := 0.4766231609253975959e-1;
  testrel(7, NE, y, f, cnt,failed);

  y := expreln(10,0.5);
  f := 1.0474240197472091943;
  testrel(8, NE, y, f, cnt,failed);

  y := expreln(100,-10000);
  f := 0.9901960878554739736060567e-2;
  testrel(9, NE, y, f, cnt,failed);

  y := expreln(100,-1000);
  f := 0.9098434063237174738539731e-1;
  testrel(10, NE, y, f, cnt,failed);

  y := expreln(1000000, 1000000.0);
  f := 1253.647575121312339;
  testrel(11, NE1, y, f, cnt,failed);

  y := expreln(1000000, 1000001.0);
  f := 1254.648201778172082;
  testrel(12, NE1, y, f, cnt,failed);

  {Tests from GSL}
  y := expreln(2, -10.0);
  f := 0.18000090799859524970;
  testrel(13, NE, y, f, cnt,failed);

  y := expreln(2, -1.0e-8);
  f := 0.9999999966666666750;
  testrel(14, NE, y, f, cnt,failed);

  y := expreln(2, 1.0e-8);
  f := 1.000000003333333342;
  testrel(15, NE, y, f, cnt,failed);

  y := expreln(2, 10.0);
  f := 440.3093158961343303;
  testrel(16, NE, y, f, cnt,failed);

  y := expreln(3, -1000.0);
  f := 0.00299400600000000000;
  testrel(17, NE, y, f, cnt,failed);

  y := expreln(3, -100.0);
  f := 0.02940600000000000000;
  testrel(18, NE, y, f, cnt,failed);

  y := expreln(3, -10.0);
  f := 0.24599972760042142509;
  testrel(19, NE, y, f, cnt,failed);

  y := expreln(3, -3.0);
  f := 0.5444917625849191238;
  testrel(20, NE, y, f, cnt,failed);

  y := expreln(3, -0.001);
  f := 0.9997500499916678570;
  testrel(21, NE, y, f, cnt,failed);

  y := expreln(3, -1.0e-8);
  f := 0.9999999975000000050;
  testrel(22, NE, y, f, cnt,failed);

  y := expreln(3, 1.0e-8);
  f := 1.0000000025000000050;
  testrel(23, NE, y, f, cnt,failed);

  y := expreln(3, 0.001);
  f := 1.0002500500083345240;
  testrel(24, NE, y, f, cnt,failed);

  y := expreln(3, 3.0);
  f := 2.574563760708370609;
  testrel(25, NE, y, f, cnt,failed);

  y := expreln(3, 3.1);
  f := 2.677241706846020625;
  testrel(26, NE, y, f, cnt,failed);

  y := expreln(3, 10.0);
  f := 131.7927947688402991;
  testrel(27, NE, y, f, cnt,failed);

  y := expreln(3, 100.0);
  f := 1.612870285089681269e+38;
  testrel(28, NE, y, f, cnt,failed);

  y := expreln(50, -1000.0);
  f := 0.04766231609253975959;
  testrel(29, NE, y, f, cnt,failed);

  y := expreln(50, -100.0);
  f := 0.3348247572345889317;
  testrel(30, NE, y, f, cnt,failed);

  y := expreln(50, -10.0);
  f := 0.8356287051853286482;
  testrel(31, NE, y, f, cnt,failed);

  y := expreln(50, -3.0);
  f := 0.9443881609152163615;
  testrel(32, NE, y, f, cnt,failed);

  y := expreln(50, -1.0);
  f := 0.980762245565660617;
  testrel(33, NE, y, f, cnt,failed);

  y := expreln(50, -1.0e-8);
  f := 1.0 -1.0e-8/51.0;
  testrel(34, NE, y, f, cnt,failed);

  y := expreln(50, 1.0e-8);
  f := 1.0 +1.0e-8/51.0;
  testrel(35, NE, y, f, cnt,failed);

  y := expreln(50,  1.0);
  f := 1.019992165836667903;
  testrel(36, NE, y, f, cnt,failed);

  y := expreln(50, 3.0);
  f := 1.0 + 0.062420575746036830695;
  testrel(37, NE, y, f, cnt,failed);

  y := expreln(50, 48.0);
  f := 7.499573876877194416;
  testrel(38, NE1, y, f, cnt,failed);

  y := expreln(50, 50.125);
  f := 9.337270531835348419;
  testrel(39, NE1, y, f, cnt,failed);

  y := expreln(50, 100.0);
  f := 8.175664432485807634e+07;
  testrel(40, NE1, y, f, cnt,failed);

  y := expreln(50, 500.0);
  f := 4.80635237066318533e+146;
  testrel(41, NE2, y, f, cnt,failed);

  y := expreln(500, -1000.0);
  f := 0.3334815803127619256;
  testrel(42, NE, y, f, cnt,failed);

  y := expreln(500, -100.0);
  f := 0.8335646217536183909;
  testrel(43, NE, y, f, cnt,failed);

  y := expreln(500, -10.0);
  f := 0.9804297803131823066;
  testrel(44, NE, y, f, cnt,failed);

  y := expreln(500, -1.0);
  f := 0.9980079602383488808;
  testrel(45, NE, y, f, cnt,failed);

  y := expreln(500, -1.0e-8);
  f := 1.0 -1.0e-8/501.0;
  testrel(46, NE, y, f, cnt,failed);

  y := expreln(500, 1.0e-8);
  f := 1.0 +1.0e-8/501.0;
  testrel(47, NE, y, f, cnt,failed);

  y := expreln(500, 1.0);
  f := 1.0019999920160634252;
  testrel(48, NE, y, f, cnt,failed);

  y := expreln(500, 100.0);
  f := 1.2492221464878287204;
  testrel(49, NE, y, f, cnt,failed);

  y := expreln(500, 500.0);
  f := 28.363019877927630858;
  testrel(50, NE2, y, f, cnt,failed);

  y := expreln(500, 1000.0);
  f := 2.403756316033530032e+68;
  testrel(51, NE3, y, f, cnt,failed);

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


end.
