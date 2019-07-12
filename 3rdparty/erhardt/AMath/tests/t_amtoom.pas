{Main test unit AMTools  (c) W.Ehrhardt 2010-2018}
{Root finding, function minimization, numerical integration}

unit t_amtoom;

interface

{$i STD.INC}

{$ifdef BIT16}
{$N+,F+}
{$endif}


procedure test_amtools;
  {-AMTools regression test via single procedure call}

procedure test_amint;
  {-Test numerical integration}

procedure test_amconvacc;
  {-Convergence acceleration methods}

implementation

uses
  amath,specfunx,amtools;

var
  totalfail: integer;
  totalcnt : integer;

var
  n: integer;
  a: extended;


{Brent's example, [28] Table 6.1}

const
  xbtf: array[1..19] of extended = (
           3.0229153,   6.6837536,  11.2387017,  19.6760001,
          29.8282273,  41.9061162,  55.9535958,  71.9856656,
          90.0088685, 110.0265327, 132.0405517, 156.0521144,
         182.0620604, 210.0711010, 240.0800483, 272.0902669,
         306.1051233, 342.1369454, 380.2687097);

  ybtf: array[1..19] of extended = (
          3.6766990169,  1.1118500100,  1.2182217637,  2.1621103109,
          3.0322905193,  3.7583856477,  4.3554103836,  4.8482959563,
          5.2587585400,  5.6036524295,  5.8956037976,  6.1438861542,
          6.3550764593,  6.5333662003,  6.6803639849,  6.7938538365,
          6.8634981053,  6.8539024631,  6.6008470481);


{---------------------------------------------------------------------------}
function btf(x: extended): extended;
var
  i: integer;
  s,t: extended;
begin
  s := 0.0;
  for i:=1 to 20 do begin
    t := (2*i-5)/(x-sqr(i));
    s := s + sqr(t);
  end;
  btf := s;
end;


{---------------------------------------------------------------------------}
function dbtf(x: extended): extended;
var
  i: integer;
  s,t: extended;
begin
  s := 0.0;
  for i:=1 to 20 do begin
    t := (x-sqr(i));
    s := s + sqr(2*i-5)/(t*t*t);
  end;
  dbtf := -2.0*s;
end;


{---------------------------------------------------------------------------}
function f4(t: extended): extended;
const
  F=2.25; s=200.0; P=10000.0; r=0.00949;
begin
  f4 := t*s*(1-F) + F*P*power(1+r,t);
end;


{---------------------------------------------------------------------------}
function f5(t: extended): extended;
const
  F=2.25; s=200.0; P=10000.0; r=0.00949;
begin
  f5 := t*s*(1-F) + F*P*exp(t*ln1p(r));
end;


{---------------------------------------------------------------------------}
function fz1(x: extended): extended;
begin
  fz1 := power(x,n)-a;
end;


{---------------------------------------------------------------------------}
function pown(x: extended): extended;
begin
  pown := power(x,n);
end;


{---------------------------------------------------------------------------}
function fz2(x: extended): extended;
begin
  fz2 := power(x,1.0/n)-a;
end;


{---------------------------------------------------------------------------}
function fz3(x: extended): extended;
begin
  fz3 := (n*x-1.0)/(n*x-x);
end;


{---------------------------------------------------------------------------}
function fz5(x: extended): extended;
begin
  fz5 := sin(x)-0.5*x;
end;


{---------------------------------------------------------------------------}
function fz6(x: extended): extended;
begin
  fz6 := cosh(x)*cos(x)+1.0;
end;


{---------------------------------------------------------------------------}
function fconst1(x: extended): extended;
begin
  fconst1 := 1.0;
end;


{---------------------------------------------------------------------------}
procedure test_amtools;
  {-AMTools regression test via single procedure call}
var
  x,y,t,z,e: extended;
  i,nm,nz,err: integer;
  sqrteps: extended;
begin

  sqrteps   := sqrt(eps_x);
  totalfail := 0;
  totalcnt  := 0;

  writeln('Test program for AMTools unit    (c) W.Ehrhardt 2010-2018');
  writeln('Testing:  zbrenty/zbrent/zeroin/zridders and localmin/mbrent/fmin');

  for i:=1 to 19 do begin
    x := fmin({$ifdef FPC_ProcVar}@{$endif}btf, sqr(i)+1e-5, sqr(i+1)-1e-5, eps_x);
    y := btf(x);
    z := xbtf[i];
    t := ybtf[i];
    if (abs(x-z) > 5e-8) or (abs(y-t) > 5e-11) then begin
      writeln('mbrent btf:', 1-x/z:30, 1-y/t:30);
      inc(totalfail);
    end;
  end;

  mbrent({$ifdef FPC_ProcVar}@{$endif}fz6, 8, 11, eps_x, x, y, i);
  z := 10.2101761228130305454682055947;
  t := -9607.99919672068567552199371256;
  inc(totalcnt);
  if (i<=0) or (abs(x-z) > sqrteps*z) or (abs(y-t) > 10*eps_x*abs(t)) then begin
    writeln('mbrent fz6:',i:5, 1-x/z:30, 1-y/t:30);
    inc(totalfail);
  end;

  mbrent({$ifdef FPC_ProcVar}@{$endif}fz6, 0, 7, eps_x, x, y, i);
  z := 3.92660231204791877823853334363;
  t := -16.9512244081414528745287948990;
  inc(totalcnt);
  if (i<=0) or (abs(x-z) > sqrteps*z) or (abs(y-t) > 10*eps_x*abs(t)) then begin
    writeln('mbrent fz6:',i:5, 1-x/z:30, 1-y/t:30);
    inc(totalfail);
  end;

  mbrent({$ifdef FPC_ProcVar}@{$endif}fz6, 14, 19, eps_x, x, y, i);
  z := 16.4933614313464097807735984324;
  t := -5145537.88082210121125558232885;
  inc(totalcnt);
  if (i<=0) or (abs(x-z) > sqrteps*z) or (abs(y-t) > 10*eps_x*abs(t)) then begin
    writeln('mbrent fz6:',i:5, 1-x/z:30, 1-y/t:30);
    inc(totalfail);
  end;

  mbrent({$ifdef FPC_ProcVar}@{$endif}fz6, 20, 24, eps_x, x, y, i);
  z := 22.7765467385260009788377002440;
  t := -2755393132.85845908456008673107;
  inc(totalcnt);
  if (i<=0) or (abs(x-z) > sqrteps*z) or (abs(y-t) > 10*eps_x*abs(t)) then begin
    writeln('mbrent fz6:',i:5, 1-x/z:30, 1-y/t:30);
    inc(totalfail);
  end;

  x := zeroin({$ifdef FPC_ProcVar}@{$endif}fz5, pi/2, pi, eps_x);
  y := 1.895494267033980947144;
  inc(totalcnt);
  if abs(x-y)>6*eps_x*y then begin
    inc(totalfail);
    writeln('zeroin fz5: ',1-x/y:30);
  end;

  x := zeroin({$ifdef FPC_ProcVar}@{$endif}fz6, 0, 4, eps_x);
  y := 1.875104068711961166445;
  inc(totalcnt);
  if abs(x-y)>6*eps_x*y then begin
    inc(totalfail);
    writeln('zeroin fz6: ',1-x/y:30);
  end;

  x := zeroin({$ifdef FPC_ProcVar}@{$endif}fz6, 2, 7, eps_x);
  y := 4.694091132974174576436;
  inc(totalcnt);
  if abs(x-y)>6*eps_x*y then begin
    inc(totalfail);
    writeln('zeroin fz6: ',1-x/y:30);
  end;

  x := zeroin({$ifdef FPC_ProcVar}@{$endif}fz6, 5, 10, eps_x);
  y := 7.854757438237612564861;
  inc(totalcnt);
  if abs(x-y)>6*eps_x*y then begin
    inc(totalfail);
    writeln('zeroin fz6: ',1-x/y:30);
  end;

  for i:=2 to 50 do begin
    {f(x) = (n*x-1)/(n*x-x)}
    n := i;
    x := zeroin({$ifdef FPC_ProcVar}@{$endif}fz3, 0.01, 100, eps_x);
    y := 1.0/n;
    inc(totalcnt);
    if abs(x-y)>10*eps_x*y then begin
      inc(totalfail);
      writeln('zeroin fz3: ',i:3,1-x/y:30);
    end;
  end;

  for i:=1 to 20 do begin
    {f(x) = x^(1/n) - n^(1/n)}
    n := i;
    a := power(n,1.0/n);
    x := zeroin({$ifdef FPC_ProcVar}@{$endif}fz2, 1, 100, eps_x);
    y := n;
    inc(totalcnt);
    if abs(x-y)>10*eps_x*y then begin
      inc(totalfail);
      writeln('zeroin fz2: ',i:3,1-x/y:30);
    end;
  end;

  for i:=1 to 50 do begin
    {f(x) = x^n - a}
    n := i;
    a := 0.25;
    x := zbrenty({$ifdef FPC_ProcVar}@{$endif}pown, a, 0, 5, eps_x, nz, err);
    (*
    {x := zeroin({$ifdef FPC_ProcVar}@{$endif}fz1, 0, 5, eps_x);}
    *)
    y := power(a,1.0/n);
    inc(totalcnt);
    if abs(x-y)>6*eps_x*y then begin
      inc(totalfail);
      writeln('zbrenty pown/a=0.25: ',i:3,1-x/y:30);
    end;
    a := 1;
    (*
    x := zeroin({$ifdef FPC_ProcVar}@{$endif}fz1, -0.999, 6, eps_x);
    *)
    x := zbrenty({$ifdef FPC_ProcVar}@{$endif}pown, a, -0.999, 6, eps_x, nz, err);
    y := power(a,1.0/n);
    inc(totalcnt);
    if abs(x-y)>6*eps_x*y then begin
      inc(totalfail);
      writeln('zbrenty pown/a=1: ',i:3,1-x/y:30);
    end;
  end;

  for i:=1 to 19 do begin
    {Derivative of Brent's sample function, given to 7 digits after .}
    x := zeroin({$ifdef FPC_ProcVar}@{$endif}dbtf, sqr(i)+1e-5, sqr(i+1)-1e-5, eps_x);
    y := xbtf[i];
    inc(totalcnt);
    if abs(x-y) > 5e-8 then begin
      inc(totalfail);
      writeln('zeroin dbft: ',i:3,1-x/y:30);
    end;
  end;

  {-------------  Ridders root finding --------------------------}
  x := zridders({$ifdef FPC_ProcVar}@{$endif}fz5, pi/2, pi, eps_x);
  y := 1.895494267033980947144;
  inc(totalcnt);
  if abs(x-y)>6*eps_x*y then begin
    inc(totalfail);
    writeln('zridders fz5: ',1-x/y:30);
  end;

  x := zridders({$ifdef FPC_ProcVar}@{$endif}fz6, 0, 4, eps_x);
  y := 1.875104068711961166445;
  inc(totalcnt);
  if abs(x-y)>6*eps_x*y then begin
    inc(totalfail);
    writeln('zridders fz6: ',1-x/y:30);
  end;

  x := zridders({$ifdef FPC_ProcVar}@{$endif}fz6, 2, 7, eps_x);
  y := 4.694091132974174576436;
  inc(totalcnt);
  if abs(x-y)>6*eps_x*y then begin
    inc(totalfail);
    writeln('zridders fz6: ',1-x/y:30);
  end;

  x := zridders({$ifdef FPC_ProcVar}@{$endif}fz6, 5, 10, eps_x);
  y := 7.854757438237612564861;
  inc(totalcnt);
  if abs(x-y)>6*eps_x*y then begin
    inc(totalfail);
    writeln('zridders fz6: ',1-x/y:30);
  end;

  for i:=2 to 50 do begin
    {f(x) = (n*x-1)/(n*x-x)}
    n := i;
    x := zridders({$ifdef FPC_ProcVar}@{$endif}fz3, 0.01, 100, eps_x);
    y := 1.0/n;
    inc(totalcnt);
    if abs(x-y)>10*eps_x*y then begin
      inc(totalfail);
      writeln('zridders fz3: ',i:3,1-x/y:30);
    end;
  end;

  for i:=1 to 20 do begin
    {f(x) = x^(1/n) - n^(1/n)}
    n := i;
    a := power(n,1.0/n);
    x := zridders({$ifdef FPC_ProcVar}@{$endif}fz2, 1, 100, eps_x);
    y := n;
    inc(totalcnt);
    if abs(x-y)>10*eps_x*y then begin
      inc(totalfail);
      writeln('zridders fz2: ',i:3,1-x/y:30);
    end;
  end;

  for i:=1 to 50 do begin
    {f(x) = x^n - a}
    n := i;
    a := 0.25;
    x := zridders({$ifdef FPC_ProcVar}@{$endif}fz1, 0, 5, eps_x);
    y := power(a,1.0/n);
    inc(totalcnt);
    if abs(x-y)>6*eps_x*y then begin
      inc(totalfail);
      writeln('zridders fz1/a=0.25: ',i:3,1-x/y:30);
    end;
    a := 1;
    x := zridders({$ifdef FPC_ProcVar}@{$endif}fz1, -0.999, 6, eps_x);
    y := power(a,1.0/n);
    inc(totalcnt);
    if abs(x-y)>6*eps_x*y then begin
      inc(totalfail);
      writeln('zridders fz1/a=1: ',i:3,1-x/y:30);
    end;
  end;

  for i:=1 to 19 do begin
    {Derivative of Brent's sample function, given to 7 digits after .}
    x := zridders({$ifdef FPC_ProcVar}@{$endif}dbtf, sqr(i)+1e-5, sqr(i+1)-1e-5, eps_x);
    y := xbtf[i];
    inc(totalcnt);
    if abs(x-y) > 5e-8 then begin
      inc(totalfail);
      writeln('zridders dbft: ',i:3,1-x/y:30);
    end;
  end;

  {---------------------------------------------------------}

  t := 1.46163214496836234126;
  localmin({$ifdef FPC_ProcVar}@{$endif}gammax,1,2,sqrt(eps_x), eps_x, x,y,i);
  inc(totalcnt);
  if abs(x-t) > 5e-8 then begin
    inc(totalfail);
    writeln('localmin gammax: ',i:3,1-x/y:30);
  end;

  t := 17.197352215395947599;
  localmin({$ifdef FPC_ProcVar}@{$endif}f4,14,22,sqrt(eps_x), eps_x, x,y,i);
  inc(totalcnt);
  if abs(x-t) > 5e-8 then begin
    inc(totalfail);
    writeln('localmin f4: ',i:3,1-x/y:30);
  end;

  localmin({$ifdef FPC_ProcVar}@{$endif}f5,14,22,sqrt(eps_x), eps_x, x,y,i);
  inc(totalcnt);
  if abs(x-t) > 5e-8 then begin
    inc(totalfail);
    writeln('localmin f5: ',i:3,1-x/y:30);
  end;

  for i:=1 to 19 do begin
    localmin({$ifdef FPC_ProcVar}@{$endif}btf, sqr(i)+1e-5, sqr(i+1)-1e-5, 0.5*sqrt(eps_x), eps_x, x,y,nm);
    z := zbrent({$ifdef FPC_ProcVar}@{$endif}dbtf, sqr(i)+1e-5, sqr(i+1)-1e-5, 0.5*sqrt(eps_x), nz, err);
    inc(totalcnt);
    if abs(x-xbtf[i]) > 5e-8  then begin
      inc(totalfail);
      writeln('localmin btf i=',i:2, ':  x-xbtf[i] = ',x-xbtf[i]:30);
    end;
    inc(totalcnt);
    if abs(y-ybtf[i]) > 5e-11 then begin
      inc(totalfail);
      writeln('localmin btf i=',i:2, ':  y-ybtf[i] = ',y-ybtf[i]:30);
    end;
    inc(totalcnt);
    if abs(x-z) > 5e-10*z then begin
      inc(totalfail);
      writeln('zbrent i=',i:2, ':  x-z = ',x-z:30);
    end;
  end;

  for i:=1 to 19 do begin
    e := 5e-8;
    x := fmin({$ifdef FPC_ProcVar}@{$endif}btf, sqr(i)+1e-5, sqr(i+1)-1e-5, e);
    y := btf(x);
    inc(totalcnt);
    if abs(x-xbtf[i]) > e*x   then begin
      inc(totalfail);
      writeln('Fmin Brent i=',i:2, ':  x-xbtf[i] = ',x-xbtf[i]:30);
    end;
    inc(totalcnt);
    if abs(y-ybtf[i]) > 5e-11 then begin
      inc(totalfail);
      writeln('Localmin Brent i=',i:2, ':  y-ybtf[i] = ',y-ybtf[i]:30);
    end;
  end;
  if totalfail=0 then writeln('All ',totalcnt, ' tests passed.')
  else writeln('*** ',totalfail,' test of ',totalcnt, ' failed.')
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}
var
  alpha: extended;   {parameter in integrands}


{---------------------------------------------------------------------------}
function qt1(x: extended): extended;
begin
  qt1 := sinc(x);
end;


{---------------------------------------------------------------------------}
function qt2(x: extended): extended;
begin
  if abs(x)<1e-10 then qt2 := 1.0
  else qt2 := tan(x)/x;
end;


{---------------------------------------------------------------------------}
function qt3(x: extended): extended;
var
  a,b: extended;
begin
  a := sqr(x-0.3)+0.01;
  b := sqr(x-0.9)+0.04;
  qt3 := 1/a + 1/b - 6;
end;


{---------------------------------------------------------------------------}
function qt4(x: extended): extended;
begin
  qt4 := sqrt(x);
end;


{---------------------------------------------------------------------------}
function qt5(x: extended): extended;
begin
  qt5 := 4.0/(1.0+sqr(x));
end;


{---------------------------------------------------------------------------}
function qt6(x: extended): extended;
begin
  qt6 := sin(10.0*x);
end;


{---------------------------------------------------------------------------}
function qt7(x: extended): extended;
begin
  qt7 := 3/(5-4*cos(x));
end;


{---------------------------------------------------------------------------}
function qt8(x: extended): extended;
begin
  {int(q10,0,1) = arctan(alpha)/alpha}
  qt8 := 1/(1+sqr(alpha*x));
end;


{---------------------------------------------------------------------------}
function qt9(x: extended): extended;
begin
  {int(q9,0,1) = ln1p(alpha)/alpha;}
  qt9 := 1/(1+alpha*x);
end;


{---------------------------------------------------------------------------}
function qt10(x: extended): extended;
begin
  {int(q10,0,1) = exp2m1(alpha + 1)/(alpha + 1)}
  qt10 := power(1+x,alpha);
end;


{---------------------------------------------------------------------------}
function qt11(x: extended): extended;
begin
  {int(q11,0,1) = 2-alpha }
  if x<= alpha then qt11 := 1
  else qt11 := 2;
end;


{---------------------------------------------------------------------------}
function qt12(x: extended): extended;
begin
  {int(q12,0,inf) = -1/4*Pi/sqrt(alpha)*ln(alpha)}
  qt12 := ln(x)/(1+alpha*x*x);
end;


{---------------------------------------------------------------------------}
function qt13(x: extended): extended;
begin
  {int(q13,-inf,inf) = sqrt(Pi/alpha)*exp(0.25/alpha)}
  qt13 := exp(-x - alpha*x*x);
end;


{---------------------------------------------------------------------------}
function qt14(x: extended): extended;
begin
  {int(q14,-inf,inf) = sqrtPi*exp(0.25*sqr(alpha))}
  qt14 := exp(-alpha*x - x*x);
end;


{---------------------------------------------------------------------------}
function qt15(x: extended): extended;
begin
  {int(q15,0,inf) = Gamma(alpha)}
  qt15 := exp(-x)*power(x,alpha-1);
end;


{---------------------------------------------------------------------------}
function qt16(x: extended): extended;
begin
  {int(q16,0,1) = -1/(sqr(alpha+1))}
  qt16 := ln(x)*power(x,alpha);
end;


{---------------------------------------------------------------------------}
function qt17(x: extended): extended;
begin
  {int(q17,0,1) = -2*(alpha-7)/sqrt(7*abs(alpha - 7)) + 2*sqrt(alpha/7)}
  qt17 := 1/sqrt(abs(x-alpha/7));
end;


{---------------------------------------------------------------------------}
function qt18(x: extended): extended;
begin
  {int(q18,0,1) = Pi/2 - Pi*ln(2)}
  qt18 := ln(x)*sqrt(x/(1-x));
end;


{---------------------------------------------------------------------------}
function qt19(x: extended): extended;
begin
  {int(q19,0,int) = Pi}
  qt19 := 1/((x+1)*sqrt(x));
end;


{---------------------------------------------------------------------------}
function qt20(x: extended): extended;
begin
  {int(q20,0,2Pi) = 2*Pi^3*J_1(60*Pi) = -2.54325961889353148987086198257}
  qt20 := x*sin(x*30.0)/sqrt(1.0-x*x/(4*pi*pi));
end;


{---------------------------------------------------------------------------}
function qt21(x: extended): extended;
begin
  {integral(1/((x^2+1/100^2)*(x-c)),x=-1..1)
    = -10000*(200*c*arctan(100)-ln(1-c)+ln(c+1))/(10000*c*c+1) for c<1
    = 10000*(-200*c*arctan(100)+ln(c-1)-ln(c+1))/(10000*c*c+1) for c>1 }
  qt21 := 1.0/(x*x + 1e-4);
end;


{---------------------------------------------------------------------------}
function qt22(x: extended): extended;
begin
  {for 0<c<2: integral(1/((x^2+1/100^2)*(x-c)),x=0..2)
   = 8/3 - ln(c) + 2*c - c^3*ln(c) + 2*c^2 + ln(2-c) + ln(2-c)*c^3}
  qt22 := x*x*x + 1.0;
end;


{---------------------------------------------------------------------------}
function qt23(x: extended): extended;
begin
  qt23 := 1.0/(5.0*x*x*x + 6.0); {GSL example with c=0}
end;


{---------------------------------------------------------------------------}
function qt24(x: extended): extended;
begin
  qt24 := exp(x);
end;


{---------------------------------------------------------------------------}
function qt25(x: extended): extended;
begin
  qt25 := power(x,alpha);
end;


{---------------------------------------------------------------------------}
function qt26(x: extended): extended;
begin
  {int(qt26,0,inf)=Pi/sin(Pi*alpha)}
  qt26 := power(x,-alpha)/(1+x);
end;


{---------------------------------------------------------------------------}
function qt27(x: extended): extended;
var
  y: extended;
begin
  {int(qt27,0,inf)=Pi_2/alpha}
  y := power(x,alpha);
  y := x*(y + 1/y);
  qt27 := 1/y;
end;


{---------------------------------------------------------------------------}
function qt28(x: extended): extended;
var
  y: extended;
begin
  {int(qt28,0,inf)=sqrt(Pi_2/alpha)}
  y := sin(x*alpha);
  qt28 := y/sqrt(x);
end;


{---------------------------------------------------------------------------}
function qt29(x: extended): extended;
var
  y: extended;
begin
  {int(qt29,0,inf)=Pi_2*exp(-alpha)}
  y := cos(x*alpha);
  qt29 := y/(1+x*x);
end;


{---------------------------------------------------------------------------}
function qt30(x: extended): extended;
var
  y: extended;
begin
  {int(qt30,0,inf)=sqrt(2*Pi*alpha)}
  y := sin(x*alpha);
  qt30 := y/sqrt(x)/x;
end;



{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
procedure test_quanc8;
  {-Test numerical integration with quanc8}
var
  ans,y,err,est,flag,d: extended;
  nf: longint;
  i,fail, cnt: integer;
const
  epsr = 1e-16;
begin

  fail := 0;
  cnt  := 0;
  writeln('Testing quanc8');

  Quanc8({$ifdef FPC_ProcVar}@{$endif}qt1,0,2,0,epsr,y,est,flag,nf);
  ans := six(2.0);
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quanc8 qt1: ',err:28, est:28);
  end;

  Quanc8({$ifdef FPC_ProcVar}@{$endif}qt2,0,2,epsr,0,y,est,flag,nf);
  inc(cnt);
  {integral(tan(x)/x, x=0..2) does not converge, the flag should be xxx.yyy }
  {with yyy about 215. For then parameters used, BP7 gives flag=19.21460253}
  {This indicates the expected problem at 2-2*0.21460253 ~= Pi/2}
  if (flag=0.0) or (frac(flag)=0.0) then begin
    inc(fail);
    writeln('quanc8 qt2: ',flag:28, est:28);
  end;

  Quanc8({$ifdef FPC_ProcVar}@{$endif}qt3,0,1,0,epsr,y,est,flag,nf);
  ans := 29.85832539549867509;
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quanc8 qt3: ',err:28, est:28);
  end;

  Quanc8({$ifdef FPC_ProcVar}@{$endif}qt4,0,4,0,epsr,y,est,flag,nf);
  ans := 16/3;  {Fix311}
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quanc8 qt4: ',err:28, est:28);
  end;

  Quanc8({$ifdef FPC_ProcVar}@{$endif}qt5,0,1,0,epsr,y,est,flag,nf);
  ans := Pi;
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quanc8 qt5: ',err:28, est:28);
  end;

  Quanc8({$ifdef FPC_ProcVar}@{$endif}qt6,0,2,0,epsr,y,est,flag,nf);
  ans := vers(20.0)/10.0;
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quanc8 qt6: ',err:28, est:28);
  end;

  Quanc8({$ifdef FPC_ProcVar}@{$endif}qt7,-10,10,0,epsr,y,est,flag,nf);
  ans := 4*(2*pi + arctan(3*tan(5))); {19.24270226683792157}
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quanc8 qt7: ',err:28, est:28);
  end;

  for i:=0 to 16 do begin
    alpha := i/8.0;
    Quanc8({$ifdef FPC_ProcVar}@{$endif}qt8,0,1,0,epsr,y,est,flag,nf);
    if i=0 then ans := 1
    else ans := arctan(alpha)/alpha;
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quanc8 qt8/',i,': ',err:28, est:28);
    end;
  end;

  for i:=-7 to 16 do begin
    alpha := i/8.0;
    Quanc8({$ifdef FPC_ProcVar}@{$endif}qt9,0,1,0,epsr,y,est,flag,nf);
    if i=0 then ans := 1
    else ans := ln1p(alpha)/alpha;
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quanc8 qt9/',i,': ',err:28, est:28);
    end;
  end;

  d := 5.0;
  for i:=-10 to 10 do begin
    alpha := i/d; {Avoid 'optimized' code}
    if i=-5 then ans := ln(2)
    else ans := exp2m1(alpha + 1)/(alpha + 1);
    Quanc8({$ifdef FPC_ProcVar}@{$endif}qt10,0,1,0,epsr,y,est,flag,nf);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quanc8 qt10/',i,': ',err:28, est:28);
    end;
  end;

  ans := 10;
  Quanc8({$ifdef FPC_ProcVar}@{$endif}fconst1,0,10,0,epsr,y,est,flag,nf);
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quanc8 fconst1: ',err:28, est:28);
  end;


  if fail=0 then writeln(' ',cnt, ' tests OK.')
  else writeln('*** ',fail,' test of ',cnt, ' failed.');

  inc(totalfail, fail);
  inc(totalcnt,cnt);

end;


{---------------------------------------------------------------------------}
procedure test_quagk;
  {-Test numerical integration with quagk}
var
  ans,y,err,est,d: extended;
  i,fail, cnt, ierr: integer;
const
  eps  = 1e-15;
  epsr = 1e-16;
begin

  fail := 0;
  cnt  := 0;
  writeln('Testing quagk');

  for i:=0 to 16 do begin
    alpha := i/8.0;
    quagk({$ifdef FPC_ProcVar}@{$endif}qt8,0,1,eps,y,est,ierr);
    if i=0 then ans := 1
    else ans := arctan(alpha)/alpha;
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt8/',i,': ',err:28, est:28);
    end;
  end;

  for i:=-7 to 16 do begin
    alpha := i/8.0;
    quagk({$ifdef FPC_ProcVar}@{$endif}qt9,0,1,eps,y,est,ierr);
    if i=0 then ans := 1
    else ans := ln1p(alpha)/alpha;
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt9/',i,': ',err:28, est:28);
    end;
  end;

  d := 5.0;
  for i:=-10 to 10 do begin
    alpha := i/d; {Avoid 'optimized' code}
    if i=-5 then ans := ln(2)
    else ans := exp2m1(alpha + 1)/(alpha + 1);
    quagk({$ifdef FPC_ProcVar}@{$endif}qt10,0,1,eps,y,est,ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt10/',i,': ',err:28, est:28);
    end;
  end;

  d := 10.0;
  for i:=0 to 10 do begin
    {jump at alpha}
    alpha := i/d; {Avoid 'optimized' code}
    ans := 2-alpha;
    quagk({$ifdef FPC_ProcVar}@{$endif}qt11,0,1,eps,y,est,ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt11/',i,': ',err:28, est:28);
    end;
  end;

  alpha := 1/128;
  while alpha <= 128 do begin
    ans := -0.25*Pi/sqrt(alpha)*ln(alpha);
    quagk({$ifdef FPC_ProcVar}@{$endif}qt12, 0, PosInf_x,eps,y,est,ierr);
    if ans=0 then err := y-ans
    else err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt12/',alpha:1:3,': ',err:28, est:28);
    end;
    alpha := 2*alpha;
  end;

  alpha := 1/8;
  while alpha <= 8 do begin
    ans := sqrt(Pi/alpha)*exp(0.25/alpha);
    quagk({$ifdef FPC_ProcVar}@{$endif}qt13, NegInf_x, PosInf_x,eps,y,est,ierr);
    if ans=0 then err := y-ans
    else err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt13/',alpha:1:3,': ',err:28, est:28);
    end;
    alpha := 2*alpha;
  end;

  alpha := 1/128;
  while alpha <= 32 do begin
    ans := sqrtPi*exp(0.25*sqr(alpha));
    quagk({$ifdef FPC_ProcVar}@{$endif}qt14, NegInf_x, PosInf_x,eps,y,est,ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt14/',alpha:1:3,': ',err:28, est:28);
    end;
    alpha := 2*alpha;
  end;

  alpha := 1.0;
  while alpha<=8 do begin
    ans := GammaX(alpha);
    quagk({$ifdef FPC_ProcVar}@{$endif}qt15, 0, PosInf_x,eps,y,est,ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt15/',alpha:1:3,': ',err:28, est:28);
    end;
    alpha := alpha+0.25;
  end;

  alpha := -3/4;
  while alpha<=2 do begin
    ans := -1/(sqr(alpha+1));
    quagk({$ifdef FPC_ProcVar}@{$endif}qt16, 0, 1,eps,y,est,ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt16/',alpha:1:3,': ',err:28, est:28);
    end;
    alpha := alpha+0.25;
  end;

  for i:=0 to 10 do begin
    alpha := i;
    if i=7 then ans:=2
    else ans := -2*(alpha-7)/sqrt(7*abs(alpha - 7)) + 2*sqrt(alpha/7);
    quagk({$ifdef FPC_ProcVar}@{$endif}qt17, 0, 1,eps,y,est,ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('quagk qt17/',alpha:1:3,': ',err:28, est:28);
    end;
  end;

  ans := Pi/2 - Pi*ln(2);
  quagk({$ifdef FPC_ProcVar}@{$endif}qt18, 0, 1,eps,y,est,ierr);
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quagk qt18',err:28, est:28);
  end;

  ans := Pi;
  quagk({$ifdef FPC_ProcVar}@{$endif}qt19, 0, PosInf_x, eps,y,est,ierr);
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quagk qt19',err:28, est:28);
  end;

  ans := -2.5432596188935314899; {ans := 2*Pi*Pi*Pi*Bessel_J1x(60*Pi);}
  quagk({$ifdef FPC_ProcVar}@{$endif}qt20, 0, 2*Pi, eps,y,est,ierr);
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quagk qt20',err:28, est:28);
  end;

  ans := 10;
  quagk({$ifdef FPC_ProcVar}@{$endif}fconst1, 0, 10, eps,y,est,ierr);
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('quagk fconst1: ',err:28, est:28);
  end;

  if fail=0 then writeln(' ',cnt, ' tests OK.')
  else writeln('*** ',fail,' test of ',cnt, ' failed.');

  inc(totalfail, fail);
  inc(totalcnt,cnt);
end;


{---------------------------------------------------------------------------}
procedure test_qawc;
  {-Test numerical integration with qawc}
var
  b,c,ans,y,err,est,t: extended;
  i,fail,cnt,ierr: integer;
  nev: longint;
const
  eps  = 1e-15;
  epsr = 1e-16;
const
  a23  : array[0..4] of extended = (
           -0.89944006957717335193136668554e-1,
           -0.356084474469606222605493250105,
           -0.253657059420611526575851472251,
           -0.170331164082325525449725822827,
           -0.126335831218586456448445523702);
begin

  fail := 0;
  cnt  := 0;
  writeln('Testing qawc');

  for i:=-15 to 32 do begin
    c := i/16;
    if (c<>0) and (c<>1) then begin
      qawc({$ifdef FPC_ProcVar}@{$endif}qt21, -1, 1, c, 0, eps, 0, y, est, nev, ierr);
      if c<1 then ans := -10000*(200*c*arctan(100)-ln(1-c)+ln(c+1))/(10000*c*c+1)
      else ans := 10000*(-200*c*arctan(100)+ln(c-1)-ln(c+1))/(10000*c*c+1);
      err := abs(1.0 - y/ans);
      inc(cnt);
      if err>epsr then begin
        inc(fail);
        writeln('qawc qt21/',i,': ',err:28, est:28);
      end;
    end;
  end;

  for i:=1 to 31 do begin
    c := i/16;
    qawc({$ifdef FPC_ProcVar}@{$endif}qt22, 0, 2, c, 0, eps, 0, y, est, nev, ierr);
    t := c*c*c;
    ans := 8/3 - ln(c) + 2*c - t*ln(c) + 2*c*c + ln(2-c) + ln(2-c)*t;
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('qawc qt22/',i,': ',err:28, est:28);
    end;
  end;

  for i:=1 to 4 do begin
    c := i;
    qawc({$ifdef FPC_ProcVar}@{$endif}qt23, -1, 5, c, 0, eps, 0, y, est, nev, ierr);
    ans := a23[i];
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('qawc qt23/',i,': ',err:28, est:28);
    end;
  end;

  c := Eix(-2.0);
  for i:=-4 to 17 do begin
    if i<>0 then begin
      b := i/16;
      qawc({$ifdef FPC_ProcVar}@{$endif}qt24, -2, b, 0, 0, eps, 0, y, est, nev, ierr);
      ans := Eix(b)-c;
      err := abs(1.0 - y/ans);
      inc(cnt);
      if err>epsr then begin
        inc(fail);
        writeln('qawc qt24/',i,': ',err:28, est:28);
      end;
    end;
  end;

  for i:=0 to 25 do begin
    c := 0.75;
    alpha := i;
    qawc({$ifdef FPC_ProcVar}@{$endif}qt25, 0, 1, c, 0, eps, 0, y, est, nev, ierr);
    if i=0 then ans := ln((1.0-c)/c)
    else ans := c*ans + 1.0/alpha;
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('qawc qt25/',i,': ',err:28, est:28);
    end;
  end;

  if fail=0 then writeln(' ',cnt, ' tests OK.')
  else writeln('*** ',fail,' test of ',cnt, ' failed.');

  inc(totalfail, fail);
  inc(totalcnt,cnt);
end;


{---------------------------------------------------------------------------}
procedure test_de;
  {-Test numerical integration with intdei/intdeo}
var
  ans,y,err,est,d: extended;
  i,fail,cnt,ierr: integer;
  nev: longint;
const
  eps  = 1e-18;
  epsr = 1e-17;
begin
  fail := 0;
  cnt  := 0;
  writeln('Testing intde/intdei/intdeo');

  for i:=0 to 16 do begin
    alpha := i/8.0;
    intde({$ifdef FPC_ProcVar}@{$endif}qt8,0,1,epsr,y,est,nev,ierr);
    if i=0 then ans := 1
    else ans := arctan(alpha)/alpha;
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('intde qt8/',i,': ',err:28, est:28);
    end;
  end;

  for i:=-7 to 16 do begin
    alpha := i/8.0;
    intde({$ifdef FPC_ProcVar}@{$endif}qt9,0,1,epsr,y,est,nev,ierr);
    if i=0 then ans := 1
    else ans := ln1p(alpha)/alpha;
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('intde qt9/',i,': ',err:28, est:28);
    end;
  end;

  d := 5.0;
  for i:=-10 to 10 do begin
    alpha := i/d; {Avoid 'optimized' code}
    if i=-5 then ans := ln(2)
    else ans := exp2m1(alpha + 1)/(alpha + 1);
    intde({$ifdef FPC_ProcVar}@{$endif}qt10,0,1,epsr,y,est,nev,ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('intde qt10/',i,': ',err:28, est:28);
    end;
  end;

  for i:=1 to 15 do begin
    alpha := i/16;
    ans := Pi/sinPi(alpha);
    intdei({$ifdef FPC_ProcVar}@{$endif}qt26, 0, eps, y, est, nev, ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('intdei qt26/',i,': ',err:28, est:28);
    end;
  end;

  for i:=1 to 16 do begin
    alpha := i/8;
    ans := Pi_2/alpha;
    intdei({$ifdef FPC_ProcVar}@{$endif}qt27, 0, eps, y, est, nev, ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('intdei qt27/',i,': ',err:28, est:28);
    end;
  end;

  for i:=1 to 16 do begin
    alpha := i/2;
    ans := sqrt(Pi_2/alpha);
    intdeo({$ifdef FPC_ProcVar}@{$endif}qt28, 0, alpha, eps, y, est, nev, ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('intdeo qt28/',i,': ',err:28, est:28);
    end;
  end;

  for i:=1 to 16 do begin
    alpha := i/4;
    ans := Pi_2*exp(-alpha);
    intdeo({$ifdef FPC_ProcVar}@{$endif}qt29, 0, alpha, eps, y, est, nev, ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('intdeo qt29/',i,': ',err:28, est:28);
    end;
  end;

  d := 7.0;
  for i:=1 to 16 do begin
    alpha := (4*i)/d; {Avoid 'optimized' code}
    ans := sqrt(TwoPi*alpha);
    intdeo({$ifdef FPC_ProcVar}@{$endif}qt30, 0, alpha, eps, y, est, nev, ierr);
    err := abs(1.0 - y/ans);
    inc(cnt);
    if err>epsr then begin
      inc(fail);
      writeln('intdeo qt30/',i,': ',err:28, est:28);
    end;
  end;

  ans := 10;
  intde({$ifdef FPC_ProcVar}@{$endif}fconst1,0,10,epsr,y,est,nev,ierr);
  err := abs(1.0 - y/ans);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('intde fconst1: ',err:28, est:28);
  end;

  if fail=0 then writeln(' ',cnt, ' tests OK.')
  else writeln('*** ',fail,' test of ',cnt, ' failed.');

  inc(totalfail, fail);
  inc(totalcnt,cnt);
end;


{---------------------------------------------------------------------------}
type
  qt_para = record a,b: extended; end;
  pqt_para = ^qt_para;


{---------------------------------------------------------------------------}
function qt31(x: extended; p: pointer): extended;
begin
  with pqt_para(p)^ do begin
    qt31 := (x-a)/power(x+b,3);
  end;
end;


{---------------------------------------------------------------------------}
function qt32(x: extended; p: pointer): extended;
begin
  with pqt_para(p)^ do begin
    qt32 := exp(-a*x)/(x+b);
  end;
end;


{---------------------------------------------------------------------------}
function qt33(x: extended; p: pointer): extended;
begin
  with pqt_para(p)^ do begin
    qt33 := exp(-a*x)/sqrt(x+b);
  end;
end;


{---------------------------------------------------------------------------}
function fconst1p(x: extended; p: pointer): extended;
begin
  fconst1p := 1;
end;


{---------------------------------------------------------------------------}
procedure test_dep;
  {-Test numerical integration with intde_p/intdei_p}
var
  x,y,f,est,err: extended;
  ierr,ia,ib,fail,cnt: integer;
  qpara: qt_para;
  n: longint;
const
  eps = 1e-18;
  epsr = 1e-15;
begin
  fail := 0;
  cnt  := 0;
  writeln('Testing intde_p/intdei_p');

  for ia:=1 to 10 do begin
    qpara.a := ia;
    for ib:=1 to 10 do begin
      qpara.b := ib;
      intdei_p({$ifdef FPC_ProcVar}@{$endif}qt31, @qpara, 0, eps, y, est, n, ierr);
      f := 0.5*(ib-ia)/sqr(ib);
      if abs(f)<0.5 then err:=y-f
      else err := 1-y/f;
      inc(cnt);
      if abs(err) > eps then begin
        inc(fail);
        writeln('dei_p/qt31', ia:3, ib:3, ' ', y:18:15, ' ',f:18:15, ' ', err:23);
      end;
    end;
  end;

  for ia:=2 to 10 do begin
    qpara.a := 0.5*ia;
    for ib:=2 to 10 do begin
      qpara.b := 0.5*ib;
      intdei_p({$ifdef FPC_ProcVar}@{$endif}qt32, @qpara, 0, eps, y, est, n, ierr);
      with qpara do begin
        f := exp(b*a)*E1x(b*a);
      end;
      if abs(f)<0.5 then err:=y-f
      else err := 1-y/f;
      inc(cnt);
      if abs(err) > eps then begin
        inc(fail);
        writeln('dei_p/qt32', ia:3, ib:3, ' ', y:18:15, ' ',f:18:15, ' ', err:23);
      end;
    end;
  end;

  for ia:=1 to 10 do begin
    qpara.a := 0.25*ia;
    for ib:=1 to 10 do begin
      qpara.b := 0.25*ib;
      intde_p({$ifdef FPC_ProcVar}@{$endif}qt32, @qpara, 1, 10, epsr, y, est, n, ierr);
      with qpara do begin
        x := a*b;
        f := exp(x)*(E1x(x+a)-E1x(x+10*a));
      end;
      if abs(f)<0.5 then err:=y-f
      else err := 1-y/f;
      inc(cnt);
      if abs(err) > epsr then begin
        inc(fail);
        writeln(' de_p/qt32', ia:3, ib:3, ' ', y:18:15, ' ',f:18:15, ' ', err:23);
      end;
    end;
  end;

  for ia:=1 to 10 do begin
    qpara.a := 0.25*ia;
    for ib:=1 to 10 do begin
      qpara.b := 0.25*ib;
      intdei_p({$ifdef FPC_ProcVar}@{$endif}qt33, @qpara, 0, eps, y, est, n, ierr);
      with qpara do begin
        f := exp(a*b)*sqrtPi*erfcx(sqrt(a*b))/sqrt(a);
      end;
      if abs(f)<0.5 then err:=y-f
      else err := 1-y/f;
      inc(cnt);
      if abs(err) > eps then begin
        inc(fail);
        writeln('dei_p/qt33', ia:3, ib:3, ' ', y:18:15, ' ',f:18:15, ' ', err:23);
      end;
    end;
  end;

  intde_p({$ifdef FPC_ProcVar}@{$endif}fconst1p,@qpara,0,10,epsr,y,est,n,ierr);
  err := abs(1.0 - y/10.0);
  inc(cnt);
  if err>epsr then begin
    inc(fail);
    writeln('intde_p fconst1_p: ',err:28, est:28);
  end;

  if fail=0 then writeln(' ',cnt, ' tests OK.')
  else writeln('*** ',fail,' test of ',cnt, ' failed.');

  inc(totalfail, fail);
  inc(totalcnt,cnt);

end;


{---------------------------------------------------------------------------}
procedure test_amint;
  {-Test numerical integration}
begin
  writeln('Test AMTools/numerical integration    (c) W.Ehrhardt 2011-2018');
  totalfail := 0;
  totalcnt  := 0;
  test_quanc8;
  test_quagk;
  test_qawc;
  test_de;
  test_dep;
  if totalfail=0 then writeln('All ',totalcnt, ' tests passed.')
  else writeln('*** ',totalfail,' test of ',totalcnt, ' failed.');
end;


{---------------------------------------------------------------------------}
{----------------------- Convergence acceleration --------------------------}
{---------------------------------------------------------------------------}

procedure test_amconvacc;
  {-Convergence acceleration methods}
const
  NMAX=100;
  relerr = 1e-18;
  reltol = 2*relerr;
var
  qnum, qden: array[0..NMAX] of extended;
  j,jc,ierr,numfail: integer;
  ref,term,sum: extended;
  res,resold,err,errold,rel: extended;
  fail,done: boolean;
type
  str8= string[8];
const
  status: array[false..true] of str8 = (' OK', ' failed!');
begin
  writeln('Test AMTools/Convergence acceleration methods    (c) W.Ehrhardt 2013');
  res := 0;
  err := 0;
  jc  := 0;
  numfail := 0;

{---------------------------------------------------------------------------}
  ref := 0.604898643421630370247265914236;  {Zeta(1/2)*(1-sqrt(2))}

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1/sqrt(j+1);
    if odd(j) then term := -term;
    resold := res;
    errold := err;
    levinu1(term, j, NMAX, qnum, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 1  Levin  u: ',jc:3, res:22:18, rel:28, status[fail]);

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1/sqrt(j+1);
    if odd(j) then term := -term;
    resold := res;
    errold := err;
    wynneps1(term, j, NMAX, true, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 2  Wynn eps: ',jc:3, res:22:18, rel:28, status[fail]);

{---------------------------------------------------------------------------}
  ref := 0.380104812609684016777542156552;  {Zeta(-1/2)*(1-2*sqrt(2))}

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := sqrt(j+1);
    if j and 1 =1  then term := -term;
    resold := res;
    errold := err;
    levinu1(term, j, NMAX, qnum, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 3  Levin  u: ',jc:3, res:22:18, rel:28, status[fail]);

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := sqrt(j+1);
    if odd(j) then term := -term;
    resold := res;
    errold := err;
    wynneps1(term, j, NMAX, true, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 4  Wynn eps: ',jc:3, res:22:18, rel:28, status[fail]);

{---------------------------------------------------------------------------}

  ref := 0.225791352644727432363097614947;

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := ln(j+2);
    if odd(j) then term := -term;
    resold := res;
    errold := err;
    levinu1(term, j, NMAX, qnum, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 5  Levin  u: ',jc:3, res:22:18, rel:28, status[fail]);

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := ln(j+2);
    if odd(j) then term := -term;
    resold := res;
    errold := err;
    wynneps1(term, j, NMAX, true,qnum, res, sum, ierr);
    err := abs(res-resold);
    if (j>1) then begin
      if (err>1.1*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 6  Wynn eps: ',jc:3, res:22:18, rel:28, status[fail]);

{---------------------------------------------------------------------------}

  ref := 0.924299897222938855959570181360;

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1/ln(j+2);
    if odd(j) then term := -term;
    resold := res;
    errold := err;
    levinu1(term, j, NMAX, qnum, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 7  Levin  u: ',jc:3, res:22:18, rel:28, status[fail]);

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1/ln(j+2);
    if odd(j) then term := -term;
    resold := res;
    errold := err;
    wynneps1(term, j, NMAX, true,qnum, res, sum, ierr);
    err := abs(res-resold);
    if (j>1) then begin
      if (err>1.1*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 8  Wynn eps: ',jc:3, res:22:18, rel:28, status[fail]);

{---------------------------------------------------------------------------}
  ref := 0.567143290409783872999968662210; {LambertW(1)}
  {Iteration sequence x(n) = exp(-x(n-1)), x(-1)=0.5}
  term := 0.5;
  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := exp(-term);
    resold := res;
    errold := err;
    wynneps1(term, j, NMAX, false, qnum, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln(' 9  Wynn eps: ',jc:3, res:22:18, rel:28, status[fail]);

{---------------------------------------------------------------------------}
  ref := 0.618033988749894848204586834360; {(sqrt(5)-1)/2}
  {Iteration sequence x(n) = 1/(1+x(n-1)), x(-1)=0.75}
  term := 0.75;
  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1/(term + 1);
    resold := res;
    errold := err;
    wynneps1(term, j, NMAX, false, qnum, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln('10  Wynn eps: ',jc:3, res:22:18, rel:28, status[fail]);

{---------------------------------------------------------------------------}
  ref := 2.07112132996782316957157299268;

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1/cosh(j);
    resold := res;
    errold := err;
    levinu1(term, j, NMAX, qnum, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>2*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > 20*reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln('11  Levin  u: ',jc:3, res:22:18, rel:28, status[fail]);

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1/cosh(j);
    resold := res;
    errold := err;
    wynneps1(term, j, NMAX, true, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  rel := 1.0-res/ref;
  fail := (abs(rel) > reltol) or (ierr<>0);
  if fail then inc(numfail);
  writeln('12  Wynn eps: ',jc:3, res:22:18, rel:28, status[fail]);

{---------------------------------------------------------------------------}
  ref := sqr(Pi)/6;

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1.0/sqr(j+1.0);
    resold := res;
    errold := err;
    levinu1(term, j, NMAX, qnum, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>2*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  {Logarithm convergence: only moderate accuracy with Levin u}
  rel := 1.0-res/ref;
  fail := (abs(rel) > 1e-12) or (ierr<>0);
  if fail then inc(numfail);
  writeln('13  Levin  u: ',jc:3, res:22:18, rel:28, status[fail], ' **');

  j := 0;
  done := false;
  while (j<=NMAX) and not done do begin
    jc := j;
    term := 1.0/sqr(j+1.0);
    resold := res;
    errold := err;
    wynneps1(term, j, NMAX, true, qden, res, sum, ierr);
    err := abs(res-resold);
    if j>1 then begin
      if (err>1.25*errold) or (err<=abs(res)*relerr) then done := true;
    end;
    inc(j);
  end;
  {Logarithmic convergence: Wynn epsilon will not accelerate!!}
  rel := 1.0-res/ref;
  fail := (abs(rel) > 5e-3) or (ierr<>0);
  if fail then inc(numfail);
  writeln('14  Wynn eps: ',jc:3, res:22:18, rel:28, status[fail], ' **');
  writeln('**  special test with logarithmic convergence');

  if numfail > 0 then writeln('** ', numfail, ' of 14 tests failed!')
  else writeln('All 14 tests passed.');
end;


end.
