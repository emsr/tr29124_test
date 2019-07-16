{Part 0 of regression test for DAMath unit,   (c) 2013  W.Ehrhardt}
{Tests for basic/statistic functions}

unit t_damat0;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$ifndef Windows}
    {$O+}
  {$endif}
{$endif}


interface


const
  sqrt2  = 1.4142135623730950488;
  sqrt3  = 1.7320508075688772935;
  sqrt_5 = 0.7071067811865475244;

type
  TPair = record
            tx,ty: double;
          end;

var
  EPS, SEPS, EPSD: double;

var
  total_cnt, total_failed: longint;

function reldev(a,b: double): double;
  {-Return abs((a-b)/b) if b<>0,  abs(a-b) if b=0, or a-b=0}

procedure Test_Consts;
procedure test_frexp_ldexp_ilogb;
procedure Test_IsInfNan;
procedure Test_dminmax;
procedure Test_sminmax;
procedure test_succ_pred;
procedure test_ceil;
procedure test_floor;
procedure test_ulp;
procedure test_rint;

procedure test_meansdev;
procedure test_stat1;
procedure test_stat_moment;

procedure test_angle2;
procedure test_geometry2;

procedure test_hex2float;

implementation

uses
  damath;


{---------------------------------------------------------------------------}
function reldev(a,b: double): double;
  {-Return abs((a-b)/b) if b<>0,  abs(a-b) if b=0, or a-b=0}
var
  c: double;
begin
  c := a-b;
  if c=0.0 then reldev := 0.0
  else begin
    if b<>0 then c := c/b;
    reldev := abs(c);
  end;
end;


{---------------------------------------------------------------------------}
function dbldiff(x,y: double): double;
  {-Diff x-y in as double}
begin
  x := x - y;
  dbldiff := x;
end;


{---------------------------------------------------------------------------}
procedure Test_Consts;
var
  e: double;
begin
  writeln('Test constants');
  {Pi := system.Pi;}
  writeln('MinSingle   : ', MinSingle   );
  writeln('MaxSingle   : ', MaxSingle   );
  writeln('PosInf_s    : ', PosInf_s    );
  writeln('NegInf_s    : ', NegInf_s    );
  writeln('NaN_s       : ', NaN_s       );
  writeln('eps_s       : ', eps_s, '   ', Sgl2Hex(eps_s));

  writeln('MinDouble   : ', MinDouble   );
  writeln('MaxDouble   : ', MaxDouble   );
  writeln('PosInf_d    : ', PosInf_d    );
  writeln('NegInf_d    : ', NegInf_d    );
  writeln('NaN_d       : ', NaN_d       );
  writeln('eps_d       : ', eps_d, '   ', Dbl2Hex(eps_d));

  e := exp(1.0);
  writeln('log2e       :', log2e :24,    ',  diff =', dbldiff(log2e  , log2(e)) :24);
  writeln('log10e      :', log10e:24,    ',  diff =', dbldiff(log10e , log10(e)):24);
  writeln('TwoPi       :', TwoPi:24,     ',  diff =', dbldiff(TwoPi  , 2.0*Pi)  :24);
  writeln('Pi_2        :', Pi_2:24,      ',  diff =', dbldiff(Pi_2   , Pi/2.0)  :24);
  writeln('Pi_4        :', Pi_4:24,      ',  diff =', dbldiff(Pi_4   , Pi/4.0)  :24);
  writeln('Sqrt_TwoPi  :', Sqrt_TwoPi:24,',  diff =', dbldiff(Sqrt_TwoPi , sqrt(2.0*Pi)):24);
  writeln('LnSqrt2Pi   :', LnSqrt2Pi:24, ',  diff =', dbldiff(LnSqrt2Pi  , ln(sqrt(2.0*Pi))):24);
  writeln('succd0      :', succd0:24,    ',  diff =', dbldiff(succd0,succd(0)):24);
  writeln('ln_succd0   :', ln_succd0:24, ', rdiff =', 1 - ln_succd0 /ln(succd(0)):24);
  writeln;

  e := Pi;
  if e <> 2.0*Pi_2  then begin
    writeln('Fatal error: 2*Pi_2 <> Pi !');
    halt;
  end;
end;


{---------------------------------------------------------------------------}
procedure test_frexp_ldexp_ilogb;
{$ifndef BIT16}
const NTST = 400000;
{$else}
const NTST = 100000;
{$endif}
  procedure remtest(a,b: double; n: integer);
  begin
    if a<>b then begin
      writeln('Error remainder test ', n, ' ', Dbl2Hex(a), ' ', Dbl2Hex(b));
    end;
  end;
var
  x,y,m,r,z: double;
  s,sm,sy: single;
  e,l,i: longint;
begin
  writeln('Test frexpd/ldexpd/scalbnd/ilogb');

  x := succd(0);
  frexpd(x,m,e);
  y := log2(x);
  if (e<>-1073) or (m<>0.5) or (round(y+1)<>e) then begin
    writeln('frexpd(succd(0) error m,e',m:22,e:8);
    halt;
  end;
  l := ilogb(x);
  if l<>-1074 then begin
    writeln('ilogb(succd(0) error: ',l);
    halt;
  end;

  s := succs(0);
  frexps(s,sm,e);
  y := log2(s);
  if (e<>-148) or (sm<>0.5) or (round(y+1)<>e) then begin
    writeln('frexpd(succs(0) error sm,e',sm:22,e:8);
    halt;
  end;

  for i:=1 to NTST do begin
    if odd(i) then x := (random-0.5)*MaxDouble
    else x := (random-0.5);
    frexpd(x,m,e);
    y := ldexpd(m,e);
    if x<>y then begin
      writeln('frexpd/ldexpd error. x,m,e,y',x:22,m:22,e:8,y:22);
      halt;
    end;
    y := scalbn(m,e);
    if x<>y then begin
      writeln('frexpd/ldexpd error. x,m,e,y',x:22,m:22,e:8,y:22);
      halt;
    end;
    l := ilogb(x);
    if abs(l)<>MaxLongint then begin
      if succ(l)<>e then begin
        writeln('frexpd/ilogb error. x,m,e,l',x:22,m:22,e:8,l:8);
        halt;
      end;
    end;
    if odd(i) then s := (random-0.5)*MaxSingle
    else s := (random-0.5);
    frexps(s,sm,e);
    sy := ldexps(sm,e);
    if s<>sy then begin
      writeln('frexps/ldexs error. s,sm,e,sy',s:22,sm:22,e:8,sy:22);
      halt;
    end;
  end;

  x := 81*succd(0);
  y := fmod(MaxDouble,x);
  y := y/succd(0);
  {(2^1024-2^(1024-53))*2^1074 mod 81 = 56}
  if round(y)<>56 then begin
    writeln('Error fmod(MaxDouble, 81*succd(0))');
  end;

  x := 81*MinDouble;
  y := fmod(MaxDouble,x);
  {(2^1024-2^(1024-53))*2^1022 mod 81 = 62}
  if round(y/MinDouble)<>62 then begin
    writeln('Error fmod(MaxDouble, 81*MinDouble)');
  end;

  x := 81;
  y := fmod(MaxDouble,x);
  {2^1024-2^(1024-53) mod 81 = 20}
  if round(y)<>20 then begin
    writeln('Error fmod(MaxDouble,81)');
  end;

  y := -0.75;
  x := 0.5;
  z := remainder(x,y);
  remtest(z, -0.25, 1);

  {The following test cases look somewhat surprising, but they}
  {give the same results with MPArith, 387 FPREM1, PurePascal.}
  r := predd(0.25*MaxDouble);

  y := -0.75*MaxDouble;
  x := 0.5*MaxDouble;
  z := remainder(x,y);
  remtest(z,-r,2);

  y := 0.75*MaxDouble;
  x := -0.5*MaxDouble;
  z := remainder(x,y);
  remtest(z,r,3);

  y := 0.75*MaxDouble;
  x := 0.5*MaxDouble;
  z := remainder(x,y);
  remtest(z,-r,4);

  y := -0.75*MaxDouble;
  x := -0.5*MaxDouble;
  z := remainder(x,y);
  remtest(z,r,5);

  y := MinDouble;
  x := -0.5*MinDouble;
  z := remainder(x,y);
  remtest(z,x,6);

  y := -MinDouble;
  x := -0.5*MinDouble;
  z := remainder(x,y);
  remtest(z,x,7);

  y := -MinDouble;
  x := 0.5*MinDouble;
  z := remainder(x,y);
  remtest(z,x,8);

  y := MinDouble;
  x := 0.5*MinDouble;
  z := remainder(x,y);
  remtest(z,x,9);

  y := MinDouble;
  x := 1.5*MinDouble;
  z := remainder(x,y);
  r := -0.5*MinDouble;
  remtest(z,r,10);

  y := MinDouble;
  x := 2.5*MinDouble;
  z := remainder(x,y);
  r := 0.5*MinDouble;
  remtest(z,r,11);

  y := 1.0;
  x := Pi;
  z := remainder(x,y);
  r := Pi;
  r := r-3.0;
  remtest(z,r,12);

  y := 1.0;
  x := Pi_2;
  z := remainder(x,y);
  r := Pi_2;
  r := r-2.0;
  remtest(z,r,13);

end;


{---------------------------------------------------------------------------}
procedure test_hex2float;
var
  d1,d2: double;
  s1,s2: single;
  k: longint;
  code: integer;
  hex: string;
const
  NTST = 20000;
begin
  writeln('Test Hex2Float/Float2Hex');
  for k:=1 to NTST do begin
    if odd(k) then d1 := (random-0.5)*MaxDouble
    else d1 := (random-0.5);
    hex := Dbl2Hex(d1);
    if k and 2 <>0 then hex := '$'+hex;
    Hex2Dbl(hex,d2,code);
    if (d1<>d2) or (code<>0) then begin
      writeln('Hex2Dbl error: ',d1:28,d2:28,code:5);
      halt;
    end;
    if odd(k) then s1 := (random-0.5)*MaxSingle
    else s1 := (random-0.5);
    hex := Sgl2Hex(s1);
    if k and 2 <>0 then hex := '$'+hex;
    Hex2Sgl(hex,s2,code);
    if (s1<>s2) or (code<>0) then begin
      writeln('Hex2Sgl error: ',s1:28,s2:28,code:5);
      halt;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure Test_IsInfNan;
begin
  writeln('Test IsInf/IsNAN/IsNaNorInf');
  {$ifdef FPC}
   {$ifdef VER1}
    writeln('**** Note: FPC 1 prints INFs as NANs!');
   {$endif}
  {$endif}

  (* double *)
  if not IsInfD(PosInf_d)     then writeln('Error: IsInfD(PosInf_d)   ');
  if not IsInfD(NegInf_d)     then writeln('Error: IsInfD(NegInf_d)   ');
  if     IsInfD(NaN_d)        then writeln('Error: IsInfD(NaN_d)      ');
  if not IsInfD(PosInf_s)     then writeln('Error: IsInfD(PosInf_s)   ');
  if not IsInfD(NegInf_s)     then writeln('Error: IsInfD(NegInf_s)   ');
  if     IsInfD(NaN_s)        then writeln('Error: IsInfD(NaN_s)      ');
  if     IsInfD(MaxDouble)    then writeln('Error: IsInfD(MaxDouble)');
  if     IsInfD(MinDouble)    then writeln('Error: IsInfD(MinDouble)');

  if     IsNaND(PosInf_d)     then writeln('Error: IsNaND(PosInf_d)   ');
  if     IsNaND(NegInf_d)     then writeln('Error: IsNaND(NegInf_d)   ');
  if not IsNaND(NaN_d)        then writeln('Error: IsNaND(NaN_d)      ');
  if     IsNaND(PosInf_s)     then writeln('Error: IsNaND(PosInf_s)   ');
  if     IsNaND(NegInf_s)     then writeln('Error: IsNaND(NegInf_s)   ');
  if not IsNaND(NaN_s)        then writeln('Error: IsNaND(NaN_s)      ');
  if     IsNaND(MaxDouble)    then writeln('Error: IsNaND(MaxDouble)');
  if     IsNaND(MinDouble)    then writeln('Error: IsNaND(MinDouble)');

  if not IsNaNorInfD(PosInf_d)   then writeln('Error: IsNaNorInfD(PosInf_d)   ');
  if not IsNaNorInfD(NegInf_d)   then writeln('Error: IsNaNorInfD(NegInf_d)   ');
  if not IsNaNorInfD(NaN_d)      then writeln('Error: IsNaNorInfD(NaN_d)      ');
  if not IsNaNorInfD(PosInf_s)   then writeln('Error: IsNaNorInfD(PosInf_s)   ');
  if not IsNaNorInfD(NegInf_s)   then writeln('Error: IsNaNorInfD(NegInf_s)   ');
  if not IsNaNorInfD(NaN_s)      then writeln('Error: IsNaNorInfD(NaN_s)      ');
  if     IsNaNorInfD(MaxDouble)  then writeln('Error: IsNaNorInfD(MaxDouble)');
  if     IsNaNorInfD(MinDouble)  then writeln('Error: IsNaNorInfD(MinDouble)');

  (* single *)
  if not IsInfS(PosInf_s)     then writeln('Error: IsInfS(PosInf_s)   ');
  if not IsInfS(NegInf_s)     then writeln('Error: IsInfS(NegInf_s)   ');
  if     IsInfS(NaN_s)        then writeln('Error: IsInfS(NaN_s)      ');
  if     IsInfS(MaxSingle)    then writeln('Error: IsInfS(MaxSingle)');
  if     IsInfS(MinSingle)    then writeln('Error: IsInfS(MinSingle)');

  if     IsNaNS(PosInf_s)     then writeln('Error: IsNaNS(PosInf_s)   ');
  if     IsNaNS(NegInf_s)     then writeln('Error: IsNaNS(NegInf_s)   ');
  if not IsNaNS(NaN_s)        then writeln('Error: IsNaNS(NaN_s)      ');
  if     IsNaNS(MaxSingle)    then writeln('Error: IsNaNS(MaxSingle)');
  if     IsNaNS(MinSingle)    then writeln('Error: IsNaNS(MinSingle)');

  if not IsNaNorInfS(PosInf_s)   then writeln('Error: IsNaNorInfS(PosInf_s)   ');
  if not IsNaNorInfS(NegInf_s)   then writeln('Error: IsNaNorInfS(NegInf_s)   ');
  if not IsNaNorInfS(NaN_s)      then writeln('Error: IsNaNorInfS(NaN_s)      ');
  if     IsNaNorInfS(MaxSingle)  then writeln('Error: IsNaNorInfS(MaxSingle)');
  if     IsNaNorInfS(MinSingle)  then writeln('Error: IsNaNorInfS(MinSingle)');
end;


{---------------------------------------------------------------------------}
procedure Test_dminmax;
var
  x,y: double;
begin
  writeln('Test mind/maxd');

  x := PosInf_d;
  y := 100;
  if mind(x,y)<>y  then writeln('Error: mind(Inf,100)');
  if maxd(x,y)<>x  then writeln('Error: maxd(Inf,100)');

  x := NegInf_d;
  y := 100;
  if mind(x,y)<>x  then writeln('Error: mind(-Inf,100)');
  if maxd(x,y)<>y  then writeln('Error: maxd(-Inf,100)');

  x := PosInf_d;
  y := predd(x);
  if mind(x,y)<>y  then writeln('Error: mind(Inf, predd(Inf))');
  if maxd(x,y)<>x  then writeln('Error: maxd(Inf, predd(Inf))');

  x := NegInf_d;
  y := succd(x);
  if mind(x,y)<>x  then writeln('Error: mind(-Inf, succd(-Inf))');
  if maxd(x,y)<>y  then writeln('Error: maxd(-Inf, succd(-Inf))');

  x := 0.0;
  y := succd(x);
  if mind(x,y)<>x  then writeln('Error: mind(0, succd(0))');
  if maxd(x,y)<>y  then writeln('Error: maxd(0, succd(0))');

  x := 0.0;
  y := predd(x);
  if mind(x,y)<>y  then writeln('Error: mind(0, predd(0))');
  if maxd(x,y)<>x  then writeln('Error: maxd(0, predd(0))');

  x := 0.0;
  y := -x;
  if mind(x,y)<>y  then writeln('Error: mind(0, -0)');
  if maxd(x,y)<>x  then writeln('Error: maxd(0, -0)');

  x := 0.0;
  y := 1.0;
  if mind(x,y)<>x  then writeln('Error: mind(0, 1)');
  if maxd(x,y)<>y  then writeln('Error: maxd(0, 1)');

  x := 0.0;
  y := -1.0;
  if mind(x,y)<>y  then writeln('Error: mind(0, -1)');
  if maxd(x,y)<>x  then writeln('Error: maxd(0, -1)');

  x := 1.0;
  y := MaxDouble;
  if mind(x,y)<>x  then writeln('Error: mind(1, MaxDouble)');
  if maxd(x,y)<>y  then writeln('Error: maxd(1, MaxDouble)');

  x := 1.0;
  y := -MaxDouble;
  if mind(x,y)<>y  then writeln('Error: mind(1, -MaxDouble)');
  if maxd(x,y)<>x  then writeln('Error: maxd(1, -MaxDouble)');

end;


{---------------------------------------------------------------------------}
procedure Test_sminmax;
var
  x,y: single;
begin
  writeln('Test mins/maxs');

  x := PosInf_d;
  y := 100;
  if mins(x,y)<>y  then writeln('Error: mins(Inf,100)');
  if maxs(x,y)<>x  then writeln('Error: maxs(Inf,100)');

  x := NegInf_d;
  y := 100;
  if mins(x,y)<>x  then writeln('Error: mins(-Inf,100)');
  if maxs(x,y)<>y  then writeln('Error: maxs(-Inf,100)');

  x := PosInf_d;
  y := preds(x);
  if mins(x,y)<>y  then writeln('Error: mins(Inf, preds(Inf))');
  if maxs(x,y)<>x  then writeln('Error: maxs(Inf, preds(Inf))');

  x := NegInf_d;
  y := succs(x);
  if mins(x,y)<>x  then writeln('Error: mins(-Inf, succs(-Inf))');
  if maxs(x,y)<>y  then writeln('Error: maxs(-Inf, succs(-Inf))');

  x := 0.0;
  y := succs(x);
  if mins(x,y)<>x  then writeln('Error: mins(0, succs(0))');
  if maxs(x,y)<>y  then writeln('Error: maxs(0, succs(0))');

  x := 0.0;
  y := preds(x);
  if mins(x,y)<>y  then writeln('Error: mins(0, preds(0))');
  if maxs(x,y)<>x  then writeln('Error: maxs(0, preds(0))');

  x := 0.0;
  y := -x;
  if mins(x,y)<>y  then writeln('Error: mins(0, -0)');
  if maxs(x,y)<>x  then writeln('Error: maxs(0, -0)');

  x := 0.0;
  y := 1.0;
  if mins(x,y)<>x  then writeln('Error: mins(0, 1)');
  if maxs(x,y)<>y  then writeln('Error: maxs(0, 1)');

  x := 0.0;
  y := -1.0;
  if mins(x,y)<>y  then writeln('Error: mins(0, -1)');
  if maxs(x,y)<>x  then writeln('Error: maxs(0, -1)');

  x := 1.0;
  y := Maxsingle;
  if mins(x,y)<>x  then writeln('Error: mins(1, Maxsingle)');
  if maxs(x,y)<>y  then writeln('Error: maxs(1, Maxsingle)');

  x := 1.0;
  y := -Maxsingle;
  if mins(x,y)<>y  then writeln('Error: mins(1, -Maxsingle)');
  if maxs(x,y)<>x  then writeln('Error: maxs(1, -Maxsingle)');

end;


{---------------------------------------------------------------------------}
procedure test_succ_pred;
type
  ThreeD = array[0..2] of THexDblW;
  ThreeS = array[0..2] of THexSglW;
const
  NumS = 20;
  NumD = 22;
const
  TSS : array[1..NumS] of ThreeS = (
                          {pred}         {succ}
          (($CCCD,$C05C), ($CCCE,$C05C), ($CCCC,$C05C)),   {-3.45}
          (($1198,$531E), ($1197,$531E), ($1199,$531E)),   {67.89Ee+10}
          (($0000,$0000), ($0001,$8000), ($0001,$0000)),   {0}
          (($0000,$8000), ($0001,$8000), ($0001,$0000)),   {-0}
          (($0000,$7F80), ($FFFF,$7F7F), ($0000,$7F80)),   {inf}
          (($0000,$FF80), ($0000,$FF80), ($FFFF,$FF7F)),   {-inf}
          (($FFFF,$7F7F), ($FFFE,$7F7F), ($0000,$7F80)),   {MaxSingle}
          (($FFFF,$FF7F), ($0000,$FF80), ($FFFE,$FF7F)),   {-MaxSingle}
          (($0000,$0080), ($FFFF,$007F), ($0001,$0080)),   {MinSingle}
          (($0000,$8080), ($0001,$8080), ($FFFF,$807F)),   {-MinSingle}
          (($0000,$8040), ($0001,$8040), ($FFFF,$803F)),   {}
          (($0000,$0040), ($FFFF,$003F), ($0001,$0040)),   {}
          (($0000,$FFC0), ($0000,$FFC0), ($0000,$FFC0)),   {Nan}
          (($FFFF,$808F), ($0000,$8090), ($FFFE,$808F)),
          (($FFFF,$8040), ($0000,$8041), ($FFFE,$8040)),
          (($FFFF,$8080), ($0000,$8081), ($FFFE,$8080)),
          (($FFFF,$007F), ($FFFE,$007F), ($0000,$0080)),
          (($FFFF,$008F), ($FFFE,$008F), ($0000,$0090)),
          (($FFFF,$FF4F), ($0000,$FF50), ($FFFE,$FF4F)),
          (($FFFF,$FFCF), ($FFFF,$FFCF), ($FFFF,$FFCF))    {Nan}
        );

const
  TSD : array[1..NumD] of ThreeD = (
          (($0000,$0000,$0000,$0000), ($0001,$0000,$0000,$8000), ($0001,$0000,$0000,$0000)),
          (($0000,$0000,$0000,$0010), ($FFFF,$FFFF,$FFFF,$000F), ($0001,$0000,$0000,$0010)),
          (($0000,$0000,$0000,$3FF0), ($FFFF,$FFFF,$FFFF,$3FEF), ($0001,$0000,$0000,$3FF0)),
          (($0000,$0000,$0000,$7FF0), ($FFFF,$FFFF,$FFFF,$7FEF), ($0000,$0000,$0000,$7FF0)),
          (($0000,$0000,$0000,$8000), ($0001,$0000,$0000,$8000), ($0001,$0000,$0000,$0000)),
          (($0000,$0000,$0000,$8010), ($0001,$0000,$0000,$8010), ($FFFF,$FFFF,$FFFF,$800F)),
          (($0000,$0000,$0000,$BFF0), ($0001,$0000,$0000,$BFF0), ($FFFF,$FFFF,$FFFF,$BFEF)),
          (($0000,$0000,$0000,$FFF0), ($0000,$0000,$0000,$FFF0), ($FFFF,$FFFF,$FFFF,$FFEF)),
          (($0000,$0000,$0001,$0000), ($FFFF,$FFFF,$0000,$0000), ($0001,$0000,$0001,$0000)),
          (($0000,$0000,$FFFF,$800F), ($0001,$0000,$FFFF,$800F), ($FFFF,$FFFF,$FFFE,$800F)),
          (($0001,$0000,$0000,$0000), ($0000,$0000,$0000,$0000), ($0002,$0000,$0000,$0000)),
          (($0001,$0000,$0000,$3FF0), ($0000,$0000,$0000,$3FF0), ($0002,$0000,$0000,$3FF0)),
          (($0001,$0000,$0000,$8000), ($0002,$0000,$0000,$8000), ($0000,$0000,$0000,$8000)),
          (($0001,$0000,$0000,$8010), ($0002,$0000,$0000,$8010), ($0000,$0000,$0000,$8010)),
          (($0001,$0000,$0000,$BFF0), ($0002,$0000,$0000,$BFF0), ($0000,$0000,$0000,$BFF0)),
          (($FFFF,$FFFF,$0000,$0000), ($FFFE,$FFFF,$0000,$0000), ($0000,$0000,$0001,$0000)),
          (($FFFF,$FFFF,$FFFF,$000F), ($FFFE,$FFFF,$FFFF,$000F), ($0000,$0000,$0000,$0010)),
          (($FFFF,$FFFF,$FFFF,$3FEF), ($FFFE,$FFFF,$FFFF,$3FEF), ($0000,$0000,$0000,$3FF0)),
          (($FFFF,$FFFF,$FFFF,$7FEF), ($FFFE,$FFFF,$FFFF,$7FEF), ($0000,$0000,$0000,$7FF0)),
          (($FFFF,$FFFF,$FFFF,$7FFF), ($FFFF,$FFFF,$FFFF,$7FFF), ($FFFF,$FFFF,$FFFF,$7FFF)),
          (($FFFF,$FFFF,$FFFF,$BFEF), ($0000,$0000,$0000,$BFF0), ($FFFE,$FFFF,$FFFF,$BFEF)),
          (($FFFF,$FFFF,$FFFF,$FFEF), ($0000,$0000,$0000,$FFF0), ($FFFE,$FFFF,$FFFF,$FFEF))
        );
var
  i, cnt, failed: integer;
  d,pd,sd: double;
  s,ps,ss: single;
begin
  writeln('Test predd/s, succd/s');
  cnt :=0;
  failed := 0;

  for i:=1 to Nums do begin
    s  := single(TSS[i][0]);
    ps := preds(s);
    ss := succs(s);
    inc(cnt,2);
    if fisNEs(ps,single(TSS[i][1])) then begin
      inc(failed);
      writeln('preds failed: ',Sgl2Hex(ps), ' <> ', Sgl2Hex(single(TSS[i][1])));
    end;
    if fisNEs(ss,single(TSS[i][2])) then begin
      inc(failed);
      writeln('succs failed: ',Sgl2Hex(ss), ' <> ', Sgl2Hex(single(TSS[i][2])));
    end;
  end;

  for i:=1 to NumD do begin
    d  := double(TSD[i][0]);
    pd := predd(d);
    sd := succd(d);
    inc(cnt,2);
    if fisNEd(pd,double(TSD[i][1])) then begin
      inc(failed);
      writeln('predd failed: ',Dbl2Hex(pd), ' <> ', Dbl2Hex(double(TSD[i][1])));
    end;
    if fisNEd(sd,double(TSD[i][2])) then begin
      inc(failed);
      writeln('succd failed: ',Dbl2Hex(sd), ' <> ', Dbl2Hex(double(TSD[i][2])));
    end;
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_ulp;
var
  da: THexDblA; d: double absolute da;
  sa: THexSglA; s: single absolute sa;
  ud: double;
  us: single;
var
  i, cnt, failed: integer;
begin
  writeln('Test ulpd/s');
  failed := 0;
  cnt := 6;
  if ulpd(PosInf_d)<>PosInf_d then begin
    writeln('error ulp(PosInf_d)');
    inc(failed);
  end;
  if ulpd(NegInf_d)<>PosInf_d then begin
    writeln('error ulp(NegInf_d)');
    inc(failed);
  end;
  if ulps(PosInf_s)<>PosInf_s then begin
    writeln('error ulp(PosInf_s)');
    inc(failed);
  end;
  if ulpd(NegInf_s)<>PosInf_s then begin
    writeln('error ulp(NegInf_s)');
    inc(failed);
  end;
  while cnt<3000 do begin
    inc(cnt,3);
    for i:=0 to sizeof(sa)-1 do sa[i] := random(256);
    for i:=0 to sizeof(da)-1 do da[i] := random(256);
    if ISNanS(s) then sa[2] := sa[2] or $40; {Clear signal}
    if ISNanD(d) then da[6] := da[6] or $08; {Clear signal}
    us := ulps(s);
    ud := ulpd(d);
    if IsNand(d) then begin
      if fisNEd(ud,d) then begin
        inc(failed);
        writeln('ulpd failed for ',Dbl2Hex(d));
      end;
    end
    else if IsInfd(d) then begin
      if ud<>PosInf_d then begin
        inc(failed);
        writeln('ulpd failed for ',Dbl2Hex(d));
      end;
    end
    else begin
      if ud <> abs(d)-predd(abs(d)) then begin
        inc(failed);
        writeln('ulpd failed for ',Dbl2Hex(d));
      end;
    end;
    if IsNanS(s) then begin
      if fisNEs(us,s) then begin
        inc(failed);
        writeln('ulps failed for ',Sgl2Hex(s));
      end;
    end
    else if IsInfs(s) then begin
      if ud<>PosInf_s then begin
        inc(failed);
        writeln('ulps failed for ',Sgl2Hex(s));
      end;
    end
    else begin
      if us <> abs(s)-preds(abs(s)) then begin
        inc(failed);
        writeln('ulps failed for ',Sgl2Hex(s));
      end;
    end;
  end;
  if failed>0 then begin
    writeln('*** failed');
    inc(total_failed);
  end
  else writeln(' tests OK');
  inc(total_cnt);
end;


{---------------------------------------------------------------------------}
procedure test_ceil;
const
  NT=8;
  TData: array[1..NT] of TPair = (
           (tx:  1.1;  ty: 2.0),
           (tx: -1.1;  ty: -1.0),
           (tx:  0.1;  ty: 1.0),
           (tx: -0.1;  ty: 0.0),
           (tx:  1.00000000000001; ty: 2.0),
           (tx: -1.00000000000001; ty: -1.0),
           (tx:  1e20;  ty: 1e20),
           (tx: -1e20;  ty: -1e20)
         );
const
  name = 'ceil/ceild';
var
  i, cnt, failed: integer;
  zl: longint;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := ceild(x);
    inc(cnt);
    if y<>z then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ', z-y:18);
    end;
    if abs(y)<MaxLongint then begin
      x := tx;
      y := ty;
      zl:= ceil(x);
      inc(cnt);
      if round(y)<>zl then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19);
      end;
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_floor;
const
  NT=8;
  TData: array[1..NT] of TPair = (
           (tx:  1.1;  ty:  1.0),
           (tx: -1.1;  ty: -2.0),
           (tx:  0.1;  ty:  0.0),
           (tx: -0.1;  ty: -1.0),
           (tx:  1.00000000000001; ty:  1.0),
           (tx: -1.00000000000001; ty: -2.0),
           (tx:  1e20;  ty: 1e20),
           (tx: -1e20;  ty: -1e20)
         );
const
  name = 'floor/floord';
var
  i, cnt, failed: integer;
  zl: longint;
  x,y,z: double;
begin
  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do with TData[i] do begin
    x := tx;
    y := ty;
    z := floord(x);
    inc(cnt);
    if y<>z then begin
      inc(failed);
      writeln(x:19,' ',y:19, ' ',z:19,' ', z-y:18);
    end;
    if abs(y)<MaxLongint then begin
      x := tx;
      y := ty;
      zl:= floor(x);
      inc(cnt);
      if round(y)<>zl then begin
        inc(failed);
        writeln(x:19,' ',y:19, ' ',z:19);
      end;
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_rint;
const
  NT=8;
  TData: array[1..NT] of double = (0.1, 0.5, 1.5, 2.5, 3.1, 3.5, 3.6, 4.0 );
  TRx  : array[1..4]  of double = (1000000000000.5, -1000000000000.5, 1000000000001.5, -10000000000000.1);
  TR_yn: array[1..4]  of double = (1000000000000.0, -1000000000000.0, 1000000000002.0, -10000000000000.0);
const
  name = 'rint';
var
  i, cnt, failed: integer;
  x,y,z: double;
begin

  writeln('Function ',name);
  cnt :=0;
  failed := 0;
  for i:=1 to NT do begin
    x := TData[i];
    y := round(x);
    z := rint(x);
    inc(cnt);
    if y<>z then begin
      inc(failed);
      writeln(x:19:2,' ',y:19:2, ' ',z:19:2,' ', z-y:18);
    end;
    x := -TData[i];
    y := round(x);
    z := rint(x);
    inc(cnt);
    if y<>z then begin
      inc(failed);
      writeln(x:19:2,' ',y:19:2, ' ',z:19:2,' ', z-y:18);
    end;
  end;
  {assume round to nearest}
  for i:=1 to 4 do begin
    x := TRx[i];
    z := rint(x);
    y := TR_yn[i];
    inc(cnt);
    if y<>z then begin
      inc(failed);
      writeln(x:19:2,' ',y:19:2, ' ',z:19:2,' ', z-y:18);
    end;
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);

end;


{---------------------------------------------------------------------------}
procedure test_meansdev;
const
  NMAX = 20;
var
  a: array[0..NMAX] of double;
  i,j,f,cnt,failed: integer;
  m,sd,m1,sd1: double;
  q: double;
begin
  writeln('Test MeanAndStdDev/mean');
  failed := 0;
  q := 1.0;
  f := 0;
  cnt := 2;
  while (q <= 1e15) and (f=0) do begin
    a[0] := q;
    a[1] := q+1.0;
    a[2] := q+2.0;
    MeanAndStdDev(a,3,m,sd);
    if abs(m-a[1])>EPSD*a[1] then begin
      inc(failed);
      f := 1;
      writeln('Special test failed. a[1]=',a[1]:16:0,'   mean=', m);
    end;
    if abs(sd-1.0)>EPSD then begin
      inc(failed);
      f := 1;
      writeln('Special test failed. a[1]=',a[1]:16:0,'   sdev=', sd);
    end;
    m := mean(a,3);
    if abs(m-a[1])>EPSD*a[1] then begin
      inc(failed);
      f := 1;
      writeln('Special test failed. a[1]=',a[1]:16:0,'   mean=', m);
    end;
    q := q*10.0;
  end;
  q := 0.5;
  f := 0;
  inc(cnt,2);
  while (q >= 1e-12) and (f=0) do begin
    a[0] := q;
    a[1] := q+1.0;
    a[2] := q+2.0;
    MeanAndStdDev(a,3,m,sd);
    if abs(m-a[1])>EPSD*a[1] then begin
      inc(failed);
      f := 1;
      writeln('Special test failed. a[1]=',a[1]:16:0,'   mean=', m);
    end;
    if abs(sd-1.0)>EPSD then begin
      inc(failed);
      f := 1;
      writeln('Special test failed. a[1]=',a[1]:16:0,'   sdev=', sd);
    end;
    m := mean(a,3);
    if abs(m-a[1])>EPSD*a[1] then begin
      inc(failed);
      f := 1;
      writeln('Special test failed. a[1]=',a[1]:16:0,'   mean=', m);
    end;
    q := q/10.0;
  end;

  randseed := 0;
  f := 0;
  i := 1;
  inc(cnt,2);
  while (i<=1000) and (f=0) do begin
    for j:=0 to NMAX do a[j] := random - 0.5;
    m1 := sum2(a,NMAX+1)/(NMAX+1);
    sd1 := 0.0;
    for j:=0 to NMAX do sd1 := sd1 + sqr(m1-a[j]);
    sd1 := sqrt(sd1/NMAX);
    MeanAndStdDev(a,NMAX+1,m,sd);
    if abs(m-m1)>1.1*EPSD*abs(m1) then begin
      inc(failed);
      f := 1;
      writeln('Random test failed.  m1=', m1:24,'  m=', m:24);
    end;
    if abs(sd-sd1)>1.1*EPSD*sd1 then begin
      inc(failed);
      f := 1;
      writeln('Random test failed. sd1=',sd1:24,'  sd=', sd:24);
    end;
    m := mean(a,NMAX+1);
    if abs(m-m1)>1.1*EPSD*abs(m1) then begin
      inc(failed);
      f := 1;
      writeln('Random test failed.  m1=', m1:24,'  mean=', m:24);
    end;
    inc(i);
  end;
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_stat1;
const
  ad: array[1..8] of double   = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);
  yd: array[1..8] of double   = (2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0);
var
  x: double;
  cnt, failed: integer;

  function relfail(a,b: double): boolean;
  begin
    inc(cnt);
    if abs(a-b)>EPSD*abs(b) then begin
      relfail := true;
      failed  := failed;
    end
    else relfail := false;
  end;

begin
  cnt :=0;
  failed := 0;
  writeln('Functions ssdev, psdev, svar, pvar, pcov, tvar ...');
  {mean=9/2, pvar=21/4, svar=6, sdev=sqrt(6), tvar=42 (!!) }
  x := ssdev(ad,8);    if relfail(x, sqrt(6.00)) then writeln('ssdev  failed');
  x := psdev(ad,8);    if relfail(x, sqrt(5.25)) then writeln('psdev  failed');
  x := svar(ad,8);     if relfail(x, 6.00) then writeln('svar failed');
  x := pvar(ad,8);     if relfail(x, 5.25) then writeln('pvar failed');
  x := tvar(ad,8);     if relfail(x, 42)   then writeln('tvar failed');
  x := pcov(ad,yd,8);  if relfail(x, 217/16) then writeln('pcov failed');
  x := CoeffVar(ad,8); if relfail(x, sqrt(6.00)/4.5) then writeln('CoeffVar failed');
  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_stat_moment;
const
  ad: array[1..7] of double   = (55, 127, 167, 175, 151, 95, 7);
var
  skew, kurt, m1, m2, m3, m4, x, tol: double;
  k,cnt,failed: integer;
  vd: array[1..4] of double;
  s: string[20];

  function relfail(a,b: double): boolean;
  begin
    inc(cnt);
    if abs(a-b)>tol*abs(b) then begin
      relfail := true;
      failed  := failed+1;
    end
    else relfail := false;
  end;
begin
  cnt :=0;
  failed := 0;
  writeln('Functions moment/x ...');

  tol := 2*EPSD;
  moment(ad, 7, m1, m2, m3, m4, skew, kurt);
  s := 'Moment  ad ';
  if relfail(m1, 111.0)                    then writeln(s,'m1   failed');
  if relfail(m2, 3328.0)                   then writeln(s,'m2   failed');
  if relfail(m3, -114102.8571428571429)    then writeln(s,'m3   failed');
  if relfail(m4, 22303305.14285714286)     then writeln(s,'m4   failed');
  if relfail(skew, -0.5943216388127454879) then writeln(s,'skew failed');
  if relfail(kurt, 2.013736263736263736)   then writeln(s,'kurt failed');

  for k:=1 to 15 do begin
    x := exp10(k);
    str(k,s);
    if k<10 then s := '0'+s;
    s := 'Moment 10^'+s+' ';
    s := s+'vd ';
    vd[1] := x;
    vd[2] := x+1;
    vd[3] := x+2;
    vd[4] := x+2;
    moment(vd, 4, m1, m2, m3, m4, skew, kurt);
    if relfail(m1, x+1.25)           then writeln(s,'m1   failed');
    if relfail(m2, 11/16)            then writeln(s,'m2   failed');
    if relfail(m3, -9/32)            then writeln(s,'m3   failed');
    if relfail(m4,  197/256)         then writeln(s,'m4   failed');
    if relfail(skew, -18/sqrt(1331)) then writeln(s,'skew failed');
    if relfail(kurt, 197/121)        then writeln(s,'kurt failed');
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt, cnt);
  inc(total_failed, failed);
end;



{---------------------------------------------------------------------------}
procedure test_angle2;
var
  cnt, failed: integer;
  f,x,y,x1,y1,x2,y2: double;
begin
  writeln('Function angle2');
  cnt :=0;
  failed := 0;

  x := ldexpd(1,-40);
  y := angle2(x,1,-x,1);
  (* Wolfram Alpha: VectorAngle[{-(1/2)^40, 1}, {(1/2)^40, 1}] with 30 digits *)
  {f = 1.818989403545856475830077623456e-12}
  f := 1.818989403545856476e-12;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x := ldexpd(1,-30);
  y := angle2(x,1,-x,1);
  {Wolfram Alpha}
  {f = 1.862645149230957030711470955369e-9;}
  f := 1.862645149230957031e-9;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  {maple angle(vector([x1,x2]),vector([y1,y2]))}
  x1 := 453558772208.0;
  x2 := 60479044456.0;
  y1 := 2956929784.0;
  y2 := 786935157943.0;
  f  := 1.434477447476657005;
  y  := angle2(x1,x2,y1,y2);
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x1 := 914825136973.0;
  x2 := 281314862289.0;
  y1 := 454100745237.0;
  y2 := 769034410057.0;
  f  := 0.7390756038472660218;
  y  := angle2(x1,x2,y1,y2);
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x1 := 620324828055.0;
  x2 := 912854227507.0;
  y1 := 525717981351.0;
  y2 := 447473570262.0;
  f  := 0.2687604622868112296;
  y  := angle2(x1,x2,y1,y2);
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x1 := 574157513825.0;
  x2 := 150650905007.0;
  y1 := 553408748182.0;
  y2 := 82815984929.0;
  f  := 0.1080569518453247518;
  y  := angle2(x1,x2,y1,y2);
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x1 := 232408486720.0;
  x2 := 571078647666.0;
  y1 := 311705727623.0;
  y2 := 740501472935.0;
  f  := 0.1193024142403457106e-1;
  y  := angle2(x1,x2,y1,y2);
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x := ldexpd(1,30);
  y := angle2(x,1,-x,1);
  f := 3.141592651727148089;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt);
  inc(total_failed, failed);
end;


{---------------------------------------------------------------------------}
procedure test_geometry2;
var
  cnt, failed: integer;
  e,f,y,x1,y1,x2,y2,x3,y3: double;
begin
  writeln('Functions orient2d, area_triangle, in_triangle/_ex');
  cnt :=0;
  failed := 0;

  e  := ldexpd(1,-35);
  x1 := 0;
  y1 := 2-e;
  x2 := 2+e;
  y2 := 0;
  x3 := 1;
  y3 := 1;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  f  := 0.4235164736271501695e-21;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  e  := ldexpd(1,-51);
  x1 := 0;
  y1 := 2-e;
  x2 := 2+e;
  y2 := 0;
  x3 := 1;
  y3 := 1;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  f  := 0.9860761315262647567646607066e-31;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  e  := ldexpd(1,-35);
  x1 := e;
  y1 := 2+e;
  x2 := 2-e;
  y2 := -e;
  x3 := 1+e;
  y3 := 1-e;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  f  := 0.1694065894508600678e-20;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  {Maple: e:=0.5^51; x1:=e; y1:=2+e; x2:=2-e; y2:=-e; x3:=1+e; y3:=1-e; }
  {triangle(ABC, [point(A,x1,y1), point(B,x2,y2), point(C,x3,y3)]): area(ABC);}
  e  := ldexpd(1,-51);
  x1 := e;
  y1 := 2+e;
  x2 := 2-e;
  y2 := -e;
  x3 := 1+e;
  y3 := 1-e;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  f  := 0.39443045261050590270586428264e-30;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  e  := ldexpd(1,-35);
  x1 := -e;
  y1 := 2;
  x2 := 2;
  y2 := +e;
  x3 := 1;
  y3 := 1;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  f  := 0.4235164736271501695e-21;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  e  := ldexpd(1,-51);
  x1 := -e;
  y1 := 2;
  x2 := 2;
  y2 := +e;
  x3 := 1;
  y3 := 1;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  f  := 0.98607613152626475676466e-31;
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x1 := 39169594160.0;
  y1 := 88430571674.0;
  x2 := 960498834085.0;
  y2 := 812920457916.0;
  x3 := 453747019461.0;
  y3 := 644031395307.0;
  f  := 105767066422336562269341.5;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x1 := 920624947349.0;
  y1 := 951053530086.0;
  x2 := 146486307198.0;
  y2 := 155590763466.0;
  x3 := 429392673709.0;
  y3 := 525428510973.0;
  f  := 30632104966163046845368.5;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  x1 := 643842443844.0;
  y1 := 131905754565.0;
  x2 := 672075358391.0;
  y2 := 135537108795.0;
  x3 := 991638155474.0;
  y3 := 452610874309.0;
  f  := 3895735405957146910534.0;
  y  := area_triangle(x1,y1,x2,y2,x3,y3);
  inc(cnt);
  if reldev(y,f) > EPS then begin
    inc(failed);
    writeln(y:25, ' ',f:25,' ',reldev(y,f):25);
  end;

  inc(cnt,2);
  e  := ldexpd(1,-33);
  if 1 <> in_triangle_ex(1,1,0,0,e,2,2,-e)   then begin inc(failed); writeln('Error 1x');  end;
  if not in_triangle(1,1,0,0,e,2,2,-e)       then begin inc(failed); writeln('Error 1');   end;

  inc(cnt,2);
  if -1 <> in_triangle_ex(1,1,0,0,-e,2,2,0)  then begin inc(failed); writeln('Error 2x');  end;
  if in_triangle(1,1,0,0,-e,2,2,0)           then begin inc(failed); writeln('Error 2');   end;

  inc(cnt,2);
  if -1<>in_triangle_ex(1,1,0,0,0,2+e,2-e,0) then begin inc(failed); writeln('Error 3x');  end;
  if in_triangle(1,1,0,0,0,2+e,2-e,0)        then begin inc(failed); writeln('Error 3');   end;

  inc(cnt,2);
  e := predd(2.0);
  if -1 <> in_triangle_ex(1,1,0,0,0,e,e,0)   then begin inc(failed); writeln('Error 4x');  end;
  if in_triangle(1,1,0,0,0,e,e,0)            then begin inc(failed); writeln('Error 4');   end;

  inc(cnt,2);
  e := succd(2.0);
  if 1 <> in_triangle_ex(1,1,0,0,0,e,e,0)    then begin inc(failed); writeln('Error 5x');  end;
  if not in_triangle(1,1,0,0,0,e,e,0)        then begin inc(failed); writeln('Error 5');   end;

  inc(cnt,2);
  if 0 <> in_triangle_ex(1,1,0,0,0,2,2,0)    then begin inc(failed); writeln('Error 6x');  end;
  if in_triangle(1,1,0,0,0,2,2,0)            then begin inc(failed); writeln('Error 6');   end;

  inc(cnt,2);
  e := predd(1.0);
  if 1 <> in_triangle_ex(1,e,0,0,0,2,2,0)    then begin inc(failed); writeln('Error 7x');  end;
  if not in_triangle(1,e,0,0,0,2,2,0)        then begin inc(failed); writeln('Error 7');   end;

  inc(cnt,2);
  if -1 <> in_triangle_ex(0,3,0,0,0,2,2,0)   then begin inc(failed); writeln('Error 8x');  end;
  if in_triangle(0,3,0,0,0,2,2,0)            then begin inc(failed); writeln('Error 8');   end;

  inc(cnt,2);
  if -1 <> in_triangle_ex(0,3,0,0,0,2,2,0)   then begin inc(failed); writeln('Error 8x');  end;
  if in_triangle(0,3,0,0,0,2,2,0)            then begin inc(failed); writeln('Error 8');   end;

  if failed>0 then writeln('*** failed: ',failed,' of ',cnt)
  else writeln(cnt:4, ' tests OK');
  inc(total_cnt);
  inc(total_failed, failed);
end;

end.
