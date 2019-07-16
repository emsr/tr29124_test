{Test program to compare precision of DAMath/math functions with mpf routines}

{Mar. 2017: Adjusted some constants to 1e9 because FPC/ARM}
{seems to return 0 for sin/cos(x) if abs(x) > 2^30, and as}
{a consequence Nan/Inf for tan/cot/csc/sec !}


program t_xdamat;

{$i STD.INC}

{$ifdef AppCons}
  {$apptype console}
{$endif}

{.$define BIGTEST}

{$ifdef BIT16}
{$N+}
{$F+}
{$endif}

{$x+}

{$ifdef WIN32or64}
  {$define USE_MATH}
{$else}
  {$ifdef CPUARM}
    {$define USE_MATH}
  {$endif}
{$endif}

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  {$ifdef UNIT_SCOPE}
      system.math,
  {$else}
    {$ifdef USE_MATH}
      math,
    {$endif}
  {$endif}

  {$ifdef MPC_Diagnostic}
    mp_supp,
  {$endif}
  DAMath, BTypes, taus88,
  mp_types, mp_base, mp_real, FDAMSupp;


const
  PREX = 160;
  NDD  = 48;

const
  DefSeed : longint = 0;

var
  tctx: taus88_ctx;


{note: math functions with x.. because most are declared as f(const x: double)}

{$ifdef CONDITIONALEXPRESSIONS} {D6+}
{$define USE_XDMATH}
{$endif}

{$ifdef DAMOnly}
  {$undef USE_MATH}
  {$undef USE_XDMATH}
{$endif}

{$ifdef USE_MATH}
{$ifdef FPC}
  {$ifdef HAS_DENORM_LIT} {2.1+}
    {$define FPCMATH}
  {$endif}
{$endif}
{$endif}

{$ifdef VER240}
  {$define USE_XXMATH}
{$endif}

{$ifdef VER250}
  {$define USE_XXMATH}
{$endif}

function xsin(x: double): double; begin xsin := system.sin(x); end;
function xcos(x: double): double; begin xcos := system.cos(x); end;
function xexp(x: double): double; begin xexp := system.exp(x); end;
function xln(x: double): double;  begin xln  := system.ln(x);  end;
function xarctan(x: double): double; begin xarctan := system.arctan(x) end;

{$ifdef USE_XXMATH}
function xexpm1(x: double): double; begin xexpm1 := system.ExpMinus1(x); end;
{$endif}

{$ifdef UNIT_SCOPE}
  function xarccos (x: double): double; begin xarccos  := system.math.arccos (x) end;
  function xarccosh(x: double): double; begin xarccosh := system.math.arccosh(x) end;
  function xarcsin (x: double): double; begin xarcsin  := system.math.arcsin (x) end;
  function xarcsinh(x: double): double; begin xarcsinh := system.math.arcsinh(x) end;
  function xarctanh(x: double): double; begin xarctanh := system.math.arctanh(x) end;
  function xcosh   (x: double): double; begin xcosh    := system.math.cosh   (x) end;
  function xln1p   (x: double): double; begin xln1p    := system.math.lnxp1  (x) end;
  function xlog10  (x: double): double; begin xlog10   := system.math.log10  (x) end;
  function xlog2   (x: double): double; begin xlog2    := system.math.log2   (x) end;
  function xsinh   (x: double): double; begin xsinh    := system.math.sinh   (x) end;
  function xtan    (x: double): double; begin xtan     := system.math.tan    (x) end;
  function xtanh   (x: double): double; begin xtanh    := system.math.tanh   (x) end;
  function xlogn (n,x: double): double; begin xlogn    := system.math.logn(n,x)  end;
  function xpower(x,y: double): double; begin xpower   := system.math.power(x,y) end;
  function xarccot (x: double): double; begin xarccot  := system.math.arccot (x) end;
  function xarccoth(x: double): double; begin xarccoth := system.math.arccoth(x) end;
  function xarccsc (x: double): double; begin xarccsc  := system.math.arccsc (x) end;
  function xarccsch(x: double): double; begin xarccsch := system.math.arccsch(x) end;
  function xarcsec (x: double): double; begin xarcsec  := system.math.arcsec (x) end;
  function xarcsech(x: double): double; begin xarcsech := system.math.arcsech(x) end;
  function xcot    (x: double): double; begin xcot     := system.math.cot    (x) end;
  function xcoth   (x: double): double; begin xcoth    := system.math.coth   (x) end;
  function xcsc    (x: double): double; begin xcsc     := system.math.csc    (x) end;
  function xcsch   (x: double): double; begin xcsch    := system.math.csch   (x) end;
  function xsec    (x: double): double; begin xsec     := system.math.sec    (x) end;
  function xsech   (x: double): double; begin xsech    := system.math.sech   (x) end;
  function xcompound(x,n: double): double;
  var z: extended; {will be a double for 64-bit}
  begin
    z := x;
    xcompound := system.math.intpower(z+1.0, round(n));
  end;
{$else}

 {$ifdef USE_MATH}
  function xarccos (x: double): double; begin xarccos  := math.arccos (x) end;
  function xarccosh(x: double): double; begin xarccosh := math.arccosh(x) end;
  function xarcsin (x: double): double; begin xarcsin  := math.arcsin (x) end;
  function xarcsinh(x: double): double; begin xarcsinh := math.arcsinh(x) end;
  function xarctanh(x: double): double; begin xarctanh := math.arctanh(x) end;
  function xcosh   (x: double): double; begin xcosh    := math.cosh   (x) end;
  function xln1p   (x: double): double; begin xln1p    := math.lnxp1  (x) end;
  function xlog10  (x: double): double; begin xlog10   := math.log10  (x) end;
  function xlog2   (x: double): double; begin xlog2    := math.log2   (x) end;
  function xsinh   (x: double): double; begin xsinh    := math.sinh   (x) end;
  function xtan    (x: double): double; begin xtan     := math.tan    (x) end;
  function xtanh   (x: double): double; begin xtanh    := math.tanh   (x) end;
  function xlogn (n,x: double): double; begin xlogn    := math.logn(n,x) end;
  function xpower(x,y: double): double; begin xpower   := math.power(x,y) end;
 {$ifdef FPCMATH}
  function xcot    (x: double): double; begin xcot     := math.cot    (x) end;
  function xcsc    (x: double): double; begin xcsc     := math.csc    (x) end;
  function xsec    (x: double): double; begin xsec     := math.sec    (x) end;
 {$endif}
 {$endif}

 {$ifdef USE_XDMATH}
  function xarccot (x: double): double; begin xarccot  := math.arccot (x) end;
  function xarccoth(x: double): double; begin xarccoth := math.arccoth(x) end;
  function xarccsc (x: double): double; begin xarccsc  := math.arccsc (x) end;
  function xarccsch(x: double): double; begin xarccsch := math.arccsch(x) end;
  function xarcsec (x: double): double; begin xarcsec  := math.arcsec (x) end;
  function xarcsech(x: double): double; begin xarcsech := math.arcsech(x) end;
  function xcot    (x: double): double; begin xcot     := math.cot    (x) end;
  function xcoth   (x: double): double; begin xcoth    := math.coth   (x) end;
  function xcsc    (x: double): double; begin xcsc     := math.csc    (x) end;
  function xcsch   (x: double): double; begin xcsch    := math.csch   (x) end;
  function xsec    (x: double): double; begin xsec     := math.sec    (x) end;
  function xsech   (x: double): double; begin xsech    := math.sech   (x) end;
  function xcompound(x,n: double): double;
  var z: extended; {will be a double for 64-bit}
  begin
    z := x;
    xcompound := math.intpower(z+1.0, round(n));
  end;
 {$endif}


{$endif}


{$ifdef USE_MATH}
{---------------------------------------------------------------------------}
function math_powm1(x,y: double): double;
begin
  math_powm1 := xpower(x,y)-1.0;
end;

{---------------------------------------------------------------------------}
function math_pow1pm1(x,y: double): double;
begin
  math_pow1pm1 := xpower(1.0+x,y)-1.0;
end;
{$endif}


{---------------------------------------------------------------------------}
procedure sep_line;
begin
  writeln;
  writeln('-----------------------------------------------------------------');
end;


{---------------------------------------------------------------------------}
function reldev(const a,b: mp_float): double;
  {-Return abs((a-b)/b) if b<>0,  abs(a-b) if b=0, or a-b=0}
var
  c: mp_float;
begin
  reldev := MaxDouble;
  mpf_initp(c,b.bitprec);
  if mp_error<>MP_OKAY then exit;
  mpf_sub(a,b,c);
  if (mpf_todouble(a)=0) and (mpf_todouble(b)=0) then reldev := 0.0
  else if s_mpf_is0(c) then reldev := 0.0
  else begin
    if not s_mpf_is0(b) then mpf_div(c,b,c);
    reldev := abs(mpf_todouble(c));
  end;
  mpf_clear(c);
end;


type
  TXFun  = function(x: double): double;
  TMProc = procedure(const x: mp_float; var y: mp_float);

type
  TXFun2  = function(x,y: double): double;
  TMProc2 = procedure(const x,y: mp_float; var z: mp_float);

var
  a,b,c,mx,mf,mmax: mp_float;


{---------------------------------------------------------------------------}
function randx(x1,x2: double): double;
begin
  randx := x1 + taus88_double53(tctx)*(x2-x1);
end;


{---------------------------------------------------------------------------}
function fmt(x: double): Bstring;
var
  a: double;
  s: BString;
const
  AMAX = 1e16;
begin
  a := abs(x);
  if (a > AMAX) or ((a > 0.0) and (a < sqrt_epsh)) then str(x:20,s)
  else str(x:1:10,s);
  fmt := s;
end;


{---------------------------------------------------------------------------}
procedure testfun(const s: mp_string; xf: TXFun; mpf: TMProc; x1,x2: double);
var
  x,y,d,u: double;
  umax, xmax, ymax, rms: double;
  i: longint;
{$ifdef BIT16}
  const n = 200;
{$else}
  {$ifdef BIGTEST}
    const n = 10000;
  {$else}
    const n = 1000;
  {$endif}
{$endif}
begin
  writeln('Test of ',s);
  writeln('at ',n,' random values in [', fmt(x1),' .. ', fmt(x2),']');
  umax := -1;
  xmax := 0;
  ymax := 0;
  rms  := 0;
  taus88_init(tctx, DefSeed);
  for i:=1 to n do begin
    x := randx(x1,x2);
    y := xf(x);
    mpf_set_dbl(mx, x);
    mpf_set_dbl(a, y);       { a = xf(x)}
    mpf(mx,mf);              {mf = mpf(x)}
    mpf_sub_dbl(mf,y,b);     { b = mpf(x)-xf(x)}
    d := reldev(a,mf);
    u := d/eps_d;
    rms := rms + sqr(u);
    if u>umax then begin
      umax := u;
      xmax := x;
      ymax := y;
      mpf_copy(mf,mmax);
    end;
  end;
  writeln('RMS = ', sqrt(rms/n):1:2,', max rel = ',umax:1:2,' eps at ');
  writeln(' x(dbl) = ',xmax:26, ' = $',dbl2hex(xmax));
  writeln(' y(dbl) = ',ymax:26, ' = $',dbl2hex(ymax));
  writeln(' y(mpf) = ',mpf_decimal(mmax,NDD));
{$ifdef debug}
  mpf_set_dbl(mx, xmax);
  mpf_mul_2k(mx,64,a);
  writeln('   xmax = 2^-64 * ', mpf_decimal(a,NDD));
{$endif}
  writeln;
end;


{---------------------------------------------------------------------------}
procedure testpow(const s: mp_string; mft: integer; xf: TXFun2; x1,x2,y1,y2: double);
var
  x,f,y,d,u: double;
  umax, xmax, zmax, ymax, rms, lnmax: double;
  i: longint;
{$ifdef BIT16}
  const n = 400;
{$else}
  {$ifdef BIGTEST}
    const n = 40000;
  {$else}
    const n = 2000;
  {$endif}
{$endif}
{$define use_eps}
begin
  case mft of
      1: writeln('Test of ',s, ', z = x^y-1  at  ',n,' random values');
      2: writeln('Test of ',s, ', z = (1+x)^y-1  at  ',n,' random values');
      3: writeln('Test of ',s, ', z = (1+x)^y  at  ',n,' random values');
      4: writeln('Test of ',s, ', z = compound(x,y)  at  ',n,' random values');
      5: writeln('Test of ',s, ', z = (1+x)^y (fast)  at  ',n,' random values');
      6: writeln('Test of ',s, ', z = comprel(x,y)  ',n,' random values');
    else writeln('Test of ',s, ', z = x^y  at  ',n,' random values');
  end;
  writeln(' with x from [', fmt(x1), ' .. ', fmt(x2),']');
  writeln('  and y from [', fmt(y1), ' .. ', fmt(y2),']');
  umax := -1;
  xmax := 0;
  zmax := 0;
  ymax := 0;
  rms  := 0;
  if mft>3 then lnmax := 680 {<ln_MaxDbl: avoid overflow with *CSD splitting}
  else lnmax := ln_MaxDbl;

  taus88_init(tctx, DefSeed);
  for i:=1 to n do begin
    x := randx(x1,x2);
    if mft<2 then u := lnMax/abs(ln(x))
    else u := lnMax/abs(ln(1+x));
    repeat
       y := randx(y1,y2);
       if (mft=4) or (mft=6) then y := int(y);
    until abs(y)<u;
    mpf_set_dbl(mx, x);
    mpf_set_dbl(a, y);
    case mft of
        1: mpf_exptm1(mx,a,mf);       {mf := x^y-1}
        2: mpf_expt1pm1(mx,a,mf);     {mf := (1+x)^y-1}
        3: mpf_expt1p(mx,a,mf);       {mf := (1+x)^y}
        4: mpf_expt1p(mx,a,mf);       {mf := (1+x)^y}
        5: mpf_expt1p(mx,a,mf);       {mf := (1+x)^y}
        6: mpf_comprel(mx,a,mf);      {mf := ((1+x)^y-1)/x}
      else mpf_expt(mx,a,mf);         {mf := x^y}
    end;
    f := xf(x,y);
    mpf_set_dbl(a, f);
    d := reldev(a,mf);
    {$ifdef use_eps}
      u := d/eps_d;
    {$else}
      u := d;
    {$endif}
    rms := rms + sqr(u);
    if u>umax then begin
      umax := u;
      xmax := x;
      zmax := f;
      ymax := y;
      mpf_copy(mf,mmax);
    end;
  end;
  {$ifdef use_eps}
    writeln('RMS = ', sqrt(rms/n):1:2,', max rel = ',umax:1:2,' eps at');
  {$else}
    writeln('RMS = ', sqrt(rms/n):26,', max = ',umax:26, '  at');
  {$endif}
  writeln(' x(dbl) = ',xmax:26, ' = $',dbl2hex(xmax));
  writeln(' y(dbl) = ',ymax:26, ' = $',dbl2hex(ymax));
  writeln(' z(dbl) = ',zmax:26, ' = $',dbl2hex(zmax));
  writeln(' z(mpf) = ',mpf_decimal(mmax,NDD));
{$ifdef debug}
  mpf_set_dbl(mx, xmax);
  mpf_mul_2k(mx,64,a);
  writeln('   xmax = 2^-64 * ', mpf_decimal(a,NDD));
{$endif}
  writeln;
end;


{---------------------------------------------------------------------------}
procedure testlogxy(const s: mp_string; xf: TXFun2; x1,x2,y1,y2: double);
var
  x,f,y,d,u: double;
  umax, xmax, zmax, ymax, rms: double;
  i: longint;
{$ifdef BIT16}
  const n = 400;
{$else}
  {$ifdef BIGTEST}
    const n = 40000;
  {$else}
    const n = 2000;
  {$endif}
{$endif}
{$define use_eps}
begin
  writeln('Test of ',s, ', z = log_y(x)  at  ',n,' random values');
  writeln(' with x from [', fmt(x1),' .. ', fmt(x2),']');
  writeln('  and y from [', fmt(y1),' .. ', fmt(y2),']');
  umax := -1;
  xmax := 0;
  zmax := 0;
  ymax := 0;
  rms  := 0;
  taus88_init(tctx, DefSeed);
  for i:=1 to n do begin
    x := randx(x1,x2);
    y := randx(y1,y2);
    mpf_set_dbl(mx, x);
    mpf_set_dbl(a, y);
    mpf_logbase(a,mx,mf);
    f := xf(y,x);            { f := logn(y,x)}
    mpf_set_dbl(a, f);
    d := reldev(a,mf);
    {$ifdef use_eps}
      u := d/eps_d;
    {$else}
      u := d;
    {$endif}
    rms := rms + sqr(u);
    if u>umax then begin
      umax := u;
      xmax := x;
      zmax := f;
      ymax := y;
      mpf_copy(mf,mmax);
    end;
  end;
  {$ifdef use_eps}
    writeln('RMS = ', sqrt(rms/n):1:2,', max rel = ',umax:1:2,' eps at');
  {$else}
    writeln('RMS = ', sqrt(rms/n):26,', max = ',umax:26, '  at');
  {$endif}
  writeln(' x(dbl) = ',xmax:26, ' = $',dbl2hex(xmax));
  writeln(' y(dbl) = ',ymax:26, ' = $',dbl2hex(ymax));
  writeln(' z(dbl) = ',zmax:26, ' = $',dbl2hex(zmax));
  writeln(' z(mpf) = ',mpf_decimal(mmax,NDD));
{$ifdef debug}
  mpf_set_dbl(mx, xmax);
  mpf_mul_2k(mx,64,a);
  writeln('   xmax = 2^-64 * ', mpf_decimal(a,NDD));
{$endif}
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_arccos;
begin
  sep_line;
  testfun('DAMath.arccos',{$ifdef FPC_ProcVar}@{$endif}arccos, {$ifdef FPC_ProcVar}@{$endif}mpf_arccos, -1, 1);
{$ifdef USE_MATH}
  testfun('math.arccos',{$ifdef FPC_ProcVar}@{$endif}xarccos, {$ifdef FPC_ProcVar}@{$endif}mpf_arccos, -1, 1);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arccosh;
begin
  sep_line;
  testfun('DAMath.arccosh',{$ifdef FPC_ProcVar}@{$endif}arccosh, {$ifdef FPC_ProcVar}@{$endif}mpf_arccosh, 1, 10);
  testfun('DAMath.arccosh',{$ifdef FPC_ProcVar}@{$endif}arccosh, {$ifdef FPC_ProcVar}@{$endif}mpf_arccosh, 10, 1E5);
{$ifdef USE_MATH}
  testfun('math.arccosh',{$ifdef FPC_ProcVar}@{$endif}xarccosh, {$ifdef FPC_ProcVar}@{$endif}mpf_arccosh, 1, 10);
  testfun('math.arccosh',{$ifdef FPC_ProcVar}@{$endif}xarccosh, {$ifdef FPC_ProcVar}@{$endif}mpf_arccosh, 10, 1E5);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arccosh1p;
begin
  sep_line;
  testfun('DAMath.arccosh1p',{$ifdef FPC_ProcVar}@{$endif}arccosh1p, {$ifdef FPC_ProcVar}@{$endif}mpf_arccosh1p, 0, 1);
  testfun('DAMath.arccosh1p',{$ifdef FPC_ProcVar}@{$endif}arccosh1p, {$ifdef FPC_ProcVar}@{$endif}mpf_arccosh1p, 1, 1E5);
end;


{---------------------------------------------------------------------------}
procedure test_arccot;
begin
  sep_line;
  testfun('DAMath.arccot',{$ifdef FPC_ProcVar}@{$endif}arccot, {$ifdef FPC_ProcVar}@{$endif}mpf_arccot, -1, -1e-10);
  testfun('DAMath.arccot',{$ifdef FPC_ProcVar}@{$endif}arccot, {$ifdef FPC_ProcVar}@{$endif}mpf_arccot, 1, 1e5);
{$ifdef USE_XDMATH}
  testfun('math.arccot',{$ifdef FPC_ProcVar}@{$endif}xarccot, {$ifdef FPC_ProcVar}@{$endif}mpf_arccot, -1, -1e-10);
  testfun('math.arccot',{$ifdef FPC_ProcVar}@{$endif}xarccot, {$ifdef FPC_ProcVar}@{$endif}mpf_arccot, 1, 1e5);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arccotc;
begin
  sep_line;
  testfun('DAMath.arccotc',{$ifdef FPC_ProcVar}@{$endif}arccotc, {$ifdef FPC_ProcVar}@{$endif}mpf_arccotc, -1e5, -1);
  testfun('DAMath.arccotc',{$ifdef FPC_ProcVar}@{$endif}arccotc, {$ifdef FPC_ProcVar}@{$endif}mpf_arccotc, -1, 1);
  testfun('DAMath.arccotc',{$ifdef FPC_ProcVar}@{$endif}arccotc, {$ifdef FPC_ProcVar}@{$endif}mpf_arccotc, 1, 1e5);
end;


{---------------------------------------------------------------------------}
procedure test_arccoth;
begin
  sep_line;
  testfun('DAMath.arccoth',{$ifdef FPC_ProcVar}@{$endif}arccoth, {$ifdef FPC_ProcVar}@{$endif}mpf_arccoth, 1+1e-10, 2);
  testfun('DAMath.arccoth',{$ifdef FPC_ProcVar}@{$endif}arccoth, {$ifdef FPC_ProcVar}@{$endif}mpf_arccoth, 2, 1e5);
{$ifdef USE_XDMATH}
  testfun('math.arccoth',{$ifdef FPC_ProcVar}@{$endif}xarccoth, {$ifdef FPC_ProcVar}@{$endif}mpf_arccoth, 1+1e-10, 2);
  testfun('math.arccoth',{$ifdef FPC_ProcVar}@{$endif}xarccoth, {$ifdef FPC_ProcVar}@{$endif}mpf_arccoth, 2, 1e5);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arccsc;
begin
  sep_line;
  testfun('DAMath.arccsc',{$ifdef FPC_ProcVar}@{$endif}arccsc, {$ifdef FPC_ProcVar}@{$endif}mpf_arccsc, 1, 2);
  testfun('DAMath.arccsc',{$ifdef FPC_ProcVar}@{$endif}arccsc, {$ifdef FPC_ProcVar}@{$endif}mpf_arccsc, 2, 1e5);
{$ifdef USE_XDMATH}
  testfun('math.arccsc',{$ifdef FPC_ProcVar}@{$endif}xarccsc, {$ifdef FPC_ProcVar}@{$endif}mpf_arccsc, 1, 2);
  testfun('math.arccsc',{$ifdef FPC_ProcVar}@{$endif}xarccsc, {$ifdef FPC_ProcVar}@{$endif}mpf_arccsc, 2, 1e5);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arccsch;
begin
  sep_line;
  testfun('DAMath.arccsch',{$ifdef FPC_ProcVar}@{$endif}arccsch, {$ifdef FPC_ProcVar}@{$endif}mpf_arccsch, 1e-10, 1);
  testfun('DAMath.arccsch',{$ifdef FPC_ProcVar}@{$endif}arccsch, {$ifdef FPC_ProcVar}@{$endif}mpf_arccsch, 1, 1e4);
{$ifdef USE_XDMATH}
  testfun('math.arccsch',{$ifdef FPC_ProcVar}@{$endif}xarccsch, {$ifdef FPC_ProcVar}@{$endif}mpf_arccsch, 1e-10, 1);
  testfun('math.arccsch',{$ifdef FPC_ProcVar}@{$endif}xarccsch, {$ifdef FPC_ProcVar}@{$endif}mpf_arccsch, 1, 1e4);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arcsec;
begin
  sep_line;
  testfun('DAMath.arcsec',{$ifdef FPC_ProcVar}@{$endif}arcsec, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsec, -2, -1);
  testfun('DAMath.arcsec',{$ifdef FPC_ProcVar}@{$endif}arcsec, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsec, 2, 1e5);
{$ifdef USE_XDMATH}
  testfun('math.arcsec',{$ifdef FPC_ProcVar}@{$endif}xarcsec, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsec, -2, -1);
  testfun('math.arcsec',{$ifdef FPC_ProcVar}@{$endif}xarcsec, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsec, 2, 1e5);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arcsech;
begin
  sep_line;
  testfun('DAMath.arcsech',{$ifdef FPC_ProcVar}@{$endif}arcsech, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsech, 1e-10, 1);
{$ifdef USE_XDMATH}
  testfun('math.arcsech',{$ifdef FPC_ProcVar}@{$endif}xarcsech, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsech, 1e-10, 1);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arcsin;
begin
  sep_line;
  testfun('DAMath.arcsin',{$ifdef FPC_ProcVar}@{$endif}arcsin, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsin, -1, 1);
{$ifdef USE_MATH}
  testfun('math.arcsin',{$ifdef FPC_ProcVar}@{$endif}xarcsin, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsin,-1, 1);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arcsinh;
begin
  sep_line;
  testfun('DAMath.arcsinh',{$ifdef FPC_ProcVar}@{$endif}arcsinh, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsinh, -1, 0);
  testfun('DAMath.arcsinh',{$ifdef FPC_ProcVar}@{$endif}arcsinh, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsinh, 1, 1e4);
{$ifdef USE_MATH}
  testfun('math.arcsinh',{$ifdef FPC_ProcVar}@{$endif}xarcsinh, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsinh, -1, 0);
  testfun('math.arcsinh',{$ifdef FPC_ProcVar}@{$endif}xarcsinh, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsinh, 1, 1e4);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arctanh;
begin
  sep_line;
  testfun('DAMath.arctanh',{$ifdef FPC_ProcVar}@{$endif}arctanh, {$ifdef FPC_ProcVar}@{$endif}mpf_arctanh, -0.5, 0);
  testfun('DAMath.arctanh',{$ifdef FPC_ProcVar}@{$endif}arctanh, {$ifdef FPC_ProcVar}@{$endif}mpf_arctanh, 0.5, 0.9999999999);
{$ifdef USE_MATH}
  testfun('math.arctanh',{$ifdef FPC_ProcVar}@{$endif}xarctanh, {$ifdef FPC_ProcVar}@{$endif}mpf_arctanh, -0.5, 0);
  testfun('math.arctanh',{$ifdef FPC_ProcVar}@{$endif}xarctanh, {$ifdef FPC_ProcVar}@{$endif}mpf_arctanh,0.5,0.9999999999);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_cosh;
begin
  sep_line;
  testfun('DAMath.cosh',{$ifdef FPC_ProcVar}@{$endif}cosh, {$ifdef FPC_ProcVar}@{$endif}mpf_cosh, -1, 0);
  testfun('DAMath.cosh',{$ifdef FPC_ProcVar}@{$endif}cosh, {$ifdef FPC_ProcVar}@{$endif}mpf_cosh, 1, 700);
{$ifdef USE_MATH}
  testfun('math.cosh',{$ifdef FPC_ProcVar}@{$endif}xcosh, {$ifdef FPC_ProcVar}@{$endif}mpf_cosh, -1, 0);
  testfun('math.cosh',{$ifdef FPC_ProcVar}@{$endif}xcosh, {$ifdef FPC_ProcVar}@{$endif}mpf_cosh, 1, 700);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_cot;
begin
  sep_line;
  testfun('DAMath.cot',{$ifdef FPC_ProcVar}@{$endif}cot, {$ifdef FPC_ProcVar}@{$endif}mpf_cot, 1e-10, 10);
  testfun('DAMath.cot',{$ifdef FPC_ProcVar}@{$endif}cot, {$ifdef FPC_ProcVar}@{$endif}mpf_cot, 1e-10, 1e9);
  testfun('DAMath.cot',{$ifdef FPC_ProcVar}@{$endif}cot, {$ifdef FPC_ProcVar}@{$endif}mpf_cot, 1e9, 5e18);;
{$ifdef FPCMATH}
  testfun('math.cot',{$ifdef FPC_ProcVar}@{$endif}xcot, {$ifdef FPC_ProcVar}@{$endif}mpf_cot,1e-10, 10);
  testfun('math.cot',{$ifdef FPC_ProcVar}@{$endif}xcot, {$ifdef FPC_ProcVar}@{$endif}mpf_cot,1e-10, 1e9);
{$endif}
{$ifdef USE_XDMATH}
  testfun('math.cot',{$ifdef FPC_ProcVar}@{$endif}xcot, {$ifdef FPC_ProcVar}@{$endif}mpf_cot,1e-10, 10);
  testfun('math.cot',{$ifdef FPC_ProcVar}@{$endif}xcot, {$ifdef FPC_ProcVar}@{$endif}mpf_cot,1e-10, 1e9);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_coth;
begin
  sep_line;
  testfun('DAMath.coth',{$ifdef FPC_ProcVar}@{$endif}coth, {$ifdef FPC_ProcVar}@{$endif}mpf_coth, -1, 0);
  testfun('DAMath.coth',{$ifdef FPC_ProcVar}@{$endif}coth, {$ifdef FPC_ProcVar}@{$endif}mpf_coth, 1, 700);
{$ifdef USE_XDMATH}
  testfun('math.coth',{$ifdef FPC_ProcVar}@{$endif}xcoth, {$ifdef FPC_ProcVar}@{$endif}mpf_coth, -1, 0);
  testfun('math.coth',{$ifdef FPC_ProcVar}@{$endif}xcoth, {$ifdef FPC_ProcVar}@{$endif}mpf_coth, 1, 700);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_csc;
begin
  sep_line;
  testfun('DAMath.csc',{$ifdef FPC_ProcVar}@{$endif}csc, {$ifdef FPC_ProcVar}@{$endif}mpf_csc, 1e-10,10);
  testfun('DAMath.csc',{$ifdef FPC_ProcVar}@{$endif}csc, {$ifdef FPC_ProcVar}@{$endif}mpf_csc, 1e-10,1e9);
{$ifdef FPCMATH}
  testfun('math.csc',{$ifdef FPC_ProcVar}@{$endif}xcsc, {$ifdef FPC_ProcVar}@{$endif}mpf_csc, 1e-10,10);
  testfun('math.csc',{$ifdef FPC_ProcVar}@{$endif}xcsc, {$ifdef FPC_ProcVar}@{$endif}mpf_csc, 1e-10,1e9);
{$endif}
{$ifdef USE_XDMATH}
  testfun('math.csc',{$ifdef FPC_ProcVar}@{$endif}xcsc, {$ifdef FPC_ProcVar}@{$endif}mpf_csc, 1e-10,10);
  testfun('math.csc',{$ifdef FPC_ProcVar}@{$endif}xcsc, {$ifdef FPC_ProcVar}@{$endif}mpf_csc, 1e-10,1e9);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_csch;
begin
  sep_line;
  testfun('DAMath.csch',{$ifdef FPC_ProcVar}@{$endif}csch, {$ifdef FPC_ProcVar}@{$endif}mpf_csch, -1, -1e-10);
  testfun('DAMath.csch',{$ifdef FPC_ProcVar}@{$endif}csch, {$ifdef FPC_ProcVar}@{$endif}mpf_csch, 1, 700);
{$ifdef USE_XDMATH}
  testfun('math.csch',{$ifdef FPC_ProcVar}@{$endif}xcsch, {$ifdef FPC_ProcVar}@{$endif}mpf_csch, -1, -1e-10);
  testfun('math.csch',{$ifdef FPC_ProcVar}@{$endif}xcsch, {$ifdef FPC_ProcVar}@{$endif}mpf_csch, 1, 700);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_exp;
begin
  sep_line;
  testfun('DAMath.exp',{$ifdef FPC_ProcVar}@{$endif}exp, {$ifdef FPC_ProcVar}@{$endif}mpf_exp, -100, 100);
  testfun('DAMath.exp',{$ifdef FPC_ProcVar}@{$endif}exp, {$ifdef FPC_ProcVar}@{$endif}mpf_exp, 100, 700);
{$ifndef DAMOnly}
  testfun('system.exp',{$ifdef FPC_ProcVar}@{$endif}xexp, {$ifdef FPC_ProcVar}@{$endif}mpf_exp,-100, 100);
  testfun('system.exp',{$ifdef FPC_ProcVar}@{$endif}xexp, {$ifdef FPC_ProcVar}@{$endif}mpf_exp, 100, 700);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_arctan;
begin
  sep_line;
  testfun('DAMath.arctan',{$ifdef FPC_ProcVar}@{$endif}arctan, {$ifdef FPC_ProcVar}@{$endif}mpf_arctan, -1e-8, 0);
  testfun('DAMath.arctan',{$ifdef FPC_ProcVar}@{$endif}arctan, {$ifdef FPC_ProcVar}@{$endif}mpf_arctan, 0, 1e9);
  testfun('DAMath.arctan',{$ifdef FPC_ProcVar}@{$endif}arctan, {$ifdef FPC_ProcVar}@{$endif}mpf_arctan, 1e9, 1e17);
{$ifndef DAMOnly}
  testfun('system.arctan',{$ifdef FPC_ProcVar}@{$endif}xarctan, {$ifdef FPC_ProcVar}@{$endif}mpf_arctan, -1e-8, 0);
  testfun('system.arctan',{$ifdef FPC_ProcVar}@{$endif}xarctan, {$ifdef FPC_ProcVar}@{$endif}mpf_arctan, 0, 1e9);
  testfun('system.arctan',{$ifdef FPC_ProcVar}@{$endif}xarctan, {$ifdef FPC_ProcVar}@{$endif}mpf_arctan, 1e9, 1e17);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_exp10;
begin
  sep_line;
  testfun('DAMath.exp10',{$ifdef FPC_ProcVar}@{$endif}exp10, {$ifdef FPC_ProcVar}@{$endif}mpf_exp10, -50, 50);
  testfun('DAMath.exp10',{$ifdef FPC_ProcVar}@{$endif}exp10, {$ifdef FPC_ProcVar}@{$endif}mpf_exp10, 50, 305);
end;


{---------------------------------------------------------------------------}
procedure test_exp2;
begin
  sep_line;
  testfun('DAMath.exp2',{$ifdef FPC_ProcVar}@{$endif}exp2, {$ifdef FPC_ProcVar}@{$endif}mpf_exp2, -100, 100);
  testfun('DAMath.exp2',{$ifdef FPC_ProcVar}@{$endif}exp2, {$ifdef FPC_ProcVar}@{$endif}mpf_exp2, 100, 1000);
end;


{---------------------------------------------------------------------------}
procedure test_expm1;
begin
  sep_line;
  testfun('DAMath.expm1',{$ifdef FPC_ProcVar}@{$endif}expm1, {$ifdef FPC_ProcVar}@{$endif}mpf_expm1, -1e-4, 1e-4);
  testfun('DAMath.expm1',{$ifdef FPC_ProcVar}@{$endif}expm1, {$ifdef FPC_ProcVar}@{$endif}mpf_expm1, -700, 700);
{$ifdef USE_XXMATH}
  testfun('system.ExpMinus1',{$ifdef FPC_ProcVar}@{$endif}xexpm1, {$ifdef FPC_ProcVar}@{$endif}mpf_expm1, -1e-4, 1e-4);
  testfun('system.ExpMinus1',{$ifdef FPC_ProcVar}@{$endif}xexpm1, {$ifdef FPC_ProcVar}@{$endif}mpf_expm1, -700, 700);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_exp2m1;
begin
  sep_line;
  testfun('DAMath.exp2m1',{$ifdef FPC_ProcVar}@{$endif}exp2m1, {$ifdef FPC_ProcVar}@{$endif}mpf_exp2m1, -1e-4, 1e-4);
  testfun('DAMath.exp2m1',{$ifdef FPC_ProcVar}@{$endif}exp2m1, {$ifdef FPC_ProcVar}@{$endif}mpf_exp2m1, -1000, 1000);
end;


{---------------------------------------------------------------------------}
procedure test_ln;
begin
  sep_line;
  testfun('DAMath.ln',{$ifdef FPC_ProcVar}@{$endif}ln, {$ifdef FPC_ProcVar}@{$endif}mpf_ln, 1e-300,1);
  testfun('DAMath.ln',{$ifdef FPC_ProcVar}@{$endif}ln, {$ifdef FPC_ProcVar}@{$endif}mpf_ln, 1,1e300);
{$ifdef USE_MATH}
  testfun('math.ln',{$ifdef FPC_ProcVar}@{$endif}xln, {$ifdef FPC_ProcVar}@{$endif}mpf_ln, 1e-300,1);
  testfun('math.ln',{$ifdef FPC_ProcVar}@{$endif}xln, {$ifdef FPC_ProcVar}@{$endif}mpf_ln, 1,1e300);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_ln1p;
begin
  sep_line;
  testfun('DAMath.ln1p',{$ifdef FPC_ProcVar}@{$endif}ln1p, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1p, -1e-7, 1e-7);
  testfun('DAMath.ln1p',{$ifdef FPC_ProcVar}@{$endif}ln1p, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1p, -0.9999999, 10);
{$ifdef USE_MATH}
  testfun('math.lnxp1',{$ifdef FPC_ProcVar}@{$endif}xln1p, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1p, -1e-7, 1e-7);
  testfun('math.lnxp1',{$ifdef FPC_ProcVar}@{$endif}xln1p, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1p, -0.9999999, 10);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_log10;
begin
  sep_line;
  testfun('DAMath.log10',{$ifdef FPC_ProcVar}@{$endif}log10, {$ifdef FPC_ProcVar}@{$endif}mpf_log10, 1e-300, 0.5);
  testfun('DAMath.log10',{$ifdef FPC_ProcVar}@{$endif}log10, {$ifdef FPC_ProcVar}@{$endif}mpf_log10, 0.5, 1e300);
{$ifdef USE_MATH}
  testfun('math.log10',{$ifdef FPC_ProcVar}@{$endif}xlog10, {$ifdef FPC_ProcVar}@{$endif}mpf_log10, 1e-300, 0.5);
  testfun('math.log10',{$ifdef FPC_ProcVar}@{$endif}xlog10, {$ifdef FPC_ProcVar}@{$endif}mpf_log10, 0.5, 1e300);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_log2;
begin
  sep_line;
  testfun('DAMath.log2',{$ifdef FPC_ProcVar}@{$endif}log2, {$ifdef FPC_ProcVar}@{$endif}mpf_log2, 1e-300, 0.5);
  testfun('DAMath.log2',{$ifdef FPC_ProcVar}@{$endif}log2, {$ifdef FPC_ProcVar}@{$endif}mpf_log2, 0.5, 1e300);
{$ifdef USE_MATH}
  testfun('math.log2',{$ifdef FPC_ProcVar}@{$endif}xlog2, {$ifdef FPC_ProcVar}@{$endif}mpf_log2, 1e-300, 0.5);
  testfun('math.log2',{$ifdef FPC_ProcVar}@{$endif}xlog2, {$ifdef FPC_ProcVar}@{$endif}mpf_log2, 0.5, 1e300);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_sec;
begin
  sep_line;
  testfun('DAMath.sec',{$ifdef FPC_ProcVar}@{$endif}sec, {$ifdef FPC_ProcVar}@{$endif}mpf_sec, -10, 10);
  testfun('DAMath.sec',{$ifdef FPC_ProcVar}@{$endif}sec, {$ifdef FPC_ProcVar}@{$endif}mpf_sec, 0, 1e9);
{$ifdef FPCMATH}
  testfun('math.sec',{$ifdef FPC_ProcVar}@{$endif}xsec, {$ifdef FPC_ProcVar}@{$endif}mpf_sec, -10, 10);
  testfun('math.sec',{$ifdef FPC_ProcVar}@{$endif}xsec, {$ifdef FPC_ProcVar}@{$endif}mpf_sec, 0, 1e9);
{$endif}
{$ifdef USE_XDMATH}
  testfun('math.sec',{$ifdef FPC_ProcVar}@{$endif}xsec, {$ifdef FPC_ProcVar}@{$endif}mpf_sec, -10, 10);
  testfun('math.sec',{$ifdef FPC_ProcVar}@{$endif}xsec, {$ifdef FPC_ProcVar}@{$endif}mpf_sec, 0, 1e9);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_sech;
begin
  sep_line;
  testfun('DAMath.sech',{$ifdef FPC_ProcVar}@{$endif}sech, {$ifdef FPC_ProcVar}@{$endif}mpf_sech, -2, 0);
  testfun('DAMath.sech',{$ifdef FPC_ProcVar}@{$endif}sech, {$ifdef FPC_ProcVar}@{$endif}mpf_sech, 2, 700);
{$ifdef USE_XDMATH}
  {D17_64 crashes for x>lnMaxDbl!!!!}
  testfun('math.sech',{$ifdef FPC_ProcVar}@{$endif}xsech, {$ifdef FPC_ProcVar}@{$endif}mpf_sech, -2, 0);
  testfun('math.sech',{$ifdef FPC_ProcVar}@{$endif}xsech, {$ifdef FPC_ProcVar}@{$endif}mpf_sech, 2, 700);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_sinh;
begin
  sep_line;
  testfun('DAMath.sinh',{$ifdef FPC_ProcVar}@{$endif}sinh, {$ifdef FPC_ProcVar}@{$endif}mpf_sinh, -1, 0);
  testfun('DAMath.sinh',{$ifdef FPC_ProcVar}@{$endif}sinh, {$ifdef FPC_ProcVar}@{$endif}mpf_sinh, 1, 700);
{$ifdef USE_MATH}
  testfun('math.sinh',{$ifdef FPC_ProcVar}@{$endif}xsinh, {$ifdef FPC_ProcVar}@{$endif}mpf_sinh, 0, 1);
  testfun('math.sinh',{$ifdef FPC_ProcVar}@{$endif}xsinh, {$ifdef FPC_ProcVar}@{$endif}mpf_sinh, 1, 700);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_sinhc;
begin
  sep_line;
  testfun('DAMath.sinhc',{$ifdef FPC_ProcVar}@{$endif}sinhc, {$ifdef FPC_ProcVar}@{$endif}mpf_sinhc, 0.0, 1e-3);
  testfun('DAMath.sinhc',{$ifdef FPC_ProcVar}@{$endif}sinhc, {$ifdef FPC_ProcVar}@{$endif}mpf_sinhc, 1e-3, 700);
end;


{---------------------------------------------------------------------------}
procedure test_sinhmx;
begin
  sep_line;
  testfun('DAMath.sinhmx',{$ifdef FPC_ProcVar}@{$endif}sinhmx, {$ifdef FPC_ProcVar}@{$endif}mpf_sinhmx, -2, 2);
  testfun('DAMath.sinhmx',{$ifdef FPC_ProcVar}@{$endif}sinhmx, {$ifdef FPC_ProcVar}@{$endif}mpf_sinhmx, 2, 200);
end;


{---------------------------------------------------------------------------}
procedure test_tan;
begin
  sep_line;
  testfun('DAMath.tan',{$ifdef FPC_ProcVar}@{$endif}tan, {$ifdef FPC_ProcVar}@{$endif}mpf_tan, -10, 10);
  testfun('DAMath.tan',{$ifdef FPC_ProcVar}@{$endif}tan, {$ifdef FPC_ProcVar}@{$endif}mpf_tan, 0, 1e9);
  testfun('DAMath.tan',{$ifdef FPC_ProcVar}@{$endif}tan, {$ifdef FPC_ProcVar}@{$endif}mpf_tan, 1e9, 5e18);;
{$ifdef USE_MATH}
  testfun('math.tan',{$ifdef FPC_ProcVar}@{$endif}xtan, {$ifdef FPC_ProcVar}@{$endif}mpf_tan, -10, 10);
  testfun('math.tan',{$ifdef FPC_ProcVar}@{$endif}xtan, {$ifdef FPC_ProcVar}@{$endif}mpf_tan, 0, 1e9);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_tanh;
begin
  sep_line;
  testfun('DAMath.tanh',{$ifdef FPC_ProcVar}@{$endif}tanh, {$ifdef FPC_ProcVar}@{$endif}mpf_tanh, 0, 2);
  testfun('DAMath.tanh',{$ifdef FPC_ProcVar}@{$endif}tanh, {$ifdef FPC_ProcVar}@{$endif}mpf_tanh, 2, 22);
{$ifdef USE_MATH}
  testfun('math.tanh',{$ifdef FPC_ProcVar}@{$endif}xtanh, {$ifdef FPC_ProcVar}@{$endif}mpf_tanh, 0, 2);
  testfun('math.tanh',{$ifdef FPC_ProcVar}@{$endif}xtanh, {$ifdef FPC_ProcVar}@{$endif}mpf_tanh, 2, 22);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_cos;
begin
  sep_line;
  testfun('DAMath.cos',{$ifdef FPC_ProcVar}@{$endif}DAMath.cos, {$ifdef FPC_ProcVar}@{$endif}mpf_cos,-10, 10);
  testfun('DAMath.cos',{$ifdef FPC_ProcVar}@{$endif}DAMath.cos, {$ifdef FPC_ProcVar}@{$endif}mpf_cos, 0, 1e9);
  testfun('DAMath.cos',{$ifdef FPC_ProcVar}@{$endif}DAMath.cos, {$ifdef FPC_ProcVar}@{$endif}mpf_cos,1e9,5e18);
{$ifndef DAMOnly}
  testfun('system.cos',{$ifdef FPC_ProcVar}@{$endif}xcos, {$ifdef FPC_ProcVar}@{$endif}mpf_cos,-10, 10);
  testfun('system.cos',{$ifdef FPC_ProcVar}@{$endif}xcos, {$ifdef FPC_ProcVar}@{$endif}mpf_cos,0, 1e9);
 {$ifndef CPUARM}
  testfun('system.cos',{$ifdef FPC_ProcVar}@{$endif}xcos, {$ifdef FPC_ProcVar}@{$endif}mpf_cos,1e9,5e18);;
 {$endif}
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_sin;
begin
  sep_line;
  testfun('DAMath.sin',{$ifdef FPC_ProcVar}@{$endif}DAMath.sin, {$ifdef FPC_ProcVar}@{$endif}mpf_sin,-10, 10);
  testfun('DAMath.sin',{$ifdef FPC_ProcVar}@{$endif}DAMath.sin, {$ifdef FPC_ProcVar}@{$endif}mpf_sin, 0, 1e9);
  testfun('DAMath.sin',{$ifdef FPC_ProcVar}@{$endif}DAMath.sin, {$ifdef FPC_ProcVar}@{$endif}mpf_sin,1e9,5e18);
{$ifndef DAMOnly}
  testfun('system.sin',{$ifdef FPC_ProcVar}@{$endif}xsin, {$ifdef FPC_ProcVar}@{$endif}mpf_sin,-10, 10);
  testfun('system.sin',{$ifdef FPC_ProcVar}@{$endif}xsin, {$ifdef FPC_ProcVar}@{$endif}mpf_sin,0, 1e9);
 {$ifndef CPUARM}
  testfun('system.sin',{$ifdef FPC_ProcVar}@{$endif}xsin, {$ifdef FPC_ProcVar}@{$endif}mpf_sin,1e9,5e18);
 {$endif}
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_power;
begin
  sep_line;
  testpow('DAMath.power',0,{$ifdef FPC_ProcVar}@{$endif}DAMath.power, 1e-3, 100, -100, 100);
  testpow('DAMath.power',0,{$ifdef FPC_ProcVar}@{$endif}DAMath.power, 10, 500, 10, 500);
{$ifdef USE_MATH}
  testpow('math.power' ,0,{$ifdef FPC_ProcVar}@{$endif}xpower, 1e-3, 100, -100, 100);
  testpow('math.power' ,0,{$ifdef FPC_ProcVar}@{$endif}xpower, 10, 500, 10, 500);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_pow1pm1;
begin
  sep_line;
  testpow('DAMath.pow1pm1',2,{$ifdef FPC_ProcVar}@{$endif}pow1pm1, -0.999, 2, -1, 1);
  testpow('DAMath.pow1pm1',2,{$ifdef FPC_ProcVar}@{$endif}pow1pm1, -0.999, 2, 1e-3, 10);
{$ifdef USE_MATH}
  testpow('math.power(1+x)-1',2,{$ifdef FPC_ProcVar}@{$endif}math_pow1pm1,  -0.999, 2, -1, 1);
  testpow('math.power(1+x)-1',2,{$ifdef FPC_ProcVar}@{$endif}math_pow1pm1, -0.999, 2, 1e-3, 10);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_powm1;
begin
  sep_line;
  testpow('DAMath.powm1'  ,1,{$ifdef FPC_ProcVar}@{$endif}powm1, 1e-3, 10, -1, 1);
  testpow('DAMath.powm1'  ,1,{$ifdef FPC_ProcVar}@{$endif}powm1, 1e-3, 10, 1e-3, 10);
{$ifdef USE_MATH}
  testpow('math.power-1',1,{$ifdef FPC_ProcVar}@{$endif}math_powm1, 1e-3, 10, -1, 1);
  testpow('math.power-1',1,{$ifdef FPC_ProcVar}@{$endif}math_powm1, 1e-3, 10, 1e-3, 10);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_logbase;
begin
  sep_line;
  testlogxy('DAMath.logbase',{$ifdef FPC_ProcVar}@{$endif}DAMath.logbase, 1e-3, 1e10, 1e-3, 1e3);
{$ifdef USE_MATH}
  testlogxy('math.logn',{$ifdef FPC_ProcVar}@{$endif}xlogn, 1e-3, 1e10, 1e-3, 1e3);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_exp3;
begin
  sep_line;
  testfun('DAMath.exp3',{$ifdef FPC_ProcVar}@{$endif}exp3, {$ifdef FPC_ProcVar}@{$endif}mpf_exp3, -50, 50);
  testfun('DAMath.exp3',{$ifdef FPC_ProcVar}@{$endif}exp3, {$ifdef FPC_ProcVar}@{$endif}mpf_exp3, 50, 640);
end;


{---------------------------------------------------------------------------}
procedure test_exp5;
begin
  sep_line;
  testfun('DAMath.exp5',{$ifdef FPC_ProcVar}@{$endif}exp5, {$ifdef FPC_ProcVar}@{$endif}mpf_exp5, -40, 40);
  testfun('DAMath.exp5',{$ifdef FPC_ProcVar}@{$endif}exp5, {$ifdef FPC_ProcVar}@{$endif}mpf_exp5, 40, 440);
end;


{---------------------------------------------------------------------------}
procedure test_exp7;
begin
  sep_line;
  testfun('DAMath.exp7',{$ifdef FPC_ProcVar}@{$endif}exp7, {$ifdef FPC_ProcVar}@{$endif}mpf_exp7, -40, 40);
  testfun('DAMath.exp7',{$ifdef FPC_ProcVar}@{$endif}exp7, {$ifdef FPC_ProcVar}@{$endif}mpf_exp7, 40, 360);
end;


{---------------------------------------------------------------------------}
procedure test_cbrt;
begin
  sep_line;
  testfun('DAMath.cbrt',{$ifdef FPC_ProcVar}@{$endif}cbrt, {$ifdef FPC_ProcVar}@{$endif}mpf_cbrt, 0, 1000);
  testfun('DAMath.cbrt',{$ifdef FPC_ProcVar}@{$endif}cbrt, {$ifdef FPC_ProcVar}@{$endif}mpf_cbrt, -1e6, -1e3);
end;


{---------------------------------------------------------------------------}
procedure test_exprel;
begin
  sep_line;
  testfun('DAMath.exprel',{$ifdef FPC_ProcVar}@{$endif}exprel, {$ifdef FPC_ProcVar}@{$endif}mpf_exprel, -700, 700);
  testfun('DAMath.exprel',{$ifdef FPC_ProcVar}@{$endif}exprel, {$ifdef FPC_ProcVar}@{$endif}mpf_exprel, -1, 1);
  testfun('DAMath.exprel',{$ifdef FPC_ProcVar}@{$endif}exprel, {$ifdef FPC_ProcVar}@{$endif}mpf_exprel, -1e-4, 1e-4);
end;


{---------------------------------------------------------------------------}
procedure test_ln1pmx;
begin
  sep_line;
  testfun('DAMath.ln1pmx',{$ifdef FPC_ProcVar}@{$endif}ln1pmx, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1pmx, -0.999, -0.5);
  testfun('DAMath.ln1pmx',{$ifdef FPC_ProcVar}@{$endif}ln1pmx, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1pmx, -0.5, 0.5);
  testfun('DAMath.ln1pmx',{$ifdef FPC_ProcVar}@{$endif}ln1pmx, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1pmx, 0.5, 5);
end;


{---------------------------------------------------------------------------}
procedure test_covers;
begin
  sep_line;
  testfun('DAMath.covers',{$ifdef FPC_ProcVar}@{$endif}covers, {$ifdef FPC_ProcVar}@{$endif}mpf_covers, -Pi_4, Pi_4);
  testfun('DAMath.covers',{$ifdef FPC_ProcVar}@{$endif}covers, {$ifdef FPC_ProcVar}@{$endif}mpf_covers, Pi_2-0.1, Pi_2+0.1);
  testfun('DAMath.covers',{$ifdef FPC_ProcVar}@{$endif}covers, {$ifdef FPC_ProcVar}@{$endif}mpf_covers, -1e6, 1e6);
end;


{---------------------------------------------------------------------------}
procedure test_vers;
begin
  sep_line;
  testfun('DAMath.vers',{$ifdef FPC_ProcVar}@{$endif}vers, {$ifdef FPC_ProcVar}@{$endif}mpf_vers, -1e-3, 1e-3);
  testfun('DAMath.vers',{$ifdef FPC_ProcVar}@{$endif}vers, {$ifdef FPC_ProcVar}@{$endif}mpf_vers, Pi-1e-3, Pi+1e-3);
  testfun('DAMath.vers',{$ifdef FPC_ProcVar}@{$endif}vers, {$ifdef FPC_ProcVar}@{$endif}mpf_vers, -1e6, 1e6);
end;


{---------------------------------------------------------------------------}
procedure test_versint;
begin
  sep_line;
  testfun('DAMath.versint',{$ifdef FPC_ProcVar}@{$endif}versint, {$ifdef FPC_ProcVar}@{$endif}mpf_versint, -1e-3, 1e-3);
  testfun('DAMath.versint',{$ifdef FPC_ProcVar}@{$endif}versint, {$ifdef FPC_ProcVar}@{$endif}mpf_versint, -2, 2);
  testfun('DAMath.versint',{$ifdef FPC_ProcVar}@{$endif}versint, {$ifdef FPC_ProcVar}@{$endif}mpf_versint, -1e6, 1e6);
end;


{---------------------------------------------------------------------------}
function sqrt1pm1_naive(x: double): double;
  {-Return sqrt(1+x)-1, accurate even for x near 0, x>=-1}
begin
  sqrt1pm1_naive := sqrt(1.0+x)-1.0;
end;


{---------------------------------------------------------------------------}
procedure test_sqrt1pm1;
begin
  sep_line;
  testfun('DAMath.sqrt1pm1',{$ifdef FPC_ProcVar}@{$endif}sqrt1pm1, {$ifdef FPC_ProcVar}@{$endif}mpf_sqrt1pm1, -1, 7);
{$ifndef DAMOnly}
  testfun('sqrt(1+x)-1',{$ifdef FPC_ProcVar}@{$endif}sqrt1pm1_naive, {$ifdef FPC_ProcVar}@{$endif}mpf_sqrt1pm1, -1, 7);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_sinc;
begin
  sep_line;
  testfun('DAMath.sinc',{$ifdef FPC_ProcVar}@{$endif}sinc, {$ifdef FPC_ProcVar}@{$endif}mpf_sinc, 0.0, 1e-3);
  testfun('DAMath.sinc',{$ifdef FPC_ProcVar}@{$endif}sinc, {$ifdef FPC_ProcVar}@{$endif}mpf_sinc, 1e-3, 1e10);;
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function expx2a(x: double): double;
begin
  expx2a := system.exp(x*abs(x));
end;


{---------------------------------------------------------------------------}
procedure test_expx2;
begin
  sep_line;
  testfun('DAMath.expx2',{$ifdef FPC_ProcVar}@{$endif}expx2, {$ifdef FPC_ProcVar}@{$endif}mpf_expx2, -2, 2);
  testfun('DAMath.expx2',{$ifdef FPC_ProcVar}@{$endif}expx2, {$ifdef FPC_ProcVar}@{$endif}mpf_expx2, -25, 0);
  testfun('DAMath.expx2',{$ifdef FPC_ProcVar}@{$endif}expx2, {$ifdef FPC_ProcVar}@{$endif}mpf_expx2, 0, 25);
{$ifndef DAMOnly}
  testfun('system.exp(x*|x|)',{$ifdef FPC_ProcVar}@{$endif}expx2a, {$ifdef FPC_ProcVar}@{$endif}mpf_expx2, -2, 2);
  testfun('system.exp(x*|x|)',{$ifdef FPC_ProcVar}@{$endif}expx2a, {$ifdef FPC_ProcVar}@{$endif}mpf_expx2, -25, 0);
  testfun('system.exp(x*|x|)',{$ifdef FPC_ProcVar}@{$endif}expx2a, {$ifdef FPC_ProcVar}@{$endif}mpf_expx2, 0, 25);
{$endif}
end;


{---------------------------------------------------------------------------}
function simple_expmx2h(x: double): double;
  {-Return exp(-0.5*x^2)}
begin
  simple_expmx2h := exp(-0.5*x*x);
end;


{---------------------------------------------------------------------------}
procedure test_expmx2h;
begin
  sep_line;
  testfun('DAMath.expmx2h',{$ifdef FPC_ProcVar}@{$endif}expmx2h, {$ifdef FPC_ProcVar}@{$endif}mpf_expmx2h, 0,4);
  testfun('DAMath.expmx2h',{$ifdef FPC_ProcVar}@{$endif}expmx2h, {$ifdef FPC_ProcVar}@{$endif}mpf_expmx2h, 0,37.5);
{$ifndef DAMOnly}
  testfun('exp(-0.5*x^2)',{$ifdef FPC_ProcVar}@{$endif}simple_expmx2h, {$ifdef FPC_ProcVar}@{$endif}mpf_expmx2h, 0,4);
  testfun('exp(-0.5*x^2)',{$ifdef FPC_ProcVar}@{$endif}simple_expmx2h, {$ifdef FPC_ProcVar}@{$endif}mpf_expmx2h, 0,37.5);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_sincPi;
begin
  sep_line;
  testfun('DAMath.sincPi',{$ifdef FPC_ProcVar}@{$endif}sincPi, {$ifdef FPC_ProcVar}@{$endif}mpf_sincPi, 0.0, 1e-3);
  testfun('DAMath.sincPi',{$ifdef FPC_ProcVar}@{$endif}sincPi, {$ifdef FPC_ProcVar}@{$endif}mpf_sincPi, 1e-3, 1e10);;
end;


{---------------------------------------------------------------------------}
procedure test_coshm1;
begin
  sep_line;
  testfun('DAMath.coshm1',{$ifdef FPC_ProcVar}@{$endif}coshm1, {$ifdef FPC_ProcVar}@{$endif}mpf_coshm1, 0, 3e-8);
  testfun('DAMath.coshm1',{$ifdef FPC_ProcVar}@{$endif}coshm1, {$ifdef FPC_ProcVar}@{$endif}mpf_coshm1, 0, 1);
  testfun('DAMath.coshm1',{$ifdef FPC_ProcVar}@{$endif}coshm1, {$ifdef FPC_ProcVar}@{$endif}mpf_coshm1, 1, 700);
end;


{---------------------------------------------------------------------------}
procedure test_archav;
begin
  sep_line;
  testfun('DAMath.archav',{$ifdef FPC_ProcVar}@{$endif}archav, {$ifdef FPC_ProcVar}@{$endif}mpf_archav, 0, 0.001);
  testfun('DAMath.archav',{$ifdef FPC_ProcVar}@{$endif}archav, {$ifdef FPC_ProcVar}@{$endif}mpf_archav, 0,1);
  testfun('DAMath.archav',{$ifdef FPC_ProcVar}@{$endif}archav, {$ifdef FPC_ProcVar}@{$endif}mpf_archav, 0.999,1);
end;


{---------------------------------------------------------------------------}
procedure test_gd;
begin
  sep_line;
  testfun('DAMath.gd',{$ifdef FPC_ProcVar}@{$endif}gd, {$ifdef FPC_ProcVar}@{$endif}mpf_gd, -2e-2, 2e-2);
  testfun('DAMath.gd',{$ifdef FPC_ProcVar}@{$endif}gd, {$ifdef FPC_ProcVar}@{$endif}mpf_gd, -50, 50);
end;


{---------------------------------------------------------------------------}
procedure test_arcgd;
begin
  sep_line;
  testfun('DAMath.arcgd',{$ifdef FPC_ProcVar}@{$endif}arcgd, {$ifdef FPC_ProcVar}@{$endif}mpf_arcgd, -2e-2, 2e-2);
  testfun('DAMath.arcgd',{$ifdef FPC_ProcVar}@{$endif}arcgd, {$ifdef FPC_ProcVar}@{$endif}mpf_arcgd, -1.57, 1.57);
end;


{---------------------------------------------------------------------------}
procedure test_degrad;
begin
  sep_line;
  testfun('DAMath.DegToRad',{$ifdef FPC_ProcVar}@{$endif}DegToRad, {$ifdef FPC_ProcVar}@{$endif}mpf_deg2rad,-10000,10000);
  testfun('DAMath.RadToDeg',{$ifdef FPC_ProcVar}@{$endif}RadToDeg, {$ifdef FPC_ProcVar}@{$endif}mpf_rad2deg,-10000,10000);
end;


{---------------------------------------------------------------------------}
procedure test_rem_2pi;
begin
  sep_line;
  testfun('DAMath.rem_2pi',{$ifdef FPC_ProcVar}@{$endif}rem_2pi, {$ifdef FPC_ProcVar}@{$endif}mpf_rem_2pi, 0, 1.125e9);
  testfun('DAMath.rem_2pi',{$ifdef FPC_ProcVar}@{$endif}rem_2pi, {$ifdef FPC_ProcVar}@{$endif}mpf_rem_2pi, 1.125e9, 1e300);
end;


{---------------------------------------------------------------------------}
procedure test_log2p1;
begin
  sep_line;
  testfun('DAMath.log2p1',{$ifdef FPC_ProcVar}@{$endif}log2p1, {$ifdef FPC_ProcVar}@{$endif}mpf_log2p1, -1e-4, 1e-4);
  testfun('DAMath.log2p1',{$ifdef FPC_ProcVar}@{$endif}log2p1, {$ifdef FPC_ProcVar}@{$endif}mpf_log2p1, -0.9999999, 1e50);
end;


{---------------------------------------------------------------------------}
procedure test_log10p1;
begin
  sep_line;
  testfun('DAMath.log10p1',{$ifdef FPC_ProcVar}@{$endif}log10p1, {$ifdef FPC_ProcVar}@{$endif}mpf_log10p1, -1e-4, 1e-4);
  testfun('DAMath.log10p1',{$ifdef FPC_ProcVar}@{$endif}log10p1, {$ifdef FPC_ProcVar}@{$endif}mpf_log10p1, -0.9999999, 1e50);
end;


{---------------------------------------------------------------------------}
procedure test_exp10m1;
begin
  sep_line;
  testfun('DAMath.exp10m1',{$ifdef FPC_ProcVar}@{$endif}exp10m1, {$ifdef FPC_ProcVar}@{$endif}mpf_exp10m1, -1e-4, 1e-4);
  testfun('DAMath.exp10m1',{$ifdef FPC_ProcVar}@{$endif}exp10m1, {$ifdef FPC_ProcVar}@{$endif}mpf_exp10m1, -2, 2);
  testfun('DAMath.exp10m1',{$ifdef FPC_ProcVar}@{$endif}exp10m1, {$ifdef FPC_ProcVar}@{$endif}mpf_exp10m1, -100, 100);
end;


{---------------------------------------------------------------------------}
procedure test_sind;
begin
  sep_line;
  testfun('DAMath.sind',{$ifdef FPC_ProcVar}@{$endif}sind, {$ifdef FPC_ProcVar}@{$endif}mpf_sind, -1e10, 1e10);
end;


{---------------------------------------------------------------------------}
procedure test_cosd;
begin
  sep_line;
  testfun('DAMath.cosd',{$ifdef FPC_ProcVar}@{$endif}cosd, {$ifdef FPC_ProcVar}@{$endif}mpf_cosd, -1e10, 1e10);
end;


{---------------------------------------------------------------------------}
procedure test_tand;
begin
  sep_line;
  testfun('DAMath.tand',{$ifdef FPC_ProcVar}@{$endif}tand, {$ifdef FPC_ProcVar}@{$endif}mpf_tand, -1e10, 1e10);
end;


{---------------------------------------------------------------------------}
procedure test_cotd;
begin
  sep_line;
  testfun('DAMath.cotd',{$ifdef FPC_ProcVar}@{$endif}cotd, {$ifdef FPC_ProcVar}@{$endif}mpf_cotd, -1e10, 1e10);
end;


{---------------------------------------------------------------------------}
procedure test_arcsind;
begin
  sep_line;
  testfun('DAMath.arcsind',{$ifdef FPC_ProcVar}@{$endif}arcsind, {$ifdef FPC_ProcVar}@{$endif}mpf_arcsind, -1, 1);
end;


{---------------------------------------------------------------------------}
procedure test_arccosd;
begin
  sep_line;
  testfun('DAMath.arccosd',{$ifdef FPC_ProcVar}@{$endif}arccosd, {$ifdef FPC_ProcVar}@{$endif}mpf_arccosd, -1, 1);
end;


{---------------------------------------------------------------------------}
procedure test_arctand;
begin
  sep_line;
  testfun('DAMath.arctand',{$ifdef FPC_ProcVar}@{$endif}arctand, {$ifdef FPC_ProcVar}@{$endif}mpf_arctand, -1, -1e-10);
  testfun('DAMath.arctand',{$ifdef FPC_ProcVar}@{$endif}arctand, {$ifdef FPC_ProcVar}@{$endif}mpf_arctand, 1, 1e5);
end;


{---------------------------------------------------------------------------}
procedure test_arccotd;
begin
  sep_line;
  testfun('DAMath.arccotd',{$ifdef FPC_ProcVar}@{$endif}arccotd, {$ifdef FPC_ProcVar}@{$endif}mpf_arccotd, -1, -1e-10);
  testfun('DAMath.arccotd',{$ifdef FPC_ProcVar}@{$endif}arccotd, {$ifdef FPC_ProcVar}@{$endif}mpf_arccotd, 1, 1e5);
end;


{---------------------------------------------------------------------------}
procedure test_arccotcd;
begin
  sep_line;
  testfun('DAMath.arccotcd',{$ifdef FPC_ProcVar}@{$endif}arccotcd, {$ifdef FPC_ProcVar}@{$endif}mpf_arccotcd, -1e5, -1);
  testfun('DAMath.arccotcd',{$ifdef FPC_ProcVar}@{$endif}arccotcd, {$ifdef FPC_ProcVar}@{$endif}mpf_arccotcd, -1, 1);
  testfun('DAMath.arccotcd',{$ifdef FPC_ProcVar}@{$endif}arccotcd, {$ifdef FPC_ProcVar}@{$endif}mpf_arccotcd, 1, 1e5);
end;


{---------------------------------------------------------------------------}
procedure test_ln1mexp;
begin
  sep_line;
  testfun('DAMath.ln1mexp',{$ifdef FPC_ProcVar}@{$endif}ln1mexp, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1mexp, -1, -1e-40);
  testfun('DAMath.ln1mexp',{$ifdef FPC_ProcVar}@{$endif}ln1mexp, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1mexp, -100,-1);
end;


{---------------------------------------------------------------------------}
procedure test_logit;
begin
  sep_line;
  testfun('DAMath.logit',{$ifdef FPC_ProcVar}@{$endif}logit, {$ifdef FPC_ProcVar}@{$endif}mpf_logit, 0,1);
end;


{---------------------------------------------------------------------------}
function dam_compound(x,n: double): double;
  {-Return (1+x)^n}
begin
  dam_compound := compound(x, round(n));
end;


{---------------------------------------------------------------------------}
procedure test_pow1p;
begin
  sep_line;
  testpow('DAMath.pow1p',3,{$ifdef FPC_ProcVar}@{$endif}pow1p, -0.999, 2, -1, 1);
  testpow('DAMath.pow1p',3,{$ifdef FPC_ProcVar}@{$endif}pow1p, -0.999, 2, -200, 600);
  testpow('DAMath.pow1p',3,{$ifdef FPC_ProcVar}@{$endif}pow1p, -0.5, 0.5, -250, 750);
end;


{---------------------------------------------------------------------------}
procedure test_pow1pf;
begin
  sep_line;
  testpow('DAMath.pow1pf',5,{$ifdef FPC_ProcVar}@{$endif}pow1pf, -0.999, 2, -1, 1);
  testpow('DAMath.pow1pf',5,{$ifdef FPC_ProcVar}@{$endif}pow1pf, -0.999, 2, -200, 600);
  testpow('DAMath.pow1pf',5,{$ifdef FPC_ProcVar}@{$endif}pow1pf, -0.5, 0.5, -250, 750);
end;



{---------------------------------------------------------------------------}
procedure test_compound;
begin
  sep_line;
  testpow('DAMath.compound',4,{$ifdef FPC_ProcVar}@{$endif}dam_compound, -0.999, 2, -1000, 10000);
  testpow('DAMath.compound',4,{$ifdef FPC_ProcVar}@{$endif}dam_compound, -0.5, 0.5, -1000, 30000);
{$ifdef USE_XDMATH}
  testpow('math.compound',4,{$ifdef FPC_ProcVar}@{$endif}xcompound, -0.999, 2, -1000, 10000);
  testpow('math.compound',4,{$ifdef FPC_ProcVar}@{$endif}xcompound, -0.5, 0.5, -1000, 30000);
{$endif}
end;


{---------------------------------------------------------------------------}
procedure test_logistic;
begin
  sep_line;
  testfun('DAMath.logistic',{$ifdef FPC_ProcVar}@{$endif}logistic, {$ifdef FPC_ProcVar}@{$endif}mpf_logistic, -0.01, 0.01);
  testfun('DAMath.logistic',{$ifdef FPC_ProcVar}@{$endif}logistic, {$ifdef FPC_ProcVar}@{$endif}mpf_logistic, -40,40);
  testfun('DAMath.logistic',{$ifdef FPC_ProcVar}@{$endif}logistic, {$ifdef FPC_ProcVar}@{$endif}mpf_logistic, -700, -40);
end;


{---------------------------------------------------------------------------}
procedure test_powpi;
var
  r,x,y,t,rms: double;
  n,m,ne: longint;
const
  nmax = 3100;
begin
  mpf_set_pi2k(a,-2);
  r := 0.0;
  m := 0;
  ne  := -Maxlongint;
  rms := 0.0;

  sep_line;
  writeln('Test of damath.powpi2k(-2,n), n=',-nmax,' .. +',nmax);
  for n:=-nmax to nmax do begin
    x := powpi2k(-2,n);
    if (not IsInfd(x)) and (x<>0.0) then begin
      mpf_expt_int(a,n,b);
      y := mpf_todouble(b);
      if y=0 then t := abs(x)
      else if x=y then t:=0
      else begin
        mpf_set_dbl(c,x);
        mpf_sub(b,c,c);
        t := abs(mpf_todouble(c)/y);
      end;
      y := t/eps_d;
      rms := rms + sqr(y);
      inc(m);
      if t>r then begin
        r := t;
        ne := n;
      end;
    end;
  end;
  writeln('RMS = ', sqrt(rms/m):1:2,', max rel = ',r/eps_d:1:2,' eps for n=',ne);
end;


{---------------------------------------------------------------------------}
procedure test_arccos1m;
begin
  sep_line;
  testfun('DAMath.arccos1m',{$ifdef FPC_ProcVar}@{$endif}arccos1m, {$ifdef FPC_ProcVar}@{$endif}mpf_arccos1m, 0, 2);
  testfun('DAMath.arccos1m',{$ifdef FPC_ProcVar}@{$endif}arccos1m, {$ifdef FPC_ProcVar}@{$endif}mpf_arccos1m, 0, 1/128);
end;


{---------------------------------------------------------------------------}
procedure test_sqrt1pmx;
begin
  sep_line;
  testfun('DAMath.sqrt1pmx',{$ifdef FPC_ProcVar}@{$endif}sqrt1pmx, {$ifdef FPC_ProcVar}@{$endif}mpf_sqrt1pmx, -1e-5, 1e-5);
  testfun('DAMath.sqrt1pmx',{$ifdef FPC_ProcVar}@{$endif}sqrt1pmx, {$ifdef FPC_ProcVar}@{$endif}mpf_sqrt1pmx, -10, 10);
  testfun('DAMath.sqrt1pmx',{$ifdef FPC_ProcVar}@{$endif}sqrt1pmx, {$ifdef FPC_ProcVar}@{$endif}mpf_sqrt1pmx, -1e10, 1e10);
end;


{---------------------------------------------------------------------------}
function xsinpi(x: double): double;
  {-Return tan(Pi*x), result will be 0 for abs(x) >= 2^64}
begin
  xsinpi := sin(Pi*x);
end;


{---------------------------------------------------------------------------}
function xcospi(x: double): double;
  {-Return cos(Pi*x), result will be 1 for abs(x) >= 2^64}
begin
  xcospi := cos(Pi*x);
end;


{---------------------------------------------------------------------------}
function xtanpi(x: double): double;
  {-Return tan(Pi*x), result will be 0 for abs(x) >= 2^64}
begin
  xtanpi := tan(Pi*x);
end;


{---------------------------------------------------------------------------}
procedure test_sinpi;
begin
  sep_line;
  testfun('DAMath.sinpi',{$ifdef FPC_ProcVar}@{$endif}sinpi,  {$ifdef FPC_ProcVar}@{$endif}mpf_sinpi, -100, 100);
  testfun('sin(pi*x)',  {$ifdef FPC_ProcVar}@{$endif}xsinpi, {$ifdef FPC_ProcVar}@{$endif}mpf_sinpi, -100, 100);
end;

{---------------------------------------------------------------------------}
procedure test_cospi;
begin
  sep_line;
  testfun('DAMath.cospi',{$ifdef FPC_ProcVar}@{$endif}cospi,  {$ifdef FPC_ProcVar}@{$endif}mpf_cospi, -100, 100);
  testfun('cos(pi*x)',  {$ifdef FPC_ProcVar}@{$endif}xcospi, {$ifdef FPC_ProcVar}@{$endif}mpf_cospi, -100, 100);
end;



{---------------------------------------------------------------------------}
procedure test_tanpi;
begin
  sep_line;
  testfun('DAMath.tanpi',{$ifdef FPC_ProcVar}@{$endif}tanpi,  {$ifdef FPC_ProcVar}@{$endif}mpf_tanpi, -10, 10);
  testfun('tan(pi*x)',   {$ifdef FPC_ProcVar}@{$endif}xtanpi, {$ifdef FPC_ProcVar}@{$endif}mpf_tanpi, -10, 10);
end;


{---------------------------------------------------------------------------}
function xlncosh(x: double): double;
begin
  xlncosh := ln(cosh(x));
end;


{---------------------------------------------------------------------------}
function xlnsinh(x: double): double;
begin
  xlnsinh := ln(cosh(x));
end;


{---------------------------------------------------------------------------}
procedure test_lncosh;
begin
  sep_line;
  testfun('DAMath.lncosh',{$ifdef FPC_ProcVar}@{$endif}lncosh,  {$ifdef FPC_ProcVar}@{$endif}mpf_lncosh, 0, 1/128);
  testfun('DAMath.lncosh',{$ifdef FPC_ProcVar}@{$endif}lncosh,  {$ifdef FPC_ProcVar}@{$endif}mpf_lncosh, 1, 700);
  testfun('ln(cosh)',     {$ifdef FPC_ProcVar}@{$endif}xlncosh, {$ifdef FPC_ProcVar}@{$endif}mpf_lncosh, 0, 1/128);
  testfun('ln(cosh)',     {$ifdef FPC_ProcVar}@{$endif}xlncosh, {$ifdef FPC_ProcVar}@{$endif}mpf_lncosh, 1, 700);
end;


{---------------------------------------------------------------------------}
procedure test_lnsinh;
begin
  sep_line;
  testfun('DAMath.lnsinh',{$ifdef FPC_ProcVar}@{$endif}lnsinh,  {$ifdef FPC_ProcVar}@{$endif}mpf_lnsinh, 1e-50, 1/128);
  testfun('DAMath.lnsinh',{$ifdef FPC_ProcVar}@{$endif}lnsinh,  {$ifdef FPC_ProcVar}@{$endif}mpf_lnsinh, 1, 700);
  testfun('ln(sinh)',     {$ifdef FPC_ProcVar}@{$endif}xlnsinh, {$ifdef FPC_ProcVar}@{$endif}mpf_lnsinh, 1e-50, 1/128);
  testfun('ln(sinh)',     {$ifdef FPC_ProcVar}@{$endif}xlnsinh, {$ifdef FPC_ProcVar}@{$endif}mpf_lnsinh, 1, 700);
end;


{---------------------------------------------------------------------------}
procedure test_ln1pexp;
begin
  sep_line;
  testfun('DAMath.ln1pexp',{$ifdef FPC_ProcVar}@{$endif}ln1pexp, {$ifdef FPC_ProcVar}@{$endif}mpf_ln1pexp, -50, 50);
end;



{---------------------------------------------------------------------------}
procedure test_fma1(const s: mp_string; x1,x2: double; fma: boolean);
var
  x,y,z,u,d,f: double;
  umax, rms: double;
  i: longint;
{$ifdef BIT16}
  const n = 1000;
{$else}
  {$ifdef BIGTEST}
    const n = 20000;
  {$else}
    const n = 2000;
  {$endif}
{$endif}
begin
  writeln('Test of ',s);
  writeln('at ',n,' random values in [', fmt(x1),' .. ', fmt(x2),']');
  umax := -1;
  rms  := 0;
  taus88_init(tctx, DefSeed);
  for i:=1 to n do begin
    x := randx(x1,x2);
    y := randx(x1,x2);
    z := randx(x1,x2);

    if fma then f := fma_d(x,y,z)
    else f := x*y + z;

    {mf = x*y + z}
    mpf_set_dbl(mx, x);
    mpf_mul_dbl(mx, y, mf);
    mpf_add_dbl(mf, z, mf);

    mpf_set_dbl(a, f);       { a = xf(x)}
    mpf_sub_dbl(mf,f,b);     { b = mpf(x)-xf(x)}

    d := reldev(a,mf);
    u := d/eps_d;
    rms := rms + sqr(u);
    if u>umax then begin
      umax := u;
      mpf_copy(mf,mmax);
    end;
  end;
  writeln('RMS = ', sqrt(rms/n):1:2,', max rel = ',umax:1:2,' eps at ');
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_fma;
begin
  sep_line;
  test_fma1('DAMath.fma_d', -128, 128, true);
  test_fma1('x*y + z', -128, 128, false);
end;


{$ifdef USE_XDMATH}
{---------------------------------------------------------------------------}
function xcomprel(R,N: double): double;
begin
  {Stripped-down version from math.annuity2}
  if R = 0.0 then
  begin
    xcomprel := N;
  end
  else
  begin
    if Abs(R) < 6.1E-5 then
    begin
      xcomprel := N*(1+(N-1)*R/2);
    end
    else
    begin
      xcomprel := (IntPower(1.0 + R, round(N))-1) / R
    end;
  end;
end;
{$endif}


{---------------------------------------------------------------------------}
procedure test_comprel;
begin
  sep_line;
  testpow('DAMath.comprel',6,{$ifdef FPC_ProcVar}@{$endif}comprel, -0.999, 2, 1, 100);
  testpow('DAMath.comprel',6,{$ifdef FPC_ProcVar}@{$endif}comprel, -0.5, 0.5, 1, 10000);
{$ifdef USE_XDMATH}
  testpow('math."comprel"',6,{$ifdef FPC_ProcVar}@{$endif}xcomprel, -0.999, 2, 1, 100);
  testpow('math."comprel"',6,{$ifdef FPC_ProcVar}@{$endif}xcomprel, -0.5, 0.5, 1, 10000);
{$endif}
end;


{---------------------------------------------------------------------------}
label
  done;

begin
  mpf_set_default_prec(PREX);
  writeln('Test DAMath V', DAMath_Version,' with MP_Arith V', MP_VERSION, '   (c) 2013-2018 W.Ehrhardt');
  writeln('Karatsuba  cutoffs: mul/sqr = ',mp_mul_cutoff,'/',mp_sqr_cutoff);
  writeln('Toom-3, BZ cutoffs: mul/sqr = ',mp_t3m_cutoff,'/',mp_t3s_cutoff,  ',  div = ',mp_bz_cutoff);
  writeln('Current mp_float default bit precision = ', mpf_get_default_prec,
          ',  decimal precision = ', mpf_get_default_prec*ln(2)/ln(10):1:1);
  writeln('Machine eps for double =', eps_d:20);

  mp_show_plus := true;

  mpf_initp3(a,b,c,mpf_get_default_prec);
  mpf_initp3(mx,mf,mmax,mpf_get_default_prec);
{
  test_comprel;
  goto done;
}
  test_arccos;
  test_arccos1m;
  test_arccosd;
  test_arccosh1p;
  test_arccosh;
  test_arccot;
  test_arccotc;
  test_arccotcd;
  test_arccotd;
  test_arccoth;
  test_arccsc;
  test_arccsch;
  test_arcgd;
  test_archav;
  test_arcsec;
  test_arcsech;
  test_arcsin;
  test_arcsind;
  test_arcsinh;
  test_arctan;
  test_arctand;
  test_arctanh;
  test_cbrt;
  test_compound;
  test_comprel;
  test_cos;
  test_cosd;
  test_cosh;
  test_coshm1;
  test_cospi;
  test_cot;
  test_cotd;
  test_coth;
  test_covers;
  test_csc;
  test_csch;
  test_degrad;
  test_exp10;
  test_exp10m1;
  test_exp2;
  test_exp2m1;
  test_exp3;
  test_exp5;
  test_exp7;
  test_exp;
  test_expm1;
  test_expmx2h;
  test_exprel;
  test_expx2;
  test_gd;
  test_ln1mexp;
  test_ln1pexp;
  test_ln1p;
  test_ln1pmx;
  test_ln;
  test_lncosh;
  test_lnsinh;
  test_log10;
  test_log10p1;
  test_log2;
  test_log2p1;
  test_logbase;
  test_logit;
  test_logistic;
  test_pow1p;
  test_pow1pf;
  test_pow1pm1;
  test_power;
  test_powm1;
  test_powpi;
  test_rem_2pi;
  test_sec;
  test_sech;
  test_sin;
  test_sinc;
  test_sincPi;
  test_sind;
  test_sinh;
  test_sinhc;
  test_sinhmx;
  test_sinpi;
  test_sqrt1pm1;
  test_sqrt1pmx;
  test_tan;
  test_tand;
  test_tanh;
  test_tanpi;
  test_vers;
  test_versint;
  test_fma;

done:
  mpf_clear3(a,b,c);
  mpf_clear3(mx,mf,mmax);

  {$ifdef MPC_Diagnostic}
    writeln;
    writeln;
    mp_dump_meminfo;
    mp_dump_diagctr;
  {$endif}

end.
