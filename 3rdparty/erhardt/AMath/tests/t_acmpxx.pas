{Test program to compare precision of AMCmplx functions with mpc routines}

program t_acmpxx;

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

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  {$ifdef MPC_Diagnostic}
    mp_supp,
  {$endif}
  taus88, BTypes, mp_types, mp_base, mp_real, mp_cmplx, amath, amcmplx;


{---------------------------------------------------------------------------}
procedure mpc_sqrt1mz2(const a: mp_complex; var b: mp_complex);
  {-Calculate b = sqrt(1-a^2)}
var
  z: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  mpc_initp(z,b.re.bitprec+16);
  if mp_error=MP_OKAY then begin
    mpc_sqr(a,z);
    mpc_sub_dbl(z,1,z);
    mpc_chs(z,z);
    mpc_sqrt(z,z);
    mpc_copyp(z,b);
    mpc_clear(z);
  end;
end;


const
  PREX = 160;
  DefSeed : longint = 0;

type
  TMProc = procedure(const x: mp_complex; var y: mp_complex);
  TCProc = procedure(const x: complex; var y: complex);

var
  tctx: taus88_ctx;
  d1: extended; {delta from 1}

var
  b,c: mp_float;
  a,mx,mf,mmax: mp_complex;


{---------------------------------------------------------------------------}
procedure sep_line;
begin
  writeln;
  writeln('-----------------------------------------------------------------');
end;


{---------------------------------------------------------------------------}
function reldevx(const a,b: mp_complex): extended;
  {-Return abs((a-b)/b)*2^b.re.bitprec, special if b=0, or a-b=0}
var
  c: mp_complex;
begin
  reldevx := -1;
  mpc_initp(c,b.re.bitprec);
  if mp_error<>MP_OKAY then exit;
  mpc_sub(a,b,c);
  if mpc_is0(c) then reldevx := 0.0
  else begin
    if not mpc_is0(b) then mpc_div(c,b,c);
    mpc_abs(c,c.re);
    reldevx := mpf_toextended(c.re);
  end;
  mpc_clear(c);
end;


{---------------------------------------------------------------------------}
function randx(x1,x2: extended): extended;
begin
  randx := x1 + taus88_double53(tctx)*(x2-x1);
end;


{---------------------------------------------------------------------------}
function fmt(x: extended): Bstring;
var
  a: extended;
  s: BString;
const
  AMAX = 1e18;
begin
  a := abs(x);
  if (a > AMAX) or ((a > 0.0) and (a < sqrt_epsh)) then str(x:20,s)
  else str(x:1:10,s);
  fmt := s;
end;



{---------------------------------------------------------------------------}
procedure testcfun(const s: mp_string; xf: TCProc; mpf: TMProc; x1,x2,y1,y2: extended);
var
  x,y,d,u: extended;
  f,z,zmax: complex;
  umax, xmax, ymax, rms: extended;
  i: longint;
const
  NDD = 50;
{$ifdef BIT16}
  const n = 200;
{$else}
  {$ifdef BIGTEST}
    const n = 10000;
  {$else}
    const n = 2000;
  {$endif}
{$endif}
{$define use_eps}
begin
  writeln('Test of ',s, '(x+i*y)  at ',n,' random values');
  writeln(' with x from [', fmt(x1),' .. ', fmt(x2),']');
  writeln('  and y from [', fmt(y1),' .. ', fmt(y2),']');
  umax := -1;
  xmax := 0;
  zmax := C_0;
  ymax := 0;
  rms  := 0;
  taus88_init(tctx, DefSeed);
  for i:=1 to n do begin
    x := randx(x1,x2);
    y := randx(y1,y2);
    z.re := x;
    z.im := y;
    xf(z,f);
    mpc_set_ext(a,f.re,f.im);  {a=xf(z)}
    mpc_set_ext(mx,x,y);
    mpf(mx,mf);                {mf = mpf(z)}
    d := reldevx(a,mf);
    {$ifdef use_eps}
      u := d/eps_x;
    {$else}
      u := d;
    {$endif}
    rms := rms + sqr(u);
    if u>umax then begin
      umax := u;
      xmax := x;
      zmax := f;
      ymax := y;
      mpc_copy(mf,mmax);
    end;
  end;
  {$ifdef use_eps}
    writeln('RMS = ', sqrt(rms/n):1:2,', max rel = ',umax:1:2,' eps at');
  {$else}
    writeln('RMS = ', sqrt(rms/n):26,', max = ',umax:26, '  at');
  {$endif}
  writeln(' x(ext) = ',xmax:26, ' = $',ext2hex(xmax));
  writeln(' y(ext) = ',ymax:26, ' = $',ext2hex(ymax));
  writeln(' fmr.re = ',zmax.re:26, '   fmr.im = ',zmax.im:26);
{$ifdef debug}
  mpf_set_ext(b, xmax);
  mpf_mul_2k(b,64,b);
  writeln('   xmax = 2^-64 * ', mpf_decimal(b,NDD));
  mpf_set_ext(b, ymax);
  mpf_mul_2k(b,64,b);
  writeln('   ymax = 2^-64 * ', mpf_decimal(b,NDD));
{$endif}
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_arccos;
begin
  sep_line;
  testcfun('arccos', {$ifdef FPC_ProcVar}@{$endif}carccos, {$ifdef FPC_ProcVar}@{$endif}mpc_arccos, -2, 2, -2, 2);
  testcfun('arccos', {$ifdef FPC_ProcVar}@{$endif}carccos, {$ifdef FPC_ProcVar}@{$endif}mpc_arccos, 1-d1, 1+d1, -d1, d1);
  testcfun('arccos', {$ifdef FPC_ProcVar}@{$endif}carccos, {$ifdef FPC_ProcVar}@{$endif}mpc_arccos, -1000, 1000, -1000, 1000);
end;


{---------------------------------------------------------------------------}
procedure test_arccosh;
begin
  sep_line;
  testcfun('arccosh', {$ifdef FPC_ProcVar}@{$endif}carccosh, {$ifdef FPC_ProcVar}@{$endif}mpc_arccosh, -2, 2, -2, 2);
  testcfun('arccosh', {$ifdef FPC_ProcVar}@{$endif}carccosh, {$ifdef FPC_ProcVar}@{$endif}mpc_arccosh, 1-d1, 1+d1, -d1, d1);
  testcfun('arccosh', {$ifdef FPC_ProcVar}@{$endif}carccosh, {$ifdef FPC_ProcVar}@{$endif}mpc_arccosh, -1e3, 1e3, -1e3, 1e3);
end;


{---------------------------------------------------------------------------}
procedure test_arcsin;
begin
  sep_line;
  testcfun('arcsin', {$ifdef FPC_ProcVar}@{$endif}carcsin, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsin, -2, 2, -2, 2);
  testcfun('arcsin', {$ifdef FPC_ProcVar}@{$endif}carcsin, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsin, 1-d1, 1+d1, -d1, d1);
  testcfun('arcsin', {$ifdef FPC_ProcVar}@{$endif}carcsin, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsin, -1000, 1000, -1000, 1000);
end;


{---------------------------------------------------------------------------}
procedure test_arcsinh;
begin
  sep_line;
  testcfun('arcsinh', {$ifdef FPC_ProcVar}@{$endif}carcsinh, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsinh, -2, 2, -2, 2);
  testcfun('arcsinh', {$ifdef FPC_ProcVar}@{$endif}carcsinh, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsinh, -d1, +d1, 1-d1, 1+d1);
  testcfun('arcsinh', {$ifdef FPC_ProcVar}@{$endif}carcsinh, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsinh, -1e3, 1e3, -1e3, 1e3);
end;


{---------------------------------------------------------------------------}
procedure test_arctan;
begin
  sep_line;
  testcfun('arctan', {$ifdef FPC_ProcVar}@{$endif}carctan, {$ifdef FPC_ProcVar}@{$endif}mpc_arctan, -2, 2, -2, 2);
  testcfun('arctan', {$ifdef FPC_ProcVar}@{$endif}carctan, {$ifdef FPC_ProcVar}@{$endif}mpc_arctan, -d1, +d1, 1-d1, 1+d1);
  testcfun('arctan', {$ifdef FPC_ProcVar}@{$endif}carctan, {$ifdef FPC_ProcVar}@{$endif}mpc_arctan, -1000, 1000, -1000, 1000);
end;


{---------------------------------------------------------------------------}
procedure test_arctanh;
begin
  sep_line;
  testcfun('arctanh', {$ifdef FPC_ProcVar}@{$endif}carctanh, {$ifdef FPC_ProcVar}@{$endif}mpc_arctanh, -2, 2, -2, 2);
  testcfun('arctanh', {$ifdef FPC_ProcVar}@{$endif}carctanh, {$ifdef FPC_ProcVar}@{$endif}mpc_arctanh, 1-d1, 1+d1, -d1, d1);
  testcfun('arctanh', {$ifdef FPC_ProcVar}@{$endif}carctanh, {$ifdef FPC_ProcVar}@{$endif}mpc_arctanh, -1e3, 1e3, -1e3, 1e3);
end;


{---------------------------------------------------------------------------}
procedure test_arccsc;
begin
  sep_line;
  testcfun('arccsc', {$ifdef FPC_ProcVar}@{$endif}carccsc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccsc, -2, 2, -2, 2);
  testcfun('arccsc', {$ifdef FPC_ProcVar}@{$endif}carccsc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccsc, 1-d1, 1+d1, -d1, d1);
  testcfun('arccsc', {$ifdef FPC_ProcVar}@{$endif}carccsc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccsc, -1000, 1000, -1000, 1000);
end;


{---------------------------------------------------------------------------}
procedure test_arccsch;
begin
  sep_line;
  testcfun('arccsch', {$ifdef FPC_ProcVar}@{$endif}carccsch, {$ifdef FPC_ProcVar}@{$endif}mpc_arccsch, -2, 2, -2, 2);
  testcfun('arccsch', {$ifdef FPC_ProcVar}@{$endif}carccsch, {$ifdef FPC_ProcVar}@{$endif}mpc_arccsch, -d1, d1, 1-d1, 1+d1);
  testcfun('arccsch', {$ifdef FPC_ProcVar}@{$endif}carccsch, {$ifdef FPC_ProcVar}@{$endif}mpc_arccsch, -1000,1000,-1000,1000);
end;


{---------------------------------------------------------------------------}
procedure test_arcsec;
begin
  sep_line;
  testcfun('arcsec', {$ifdef FPC_ProcVar}@{$endif}carcsec, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsec, -2, 2, -2, 2);
  testcfun('arcsec', {$ifdef FPC_ProcVar}@{$endif}carcsec, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsec, 1-d1, 1+d1, -d1, d1);
  testcfun('arcsec', {$ifdef FPC_ProcVar}@{$endif}carcsec, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsec, -1000, 1000, -1000, 1000);
end;


{---------------------------------------------------------------------------}
procedure test_arcsech;
begin
  sep_line;
  testcfun('arcsech', {$ifdef FPC_ProcVar}@{$endif}carcsech, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsech, -2, 2, -2, 2);
  testcfun('arcsech', {$ifdef FPC_ProcVar}@{$endif}carcsech, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsech, 1-d1, 1+d1, -d1, d1);
  testcfun('arcsech', {$ifdef FPC_ProcVar}@{$endif}carcsech, {$ifdef FPC_ProcVar}@{$endif}mpc_arcsech, -1000,1000,-1000,1000);
end;


{---------------------------------------------------------------------------}
procedure test_cos;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('cos', {$ifdef FPC_ProcVar}@{$endif}ccos, {$ifdef FPC_ProcVar}@{$endif}mpc_cos, -2, 2, -2, 2);
  testcfun('cos', {$ifdef FPC_ProcVar}@{$endif}ccos, {$ifdef FPC_ProcVar}@{$endif}mpc_cos, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_cosh;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('cosh', {$ifdef FPC_ProcVar}@{$endif}ccosh, {$ifdef FPC_ProcVar}@{$endif}mpc_cosh, -2, 2, -2, 2);
  testcfun('cosh', {$ifdef FPC_ProcVar}@{$endif}ccosh, {$ifdef FPC_ProcVar}@{$endif}mpc_cosh, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_cot;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('cot', {$ifdef FPC_ProcVar}@{$endif}ccot, {$ifdef FPC_ProcVar}@{$endif}mpc_cot, -2, 2, -2, 2);
  testcfun('cot', {$ifdef FPC_ProcVar}@{$endif}ccot, {$ifdef FPC_ProcVar}@{$endif}mpc_cot, -a, a, -20, 20);
end;


{---------------------------------------------------------------------------}
procedure test_coth;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('coth', {$ifdef FPC_ProcVar}@{$endif}ccoth, {$ifdef FPC_ProcVar}@{$endif}mpc_coth, -2, 2, -2, 2);
  testcfun('coth', {$ifdef FPC_ProcVar}@{$endif}ccoth, {$ifdef FPC_ProcVar}@{$endif}mpc_coth, -20, -20, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_exp;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('exp', {$ifdef FPC_ProcVar}@{$endif}cexp, {$ifdef FPC_ProcVar}@{$endif}mpc_exp, -2, 2, -2, 2);
  testcfun('exp', {$ifdef FPC_ProcVar}@{$endif}cexp, {$ifdef FPC_ProcVar}@{$endif}mpc_exp, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_exp2;
var
  a: extended;
begin
  sep_line;
  a := 1000;
  testcfun('exp2', {$ifdef FPC_ProcVar}@{$endif}cexp2, {$ifdef FPC_ProcVar}@{$endif}mpc_exp2, -2, 2, -2, 2);
  testcfun('exp2', {$ifdef FPC_ProcVar}@{$endif}cexp2, {$ifdef FPC_ProcVar}@{$endif}mpc_exp2, -a, a, -a, a);
end;

{---------------------------------------------------------------------------}
procedure test_exp10;
var
  a: extended;
begin
  sep_line;
  a := 500;
  testcfun('exp10', {$ifdef FPC_ProcVar}@{$endif}cexp10, {$ifdef FPC_ProcVar}@{$endif}mpc_exp10, -2, 2, -2, 2);
  testcfun('exp10', {$ifdef FPC_ProcVar}@{$endif}cexp10, {$ifdef FPC_ProcVar}@{$endif}mpc_exp10, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_expm1;
begin
  sep_line;
  testcfun('expm1', {$ifdef FPC_ProcVar}@{$endif}cexpm1, {$ifdef FPC_ProcVar}@{$endif}mpc_expm1, -d1, d1, -d1, d1);
  testcfun('expm1', {$ifdef FPC_ProcVar}@{$endif}cexpm1, {$ifdef FPC_ProcVar}@{$endif}mpc_expm1, -2, 2, -2, 2);
end;


{---------------------------------------------------------------------------}
procedure test_ln;
begin
  sep_line;
  testcfun('ln', {$ifdef FPC_ProcVar}@{$endif}cln, {$ifdef FPC_ProcVar}@{$endif}mpc_ln, -2, 2, -2, 2);
  testcfun('ln', {$ifdef FPC_ProcVar}@{$endif}cln, {$ifdef FPC_ProcVar}@{$endif}mpc_ln, -d1, d1, -d1, d1);
  testcfun('ln', {$ifdef FPC_ProcVar}@{$endif}cln, {$ifdef FPC_ProcVar}@{$endif}mpc_ln, -1e50, 1e50, -1e50, 1e50);
end;


{---------------------------------------------------------------------------}
procedure test_ln1p;
begin
  sep_line;
  testcfun('ln1p', {$ifdef FPC_ProcVar}@{$endif}cln1p, {$ifdef FPC_ProcVar}@{$endif}mpc_ln1p, -d1, d1, -d1, d1);
  testcfun('ln1p', {$ifdef FPC_ProcVar}@{$endif}cln1p, {$ifdef FPC_ProcVar}@{$endif}mpc_ln1p, -10, 10, -10, 10);
  testcfun('ln1p', {$ifdef FPC_ProcVar}@{$endif}cln1p, {$ifdef FPC_ProcVar}@{$endif}mpc_ln1p, 0, 0, -1e-8, 1e-8);
end;


{---------------------------------------------------------------------------}
procedure test_sec;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('sec', {$ifdef FPC_ProcVar}@{$endif}csec, {$ifdef FPC_ProcVar}@{$endif}mpc_sec, -2, 2, -2, 2);
  testcfun('sec', {$ifdef FPC_ProcVar}@{$endif}csec, {$ifdef FPC_ProcVar}@{$endif}mpc_sec, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_sech;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('sech', {$ifdef FPC_ProcVar}@{$endif}csech, {$ifdef FPC_ProcVar}@{$endif}mpc_sech, -2, 2, -2, 2);
  testcfun('sech', {$ifdef FPC_ProcVar}@{$endif}csech, {$ifdef FPC_ProcVar}@{$endif}mpc_sech, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_csc;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('csc', {$ifdef FPC_ProcVar}@{$endif}ccsc, {$ifdef FPC_ProcVar}@{$endif}mpc_csc, -2, 2, -2, 2);
  testcfun('csc', {$ifdef FPC_ProcVar}@{$endif}ccsc, {$ifdef FPC_ProcVar}@{$endif}mpc_csc, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_csch;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('csch', {$ifdef FPC_ProcVar}@{$endif}ccsch, {$ifdef FPC_ProcVar}@{$endif}mpc_csch, -2, 2, -2, 2);
  testcfun('csch', {$ifdef FPC_ProcVar}@{$endif}ccsch, {$ifdef FPC_ProcVar}@{$endif}mpc_csch, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_sin;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('sin', {$ifdef FPC_ProcVar}@{$endif}csin, {$ifdef FPC_ProcVar}@{$endif}mpc_sin, -2, 2, -2, 2);
  testcfun('sin', {$ifdef FPC_ProcVar}@{$endif}csin, {$ifdef FPC_ProcVar}@{$endif}mpc_sin, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_sinh;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('sinh', {$ifdef FPC_ProcVar}@{$endif}csinh, {$ifdef FPC_ProcVar}@{$endif}mpc_sinh, -2, 2, -2, 2);
  testcfun('sinh', {$ifdef FPC_ProcVar}@{$endif}csinh, {$ifdef FPC_ProcVar}@{$endif}mpc_sinh, -a, a, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_sqr;
begin
  sep_line;
  testcfun('sqr', {$ifdef FPC_ProcVar}@{$endif}csqr, {$ifdef FPC_ProcVar}@{$endif}mpc_sqr, -2, 2, -2, 2);
  testcfun('sqr', {$ifdef FPC_ProcVar}@{$endif}csqr, {$ifdef FPC_ProcVar}@{$endif}mpc_sqr, -1e10, 1e10, -1e10, 1e10);
end;


{---------------------------------------------------------------------------}
procedure test_sqrt;
begin
  sep_line;
  testcfun('sqrt', {$ifdef FPC_ProcVar}@{$endif}csqrt, {$ifdef FPC_ProcVar}@{$endif}mpc_sqrt, -2, 2, -2, 2);
  testcfun('sqrt', {$ifdef FPC_ProcVar}@{$endif}csqrt, {$ifdef FPC_ProcVar}@{$endif}mpc_sqrt, -1e10, 1e10, -1e10, 1e10);
end;


{---------------------------------------------------------------------------}
procedure test_tan;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('tan', {$ifdef FPC_ProcVar}@{$endif}ctan, {$ifdef FPC_ProcVar}@{$endif}mpc_tan, -2, 2, -2, 2);
  testcfun('tan', {$ifdef FPC_ProcVar}@{$endif}ctan, {$ifdef FPC_ProcVar}@{$endif}mpc_tan, -a, a, -20, 20);
end;


{---------------------------------------------------------------------------}
procedure test_tanh;
var
  a: extended;
begin
  sep_line;
  a := 0.5*ln_MaxExt;
  testcfun('tanh', {$ifdef FPC_ProcVar}@{$endif}ctanh, {$ifdef FPC_ProcVar}@{$endif}mpc_tanh, -2, 2, -2, 2);
  testcfun('tanh', {$ifdef FPC_ProcVar}@{$endif}ctanh, {$ifdef FPC_ProcVar}@{$endif}mpc_tanh, -20, -20, -a, a);
end;


{---------------------------------------------------------------------------}
procedure test_arccoth;
begin
  sep_line;
  testcfun('arccoth', {$ifdef FPC_ProcVar}@{$endif}carccoth, {$ifdef FPC_ProcVar}@{$endif}mpc_arccoth, -2, 2, -2, 2);
  testcfun('arccoth', {$ifdef FPC_ProcVar}@{$endif}carccoth, {$ifdef FPC_ProcVar}@{$endif}mpc_arccoth, 1-d1, 1+d1,-d1,d1);
  testcfun('arccoth', {$ifdef FPC_ProcVar}@{$endif}carccoth, {$ifdef FPC_ProcVar}@{$endif}mpc_arccoth, -1000,1000,-1000,1000);
end;


{---------------------------------------------------------------------------}
procedure test_arccothc;
begin
  sep_line;
  testcfun('arccothc', {$ifdef FPC_ProcVar}@{$endif}carccothc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccothc, -2, 2, -2, 2);
  testcfun('arccothc', {$ifdef FPC_ProcVar}@{$endif}carccothc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccothc, 1-d1, 1+d1,-d1,d1);
  testcfun('arccothc', {$ifdef FPC_ProcVar}@{$endif}carccothc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccothc, -1e3,1e3,-1e3,1e3);
end;


{---------------------------------------------------------------------------}
procedure test_arccot;
begin
  sep_line;
  testcfun('arccot', {$ifdef FPC_ProcVar}@{$endif}carccot, {$ifdef FPC_ProcVar}@{$endif}mpc_arccot, -2, 2, -2, 2);
  testcfun('arccot', {$ifdef FPC_ProcVar}@{$endif}carccot, {$ifdef FPC_ProcVar}@{$endif}mpc_arccot, -d1,d1, 1-d1, 1+d1);
  testcfun('arccot', {$ifdef FPC_ProcVar}@{$endif}carccot, {$ifdef FPC_ProcVar}@{$endif}mpc_arccot, -1000,1000,-1000,1000);
end;


{---------------------------------------------------------------------------}
procedure test_arccotc;
begin
  sep_line;
  testcfun('arccotc', {$ifdef FPC_ProcVar}@{$endif}carccotc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccotc, -2, 2, -2, 2);
  testcfun('arccotc', {$ifdef FPC_ProcVar}@{$endif}carccotc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccotc, -d1, +d1, 1-d1, 1+d1);
  testcfun('arccotc', {$ifdef FPC_ProcVar}@{$endif}carccotc, {$ifdef FPC_ProcVar}@{$endif}mpc_arccotc, -1000,1000,-1000,1000);
end;


{---------------------------------------------------------------------------}
procedure test_sqrt1mz2;
begin
  sep_line;
  testcfun('sqrt1mz2', {$ifdef FPC_ProcVar}@{$endif}csqrt1mz2, {$ifdef FPC_ProcVar}@{$endif}mpc_sqrt1mz2,
           -1e-7, 1e-7, -1e-7, 1e-7);
  testcfun('sqrt1mz2', {$ifdef FPC_ProcVar}@{$endif}csqrt1mz2, {$ifdef FPC_ProcVar}@{$endif}mpc_sqrt1mz2,
           1-1e-7, 1+1e-7, -1e-7, 1e-7);
  testcfun('sqrt1mz2', {$ifdef FPC_ProcVar}@{$endif}csqrt1mz2, {$ifdef FPC_ProcVar}@{$endif}mpc_sqrt1mz2,
           -1e7, 1e7, -1e7, 1e7);
  testcfun('sqrt1mz2', {$ifdef FPC_ProcVar}@{$endif}csqrt1mz2, {$ifdef FPC_ProcVar}@{$endif}mpc_sqrt1mz2,
           -1e2500, 1e2500, -1e2500, 1e2500);
  testcfun('sqrt1mz2', {$ifdef FPC_ProcVar}@{$endif}csqrt1mz2, {$ifdef FPC_ProcVar}@{$endif}mpc_sqrt1mz2,
           -1000, 1000, -1e2500, 1e2500);
  testcfun('sqrt1mz2', {$ifdef FPC_ProcVar}@{$endif}csqrt1mz2, {$ifdef FPC_ProcVar}@{$endif}mpc_sqrt1mz2,
           -1e2500, 1e2500, -1000, 1000);

end;


{---------------------------------------------------------------------------}
procedure test_agm1;
var
  a: extended;
begin
  sep_line;
  a := sqrt(Sqrt_MaxExt);
  testcfun('agm1', {$ifdef FPC_ProcVar}@{$endif}cagm1, {$ifdef FPC_ProcVar}@{$endif}mpc_agm1, -100, 100, -100, 100);
  testcfun('agm1', {$ifdef FPC_ProcVar}@{$endif}cagm1, {$ifdef FPC_ProcVar}@{$endif}mpc_agm1, -a, a, -a, a);
end;


label
  done;

begin
  mpf_set_default_prec(PREX);

  writeln('Test AMCmplx V', AMCmplx_Version,' with MPArith V', MP_VERSION, '   (c) 2014-2016 W.Ehrhardt');
  writeln('Current mp_float default bit precision = ', mpf_get_default_prec,
          ',  decimal precision = ', mpf_get_default_prec*ln(2)/ln(10):1:1);
  writeln('Machine eps for extended =', eps_x:20);
  mp_show_plus := true;

  mpf_initp2(b,c,mpf_get_default_prec);
  mpc_initp4(a,mx,mf,mmax,mpf_get_default_prec);

  d1 := ldexp(1,-16);

{
  test_ln1p;
  goto done;
}

  test_agm1;
  test_arccos;
  test_arccosh;
  test_arccot;
  test_arccotc;
  test_arccoth;
  test_arccothc;
  test_arccsc;
  test_arccsch;
  test_arcsec;
  test_arcsech;
  test_arcsin;
  test_arcsinh;
  test_arctan;
  test_arctanh;
  test_cos;
  test_cosh;
  test_cot;
  test_coth;
  test_csc;
  test_csch;
  test_exp;
  test_expm1;
  test_exp2;
  test_exp10;
  test_ln1p;
  test_ln;
  test_sec;
  test_sech;
  test_sin;
  test_sinh;
  test_sqr;
  test_sqrt;
  test_sqrt1mz2;
  test_tan;
  test_tanh;

done:
  mpf_clear2(b,c);
  mpc_clear4(a,mx,mf,mmax);

  {$ifdef MPC_Diagnostic}
    writeln;
    writeln;
    mp_dump_meminfo;
    mp_dump_diagctr;
  {$endif}

end.
