{Test program for MPArith/mp_rat, (c) W.Ehrhardt 2004-2012}
{Many test cases derived from M.J. Fromberger's IMath }
{Some test files are calculated with ARIBAS or Pari/GP}

program t_chkrat;


{$i STD.INC}
{$i mp_conf.inc}

{$x+}  {pchar I/O}
{$i+}  {RTE on I/O error}

{$ifndef FPC}
{$N+}
{$endif}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  BTypes,mp_types, mp_base,
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  {$ifdef MPC_Diagnostic}
    mp_supp,
  {$endif}
  mp_ratio;

var
  a,b,c,t: mp_rat;
  x,y: mp_int;

const
  LMAX=32000;

var
  tf: text;
  totalfailed: longint;
  line: array[0..LMAX] of char8;

const
  Prefix = 'm#';

const
  fn_add  = Prefix+'radd.dat';
  fn_div  = Prefix+'rdiv.dat';
  fn_expt = Prefix+'rexpt.dat';
  fn_inv  = Prefix+'rdiv.dat';
  fn_mul  = Prefix+'rmul.dat';
  fn_sub  = Prefix+'rsub.dat';
  fn_flt  = Prefix+'rflt.dat';


{---------------------------------------------------------------------------}
function open_file(const name: mp_string): boolean;
  {-Open data file (with error message), true if error}
begin
  {$i-}
  assign(tf,string(name));
  reset(tf);
  {$i+}
  open_file := false;
  if IOResult<>0 then begin
    writeln('File reset error: ', name);
    open_file := true;
  end;
end;


{---------------------------------------------------------------------------}
procedure readline;
  {-Read a line into buffer}
var
  c: char8;
  i: integer;
begin
  readln(tf,line);
  for i:=0 to sizeof(line)-1 do begin
    c := line[i];
    if (c=#10) or (c=#13) or (c=#0) then begin
      line[i] := #0;
      break;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure add_test;
var
  i,f: integer;
begin
  writeln('mpr_add test: ', fn_add);
  if open_file(fn_add) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    case random(3) of
      0: begin
           mpr_copy(a,t);
           mpr_add(t,b,t);
         end;
      1: begin
           mpr_copy(b,t);
           mpr_add(a,t,t);
         end;
      2: mpr_add(a,b,t);
    end;

    if not mpr_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', mpr_decimal(a));
      writeln(' b: ', mpr_decimal(b));
      writeln(' c: ', mpr_decimal(c));
      writeln(' t: ', mpr_decimal(t));
    end;
    inc(i);
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure sub_test;
var
  i,f: integer;
begin
  writeln('mpr_sub test: ', fn_sub);
  if open_file(fn_sub) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    case random(3) of
      0: begin
           mpr_copy(a,t);
           mpr_sub(t,b,t);
         end;
      1: begin
           mpr_copy(b,t);
           mpr_sub(a,t,t);
         end;
      2: mpr_sub(a,b,t);
    end;
    if not mpr_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure mul_test;
var
  i,f: integer;
begin
  writeln('mpr_mul test: ', fn_mul);
  if open_file(fn_mul) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    case random(3) of
      0: begin
           mpr_copy(a,t);
           mpr_mul(t,b,t);
         end;
      1: begin
           mpr_copy(b,t);
           mpr_mul(a,t,t);
         end;
      2: mpr_mul(a,b,t);
    end;
    if not mpr_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', mpr_decimal(a));
      writeln(' b: ', mpr_decimal(b));
      writeln(' c: ', mpr_decimal(c));
      writeln(' t: ', mpr_decimal(t));
      inc(f);
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure div_test;
var
  i,f: integer;
begin
  writeln('mpr_div test: ', fn_div);
  if open_file(fn_div) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    case random(3) of
      0: begin
           mpr_copy(a,t);
           mpr_div(t,b,t);
         end;
      1: begin
           mpr_copy(b,t);
           mpr_div(a,t,t);
         end;
      2: mpr_div(a,b,t);
    end;

    if not mpr_is_eq(c,t) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', mpr_decimal(a));
      writeln(' b: ', mpr_decimal(b));
      writeln(' c: ', mpr_decimal(c));
      writeln(' t: ', mpr_decimal(t));
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure add1_test;
var
  i,f: integer;
begin
  writeln('mpr_add1 test: ', fn_add);
  if open_file(fn_add) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    if mp_is1(b.den) then begin
      mpr_add_mpi(a,b.num,t);
      if not mpr_is_eq(t,c) then begin
        writeln('Failed! #',i);
        inc(f);
        writeln(' a: ', mpr_decimal(a));
        writeln(' b: ', mpr_decimal(b));
        writeln(' c: ', mpr_decimal(c));
        writeln(' t: ', mpr_decimal(t));
      end;
      inc(i);
    end;
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure sub1_test;
var
  i,f: integer;
begin
  writeln('mpr_sub1 test: ', fn_sub);
  if open_file(fn_sub) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    if mp_is1(b.den) then begin
      mpr_sub_mpi(a, b.num,t);
      if not mpr_is_eq(t,c) then begin
        writeln('Failed! #',i);
        inc(f);
      end;
      inc(i)
    end
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure mul1_test;
var
  i,f: integer;
begin
  writeln('mpr_mul1 test: ', fn_mul);
  if open_file(fn_mul) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    if mp_is1(b.den) then begin
      mpr_mul_mpi(a,b.num,t);
      if not mpr_is_eq(t,c) then begin
        writeln('Failed! #',i);
        inc(f);
        writeln(' a: ', mpr_decimal(a));
        writeln(' b: ', mpr_decimal(b));
        writeln(' c: ', mpr_decimal(c));
        writeln(' t: ', mpr_decimal(t));
        inc(f);
      end;
      inc(i)
    end;
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure div1_test;
var
  i,f: integer;
begin
  writeln('mpr_div1 test: ', fn_div);
  if open_file(fn_div) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    if mp_is1(b.den) then begin
      mpr_div_mpi(a,b.num,t);
      if not mpr_is_eq(c,t) then begin
        writeln('Failed! #',i);
        inc(f);
        writeln(' a: ', mpr_decimal(a));
        writeln(' b: ', mpr_decimal(b));
        writeln(' c: ', mpr_decimal(c));
        writeln(' t: ', mpr_decimal(t));
      end;
      inc(i)
    end;
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure divint_test;
var
  i,f: integer;
  k: longint;
begin
  writeln('mpr_div_int test: ', fn_div);
  if open_file(fn_div) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    if (mp_is1(b.den)) and mp_is_longint(b.num,k) then begin
      mpr_div_int(a,k,t);
      if not mpr_is_eq(c,t) then begin
        writeln('Failed! #',i);
        inc(f);
        writeln(' a: ', mpr_decimal(a));
        writeln(' k: ', k);
        writeln(' c: ', mpr_decimal(c));
        writeln(' t: ', mpr_decimal(t));
      end;
      inc(i)
    end;
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure mulint_test;
var
  i,f: integer;
  k: longint;
begin
  writeln('mpr_mul_int test: ', fn_mul);
  if open_file(fn_mul) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    if (mp_is1(b.den)) and mp_is_longint(b.num,k) then begin
      mpr_mul_int(a,k,t);
      if not mpr_is_eq(t,c) then begin
        writeln('Failed! #',i);
        inc(f);
        writeln(' a: ', mpr_decimal(a));
        writeln(' k: ', k);
        writeln(' c: ', mpr_decimal(c));
        writeln(' t: ', mpr_decimal(t));
        inc(f);
      end;
      inc(i)
    end;
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure expt_test;
var
  i,f: integer;
  x: longint;
begin
  writeln('mpr_expt test: ', fn_expt);
  if open_file(fn_expt) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    {write(i,#13);}
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readln(tf, x);
    readline;
    mpr_read_radix(c,line,10);
    case random(2) of
      0: begin
           mpr_copy(a,t);
           mpr_expt(t,x,t);
         end;
      1: mpr_expt(a,x,t);
    end;
    if not mpr_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure inv_test;
var
  i,f: integer;
begin
  {use div data for inv test:  a/b = a*t = c  with t=inv(b)}
  writeln('mpr_inv test: ', fn_inv);
  if open_file(fn_inv) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    case random(2) of
      0: begin
           mpr_copy(b,t);
           mpr_inv(t,t);
         end;
      1: begin
           mpr_inv(b,t);
         end;
    end;
    mpr_mul(a,t,t);
    if not mpr_is_eq(c,t) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', mpr_decimal(a));
      writeln(' b: ', mpr_decimal(b));
      writeln(' c: ', mpr_decimal(c));
      writeln(' t: ', mpr_decimal(t));
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure tofloat_test;
var
  i,f: integer;
  radix, prec: word;
  t1,t2,r1,r2,tt: mp_string;
begin
  writeln('mpr_tofloat test: ', fn_flt);
  if open_file(fn_flt) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readln(tf, radix);
    readln(tf, prec);
    readln(tf, t1);
    readln(tf, r1);
    mp_roundfloat := false;
    t2 := mpr_tofloat_str(a, radix, prec);
    tt := mpr_tofloat_str(a, radix, prec+10);
    mp_roundfloat := true;
    r2 := mpr_tofloat_str(a, radix, prec);
    if (t1<>t2) or (r1<>r2) then begin
      writeln('Failed! #',i);
      writeln(' a: ', mpr_decimal(a));
      writeln(' r: ', radix);
      writeln(' p: ', prec);
      writeln('t1: ', t1);
      writeln('t2: ', t2);
      writeln('r1: ', r1);
      writeln('r2: ', r2);
      writeln('tt: ', tt);
      writeln('CM:', mp_ucrmap[radix-1], '    CH:', mp_ucrmap[radix div 2]);
      inc(f);
    end;
    if i and 63 = 0 then write('.');
    inc(i)
  end;
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure testradixconv;
  {-Test input/output radix conversion}
var
  cnt,f, radix,digs: word;
  x,y: mp_int;
{$ifndef BIT16}
  const dmax=1000;
  var astr: ansistring;
{$else}
  {$ifdef IN_IDE}
    const dmax=100;
  {$else}
    const dmax=200;
  {$endif}
{$endif}
begin
  mp_init2(x,y);
  writeln('Test radix conversion (MaxRadix=',MaxRadix,'):  ');
  f := 0;
  cnt := 0;
  for radix:=2 to MaxRadix do begin
    digs := 2;
    while digs<=dmax do begin
      mp_rand_radix(x,radix,1+random(digs));
      if random(2)=0 then mp_chs(x,x);
      mp_rand_radix(y,radix,1+random(digs));
      mpr_set(a,x,y);
      mp_show_plus := random(2)=0;
      mp_uppercase := random(2)=0;
      mpr_toradix_n(a, line, radix, LMAX);
      {$ifndef BIT16}
        astr := mpr_radix_astr(a,radix);
        mpr_read_radix(b, pchar8(astr),radix);
        inc(cnt);
        if mpr_is_ne(a,b) then begin
          writeln('Diff for mpr_radix_astr: Radix=',radix, '  digits=',digs);
          write('a='); mpr_output_radix(a,radix); writeln;
          write('b='); mpr_output_radix(b,radix); writeln;
          inc(f);
        end;
      {$endif}
      inc(cnt);
      mpr_toradix_n(a, line, radix, LMAX);
      mpr_read_radix(b,line,radix);
      if mpr_is_ne(a,b) then begin
        writeln('Diff mpr_toradix_n: Radix=',radix, '  digits=',digs);
        write('a='); mpr_output_radix(a,radix); writeln;
        write('b='); mpr_output_radix(b,radix); writeln;
        inc(f);
      end;
      digs := 2 + 14*digs div 10;
    end;
    write('.');
  end;
  mp_clear2(x,y);
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure testfloatconv;
  {-Test input/output radix conversion}
var
  cnt,f, radix,digs,prec: word;
  x,y: mp_int;
{$ifndef BIT16}
  const dmax=2000;
{$else}
  {$ifdef IN_IDE}
    const dmax=200;
  {$else}
    const dmax=400;
  {$endif}
{$endif}
begin
  mp_init2(x,y);
  writeln('Test float conversion:  ');
  f := 0;
  cnt := 0;
  for radix:=2 to MaxRadix do begin
    digs := 2;
    while digs<=dmax do begin
      inc(cnt);
      mp_roundfloat := random(2)=1;
      mp_rand_radix(x,radix,1+random(digs));
      if random(2)=0 then mp_chs(x,x);
      prec := 1+random(digs);
      mp_rand_radix(y,radix,prec);
      mpr_set(a,x,y);
      mp_show_plus := random(2)=0;
      mp_uppercase := random(2)=0;
      mpr_tofloat_n(a, line, radix, prec, LMAX);
      mpr_read_float_radix(b,line,radix);
      mpr_sub(a,b,c);
      mpr_abs(c,c);
      mp_set(t.num,1);
      mp_set_pow(t.den, radix, prec);
      if mp_roundfloat then mp_shl(t.den, 1, t.den);
      if mpr_is_gt(c,t) then begin
        writeln('Diff Radix=',radix, '  digits=',digs);
        write('a='); mpr_output_radix(a,radix); writeln;
        write('b='); mpr_output_radix(b,radix); writeln;
        write('c='); writeln(mpr_tofloat_str(c,radix,prec+3));
        write('t='); writeln(mpr_tofloat_str(t,radix,prec+3));
        writeln('float=',line);
        inc(f);
      end;
      digs := 2 + 14*digs div 10;
    end;
    write('.');
  end;
  mp_clear2(x,y);
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure t_read_double;
const
  EpsDbl = 2.220446049250313E-16;   { 2^(-52) }
  MaxDbl = 1.797693134862315E+308;  { 2^1024 }
  MinDbl = 2.225073858507202E-308;  { 2^(-1022) }
var
  i,e,f,cnt: integer;
  h: double;

  procedure onetest(d: double);
  var
    g: double;
  begin
    inc(cnt);
    mpr_read_double(a,d);
    g := mpr_todouble(a);
    if abs(d-g)>EpsDbl*abs(d) then begin
      writeln('d=',d, '    g=',g);
      inc(f);
    end;
  end;
begin
  writeln('Test read_double:  ');
  f := 0;
  cnt := 0;
  onetest(0.0);
  onetest(EpsDbl);
  onetest(MaxDbl);
  onetest(MinDbl);
  onetest(-EpsDbl);
  onetest(-MaxDbl);
  onetest(-MinDbl);
  for i:=1 to 1000 do begin
    h := random-0.5;
    e := random(2040);
    h := ldexpd(h,e-1020);
    if (abs(h)>=MinDbl) and (h<=MaxDbl) then onetest(h);
  end;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_checksum;
  {-adler self test and checksums for i/2, i=-2..2}
const
  vec: array[1..3] of byte = ($11,$22,$33);
{$ifdef MP_32BIT}
  TV: array[-2..2] of longint = ($540006,$580007,$1E0003,$480006,$440005);  {4 byte mp_digits}
{$else}
  TV: array[-2..2] of longint = ($400006,$420007,$180003,$360006,$340005);  {2 byte mp_digits}
{$endif}
var
  i,f: integer;
  A1, AF: longint;
begin
  AF := 1;
  A1 := 1;
  f  := 0;
  s_mp_checksum(AF, @vec, sizeof(vec));
  for i:=1 to sizeof(vec) do s_mp_checksum(A1, @vec[i], 1);
  if (A1<>$00AD0067) or (AF<>$00AD0067) then begin
    writeln('Check of s_mp_checksum failed!');
    inc(f);
  end;
  for i:=-2 to 2 do begin
    mpr_set_int(a,i,2);
    AF := mpr_checksum(a);
    if AF<>TV[i] then inc(f);
  end;
  if f<>0 then begin
    writeln('No. of mp_checksum test failed: ',f);
    inc(totalfailed,f);
  end;
end;


{---------------------------------------------------------------------------}
procedure divmul2k_test;
var
  i,f: integer;
  k: longint;
{$ifdef BIT16}
const
  KR=15000;
{$else}
const
  KR=60000;
{$endif}
begin
  writeln('mpr_div/mul_2k test: ', fn_add);
  if open_file(fn_add) then exit;
  i:=0;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpr_read_radix(a,line,10);
    readline;
    mpr_read_radix(b,line,10);
    readline;
    mpr_read_radix(c,line,10);
    k := random(KR);
    k := k - (KR div 2);
    mp_2expt(x,abs(k));
    if i and 7 = 0 then begin
      if k<0 then write('-') else write('+');
    end;

    inc(i);
    mpr_mul_2k(c,k,a);
    if k<0 then mpr_div_mpi(c,x,b)
    else mpr_mul_mpi(c,x,b);
    if not mpr_is_eq(a,b) then begin
      writeln('mpr_mul_2k failed! #',i);
      inc(f);
      writeln(' a: ', mpr_decimal(a));
      writeln(' b: ', mpr_decimal(b));
      writeln(' c: ', mpr_decimal(c));
      writeln(' k: ', k);
    end;

    inc(i);
    mpr_div_2k(c,k,a);
    if k<0 then mpr_mul_mpi(c,x,b)
    else mpr_div_mpi(c,x,b);
    if not mpr_is_eq(a,b) then begin
      writeln('mpr_div_2k failed! #',i);
      inc(f);
      writeln(' a: ', mpr_decimal(a));
      writeln(' b: ', mpr_decimal(b));
      writeln(' c: ', mpr_decimal(c));
      writeln(' k: ', k);
    end;

    inc(i);
    mpr_mul_2(c,a);
    mpr_add(c,c,b);
    if not mpr_is_eq(a,b) then begin
      writeln('mpr_mul_2 failed! #',i);
      inc(f);
      writeln(' a: ', mpr_decimal(a));
      writeln(' b: ', mpr_decimal(b));
      writeln(' c: ', mpr_decimal(c));
    end;

    inc(i);
    mpr_div_2(c,a);
    mpr_div_int(c,2,b);
    if not mpr_is_eq(a,b) then begin
      writeln('mpr_div_2 failed! #',i);
      inc(f);
      writeln(' a: ', mpr_decimal(a));
      writeln(' b: ', mpr_decimal(b));
      writeln(' c: ', mpr_decimal(c));
    end;
  end;
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


label
  done;

begin

  writeln('Test of MPArith V', MP_VERSION, '   [mp_ratio]    (c) 2007-2012 W.Ehrhardt ');
  writeln('MAXDigits = ',MAXDigits, ',  MP_MAXBIT = ',MP_MAXBIT, ',  MAXRadix = ', MAXRadix);
  writeln('Karatsuba  cutoffs: mul/sqr = ',mp_mul_cutoff,'/',mp_sqr_cutoff);
  writeln('Toom-3, BZ cutoffs: mul/sqr = ',mp_t3m_cutoff,'/',mp_t3s_cutoff,  ',  div = ',mp_bz_cutoff);
  writeln;

  mpr_init4(a,b,c,t);
  mp_init2(x,y);

  totalfailed := 0;

  testradixconv;
  if totalfailed<>0 then begin
    writeln('Failures in testradixconv! Program aborted!');
    goto done;
  end;

  test_checksum;
  testfloatconv;

  mp_show_plus := false;
  mp_uppercase := true;
  mp_roundfloat := true;

  tofloat_test;

  mp_roundfloat := false;
  mp_uppercase := false;

  t_read_double;

  add_test;
  sub_test;
  mul_test;
  div_test;

  add1_test;
  sub1_test;
  mul1_test;
  div1_test;

  mulint_test;
  divint_test;
  divmul2k_test;

  inv_test;
  expt_test;

done:

  mp_clear2(x,y);
  mpr_clear4(a,b,c,t);

  if totalfailed<>0 then begin
    writeln('************************');
    writeln('****** Total failed: '#7#7#7, totalfailed);
    writeln('************************');
  end
  else writeln('All test OK: no failures');

  writeln;

  {$ifdef MPC_Diagnostic}
    mp_dump_meminfo;
    mp_dump_diagctr;
  {$endif}

end.
