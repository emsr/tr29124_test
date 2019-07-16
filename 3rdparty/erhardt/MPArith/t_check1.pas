{Part 1 of main check program for mp_ints}

unit t_check1;

{$x+}  {pchar I/O}

{$i STD.INC}
{$i mp_conf.inc}

{$i+}  {RTE on I/O error}
{$ifndef FPC}
{$N+}
{$endif}

interface

uses
  mp_types, mp_rsa;

var
  totalfailed: longint;

var
  a,b,c,d,e,s,t,z: mp_int;
  prk: TPrivateKey;

procedure testradixconv;
procedure test_checksum;

procedure add_test;
procedure sub_test;
procedure mul_test1;
procedure mul_test2;
procedure sqr_test;
procedure mod_test;
procedure gr_mod_test;
procedure reduce_test;
procedure div_test;
procedure expt_test;
procedure invmod_test;
procedure emod_test;
procedure emod2_test;
procedure red2k_test;
procedure egcd_test;
procedure jac_test;
procedure kron_test;
procedure lucas_test;
procedure coshmult_test;
procedure lucfib_test;
procedure probprime_test;
procedure nroottest;
procedure facttest;
procedure pfdu_test;
procedure rsa_test1;
procedure rsa_test2;
procedure rsa_test3;
procedure sqrtmod_test;
procedure sqrtmod2k_test;
procedure sqrtmodpq_test;
procedure cbrtmod_test;
procedure cbrtmodpq_test;
procedure bool_test;
procedure pell1_test;
procedure pell4_test;
procedure cornacchia_test;
procedure test_rnr;
procedure test_rqffu;
procedure powerd_test;
procedure binomial_test;
procedure swing_test;
procedure digitsum_test;
procedure reverse_check;


implementation

uses
  BTypes, mp_base, mp_prime, mp_modul, mp_numth, mp_pfu, pfdu;

const
  LMAX=4096;

var
  tf: text;
  line: array[0..LMAX] of char8;

const
  Prefix = 'm#';

const
  fn_add    = Prefix+'add.dat';
  fn_div    = Prefix+'div.dat';
  fn_emod   = Prefix+'emod.dat';
  fn_emod2  = Prefix+'emod2.dat';
  fn_emod3  = Prefix+'red2k.dat';
  fn_expt   = Prefix+'expt.dat';
  fn_invmod = Prefix+'invmod.dat';
  fn_mod    = Prefix+'mod.dat';
  fn_mul1   = Prefix+'mul1.dat';
  fn_mul2   = Prefix+'mul2.dat';
  fn_sqr    = Prefix+'sqr.dat';
  fn_sub    = Prefix+'sub.dat';
  fn_xgcd   = Prefix+'xgcd.dat';
  fn_jac    = Prefix+'jac.dat';
  fn_luc    = Prefix+'lucas.dat';
  fn_coshm  = Prefix+'coshm.dat';
  fn_lucfib = Prefix+'lucfib.dat';
  fn_kron   = Prefix+'kron.dat';
  fn_pell1  = Prefix+'pell1.dat';
  fn_pell4  = Prefix+'pell4.dat';
  fn_bool   = Prefix+'bool.dat';
  fn_rqffu  = Prefix+'rqffu.dat';
  fn_powerd = Prefix+'powerd.dat';


{---------------------------------------------------------------------------}
function decimal(const a: mp_int): pchar8;
  {-do not use mp_decimal, because a may have more than 255 digits}
begin
  mp_toradix_n(a, pchar8(@line), 10, sizeof(line));
  decimal := pchar8(@line);
end;


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
  writeln('mp_add test: ', fn_add);
  if open_file(fn_add) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    case random(3) of
      0: begin
           mp_copy(a,t);
           mp_add(t,b,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_add(a,t,t);
         end;
      2: mp_add(a,b,t);
    end;

    if not mp_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure sub_test;
var
  i,f: integer;
begin
  writeln('mp_sub test: ', fn_sub);
  if open_file(fn_sub) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    case random(3) of
      0: begin
           mp_copy(a,t);
           mp_sub(t,b,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_sub(a,t,t);
         end;
      2: mp_sub(a,b,t);
    end;
    if not mp_is_eq(t,c) then begin
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
procedure mul_test1;
var
  i,f: integer;
begin
  writeln('mp_mul test 1: ', fn_mul1);
  if open_file(fn_mul1) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    case random(3) of
      0: begin
           mp_copy(a,t);
           mp_mul(t,b,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_mul(a,t,t);
         end;
      2: mp_mul(a,b,t);
    end;
    if not mp_is_eq(t,c) then begin
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
procedure mul_test2;
var
  i,f: integer;
begin
  writeln('mp_mul test 2: ', fn_mul2);
  if open_file(fn_mul2) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    case random(3) of
      0: begin
           mp_copy(a,t);
           mp_mul(t,b,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_mul(a,t,t);
         end;
      2: mp_mul(a,b,t);
    end;
    if not mp_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
    end;
    inc(i);
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure sqr_test;
var
  i,f: integer;
begin
  writeln('mp_sqr test: ', fn_sqr);
  if open_file(fn_sqr) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(c,line,10);
    case random(2) of
      0: begin
           mp_copy(a,t);
           mp_sqr(t,t);
         end;
      1: mp_sqr(a,t);
    end;
    if not mp_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure mod_test;
var
  i,f: integer;
begin
  writeln('mp_mod test: ', fn_mod);
  if open_file(fn_mod) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    case random(3) of
      0: begin
           mp_copy(a,t);
           mp_mod(t,b,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_mod(a,t,t);
         end;
      2: mp_mod(a,b,t);
    end;
    if not mp_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' t: ', decimal(t));
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure gr_mod_test;
var
  i,f: integer;
begin
  writeln('mp_gr_mod test: ', fn_mod);
  if open_file(fn_mod) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    {skip cases b<2}
    if mp_cmp_d(b,1)=MP_GT then begin
      mp_gr_setup(s,b);
      mp_copy(a,t);
      mp_gr_mod(t,b,s);
      if not mp_is_eq(t,c) then begin
        writeln('Failed! #',i);
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
procedure reduce_test;
var
  i,f: integer;
begin
  writeln('mp_sqrmod/mp_reduce test: ', fn_mod);
  if open_file(fn_mod) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    {t = a mod b}
    mp_mod(a,b,t);

    {c = sqr(t) mod b}
    mp_sqrmod(t,b,c);

    {t = sqr(t) mod b with Barrett}
    mp_reduce_setup(z,b);
    mp_sqr(t,t);
    mp_reduce(t,b,z);
    if not mp_is_eq(t,c) then begin
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
procedure div_test;
var
  i,f: integer;
begin
  writeln('mp_divrem test: ', fn_div);
  if open_file(fn_div) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readline;
    mp_read_radix(d,line,10);
    case random(6) of
      0: begin
           mp_copy(a,s);
           mp_divrem(s,b,@s,@t);
         end;
      1: begin
           mp_copy(b,s);
           mp_divrem(a,s,@s,@t);
         end;
      2: begin
           mp_copy(a,t);
           mp_divrem(t,b,@s,@t);
         end;
      3: begin
           mp_copy(b,t);
           mp_divrem(a,t,@s,@t);
         end;
      4: begin
           mp_copy(a,s);
           mp_copy(b,t);
           mp_divrem(s,t,@s,@t);
         end;
      5: mp_divrem(a,b,@s,@t);

        {note that the case
           mp_copy(a,s);
           mp_copy(b,t);
           mp_divrem(t,s,@s,@t);
         is not allowed}
    end;

    if not (mp_is_eq(c,s) and mp_is_eq(t,d)) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' d: ', decimal(d));
      writeln(' s: ', decimal(s));
      writeln(' t: ', decimal(t));
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure expt_test;
var
  i,j,f: integer;
  x: longint;
begin
  writeln('mp_expt/mp_expt_int test: ', fn_expt);
  if open_file(fn_expt) then exit;
  i:=0;
  f:=0;
  j:=1;
  while not eof(tf) do begin
    inc(i);
    if i and 1 = 0 then write('.');
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    case random(3) of
      0: begin
           mp_copy(a,t);
           mp_expt(t,b,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_expt(a,t,t);
         end;
      2: mp_expt(a,b,t);
    end;
    if not mp_is_eq(t,c) then begin
      writeln('mp_expt failed! #',i);
      inc(f);
    end;
    inc(j);
    if mp_is_longint(b,x) then begin
      mp_expt_int(a,x,t);
      if not mp_is_eq(t,c) then begin
        writeln('mp_expt_int failed! #',i);
        inc(f);
      end;
      inc(j)
    end;
  end;
  writeln;
  close(tf);
  writeln('No. of checks: ',j,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure invmod_test;
var
  i,f: integer;
begin
  writeln('mp_invmod test: ', fn_invmod);
  if open_file(fn_invmod) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    {write(i,#13);}
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    case random(3) of
      0: begin
           mp_copy(a,t);
           mp_invmod(t,b,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_invmod(a,t,t);
         end;
      2: mp_invmod(a,b,t);
    end;
    if not mp_is_eq(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' t: ', decimal(t));
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure emod_test;
var
  i,f: integer;
begin
  writeln('mp_exptmod test 1: ', fn_emod);
  if open_file(fn_emod) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readline;
    mp_read_radix(d,line,10);
    case random(4) of
      0: begin
           mp_copy(a,t);
           mp_exptmod(t,b,c,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_exptmod(a,t,c,t);
         end;
      2: begin
           mp_copy(c,t);
           mp_exptmod(a,b,t,t);
         end;
      3: begin
           if b.used=1 then begin
             mp_exptmod_d(a,b.pdigits^[0],c,t);
             write('+');
           end
           else mp_exptmod(a,b,c,t);
         end;
    end;
    if not mp_is_eq(t,d) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' d: ', decimal(d));
      writeln(' t: ', decimal(t));
    end;
    inc(i);
    if (i and 15)=0 then write('.');
  end;
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure emod2_test;
var
  i,f: integer;
begin
  writeln('mp_exptmod test 2: ', fn_emod2);
  if open_file(fn_emod2) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readline;
    mp_read_radix(d,line,10);
    case random(4) of
      0: begin
           mp_copy(a,t);
           mp_exptmod(t,b,c,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_exptmod(a,t,c,t);
         end;
      2: begin
           mp_copy(c,t);
           mp_exptmod(a,b,t,t);
         end;
      3: mp_exptmod(a,b,c,t);
    end;
    if not mp_is_eq(t,d) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' d: ', decimal(d));
      writeln(' t: ', decimal(t));
    end;
    inc(i);
    if (i and 15)=0 then write('.');
  end;
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure red2k_test;
var
  i,f: integer;
begin
  writeln('mp_exptmod test 3: ', fn_emod3);
  if open_file(fn_emod3) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readline;
    mp_read_radix(d,line,10);
    case random(4) of
      0: begin
           mp_copy(a,t);
           mp_exptmod(t,b,c,t);
         end;
      1: begin
           mp_copy(b,t);
           mp_exptmod(a,t,c,t);
         end;
      2: begin
           mp_copy(c,t);
           mp_exptmod(a,b,t,t);
         end;
      3: mp_exptmod(a,b,c,t);
    end;
    if not mp_is_eq(t,d) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' d: ', decimal(d));
      writeln(' t: ', decimal(t));
    end;
    inc(i);
    if (i and 15)=0 then write('.');
  end;
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure egcd_test;
var
  i,f,k: integer;
  g,q: longint;
  x,y: mp_int;
begin
  writeln('(x)gcd and (x)lcm tests: ', fn_xgcd);
  if open_file(fn_xgcd) then exit;
  i:=0;
  k:=1;
  f:=0;
  mp_init2(x,y);
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readline;
    mp_read_radix(d,line,10);
    readline;
    mp_read_radix(e,line,10);

    if (mp_is_longint(b,q)) and (q<>0) then begin
      inc(i);
      g := mp_gcd_int(a,q);
      mp_sub_int(c,g,t);
      if not mp_is0(t) then begin
        writeln('mp_gcd_int failed! #',k);
        inc(f);
        writeln(' a: ', decimal(a));
        writeln(' b: ', decimal(b));
        writeln(' q: ', q);
        writeln(' c: ', decimal(c));
        writeln(' g: ', g);
        writeln(' t: ', decimal(t));
      end;
    end;

    {Note: mp_xgcd(a,b,@s,@t,@z) is clean against duplicated addrs of }
    {parameters (code review, pointer deref. only for final results)  }
    inc(i);
    mp_xgcd(a,b,@s,@t,@z);
    if not (mp_is_eq(c,z) and mp_is_eq(d,s) and mp_is_eq(e,t)) then begin
      writeln('mp_xgcd failed! #',k);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' d: ', decimal(d));
      writeln(' e: ', decimal(e));
      writeln(' s: ', decimal(s));
      writeln(' t: ', decimal(t));
      writeln(' z: ', decimal(z));
    end;

    {Check gcd=c and s*a + b*t = c, NOT d=s and e=t because s,t from}
    {binary xgcd do normally not equal those values from Euclid xgcd}
    inc(i);
    mp_xgcd_bin(a,b,@s,@t,@z);
    mp_mul(a,s,d);
    mp_mul(b,t,e);
    mp_add(e,d,e);
    if not (mp_is_eq(c,z) and mp_is_eq(e,c)) then begin
      writeln('mp_xgcd_bin failed! #',k);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' e: ', decimal(e));
      writeln(' s: ', decimal(s));
      writeln(' t: ', decimal(t));
      writeln(' z: ', decimal(z));
    end;

    case random(3) of
      0: begin
           mp_copy(a,z);
           mp_gcd(z,b,z);
         end;
      1: begin
           mp_copy(b,z);
           mp_gcd(a,z,z);
         end;
      2: mp_gcd(a,b,z);
    end;

    inc(i);
    if not mp_is_eq(c,z) then begin
      writeln('mp_gcd failed! #',k);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' z: ', decimal(z));
    end;

    case random(3) of
      0: begin
           mp_copy(a,c);
           mp_gcd_euclid(c,b,c);
         end;
      1: begin
           mp_copy(b,c);
           mp_gcd_euclid(a,c,c);
         end;
      2: mp_gcd_euclid(a,b,c);
    end;

    inc(i);
    if not mp_is_eq(c,z) then begin
      writeln('mp_gcd_euclid failed! #',k);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' z: ', decimal(z));
    end;

    case random(3) of
      0: begin
           mp_copy(a,c);
           mp_gcd_ml(c,b,c);
         end;
      1: begin
           mp_copy(b,c);
           mp_gcd_ml(a,c,c);
         end;
      2: mp_gcd_ml(a,b,c);
    end;

    inc(i);
    if not mp_is_eq(c,z) then begin
      writeln('mp_gcd_ml failed! #',k);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' z: ', decimal(z));
    end;

    {c = z = gcd(a,b)}
    case random(3) of
      0: begin
           mp_copy(a,c);
           mp_lcm(c,b,c);
         end;
      1: begin
           mp_copy(b,c);
           mp_lcm(a,c,c);
         end;
      2: mp_lcm(a,b,c);
    end;
    {gcd*lcm, both should be >=0}
    mp_mul(c,z,d);
    mp_mul(a,b,t);
    mp_abs(t,t);

    inc(i);
    if not mp_is_eq(d,t) then begin
      writeln('mp_lcm failed! #',k);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' d: ', decimal(z));
      writeln(' t: ', decimal(z));
    end
    else begin
      mp_xlcm(a,b,s,t,z);
      mp_mul(t,z,d);
      inc(i);
      if mp_is_ne(s,c) or mp_is_ne(d,c) then begin
        writeln('mp_xlcm failed! #',k);
        inc(f);
        writeln(' a: ', decimal(a));
        writeln(' b: ', decimal(b));
        writeln(' c: ', decimal(c));
        writeln(' d: ', decimal(d));
        writeln(' s: ', decimal(s));
        writeln(' t: ', decimal(t));
        writeln(' z: ', decimal(z));
      end;
      if not mp_is0(c) then begin
        {lcm<>0 -> t<>0, z<>0}
        mp_mod(a,t,d);
        mp_mod(b,z,e);
        inc(i);
        if (not mp_gcd1(t,z,s)) or (not mp_is0(d)) or (not mp_is0(e)) then begin
          writeln('mp_xlcm failed! #',k);
          inc(f);
          writeln(' a: ', decimal(a));
          writeln(' b: ', decimal(b));
          writeln(' c: ', decimal(c));
          writeln(' d: ', decimal(d));
          writeln(' e: ', decimal(e));
          writeln(' s: ', decimal(s));
          writeln(' t: ', decimal(t));
          writeln(' z: ', decimal(z));
        end;
      end;
    end;


    inc(k);
    if (k and 15)=0 then write('.');
  end;
  mp_clear2(x,y);
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure jac_test;
var
  i,f,j,k,m: integer;
  jac: longint;
begin
  writeln('mp_jacobi test: ', fn_jac);
  if open_file(fn_jac) then exit;
  i:=1;
  f:=0;
  k:=0;
  m:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readln(tf, j);
    {mp_jac returns integer, no check with duplicated addrs of parameters}
    jac := mp_jacobi(a,b);
    if j<>jac then begin
      writeln('mp_jacobi #',i);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', j);
      writeln(' t: ', jac);
      inc(f);
    end;
    if mp_is_longint(a,jac) then begin
      inc(k);
      if mp_error=MP_OKAY then begin
        jac := mp_jacobi_lm(jac,b);
        if j<>jac then begin
          writeln('mp_jacobi_lm #',i);
          writeln(' a: ', decimal(a));
          writeln(' b: ', decimal(b));
          writeln(' c: ', j);
          writeln(' t: ', jac);
          inc(f);
        end;
      end;
    end;
    if mp_is_longint(b,jac) then begin
      inc(m);
      if mp_error=MP_OKAY then begin
        jac := mp_jacobi_ml(a,jac);
        if j<>jac then begin
          writeln('mp_jacobi_ml #',i);
          writeln(' a: ', decimal(a));
          writeln(' b: ', decimal(b));
          writeln(' c: ', j);
          writeln(' t: ', jac);
          inc(f);
        end;
      end;
    end;
    inc(i);
    if (i and 15)=0 then write('.');
  end;
  writeln;
  close(tf);
  writeln('No. of checks: ',i,'+',k,'+',m,',  failed: ',f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure kron_test;
var
  i,f,j,k,m,n: integer;
  kron,L: longint;
begin
  writeln('mp_kronecker test: ', fn_kron);
  if open_file(fn_kron) then exit;
  i:=1;
  f:=0;
  k:=0;
  m:=0;
  n:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readln(tf, j);
    {mp_kronecker returns integer, no check with duplicated addrs of parameters}
    kron := mp_kronecker(a,b);
    if j<>kron then begin
      writeln('mp_kronecker #',i);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', j);
      writeln(' t: ', kron);
      inc(f);
    end;
    if mp_is_longint(a,L) and mp_isodd(b) and (mp_cmp_d(b,2)=MP_GT) then begin
      inc(k);
      if mp_error=MP_OKAY then begin
        L := mp_jacobi_lm(L,b);
        if j<>L then begin
          writeln('mp_jacobi_lm #',i);
          writeln(' a: ', decimal(a));
          writeln(' b: ', decimal(b));
          writeln(' c: ', j);
          writeln(' t: ', L);
          inc(f);
        end;
      end;
    end;
    if mp_is_longint(b, L) and mp_isodd(b) and (mp_cmp_d(b,2)=MP_GT) then begin
      inc(m);
      if mp_error=MP_OKAY then begin
        L := mp_jacobi_ml(a,L);
        if j<>L then begin
          writeln('mp_jacobi_ml #',i);
          writeln(' a: ', decimal(a));
          writeln(' b: ', decimal(b));
          writeln(' c: ', j);
          writeln(' t: ', L);
          inc(f);
        end;
      end;
    end;
    if mp_is_longint(a, kron) and mp_is_longint(b, L) then begin
      inc(n);
      if mp_error=MP_OKAY then begin
        L := kronecker32(kron,L);
        if j<>L then begin
          writeln('kronecker32 #',i);
          writeln(' a: ', decimal(a));
          writeln(' b: ', decimal(b));
          writeln(' c: ', j);
          writeln(' t: ', L);
          inc(f);
        end;
      end;
    end;
    inc(i);
    if (i and 15)=0 then write('.');
  end;
  writeln;
  close(tf);
  writeln('No. of checks: ',i,'+',k,'+',m,'+',n,',  failed: ',f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure lucas_test;
var
  i,f: integer;
  k: longint;
begin
  writeln('mp_lucasv test: ', fn_luc);
  if open_file(fn_luc) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readln(tf,k);
    readline;
    mp_read_radix(d,line,10);
    {t = v[a,b,k]}
    {check if mp_lucasv works if addr(t)=addr(a) or =addr(b)}

    case random(3) of
      0: begin
           mp_copy(a,t);
           mp_lucasv(t,b,k,t);;
         end;
      1: begin
           mp_copy(b,t);
           mp_lucasv(a,t,k,t);
         end;
      2: mp_lucasv(a,b,k,t);
    end;

    if not mp_is_eq(t,d) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' k: ', k);
      writeln(' d: ', decimal(d));
      writeln(' t: ', decimal(t));
    end;
    inc(i);
    if (i and 15)=0 then write('.');
  end;
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure coshmult_test;
var
  i,f: integer;
begin
  writeln('s_mp_coshmult test: ', fn_coshm);
  if open_file(fn_coshm) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readline;
    mp_read_radix(d,line,10);
    mp_reduce_setup(s, c);

    case random(4) of
      0: begin
           mp_copy(a,t);
           s_mp_coshmult(t,b,c,s,t);
         end;
      1: begin
           mp_copy(b,t);
           s_mp_coshmult(a,t,c,s,t);
         end;
      2: begin
           mp_copy(c,t);
           s_mp_coshmult(a,b,t,s,t);
         end;
      3: s_mp_coshmult(a,b,c,s,t);
    end;
    if not mp_is_eq(t,d) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' d: ', decimal(d));
      writeln(' t: ', decimal(t));
    end;
    inc(i);
    if (i and 15)=0 then write('.');
  end;
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure lucfib_test;
var
  i,f: integer;
  n: longint;
begin
  writeln('Lucas/Fibonacci test: ', fn_lucfib);
  if open_file(fn_lucfib) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    {write(i,#13);}
    readln(tf);
    readln(tf,n);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    mp_lucas(n,t);
    mp_fib(n,z);
    if mp_is_ne(t,a) or mp_is_ne(z,b)  then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' n: ', n);
      writeln(' a: ', decimal(a));
      writeln(' t: ', decimal(t));
      writeln(' b: ', decimal(b));
      writeln(' z: ', decimal(z));
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure probprime_test;
var
  pa,pb,pc: boolean;
  f: integer;
begin
  {Pinch:  68528663395046912244223605902738356719751082784386681071 is }
  {a strong pseudoprime for all bases up to 100, that is, for 25 primes}

  {HAC 4.7: Arnault [56] found the following 46-digit composite integer}
  {n = 1195068768795265792518361315725116351898245581 that is a strong }
  {pseudoprime to all the 11 prime bases up to 31.                     }

  {a,b: composite; c: prime}

  mp_read_decimal(a,'68528663395046912244223605902738356719751082784386681071');
  mp_read_decimal(b,'1195068768795265792518361315725116351898245581');
  mp_mersenne(521,c);
  f := 0;

  write('Probable prime test  ');

  mp_miller_rabin(a, 0, pa); write('.');
  mp_miller_rabin(b, 0, pb); write('.');
  mp_miller_rabin(c, 0, pc); write('.');
  if pa or pb or (not pc) then begin
    writeln('mp_miller_rabin failed!');
    if pa then begin
      writeln(' a: ', decimal(a));
      inc(f);
    end;
    if pa then begin
      writeln(' b: ', decimal(b));
      inc(f);
    end;
    if not pc then begin
      writeln(' c: Mersenne(521)');
      inc(f);
    end;
  end;

  pa := mp_is_pprime(a); write('.');
  pb := mp_is_pprime(b); write('.');
  pc := mp_is_pprime(c); writeln('.');
  if pa or pb or (not pc) then begin
    writeln('mp_is_pprime failed!');
    if pa then begin
      writeln(' a: ', decimal(a));
      inc(f);
    end;
    if pa then begin
      writeln(' b: ', decimal(b));
      inc(f);
    end;
    if not pc then begin
      writeln(' c: Mersenne(521)');
      inc(f);
    end;
  end;
  writeln('No. of checks: 6, failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure nroottest;
{$ifndef BIT16}
  {.$define BIGTEST}
  {$ifdef BIGTEST}
   const
     maxorder = 300; { check mp_n_root(a^n,n) up to order(a) = 10^maxorder}
     maxpower = 50;  { max n}
     ProgMask = 15;
  {$else}
  const
    maxorder = 120; { check mp_n_root(a^n,n) up to order(a) = 10^maxorder}
    maxpower = 27;  { max n}
    ProgMask = 15;
  {$endif}
{$else}
const
  maxorder = 60; { check mp_n_root(a^n,n) up to order(a) = 10^maxorder}
  maxpower = 21; { max n}
  ProgMask = 15;
{$endif}
var
  f  : integer;
  d  : word;
  Cnt: longint;

  procedure checkorder;
    {-Check order d}
  var
    n, k: integer;
    p: longint;
    pp,ip: boolean;
  begin
    for k:=1 to 8 do begin
      repeat
        mp_rand_radix(a, 10, d);
      until not mp_is1(a);
      for n:=1 to MaxPower do begin
        if odd(n) and (random(2)=0) then a.sign := MP_NEG else a.sign := MP_ZPOS;
        mp_expt_w(a,n,b);
        {b = a^n}
        inc(cnt);
        if random(2)=0 then begin
          mp_n_root(b,n,c);
        end
        else begin
          mp_copy(b,c);
          mp_n_root(c,n,c);
        end;
        if not mp_is_eq(a,c) then begin
          writeln('mp_n_root failed!');
          inc(f);
          writeln('a=',mp_decimal(a));
          writeln('b=',mp_decimal(b));
          writeln('c=',mp_decimal(c));
          writeln('n=',n);
        end;

        {test mp_is_power}
        mp_is_power(b,z,p);
        inc(cnt);
        if p=1 then begin
          if n<>1 then begin
            writeln('mp_is_power failed!');
            inc(f);
            writeln('n=',n, ' but p=1');
            writeln('a=',mp_decimal(a));
            writeln('b=',mp_decimal(b));
            writeln('z=',mp_decimal(z));
          end;
        end
        else begin
          mp_expt_int(z,p,t);
          if mp_is_ne(t,b) then begin
            writeln('mp_is_power failed! a^n = b <> z^p');
            writeln('a=',mp_decimal(a));
            writeln('b=',mp_decimal(b));
            writeln('n=',n);
            writeln('z=',mp_decimal(z));
            writeln('p=',p);
            inc(f);
          end;
        end;

        {test mp_is_power_max}
        mp_is_power_max(b,z,p);
        inc(cnt);
        if p=1 then begin
          if n<>1 then begin
            writeln('mp_is_power_max failed!');
            inc(f);
            writeln('n=',n, ' but p=1');
            writeln('a=',mp_decimal(a));
            writeln('b=',mp_decimal(b));
            writeln('z=',mp_decimal(z));
          end;
        end
        else begin
          mp_expt_int(z,p,t);
          if mp_is_ne(t,b) then begin
            writeln('mp_is_power_max failed! a^n = b <> z^p');
            writeln('a=',mp_decimal(a));
            writeln('b=',mp_decimal(b));
            writeln('n=',n);
            writeln('z=',mp_decimal(z));
            writeln('p=',p);
            inc(f);
          end;
        end;

        {Test mp_is_PrimePower}
        inc(Cnt);
        pp := mp_is_PrimePower(b,z,p);
        mp_expt_int(z,p,t);
        if mp_is_ne(t,b) then begin
          writeln('mp_is_primepower failed! a^n = b <> z^p');
          writeln('a=',mp_decimal(a));
          writeln('b=',mp_decimal(b));
          writeln('n=',n);
          writeln('z=',mp_decimal(z));
          writeln('p=',p);
          inc(f);
        end
        else begin
          ip := mp_is_pprime(z);
          if ip then begin
            if (p=1) and pp then begin
              writeln('mp_is_primepower failed!');
              inc(f);
              writeln('a=z^1, z prime and mp_is_primepower=true');
              writeln('a=',mp_decimal(a));
              writeln('z=',mp_decimal(z));
            end
            else if (p>1) and (not pp) then begin
              writeln('mp_is_primepower failed!');
              inc(f);
              writeln('a=z^p, z prime, p>1 and mp_is_primepower=false');
              writeln('a=',mp_decimal(a));
              writeln('z=',mp_decimal(z));
              writeln('p=',p);
            end;
          end
          else begin
            if pp then begin
              writeln('mp_is_primepower failed!');
              inc(f);
              writeln('a=z^p, z not prime, but mp_is_primepower=true');
              writeln('a=',mp_decimal(a));
              writeln('b=',mp_decimal(b));
              writeln('n=',n);
              writeln('z=',mp_decimal(z));
              writeln('p=',p);
            end;
          end;
        end;
        inc(Cnt);
        if odd(n) then mp_dec(b) else mp_inc(b);
        mp_n_root2(b,n,c,@t);
        mp_expt_w(c,n,s);
        mp_add(s,t,s);
        if mp_is_ne(s,b) then begin
          writeln('mp_n_root2 failed! s=c^n+t <> b');
          writeln('s=',mp_decimal(s));
          writeln('c=',mp_decimal(c));
          writeln('n=',n);
          writeln('t=',mp_decimal(t));
          writeln('b=',mp_decimal(b));
          inc(f);
        end;
        if n and ProgMask = 1 then begin
          if a.sign=MP_NEG then write('-') else write('+');
        end;
      end;
    end;
  end;

begin
  d := 1;
  f := 0;
  cnt := 0;
  writeln('mp_n_root/mp_is_(prime)power test:  ');
  while d <= MaxOrder do begin
    checkorder;
    d := 1+ 3*d div 2;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure mp_fact0(N: word; var a: mp_int);
  {-Calculate a = factorial(n)}
var
  d: mp_int;
  i: word;
begin
  if mp_error<>MP_OKAY then exit;
  mp_set(a,1);
  mp_init(d);
  if mp_error=MP_OKAY then begin
    mp_set(d, 1);
    {multiply all remaining factors}
    for i:=2 to N do begin
      mp_inc(d);
      mp_mul(a,d,a);
      if mp_error<>MP_OKAY then break;
    end;
    mp_clear(d);
  end;
end;


{---------------------------------------------------------------------------}
procedure facttest;
const
  AF = 9;
  fn : array[1..AF] of word = (0,1,2,15,99,255,512,1111,2222);
var
  f,i: integer;
  n  : word;
begin
  f := 0;
  write('mp_fact test:  ');
  for i:=1 to AF do begin
    n := fn[i];
    mp_fact0(n,a);
    mp_fact(n,b);
    if mp_is_ne(a,b) then begin
      writeln('mp_fact failed for n = ',n);
      inc(f);
    end;
    s_mp_recsplit_fact(n,b);
    if mp_is_ne(a,b) then begin
      writeln('s_mp_recsplit_fact failed for n = ',n);
      inc(f);
    end;
    s_mp_borsch_fact(n,b);
    if mp_is_ne(a,b) then begin
      writeln('s_mp_borsch_fact failed for n = ',n);
      inc(f);
    end;
    write('.');
  end;
  writeln;
  writeln('No. of checks: ',3*AF,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure pfdu_test;
var
  ctx: pfd_ctx;
  ps:  boolean;
begin
  writeln('prime factor decomposition test for 36! + 1:');
  ps := mp_show_progress;
  mp_show_progress := false;
  pfd_initialize(ctx);
  mp_fact(36,a);
  mp_inc(a);
  pfd_factor(ctx, a);
(*
  writeln('prime factor decomposition test for (2^31-1)^15');
  mp_mersenne(31,a);
  mp_expt_w(a,15,a);
  pfd_factor(ctx, a);
*)
  pfd_finalize(ctx);
  mp_show_progress := ps;
  writeln;
end;


{---------------------------------------------------------------------------}
procedure rsa_test1;
{$ifndef BIT16}
const jmax=16;
{$else}
const jmax=8;
{$endif}
var
  osize: word;
  f,i,j,cnt: integer;
  ml,plen,clen: word;
  pt: array[0..128] of char8;
  ct: array[0..128] of char8;
  OK,fail: boolean;
begin
  cnt := 0;
  f := 0;
  write('RSA test 1:  ');
  for j:=1 to jmax do begin
    write('.');
    osize := RSA_MINSIZE + 4*j;
    if odd(j) then mp_set_int(e,65537)
    else mp_set_int(e, nextprime32(10*j));
    mp_rsa_keygen1(e, osize, d, z);

    inc(cnt);
    mp_rsa_recover_pq(z,e,d,s,t,fail);
    mp_mul(s,t,b);
    if fail or mp_is_ne(z,b) or mp_is1(s) or mp_is1(t) then begin
      inc(f);
      writeln('mp_rsa_recover_pq failed for j=',j);
    end;

    ml := mp_pkcs1v15_maxlen(z);
    while ml>1 do begin
      inc(cnt);
      for i:=0 to ml-1 do line[i]:=char8(random(256));
      mp_pkcs1v15_encrypt(e, z, 0, @line, @ct, ml, sizeof(ct), clen);
      mp_pkcs1v15_decrypt(d, z, @ct, @pt, clen, sizeof(pt), plen);
      OK := plen=ml;
      if OK then begin
        for i:=0 to ml-1 do begin
          if pt[i]<>line[i] then begin
            OK := false;
            break;
          end;
        end;
      end;
      if not OK then inc(f);
      dec(ml);
    end;
  end;
  writeln;
  writeln('Max bits size = ',8*osize);
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure rsa_test2;
{$ifndef BIT16}
const jmax=16;
{$else}
const jmax=8;
{$endif}
var
  osize: word;
  f,i,j,cnt: integer;
  ml,plen,clen: word;
  pt: array[0..128] of char8;
  ct: array[0..128] of char8;
  OK,fail: boolean;
begin
  cnt := 0;
  f := 0;
  write('RSA test 2:  ');
  for j:=1 to jmax do begin
    write('.');
    osize := RSA_MINSIZE + 4*j;
    if odd(j) then mp_set_int(e,65537)
    else mp_set_int(e, nextprime32(10*j));
    {generate CRT private key}
    mp_rsa_keygen2(e, osize, z, prk);
    {Calc single decryption exponent from CRT private key}
    mp_rsa_calc_d(e,prk,d);

    inc(cnt);
    mp_rsa_recover_pq(z,e,d,s,t,fail);
    mp_mul(s,t,b);
    if fail or mp_is_ne(z,b) or mp_is1(s) or mp_is1(t) then begin
      inc(f);
      writeln('mp_rsa_recover_pq failed for j=',j);
    end;

    inc(cnt);
    mp_rsa_recover_pq2(z,e,prk.dp,s,t,fail);
    mp_mul(s,t,b);
    if fail or mp_is_ne(z,b) or mp_is1(s) or mp_is1(t) then begin
      inc(f);
      writeln('mp_rsa_recover_pq2 failed for j=',j);
    end;

    ml := mp_pkcs1v15_maxlen(z);
    while ml>1 do begin
      inc(cnt);
      for i:=0 to ml-1 do line[i]:=char8(random(256));
      mp_pkcs1v15_encrypt(e, z, 0, @line, @ct, ml, sizeof(ct), clen);
      {Test CRT decryption}
      fillchar(pt, sizeof(pt), 0);
      mp_pkcs1v15_decrypt2(prk, z, @ct, @pt, clen, sizeof(pt), plen);
      OK := plen=ml;
      if OK then begin
        for i:=0 to ml-1 do begin
          if pt[i]<>line[i] then begin
            OK := false;
            break;
          end;
        end;
      end;
      if not OK then inc(f);
      {Test standard decryption}
      fillchar(pt, sizeof(pt), 0);
      mp_pkcs1v15_decrypt(d, z, @ct, @pt, clen, sizeof(pt), plen);
      OK := plen=ml;
      inc(cnt);
      if OK then begin
        for i:=0 to ml-1 do begin
          if pt[i]<>line[i] then begin
            OK := false;
            break;
          end;
        end;
      end;
      if not OK then inc(f);
      dec(ml);
    end;
  end;
  writeln;
  writeln('Max bits size = ',8*osize);
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure rsa_test3;
var
  n,d,e,p,q,m,t: mp_int;
  prk: TPrivateKey;
{Test data from Crypto++/Botan}
const
  {MD2}
  d1p: pchar8 = '33d48445c859e52340de704bcdda065fbb4058d740bd1d67d29e9c146c11cf61';
  d1q: pchar8 = '335e8408866b0fd38dc7002d3f972c67389a65d5d8306566d5c4f2a5aa52628b';
  d1h: array[0..15] of byte = ($1d,$32,$de,$00,$9f,$9c,$56,$ea,$46,$36,$d3,$9a,$af,$fd,$ae,$a1);
  d1s: array[0..63] of byte = ($05,$fa,$6a,$81,$2f,$c7,$df,$8b,$f4,$f2,$54,$25,$09,$e0,$3e,$84,
                               $6e,$11,$b9,$c6,$20,$be,$20,$09,$ef,$b4,$40,$ef,$bc,$c6,$69,$21,
                               $69,$94,$ac,$04,$f3,$41,$b5,$7d,$05,$20,$2d,$42,$8f,$b2,$a2,$7b,
                               $5c,$77,$df,$d9,$b1,$5b,$fc,$3d,$55,$93,$53,$50,$34,$10,$c1,$e1);
  {MD5}
  d2p: pchar8 = 'D53E981F3D9AE9628B5038C8E48CBB944534522E8293145A3A98B855C4BB091E'+
                '67493348454DFDE02FF3C7148E313B917A199415937407B4856A98E50A570BDF';
  d2q: pchar8 = 'CE8671329A80756167093EEFDB10D2E0E0906BDBC58C4A1A8E8FF1CD2AD25086'+
                '8C79F360A357B7EDC1A7220CF698D0565385ECCC9FFBB89EE76EFAA6B70E8881';
  d2h: array[0..15] of byte = ($65,$37,$c3,$96,$19,$6a,$1d,$f1,$fd,$9f,$f0,$ec,$c7,$3e,$3c,$2d);
  d2s: array[0..127] of byte= ($31,$9f,$7d,$a1,$44,$8e,$5b,$a7,$ed,$ac,$7a,$5f,$b4,$22,$a4,$01,
                               $48,$1e,$89,$5e,$05,$08,$d1,$c0,$fd,$ed,$a2,$ac,$51,$de,$1d,$39,
                               $91,$3f,$0d,$41,$2e,$6e,$6d,$93,$13,$14,$19,$92,$a2,$02,$fb,$ef,
                               $f3,$bd,$33,$35,$42,$c8,$8f,$62,$64,$57,$04,$61,$90,$ae,$b1,$6f,
                               $f2,$a4,$99,$df,$58,$20,$24,$0a,$52,$48,$07,$44,$45,$b2,$d5,$4d,
                               $df,$0c,$29,$8f,$57,$b6,$1d,$89,$ee,$ab,$e7,$ab,$c7,$28,$d4,$bd,
                               $e8,$28,$34,$ba,$59,$4c,$22,$31,$f2,$75,$7a,$7f,$cd,$04,$70,$39,
                               $d3,$a3,$fe,$22,$05,$71,$f3,$0c,$41,$b5,$c2,$5f,$dd,$e4,$fe,$87);
  {RMD160}
  d3p: pchar8 = 'BB305054066BEE5B66E9C651583F6B5F4005D3CE970520CDF277EF463EF1EF1B10E9428C6BCA4254F42E40F0040C7577';
  d3q: pchar8 = 'BE54384C3A200BBB597B8A59F8165D553B24B2A9F96236521580AB79924EBC6095E0905E58A7C1BB302DFDA4FE489395';
  d3h: array[0..19] of byte = ($e3,$f6,$39,$42,$69,$26,$c4,$a7,$da,$0e,$39,$a0,$0b,$a0,$05,$7c,
                               $68,$ec,$f9,$5a);
  d3s: array[0..95] of byte = ($05,$c9,$08,$2d,$af,$39,$a7,$a7,$06,$f6,$ae,$5c,$c1,$fe,$1a,$26,
                               $f8,$20,$30,$12,$57,$dc,$dd,$b7,$b4,$dd,$91,$db,$01,$b4,$72,$e4,
                               $c8,$55,$e1,$ae,$61,$08,$c3,$29,$67,$65,$9c,$69,$16,$91,$97,$94,
                               $9b,$a5,$8a,$19,$4f,$fd,$e8,$44,$99,$11,$fd,$5a,$8e,$77,$51,$a4,
                               $e1,$5c,$ce,$d1,$18,$df,$2f,$f6,$77,$d8,$98,$27,$99,$a7,$6b,$59,
                               $40,$9b,$c7,$6a,$31,$e6,$96,$b2,$6f,$ea,$32,$4e,$c2,$cb,$65,$58);
  {SHA1}
  d4p: pchar8 = 'D7103CD676E39824E2BE50B8E6533FE7CB7484348E283802AD2B8D00C80D19DF';
  d4q: pchar8 = 'C89996DC169CEB3F227958275968804D4BE9FC4012C3219662F1A438C9950BB3';
  d4h: array[0..19] of byte = ($a9,$4a,$8f,$e5,$cc,$b1,$9b,$a6,$1c,$4c,$08,$73,$d3,$91,$e9,$87,
                               $98,$2f,$bb,$d3);
  d4s: array[0..63] of byte = ($a7,$e0,$0c,$e4,$39,$1f,$91,$4d,$82,$15,$8d,$9b,$73,$27,$59,$80,
                               $8e,$25,$a1,$c6,$38,$3f,$e8,$7a,$51,$99,$15,$76,$50,$d4,$29,$6c,
                               $f6,$12,$e9,$ff,$80,$9e,$68,$6a,$0a,$f3,$28,$23,$83,$06,$e7,$99,
                               $65,$f6,$d0,$13,$81,$38,$82,$9d,$9a,$1a,$22,$76,$43,$06,$f6,$ce);
  {SHA224}
  d5p: pchar8 = 'F75306FB8700184C998959EB1D271B3DBAB883726D270E21B91CA78C1BD148D1'+
                'D428533853A88E4AB4E4CB033CF27E2A8A6C1A12880206D6B4D74A9373F35B5F';
  d5q: pchar8 = 'E615C504FD127394B4522FF5BC0143C93A07636D184B9B40DF6461E3BAACC23A'+
                'A1937C52B06E95849692465060D092713D36B9E49FD76CA86D2389CF1B0F1A73';
  d5h: array[0..27] of byte = ($64,$e1,$c8,$e2,$3d,$63,$9a,$7a,$f3,$e4,$74,$32,$91,$06,$f9,$47,
                               $8f,$bd,$e7,$5f,$72,$38,$7c,$60,$91,$ee,$ab,$0c);
  d5s: array[0..127] of byte= ($04,$c4,$01,$28,$38,$13,$19,$1b,$75,$0c,$8f,$6c,$ee,$6c,$41,$c1,
                               $c9,$ad,$b7,$cb,$e2,$66,$85,$e5,$68,$d0,$46,$d7,$f0,$73,$77,$5f,
                               $19,$f9,$ac,$a0,$a4,$3a,$b9,$a6,$c2,$00,$c0,$71,$74,$c6,$ac,$6a,
                               $67,$e8,$9b,$63,$67,$db,$a7,$13,$5d,$c5,$c2,$23,$62,$18,$be,$42,
                               $b7,$de,$c0,$a4,$0c,$1e,$16,$5b,$f2,$4f,$9b,$33,$e6,$9a,$b1,$96,
                               $b8,$fa,$ad,$bd,$ca,$a6,$57,$51,$2b,$16,$41,$af,$7d,$0d,$bb,$b9,
                               $be,$1b,$fd,$88,$95,$c6,$b7,$9e,$a3,$0c,$10,$4a,$5d,$f1,$38,$01,
                               $82,$d3,$d2,$c9,$b5,$53,$d5,$9f,$56,$4d,$b8,$2d,$79,$40,$e1,$8a);
  {SHA256}
  d6p: pchar8 = 'B6661EB38B9C75A91E65FE5E6743C6949D6415D57CC4BCE84ED81E39EB8B4A77A867A15A34155B7E1CCCF972D6A06D73';
  d6q: pchar8 = 'B85D9FA07D4F7CE1CF05150256FCE15BF7E647CDA350DDC574E72104DBDAF1F1F87930D1FAD8C8E3094746ABE9E8C46B';
  d6h: array[0..31] of byte = ($a9,$75,$b4,$bb,$89,$42,$94,$5e,$43,$a2,$23,$9a,$ba,$9b,$fb,$a0,
                               $0c,$62,$29,$e7,$0c,$21,$72,$77,$a3,$93,$1d,$9c,$10,$4b,$b8,$58);
  d6s: array[0..95] of byte = ($27,$5a,$b9,$84,$6d,$d6,$99,$9f,$df,$4f,$c4,$6f,$0e,$45,$4b,$3c,
                               $ff,$8c,$12,$50,$96,$c0,$88,$0d,$a9,$57,$d6,$57,$66,$a6,$29,$8e,
                               $db,$41,$ae,$0f,$33,$db,$9f,$b1,$5f,$93,$69,$e9,$8b,$27,$38,$b3,
                               $7f,$f7,$16,$51,$75,$e5,$29,$b2,$34,$99,$61,$c4,$8b,$94,$59,$43,
                               $2f,$d5,$87,$95,$05,$2a,$dc,$25,$06,$d0,$b5,$f8,$37,$87,$9e,$77,
                               $1c,$66,$f5,$98,$80,$f1,$37,$1d,$78,$29,$0b,$03,$f8,$64,$54,$f9);
  {SHA384}
  d7p: pchar8 = 'D6392D7E84F9A7233FAEE9B2F386C7921BC9974393EB3581EAFA66D8E7DDFB75BC2F0AF8A6746563FF9C80A420B56BD1';
  d7q: pchar8 = 'D6D02FCFD0CEBEC15844AB927C37BAE12E8EFD552FF706B3B632A30E5FC14C39BC0ACF423C73C412BC8AF0A473276283';
  d7h: array[0..47] of byte = ($9f,$5c,$8b,$a9,$c1,$d6,$39,$05,$63,$8e,$0d,$25,$4c,$3e,$d4,$dd,
                               $ce,$0e,$bb,$c2,$8a,$39,$1d,$ba,$42,$9d,$b9,$f6,$bc,$b2,$eb,$8f,
                               $87,$23,$85,$c6,$6d,$78,$df,$f1,$cd,$9f,$9c,$cb,$a9,$4e,$65,$8f);
  d7s: array[0..95] of byte = ($8c,$f2,$44,$97,$f9,$c2,$dd,$b9,$18,$38,$22,$44,$4f,$20,$74,$ca,
                               $22,$c5,$d5,$1b,$52,$e3,$f5,$13,$91,$9e,$a6,$9b,$da,$00,$82,$b9,
                               $c4,$06,$80,$eb,$81,$3b,$e9,$14,$e9,$5e,$b9,$da,$cd,$3a,$3a,$23,
                               $8a,$aa,$67,$5f,$b7,$de,$fd,$55,$d2,$f9,$81,$d5,$b3,$d6,$b9,$46,
                               $7d,$b8,$06,$00,$81,$ae,$6f,$56,$ec,$7c,$0f,$c6,$e6,$e0,$78,$d6,
                               $14,$49,$7e,$cf,$89,$6d,$e4,$da,$a4,$f6,$79,$f6,$83,$fe,$0a,$2d);
  {SHA512}
  d8p: pchar8 = 'FF1414834AF27B94F7818DC435D95D756FCC9F4F4DD9B25FFDC23BDACB511DAA9F76E39B411DAD20D543AF198F8CC867';
  d8q: pchar8 = 'F025AE2AC9135A5E23B1D90E17DEEAA872D45E6700CF011015E949429D320142EC0A802B7863C8BDF99313FDAF2BF5CF';
  d8h: array[0..63] of byte = ($87,$45,$c6,$e4,$25,$bc,$5b,$85,$60,$31,$1e,$19,$28,$84,$fb,$76,
                               $8b,$ca,$76,$54,$da,$05,$85,$2f,$b1,$6e,$4e,$48,$ad,$c4,$44,$fe,
                               $af,$04,$93,$fc,$43,$61,$25,$b2,$6e,$84,$4e,$50,$9e,$17,$64,$d9,
                               $6d,$a1,$0e,$bf,$87,$cc,$e9,$57,$c9,$30,$45,$60,$97,$ce,$37,$86);
  d8s: array[0..95] of byte = ($0c,$98,$a0,$91,$c5,$b5,$3a,$71,$f0,$ba,$4f,$d4,$ce,$5f,$58,$b0,
                               $77,$30,$70,$c0,$5d,$b9,$16,$b4,$fc,$fe,$52,$98,$26,$7f,$28,$c4,
                               $57,$32,$2a,$3a,$b2,$62,$e4,$40,$79,$ce,$44,$0d,$b5,$69,$97,$af,
                               $cf,$30,$af,$33,$e1,$3d,$e7,$74,$d0,$66,$f0,$7b,$92,$cb,$b8,$b6,
                               $32,$1f,$2d,$d9,$67,$a9,$45,$24,$df,$05,$0a,$1c,$67,$aa,$a0,$5d,
                               $28,$c2,$32,$47,$2f,$ab,$21,$4c,$e8,$19,$79,$2d,$d1,$14,$77,$50);

var
  sm: array[0..255] of byte;
  lt,osize: word;
  f,cnt: integer;
  veri,sign1,sign2: boolean;

  function  Compare(psrc, pdest: pointer; size: integer): boolean;
    {-Compare memory block}
  var
    i: longint;
  begin
    if size>0 then begin
      Compare := false;
      if (psrc=nil) or (pdest=nil) then exit;
      for i:=1 to size do begin
        if pByte(psrc)^<>pByte(pdest)^ then exit;
        inc(Ptr2Inc(psrc));
        inc(Ptr2Inc(pdest));
      end;
    end;
    Compare := true;
  end;

  procedure test(hn: mp_string);
  begin
    inc(cnt,3);
    if not (veri and sign1 and sign2) then begin
      if not veri  then inc(f);
      if not sign1 then inc(f);
      if not sign2 then inc(f);
      writeln('Failure for ', hn:6,'  verify: ', veri:5, '    sign1: ', sign1:5, '    sign2: ', sign2:5);
    end;
  end;

begin
  cnt := 0;
  f := 0;
  mp_init7(n,d,e,p,q,m,t);
  mp_rsa_init_private(prk);
  writeln('RSA test 3 (sign/verify)');

  mp_read_radix(p,d1p,16);
  mp_read_radix(q,d1q,16);
  mp_set_int(e,$10001);
  mp_rsa_calc_nd(e,p,q,n,d);
  mp_rsa_calc_private(e,p,q,prk);
  osize := mp_unsigned_bin_size(n);
  veri  := mp_pkcs1v15_verify(e, n, SA_MD2, @d1h, @d1s, sizeof(d1h), osize);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign(d,n,SA_MD2, @d1h, @sm, sizeof(d1h), sizeof(sm), lt);
  sign1 := (lt=osize) and compare(@sm,@d1s,lt);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign2(prk,n,SA_MD2, @d1h, @sm, sizeof(d1h), sizeof(sm), lt);
  sign2 := (lt=osize) and compare(@sm,@d1s,lt);
  test('MD2');

  mp_read_radix(p,d2p,16);
  mp_read_radix(q,d2q,16);
  mp_set_int(e,$10003);
  mp_rsa_calc_nd(e,p,q,n,d);
  mp_rsa_calc_private(e,p,q,prk);
  osize := mp_unsigned_bin_size(n);
  veri  := mp_pkcs1v15_verify(e, n, SA_MD5, @d2h, @d2s, sizeof(d2h), osize);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign(d,n,SA_MD5, @d2h, @sm, sizeof(d2h), sizeof(sm), lt);
  sign1 := (lt=osize) and compare(@sm,@d2s,lt);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign2(prk,n,SA_MD5, @d2h, @sm, sizeof(d2h), sizeof(sm), lt);
  sign2 := (lt=osize) and compare(@sm,@d2s,lt);
  test('MD5');

  mp_read_radix(p,d3p,16);
  mp_read_radix(q,d3q,16);
  mp_set_int(e,$10001);
  mp_rsa_calc_nd(e,p,q,n,d);
  mp_rsa_calc_private(e,p,q,prk);
  osize := mp_unsigned_bin_size(n);
  veri  := mp_pkcs1v15_verify(e, n, SA_RMD160, @d3h, @d3s, sizeof(d3h), osize);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign(d,n, SA_RMD160, @d3h, @sm, sizeof(d3h), sizeof(sm), lt);
  sign1 := (lt=osize) and compare(@sm,@d3s,lt);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign2(prk,n,SA_RMD160, @d3h, @sm, sizeof(d3h), sizeof(sm), lt);
  sign2 := (lt=osize) and compare(@sm,@d3s,lt);
  test('RMD160');

  mp_read_radix(p,d4p,16);
  mp_read_radix(q,d4q,16);
  mp_set_int(e,$10001);
  mp_rsa_calc_nd(e,p,q,n,d);
  mp_rsa_calc_private(e,p,q,prk);
  osize := mp_unsigned_bin_size(n);
  veri  := mp_pkcs1v15_verify(e, n, SA_SHA1, @d4h, @d4s, sizeof(d4h), osize);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign(d,n,SA_SHA1, @d4h, @sm, sizeof(d4h), sizeof(sm), lt);
  sign1 := (lt=osize) and compare(@sm,@d4s,lt);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign2(prk,n,SA_SHA1, @d4h, @sm, sizeof(d4h), sizeof(sm), lt);
  sign2 := (lt=osize) and compare(@sm,@d4s,lt);
  test('SHA1');

  mp_read_radix(p,d5p,16);
  mp_read_radix(q,d5q,16);
  mp_set_int(e,$10001);
  mp_rsa_calc_nd(e,p,q,n,d);
  mp_rsa_calc_private(e,p,q,prk);
  osize := mp_unsigned_bin_size(n);
  veri  := mp_pkcs1v15_verify(e, n, SA_SHA224, @d5h, @d5s, sizeof(d5h), osize);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign(d,n,SA_SHA224, @d5h, @sm, sizeof(d5h), sizeof(sm), lt);
  sign1 := (lt=osize) and compare(@sm,@d5s,lt);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign2(prk,n,SA_SHA224, @d5h, @sm, sizeof(d5h), sizeof(sm), lt);
  sign2 := (lt=osize) and compare(@sm,@d5s,lt);
  test('SHA224');

  mp_read_radix(p,d6p,16);
  mp_read_radix(q,d6q,16);
  mp_set_int(e,$10001);
  mp_rsa_calc_nd(e,p,q,n,d);
  mp_rsa_calc_private(e,p,q,prk);
  osize := mp_unsigned_bin_size(n);
  veri  := mp_pkcs1v15_verify(e, n, SA_SHA256, @d6h, @d6s, sizeof(d6h), osize);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign(d,n,SA_SHA256, @d6h, @sm, sizeof(d6h), sizeof(sm), lt);
  sign1 := (lt=osize) and compare(@sm,@d6s,lt);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign2(prk,n,SA_SHA256, @d6h, @sm, sizeof(d6h), sizeof(sm), lt);
  sign2 := (lt=osize) and compare(@sm,@d6s,lt);
  test('SHA256');

  mp_read_radix(p,d7p,16);
  mp_read_radix(q,d7q,16);
  mp_set_int(e,$10001);
  mp_rsa_calc_nd(e,p,q,n,d);
  mp_rsa_calc_private(e,p,q,prk);
  osize := mp_unsigned_bin_size(n);
  veri  := mp_pkcs1v15_verify(e, n, SA_SHA384, @d7h, @d7s, sizeof(d7h), osize);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign(d,n,SA_SHA384, @d7h, @sm, sizeof(d7h), sizeof(sm), lt);
  sign1 := (lt=osize) and compare(@sm,@d7s,lt);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign2(prk,n,SA_SHA384, @d7h, @sm, sizeof(d7h), sizeof(sm), lt);
  sign2 := (lt=osize) and compare(@sm,@d7s,lt);
  test('SHA384');

  mp_read_radix(p,d8p,16);
  mp_read_radix(q,d8q,16);
  mp_set_int(e,$10001);
  mp_rsa_calc_nd(e,p,q,n,d);
  mp_rsa_calc_private(e,p,q,prk);
  osize := mp_unsigned_bin_size(n);
  veri := mp_pkcs1v15_verify(e, n, SA_SHA512, @d8h, @d8s, sizeof(d8h), osize);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign(d,n,SA_SHA512, @d8h, @sm, sizeof(d8h), sizeof(sm), lt);
  sign1 := (lt=osize) and compare(@sm,@d8s,lt);
  fillchar(sm, sizeof(sm), 0);
  mp_pkcs1v15_sign2(prk,n,SA_SHA512, @d8h, @sm, sizeof(d8h), sizeof(sm), lt);
  sign2 := (lt=osize) and compare(@sm,@d8s,lt);
  inc(cnt,3);
  if not (veri and sign1 and sign2) then begin
    if not veri  then inc(f);
    if not sign1 then inc(f);
    if not sign2 then inc(f);
    writeln('SHA512   Verify: ', veri:5, '    sign1: ', sign1:5, '    sign2: ', sign2:5);
  end;

  mp_rsa_clear_private(prk);
  mp_clear7(n,d,e,p,q,m,t);
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure sqrtmod_test;
  {-Check mp_sqrtmod, mp_sqrtmodpk, mp_sqrtmodp2}
var
  i,Err,f,k,cnt: integer;
{$ifndef BIT16}
const
  TCnt=50;
  BitSize=384;
  KCnt=8;
{$else}
const
  TCnt=50;
  BitSize=128;
  KCnt=5;
{$endif}
begin
  f := 0;
  cnt := 0;
  writeln('Test mp_sqrtmod and mp_sqrtmodpk: k=2..',KCnt);
  for i:=1 to TCnt do begin
    {s prime}
    mp_rand_prime(BitSize, pt_normal, s);
    mp_set(a,2);
    while mp_jacobi(a,s)<>1 do mp_inc(a);
    inc(cnt);
    mp_sqrtmod(a,s,c,Err);
    if Err<>0 then begin
      writeln('Err=',Err);
      inc(f);
    end
    else begin
      mp_sqrmod(c,s,d);
      {$ifdef MPC_UseKONG}
        k := s.pdigits^[0] and 15;
        if k=9 then write(k) else write(k and 7);
      {$else}
        write(s.pdigits^[0] and 7);
      {$endif}
      if mp_is_ne(a,d) then begin
        writeln;
        writeln('s^2: ', mp_decimal(t));
        writeln('s&7: ', s.pdigits^[0] and 7);
        writeln('  a: ', mp_decimal(a));
        writeln(' ûa: ', mp_decimal(b));
        inc(f);
      end;
      for k:=2 to KCnt do begin
        inc(cnt);
        {d = s^k, c = sqrt(a) mod s^k}
        if k=2 then begin
          mp_sqr(s,d);
          mp_sqrtmodp2(a,s,c,Err);
        end
        else begin
          mp_mul(s,d,d);
          mp_sqrtmodpk(a,s,k,c,Err);
        end;
        if Err=0 then begin
          mp_sqrmod(c,d,e);
          mp_submod(e,a,d,e);
          if not mp_iszero(e) then begin
            writeln;
            writeln('  k: ', k);
            writeln('  s: ', mp_decimal(s));
            writeln('s&7: ', s.pdigits^[0] and 7);
            writeln('  a: ', mp_decimal(a));
            writeln(' ûa: ', mp_decimal(c));
            inc(f);
          end;
        end
        else begin
          writeln('Err=',Err, ', k=',k);
          inc(f);
          break;
        end;
      end;
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure sqrtmod2k_test;
  {-Check mp_sqrtmod2k}
var
  i,Err,f,k,cnt,j: integer;
{$ifndef BIT16}
const
  MaxBit=1024;
{$else}
const
  MaxBit=256;
{$endif}
begin
  f := 0;
  cnt := 0;
  k := 6;
  writeln('Test mp_sqrtmod2k');
  while k<MaxBit do begin
    {s := 2^k}
    mp_2expt(s,k);
    for i:=1 to 10 do begin
      inc(cnt);
      {construct a which have sqrtmods}
      mp_rand_bits(a,k-3);
      case random(5) of
         0:   begin
                {odd(a) < 0, a mod 8 = 1}
                mp_chs(a,a);
                mp_shl(a,3,a);
                mp_sub_d(a,7,a);
              end;
         1:   begin
                {even a = r^(2*j) }
                j := random(k) div 2;
                mp_2expt(a,2*j);
              end;
        else  begin
                {odd(a) > 0, a mod 8  = 1}
                mp_shl(a,3,a);
                mp_inc(a);
              end;
      end;
      mp_sqrtmod2k(a,k,b,Err);
      if Err<>0 then begin
        writeln('Err=',Err);
        mp_writeln('a=', a);
        inc(f);
        break;
      end
      else begin
        {Check b^2 - a = 0 mod s}
        mp_sqrmod(b,s,d);
        mp_submod(d,a,s,t);
        if not mp_iszero(t) then begin
          writeln;
          writeln('  k: ', k);
          writeln('  s: ', mp_decimal(t));
          writeln('  a: ', mp_decimal(a));
          writeln(' ûa: ', mp_decimal(b));
          inc(f);
        end;
      end;
    end;
    inc(k);
    if k and 15 = 0 then write('.');
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure cbrtmod_test;
  {-Check mp_cbrtmod_ex and mp_cbrtmodpk}
var
  i,j,Err,f,k,cnt: integer;
  m3: mp_digit;
  isCR: boolean;
{$ifndef BIT16}
const
  PCnt=16;
  BFac=20;
  KCnt=8;
{$else}
const
  PCnt=16;
  BFac=10;
  KCnt=5;
{$endif}
begin
  f := 0;
  cnt := 0;
  writeln('Test mp_cbrtmod_ex and mp_cbrtmodpk: k=2..',KCnt);
  for i:=1 to PCnt do begin
    mp_rand_prime(PCnt*BFac+64, pt_normal, s);
    mp_set_int(a,-1);
    inc(cnt);
    if not s_mp_is_cubres(a, s) then begin
      inc(f);
      writeln('-1 not CR mod ',mp_decimal(s));
    end;
    mp_div_d(s,3,@t,m3);
    for j:=1 to 4 do begin
      mp_rand_bits(a,60);
      isCR := false;
      if m3=1 then begin
        {if s mod 3 = 1 test up to 2 values for 3rd power residue}
        for k:=1 to 2 do begin
          if s_mp_is_cubres(a, s) then begin
            isCR := true;
            break;
          end;
          mp_inc(a);
        end;
      end;
      inc(cnt);
      {b = cbrt(a) mod s}
      mp_cbrtmod_ex(a,s,b,@z,Err);
      if Err<>0 then begin
        if m3=1 then begin
          {No error if s mod 3 = 1 and non-residue}
          if isCR then begin
            writeln('Err=',Err);
            inc(f);
          end
          else write('-');
        end
        else begin
          writeln('Err=',Err);
          inc(f);
        end;
      end
      else begin
        write(m3);
        mp_sqrmod(b,s,c);
        mp_mulmod(c,b,s,c);
        Err := 0;
        if mp_is_ne(a,c) then begin
          writeln;
          writeln('        a: ', mp_decimal(a));
          writeln('        s: ', mp_decimal(s));
          writeln('b=a^(1/3): ', mp_decimal(b));
          writeln('      b^3: ', mp_decimal(c));
          inc(Err);
        end;
        if not mp_is1(z) then begin
          {check other two roots}
          for k:=1 to 2 do begin
            mp_mulmod(b,z,s,b);
            mp_sqrmod(b,s,c);
            mp_mulmod(c,b,s,c);
            if mp_is_ne(a,c) then begin
              writeln;
              writeln('            a: ', mp_decimal(a));
              writeln('            s: ', mp_decimal(s));
              writeln('            z: ', mp_decimal(z));
              writeln('b=a^(1/3)*z^',k,': ', mp_decimal(b));
              writeln('          b^3: ', mp_decimal(c));
              inc(Err);
            end;
          end;
        end;
        if Err<>0 then inc(f)
        else begin
          mp_copy(s,d);
          for k:=2 to KCnt do begin
            inc(cnt);
            {d = s^k, c = cqrt(a) mod s^k}
            mp_mul(s,d,d);
            mp_cbrtmodpk(a,s,k,c,Err);
            if Err=0 then begin
              mp_sqrmod(c,d,e);
              mp_mulmod(e,c,d,e);
              mp_submod(e,a,d,e);
              if not mp_iszero(e) then begin
                writeln;
                writeln('      k: ', k);
                writeln('      s: ', mp_decimal(s));
                writeln('s mod 3: ', m3);
                writeln('      a: ', mp_decimal(a));
                writeln('a^(1/3): ', mp_decimal(c));
                inc(f);
              end;
            end
            else begin
              writeln('Err=',Err, ', k=',k);
              inc(f);
              break;
            end;
          end;
        end;
      end;
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure testradixconv;
  {-Test input/output radix conversion}
var
  cnt,f, radix,digs: word;
{$ifndef BIT16}
 const dmax=2000;
{$else}
 const dmax=500;
{$endif}
begin
  writeln('Test radix conversion (MaxRadix=',MaxRadix,'):  ');
  f := 0;
  cnt := 0;
  for radix:=2 to MaxRadix do begin
    digs := 2;
    while digs<=dmax do begin
      mp_rand_radix(a,radix,digs);
      if random(2)=0 then mp_chs(a,a);
      mp_show_plus := random(2)=0;
      mp_uppercase := random(2)=0;
      mp_toradix_n(a, line, radix, LMAX);
      mp_read_radix(b,line,radix);
      inc(cnt);
      if mp_is_ne(a,b) then begin
        writeln('Diff Radix=',radix, '  digits=',digs);
        write('a='); mp_output_radix(a,radix); writeln;
        write('b='); mp_output_radix(b,radix); writeln;
        inc(f);
      end;
      digs := 2 + 14*digs div 10;
    end;
    write('.');
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure pell4_test;
var
  i,f: integer;
  r,q: integer;
begin
  writeln('mp_pell4 test: ', fn_pell4);
  if open_file(fn_pell4) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readln(tf, q);
    mp_pell4(a,s,t,r);
    if not (mp_is_eq(b,s) and mp_is_eq(c,t) and (q=r) ) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' r: ', r);
      writeln(' s: ', decimal(s));
      writeln(' t: ', decimal(t));
      writeln(' q: ', q);
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure pell1_test;
var
  i,f: integer;
begin
  writeln('mp_pell1 test: ', fn_pell1);
  if open_file(fn_pell1) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    mp_pell(a,s,t);
    if not (mp_is_eq(b,s) and mp_is_eq(c,t) ) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', decimal(a));
      writeln(' b: ', decimal(b));
      writeln(' c: ', decimal(c));
      writeln(' s: ', decimal(s));
      writeln(' t: ', decimal(t));
    end;
    inc(i)
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_rqffu;
var
  f,h,i,m,n,norm: integer;
  d,j: longint;
  reg, fx, fy, px, py, e, rc, rd: double;
const
  eps = 1e-15;
  sx  = '35040408137079754042430530'; {direct values for d=54321,computed with}
  sy  = '150343715133198616100976';   {Aribas Pell4, verified with Pari. Note factor 2}
begin
  {Cohen [24], Table B.2, recomputed with Pari script:}
  { for (d=2, 500,
         if (isfundamental(d)==1,
             e = quadunit(d);
             printp(d, " \t ", qfbclassno(d), " \t  ", quadregulator(d)," \t ",norm(e), " \t ", e))
         );
  }
  writeln('Test mp_rqffu and s_mp_rqffu: ', fn_rqffu);

  if open_file(fn_rqffu) then exit;
  readln(tf); {comment line}

  i := 0;
  f := 0;
  while not eof(tf) do begin
    inc(i);
    readln(tf,d,h,reg,norm,fx,fy);
    if h<0 then write('?');
    n := mp_rqffu(d,s,t);
    if n=0 then begin
      writeln('Error return from mp_rqffu for d=',d);
      inc(f);
      continue;
    end;
    px := mp_todouble(s);
    py := mp_todouble(t);
    {Fundamental unit e = (px + py*sqrt(disc))/2}
    e  := 0.5*(px+py*sqrt(d));
    {Computed regulator rc = ln(e)}
    rc := ln(e);
    rd := abs(1-rc/reg);
    {Compute representation e = px + py*w as in Cohen's table}
    if odd(d) then begin
      px := 0.5*(px - py);
    end
    else begin
      px := 0.5*px;
    end;
    {Compute b = norm of e, should be = n = norm}
    mp_sqr(s,a);
    mp_sqr(t,b);
    mp_mul_int(b,d,b);
    mp_sub(a,b,b);
    mp_shr(b,2,b);
    m := 0;
    if abs(fx-px)>eps then begin
      writeln('fx<>px for d=',d);
      inc(m);
    end;
    if abs(fy-py)>eps then begin
      writeln('fy<>py for d=',d);
      inc(m);
    end;
    if (not mp_is_longint(b,j)) or (j<>norm) or (n<>norm) then begin
      writeln('Wrong norm for d=',d);
      inc(m);
    end;
    if rd>1e-8 then begin
      writeln('Wrong regulator for d=',d);
      inc(m);
    end;
    if m=0 then begin
      mp_set_int(c,d);
      n  := s_mp_rqffu(c,a,b);
      if n=0 then begin
        writeln('Error return from s_mp_rqffu for d=',d);
        inc(m);
      end
      else begin
        if mp_is_ne(a,s) or mp_is_ne(b,t) then begin
          writeln('Different results for d=',d);
          inc(m);
        end;
      end;
    end;
    if m<>0 then inc(f);
  end;
  close(tf);
  {Bigger values}
  inc(i);
  m := 0;
  d := 54321;
  mp_set_int(c,d);
  mp_read_decimal(a,sx);
  mp_read_decimal(b,sy);
  n := mp_rqffu(d,s,t);
  if (n=0) or mp_is_ne(a,s) or mp_is_ne(b,t) then begin
    writeln('mp_rqffu failed for d=',d);
    inc(m);
  end;
  n  := s_mp_rqffu(c,s,t);
  if (n=0) or mp_is_ne(a,s) or mp_is_ne(b,t) then begin
    writeln('s_mp_rqffu failed for d=',d);
    inc(m);
  end;
  if m<>0 then inc(f);
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure sqrtmodpq_test;
  {-Check mp_sqrtmodpq}
var
  i,Err,f,k,cnt: integer;
{$ifndef BIT16}
const
  MaxBit=384;
{$else}
const
  MaxBit=128;
{$endif}
begin
  f := 0;
  k := 8;
  writeln('Test mp_sqrtmodpq');
  cnt := 0;
  {a,b,c,d,e,s,t,z: mp_int;}

  while k<=MaxBit do begin
    {s,t primes, z=s*t}
    {calc roots b,c of a mod z}
    mp_rand_prime(k, pt_normal, s);
    repeat
      mp_rand_prime(k, pt_normal, t);
    until mp_is_ne(s,t);
    mp_mul(s,t,z);

    for i:=1 to 5 do begin
      repeat
        mp_rand_bits(a,k-1);
        mp_mod(a,z,a);
      until (mp_jacobi(a,s)=1) and (mp_jacobi(a,t)=1);
      mp_sqrtmodpq(a,s,t,b,c,Err);
      inc(cnt);
      if Err<>0 then begin
        writeln('Err=',Err);
        mp_writeln('a=', a);
        inc(f);
        break;
      end
      else begin
        {check if d=b^2=a mod z, e=c^2=a mod z}
        mp_sqrmod(b,z,d);
        mp_sqrmod(c,z,e);
        if mp_is_ne(a,d) or mp_is_ne(a,e) then begin
          mp_writeln(' p=',s);
          mp_writeln(' q=',t);
          mp_writeln(' z=',z);
          mp_writeln(' a=',a);
          mp_writeln(' d=',d);
          mp_writeln(' e=',e);
          inc(f);
        end;
      end;
    end;
    mp_sqrtmodpq(a,s,s,b,c,Err);
    inc(cnt);
    {check case p=q}
    if Err<>0 then begin
      writeln('Err=',Err);
      mp_writeln('a=', a);
      inc(f);
      break;
    end
    else begin
      mp_sqr(s,z);
      {check if d=b^2=a mod z}
      mp_sqrmod(b,z,d);
      if mp_is_ne(a,d) then begin
        mp_writeln(' p=',s);
        mp_writeln(' z=',z);
        mp_writeln(' a=',a);
        mp_writeln(' d=',d);
        inc(f);
      end;
    end;
    inc(k,8);
    if k and 15 = 0 then write('.');
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure cbrtmodpq_test;
  {-Check mp_sbrtmodpq}
var
  i,Err,f,k,cnt: integer;
{$ifndef BIT16}
const
  MaxBit=384;
{$else}
const
  MaxBit=128;
{$endif}
begin
  f := 0;
  k := 8;
  writeln('Test mp_cbrtmodpq');
  cnt := 0;
  {a,b,c,d,e,s,t,z: mp_int;}

  while k<=MaxBit do begin
    {s,t primes, z=s*t}
    {calc roots b,c of a mod z}
    mp_rand_prime(k, pt_normal, s);
    repeat
      mp_rand_prime(k, pt_normal, t);
    until mp_is_ne(s,t);
    mp_mul(s,t,z);

    for i:=1 to 5 do begin
      repeat
        mp_rand_bits(a,k-1);
        mp_mod(a,z,a);
      until s_mp_is_cubres(a,s) and s_mp_is_cubres(a,t);
      mp_cbrtmodpq(a,s,t,b,Err);
      inc(cnt);
      if Err<>0 then begin
        writeln('Err=',Err);
        mp_writeln('a=', a);
        inc(f);
        break;
      end
      else begin
        {check if d=b^3=a mod z}
        mp_sqrmod(b,z,d);
        mp_mulmod(d,b,z,d);
        if mp_is_ne(a,d) then begin
          mp_writeln(' p=',s);
          mp_writeln(' q=',t);
          mp_writeln(' z=',z);
          mp_writeln(' a=',a);
          mp_writeln(' d=',d);
          inc(f);
        end;
      end;
    end;
    {check case p=q}
    mp_cbrtmodpq(a,s,s,b,Err);
    inc(cnt);
    if Err<>0 then begin
      writeln('Err=',Err);
      mp_writeln('a=', a);
      inc(f);
      break;
    end
    else begin
      mp_sqr(s,z);
      {check if d=b^3=a mod z}
      mp_sqrmod(b,z,d);
      mp_mulmod(d,b,z,d);
      if mp_is_ne(a,d) then begin
        mp_writeln(' p=',s);
        mp_writeln(' z=',z);
        mp_writeln(' a=',a);
        mp_writeln(' d=',d);
        inc(f);
      end;
    end;

    inc(k,8);
    if k and 15 = 0 then write('.');
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_checksum;
const
  vec: array[1..3] of byte = ($11,$22,$33);
{$ifdef MP_32BIT}
  TV: array[-1..1] of longint = ($1E0004,$60001,$180003);  {4 byte mp_digits}
{$else}
  TV: array[-1..1] of longint = ($160004,$60001,$120003);  {2 byte mp_digits}
{$endif}
var
  i,f: integer;
  A1,AF: longint;
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
  for i:=-1 to 1 do begin
    mp_set_int(a,i);
    AF := mp_checksum(a);
    if AF<>TV[i] then inc(f);
  end;
  if f<>0 then begin
    writeln('No. of mp_checksum test failed: ',f);
    inc(totalfailed,f);
  end;
end;


{---------------------------------------------------------------------------}
procedure bool_test;
var
  i,j,f: integer;
begin
  {test cases calculated with NX}
  writeln('mp_and/or/xor test: ', fn_bool);
  if open_file(fn_bool) then exit;
  i:=0;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readline;
    mp_read_radix(d,line,10);
    readline;
    mp_read_radix(e,line,10);
    {and}
    for j:=0 to 2 do begin
      inc(i);
      if j=0 then begin
        mp_copy(a,z);
        mp_and(z,b,z);
      end
      else if j=1 then begin
        mp_copy(b,z);
        mp_and(a,z,z);
      end
      else mp_and(a,b,z);
      if mp_is_ne(c,z) then begin
        writeln('mp_and failed! #',i);
        inc(f);
        writeln(' a: ', decimal(a));
        writeln(' b: ', decimal(b));
        writeln(' c: ', decimal(c));
        writeln(' z: ', decimal(z));
      end;
    end;
    {or}
    for j:=0 to 2 do begin
      inc(i);
      if j=0 then begin
        mp_copy(a,z);
        mp_or(z,b,z);
      end
      else if j=1 then begin
        mp_copy(b,z);
        mp_or(a,z,z);
      end
      else mp_or(a,b,z);
      if mp_is_ne(d,z) then begin
        writeln('mp_or failed! #',i);
        inc(f);
        writeln(' a: ', decimal(a));
        writeln(' b: ', decimal(b));
        writeln(' d: ', decimal(c));
        writeln(' z: ', decimal(z));
      end;
    end;
    {xor}
    for j:=0 to 2 do begin
      inc(i);
      if j=0 then begin
        mp_copy(a,z);
        mp_xor(z,b,z);
      end
      else if j=1 then begin
        mp_copy(b,z);
        mp_xor(a,z,z);
      end
      else mp_xor(a,b,z);
      if mp_is_ne(e,z) then begin
        writeln('mp_xor failed! #',i);
        inc(f);
        writeln(' a: ', decimal(a));
        writeln(' b: ', decimal(b));
        writeln(' d: ', decimal(c));
        writeln(' z: ', decimal(z));
      end;
    end;
    write('.');
  end;
  close(tf);
  writeln;
  writeln('No. of checks: ',i,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


(*
D:\MATH\Calc_KM>calc_win32.exe

                        CALC

                A NUMBER THEORY CALCULATOR
                K.R.MATTHEWS, 27th July 2011

> a=2^81+17
> b=997*a
2410598084311570574364131893

> cornacchia(2,5,b)
X=34707886933112, Y=514442142169
X=6620400218468, Y=21554297419337

> cornacchia(1,1,b)
X=46419794101458, Y=15993773785127
X=37099925472198, Y=32158880799383

> a=2^81+17
> b=60013*a
145102530425065481323284500797

> cornacchia(1,1,b)
X=367413518116066, Y=100547685854229
X=271108208947046, Y=267587125001541

> cornacchia(1,3,b)
X=311619996205463, Y=126485188052426
X=115761101081447, Y=209524778090386

> cornacchia(1,7,b)
X=9394418820753, Y=143931667973822
X=363819928338945, Y=42657423015914

> cornacchia(1,19,b)
X=216916785997911, Y=71836670879302

> cornacchia(1,37,b)
X=380188300220883, Y=3888259529422
X=230629014965133, Y=49841041343422

> cornacchia(2,5,b)
X=148116157800386, Y=142285444105831
X=196958678979526, Y=116204206442233

> cornacchia(3,37,b)
X=192786408078040, Y=30136057690741

> cornacchia(2,29,b)
X=269004043732618, Y=3601624806059
X=217555667482594, Y=41705698575925

> cornacchia(5,17,b)
X=158476830245957, Y=33892565383816
X=31379649496219, Y=90806555351776

> cornacchia(5,13,b)
X=121259578960029, Y=74205070828828
X=73455449487384, Y=95322939013447

> cornacchia(7,37,b)
X=105113802923940, Y=42794338306109

b = 60013*4294967311
> cornacchia(3,7,b)
X=8580335, Y=2295568
X=9195759, Y=762320

> cornacchia(3,1,b)
X=8829593, Y=4885564
X=5334223, Y=13129816

> a
4294967311
> b=997*a

> cornacchia(1,3,b)
X=340252, Y=1178461
X=961472, Y=1057931

> cornacchia(1,23,b)
X=944078, Y=383961
X=757778, Y=401511

> cornacchia(3,19,b)
X=1065259, Y=214936
X=571667, Y=416860
*)


{---------------------------------------------------------------------------}
procedure cornacchia_test;
var
  i,Err,f,cnt: integer;
  OK: boolean;
  bs: word;
{$ifndef BIT16}
const
  BitSize=512;
{$else}
const
  BitSize=256;
{$endif}
type
 t4tr = record
          d: integer;
          x: pchar8;
          y: pchar8;
        end;
{Calculated with: CALC - a number theory calculator, K.R.Matthews}
{available from http://www.numbertheory.org/calc/krm_calc.html}
const
  tests: array[1..5] of t4tr = (
           (d: -43; X:'45656439196037'; Y:'3016889043405'),
           (d: -67; X:'21533115172379'; Y:'5480229090897'),
           (d:-211; X:'30337247743263'; Y:'2715177176395'),
           (d:-307; X:'47866492077971'; Y:' 775603190127'),
           (d:-499; X:'28743118396413'; Y:'1818251501825')
         );
type
 TCPKR = record
          a,b,k: word;
          p: longint;
          x: pchar8;
          y: pchar8;
        end;
const
  NCXT = 21;
const
  xtests: array[1..NCXT] of TCPKR = (
           (a: 1; b: 1; k: 1; p: 2;          X:'1';  Y:'1'),
           (a: 1; b: 1; k: 1; p: 5;          X:'2';  Y:'1'),
           (a: 1; b: 1; k: 2; p: 5;          X:'4';  Y:'3'),
           (a: 5; b: 2; k: 1; p: 997;        X:'11'; Y:'14'),
           (a: 2; b:11; k: 1; p: 997;        X:'19'; Y:'5'),
           (a: 1; b: 1; k:15; p: 997;        X:'29766257562016024252397'; Y:'8360785584538660197522'),
           (a: 1; b: 6; k: 1; p: MaxLongint; X:'44029'; Y:'5901'),
           (a: 5; b: 2; k: 1; p: MaxLongint; X:'15661'; Y:'21461'),
           (a: 2; b: 5; k: 1; p: MaxLongint; X:'21461'; Y:'15661'),
           (a: 3; b: 1; k: 5; p: MaxLongint; X:'52504498159159928532381';  Y:'193395343996221355924282'),
           (a: 3; b:19; k: 5; p: MaxLongint; X:'113111620127229388498231'; Y:'19586797697526210438386'),
           (a: 5; b: 2; k: 5; p: MaxLongint; X:'86916907386052869622555';  Y:'62845771983777943024271'),
           (a: 7; b: 1; k: 5; p: MaxLongint; X:'80721804012440188212637';  Y:'7736831334134313690432'),
           (a: 7; b:15; k: 5; p: MaxLongint; X:'68213405464036673894251';  Y:'29552715910437062147580'),
           (a:10; b:23; k: 5; p: MaxLongint; X:'66707836248743849013652';  Y:'7140125106565126904373'),
           (a: 1; b: 5; k: 6; p: 100000007;  X:'268842986016149464338502'; Y:'430749084472618226730723'),
           (a: 2; b:17; k: 6; p: 100000007;  X:'91251701252348796784364';  Y:'240507633979021384760911'),
           (a: 7; b: 1; k: 6; p: 100000007;  X:'599999860000002940000';    Y:'999998950000073499999657'),
           (a: 7; b: 9; k: 6; p: 100000007;  X:'599999860000002940000';    Y:'333332983333357833333219'),
           (a: 1; b: 7; k: 7; p: 100000007;  X:'9999985300001714999975990000';  Y:'6999997550000102899999657'),
           (a: 2; b: 5; k: 7; p: 100000007;  X:'5306110869911201931560791797';  Y:'2956025164724774997324054655')
         );

type
 TCPQR = record
          p,q: byte;
          a,b: word;
          x: pchar8;
          y: pchar8;
        end;
const
  NPQ=11;
const
  FPQ: array[0..3] of pchar8 = ('997','60013','4294967311','2417851639229258349412369');
const
  pqtests: array[1..NPQ] of TCPQR = (
             (p: 0; q: 3; a: 1; b:  1; X:'37099925472198';   Y:'32158880799383'),
             (p: 0; q: 3; a: 2; b:  5; X:'34707886933112';   Y:'514442142169'),
             (p: 1; q: 3; a: 3; b: 37; X:'192786408078040';  Y:'30136057690741'),
             (p: 1; q: 3; a: 7; b: 37; X:'105113802923940';  Y:'42794338306109'),
             (p: 1; q: 3; a: 2; b:  5; X:'148116157800386';  Y:'142285444105831'),
             (p: 1; q: 3; a: 1; b: 19; X:'216916785997911';  Y:'71836670879302'),
             (p: 1; q: 2; a: 3; b:  7; X:'9195759';  Y:'762320'),
             (p: 1; q: 2; a: 3; b:  1; X:'5334223';  Y:'13129816'),
             (p: 0; q: 2; a: 1; b:  3; X:'340252';   Y:'1178461'),
             (p: 0; q: 2; a: 1; b: 23; X:'757778';   Y:'401511'),
             (p: 0; q: 2; a: 3; b: 19; X:'571667';   Y:'416860'));
begin
  f := 0;
  cnt := 0;
  writeln('Test mp_cornacchia/4/pk/pq: ');
  for i:=1 to 5 do begin
    inc(cnt);
    mp_set_int(d,tests[i].d);
    mp_mersenne(89,e);
    mp_cornacchia4(d,e, a,b, Err);
    if Err=0 then begin
      mp_read_decimal(c,tests[i].x);
      OK := mp_is_eq(a,c);
      mp_read_decimal(c,tests[i].y);
      OK := OK and mp_is_eq(b,c);
      {test a^2 - d*b2^ = 4*e}
      mp_sqr(b,b);
      mp_sqr(a,a);
      mp_mul(b,d,b);
      mp_sub(a,b,b);
      mp_shl(e,2,c);
      OK := OK and mp_is_eq(b,c);
      if not OK then begin
        writeln('Failure mp_cornacchia4 test case ',i);
        inc(f);
      end;
    end
    else begin
      writeln('Error mp_cornacchia4 test case ',i);
      inc(f);
    end;
  end;

  bs := 32;
  while bs<BitSize do begin
    inc(cnt);
    write('.');
    {check  s^2 + t^2 = z has solution for prime z iff z=1 mod 4}
    mp_rand_prime(bs, pt_normal, z);
    mp_set(d,1);
    mp_cornacchia(d,z,s,t,Err);
    if Err=0 then begin
      mp_sqr(s,a);
      mp_sqr(t,b);
      mp_add(a,b,c);
      if mp_is_ne(z,c) then begin
        inc(f);
        writeln('z=',mp_decimal(z));
      end
      else if z.pdigits^[0] and 3 <> 1 then begin
        inc(f);
        writeln('Err=',Err, ',  z=',mp_decimal(z));
      end;
    end
    else begin
      if z.pdigits^[0] and 3 = 1 then begin
        inc(f);
        writeln('Err=',Err, ',  z=',mp_decimal(z));
      end;
    end;
    {cornacchia tests prime with small d}
    for i:=2 to 5 do begin
      mp_set(d,mp_digit(i));
      mp_chs(d,a);
      if mp_jacobi(a,z)=1 then begin
        mp_cornacchia(d,z,s,t,Err);
        if Err=0 then begin
          inc(cnt);
          mp_sqr(s,a);
          mp_sqr(t,b);
          mp_mul(b,d,b);
          mp_add(a,b,c);
          if mp_is_ne(z,c) then begin
            inc(f);
            writeln('Failure mp_cornacchia:');
            writeln('z=',mp_decimal(z));
            writeln('d=',mp_decimal(d));
          end;
        end;
      end;
    end;
    for i:=2 to 5 do begin
      mp_set(d,mp_digit(i));
      mp_chs(d,d);
      if mp_jacobi(d,z)=1 then begin
        mp_cornacchia4(d,z,s,t,Err);
        if Err=0 then begin
          inc(cnt);
          mp_sqr(s,a);
          mp_sqr(t,b);
          mp_mul(b,d,b);
          mp_sub(a,b,c);
          mp_shl(z,2,a);
          if mp_is_ne(a,c) then begin
            inc(f);
            writeln('Failure mp_cornacchia4:');
            writeln('z=',mp_decimal(z));
            writeln('d=',mp_decimal(d));
          end;
        end;
      end;
    end;
    bs := round(bs*1.05);
  end;

  {mp_cornacchia_pk}
  for i:=1 to NCXT do begin
    write('.');
    inc(cnt);
    mp_set(a, xtests[i].a);
    mp_set(b, xtests[i].b);
    mp_set_int(z, xtests[i].p);
    mp_read_decimal(c,xtests[i].x);
    mp_read_decimal(d,xtests[i].y);
    mp_cornacchia_pk(a,b,z,xtests[i].k,s,t,err);
    if (err<>0) or mp_is_ne(c,s) or mp_is_ne(d,t) then begin
      writeln('Failure mp_cornacchia_pk for i=',i);
      inc(f);
    end;
  end;

  {mp_cornacchia_pq test cases: Note that mp_cornacchia_pq returns at}
  {most one solution for the test cases, while CALC may return two!  }
  for i:=1 to NPQ do begin
    write('.');
    inc(cnt);
    mp_set(a, pqtests[i].a);
    mp_set(b, pqtests[i].b);
    mp_read_decimal(e,FPQ[pqtests[i].p]);
    mp_read_decimal(z,FPQ[pqtests[i].q]);
    mp_read_decimal(d,pqtests[i].y);
    mp_read_decimal(c,pqtests[i].x);
    mp_read_decimal(d,pqtests[i].y);
    mp_cornacchia_pq(a,b,e,z,s,t,err);
    if (err<>0) or mp_is_ne(c,s) or mp_is_ne(d,t) then begin
      writeln('Failure mp_cornacchia_pq for i=',i);
      inc(f);
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure test_rnr;
type
  tr = array[1..4] of pchar8;
const
  {Maple VR4: m:=2^63; q:=isqrt(iquo(m,2)); u:= 100000000000006 mod m; iratrecon(u,m,q,q,'n','d'); n;d;}
  tests: array[1..14] of tr = (
         ('9223372036854775808', '100000000000006', '90545726',    '1401491381'),
         ('9223372036854775808', '100000000000006', '90545726',    '1401491381'),
         ('9223372036854775808', '100000000000005', '-1310945655', '1401491381'),
         ('9223372036854775808', '100000000000004', '0', '0'),
         ('9223372036854775808', '100000000000003', '0', '0'),
         ('9223372036854775808', '100000000000002', '0', '0'),
         ('9223372036854775808', '100000000000001', '0', '0'), {1492323128/1031080760}  {gcd<>1!}
         ('9223372036854775808', '100000000000000', '0', '0'), {461242368/1031080760}   {gcd<>1!}

         ('633825300114114700748351602688',  '605853070219409538428464470774',  '14', '997'),
         ('633825300114114700748351602688',  '422762110908612112334657789155',  '15', '997'),
         ('633825300114114700748351602688',  '422762110908612112334657789155',  '15', '997'),
         ('633825300114114700748351602688',  '239671151597814686240851107536',  '16', '997'),
         ('633825300114114700748351602688',  '394154148516300014507500495152', '-16', '997'),
         ('633825300114114700748351602688',  '0', '0', '1'));
const
  test2: array[1..18] of array[0..5] of shortint = (
           ( 1, 1,  1, 1,  1, 1),
           ( 2, 1,  2, 1,  2, 1),
           ( 3, 1,  0, 0,  3, 1),
           ( 0, 0,  0, 0,  4, 1),
           ( 0, 0,  1, 4,  0, 0),
           (-1, 3, -1, 3,  0, 0),
           ( 2, 3,  2, 3,  0, 0),
           (-3, 2,  0, 0, -3, 2),
           (-1, 2, -1, 2, -1, 2),
           ( 1, 2,  1, 2,  1, 2),
           ( 3, 2,  0, 0,  3, 2),
           (-2, 3, -2, 3,  0, 0),
           ( 1, 3,  1, 3,  0, 0),
           ( 0, 0, -1, 4,  0, 0),
           ( 0, 0,  0, 0, -4, 1),
           (-3, 1,  0, 0, -3, 1),
           (-2, 1, -2, 1, -2, 1),
           (-1, 1, -1, 1, -1, 1));
var
  i,j,cnt,f: integer;
begin
  f := 0;
  cnt := 0;
  write('mp_rnr/mp_rnr2 test ');
  for i:=1 to 14 do begin
    write('.');
    mp_read_decimal(b,tests[i][1]);
    mp_read_decimal(a,tests[i][2]);
    mp_read_decimal(c,tests[i][3]);
    mp_read_decimal(d,tests[i][4]);
    mp_rnr(a,b,s,t);
    inc(cnt);
    if mp_is_ne(c,t) or mp_is_ne(d,s) then begin
      inc(f);
      writeln('mp_rnr failed for ', i);
      writeln('a=',mp_decimal(a));
      writeln('b=',mp_decimal(b));
      writeln('c=',mp_decimal(c));
      writeln('d=',mp_decimal(d));
      writeln('s=',mp_decimal(s));
      writeln('t=',mp_decimal(t));
    end;
    mp_shr(b,1,e);
    mp_sqrt(e,e);
    mp_rnr2(a,b,e,e,s,t);
    inc(cnt);
    if mp_is_ne(c,s) or mp_is_ne(d,t) then begin
      inc(f);
      writeln('mp_rnr2 failed for ', i);
      writeln('a=',mp_decimal(a));
      writeln('b=',mp_decimal(b));
      writeln('c=',mp_decimal(c));
      writeln('d=',mp_decimal(d));
      writeln('s=',mp_decimal(s));
      writeln('t=',mp_decimal(t));
    end;
  end;

  {Test values from M. Monagan (MQRR.PDF), except for i=5,j=1. There seems}
  {to be a superfluous '-'. The mp_rnr2 result '1/4' is confirmed by Maple.}
  mp_set_short(b,19);
  for i:=1 to 18 do begin
    write('.');
    mp_set_short(a,i);
    for j:=0 to 2 do begin
      if j=0 then begin
        mp_set(c,3);
        mp_set(d,3);
      end
      else if j=1 then begin
        mp_set(c,2);
        mp_set(d,4);
      end
      else begin
        mp_set(c,4);
        mp_set(d,2);
      end;
      mp_set_short(e,test2[i][2*j]);
      mp_set_short(s,test2[i][2*j+1]);
      mp_rnr2(a,b,c,d,t,z);
      inc(cnt);
      if mp_is_ne(e,t) or mp_is_ne(s,z) then begin
        inc(f);
        writeln('mp_rnr2 failed for i/j = ', i, '/', j);
        writeln('a=',mp_decimal(a));
        writeln('c=',mp_decimal(c));
        writeln('d=',mp_decimal(d));
        writeln('e=',mp_decimal(e));
        writeln('s=',mp_decimal(s));
        writeln('t=',mp_decimal(t));
        writeln('z=',mp_decimal(z));
      end;
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure powerd_test;
var
  i,f: integer;
  n: longint;
begin
  writeln('mp_powerd test: ', fn_powerd);
  if open_file(fn_powerd) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    inc(i);
    readline;
    mp_read_radix(a,line,10);
    readline;
    mp_read_radix(b,line,10);
    readline;
    mp_read_radix(c,line,10);
    readline;
    mp_read_radix(d,line,10);
    readline;
    mp_read_radix(e,line,10);
    readline;
    mp_read_radix(s,line,10);
    if not mp_is_longint(d,n) then begin
      writeln('Failed! #',i,' no longint: ', mp_decimal(d));
      inc(f);
    end
    else begin
      mp_powerd(a,b,c,n,t,z);
      if mp_is_ne(t,e) or mp_is_ne(z,s) then begin
        writeln('Failed! #',i);
        inc(f);
      end;
    end;
  end;
  close(tf);
  writeln('No. of checks: ',i,', failed: ',f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure binomial_test;
var
  i,cnt,f: integer;
  n,k,r1,r2: longint;
{$ifdef BIT16}
const
  maxf = 12000;
{$else}
const
  maxf = 36000;
{$endif}

begin
  f := 0;
  cnt := 0;
  writeln('mp_binomial test: ');

  mp_binomial(100000, 70000, a);
  mp_mod_int(a, longint(1) shl 30, r1);
  mp_mod_int(a, MaxLongint, r2);
  inc(cnt);
  if (r1<>734402688) or (r2<>720191980) then begin
    inc(f);
    writeln('Failed for n=100000, k=70000');
  end;

  mp_binomial(10000000, 60, a);
  mp_mod_int(a, longint(1) shl 30, r1);
  mp_mod_int(a, MaxLongint, r2);
  inc(cnt);
  if (r1<>103482144) or (r2<>1218222532) then begin
    inc(f);
    writeln('Failed for n=10000000, k=60');
  end;

  for i:= 1 to 18 do begin
    write('.');
    repeat
      n := random(maxf);
    until n > 2;
    k := random(n);
    mp_product(a,n-k+1,n);
    mp_fact(k,b);
    mp_div(a,b,a);
    mp_binomial(n,k,c);
    inc(cnt);
    if mp_is_ne(a,c) then begin
      inc(f);
      writeln('Failed for n=',n, ', k=',k);
    end;
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure swing_test;
var
  i,cnt,f: integer;
  n: longint;
{$ifdef BIT16}
const
  maxf = 2000;
{$else}
const
  maxf = 10000;
{$endif}

begin
  f := 0;
  cnt := 0;
  writeln('mp_(fact)_swing test: ');

  for i:= 1 to 10 do begin
    repeat
      n := random(maxf);
    until n > 2;
    n := n*10;
    write('.');
    mp_fact(n,a);
    mp_fact_swing(n,b);
    if mp_is_ne(a,b) then begin
      inc(f);
      writeln('Failed: mp_fact_swing for n=',n);
    end;
    mp_swing(n,b);
    mp_fact(n div 2, c);
    mp_sqr(c,d);
    mp_mul(d,b,d);
    if mp_is_ne(a,d) then begin
      inc(f);
      writeln('Failed: mp_swing for n=',n);
    end;
    inc(cnt);
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure digitsum_test;
var
  i,cnt,f: integer;
  n: longint;
const
  ds: array[1..8] of longint =
       (274, 380, 622, 1471, {2^456, radix 3,5,10,36}
        63406, 86272, 135178, 339411); {2^100000, radix 3,5,10,36}
        {Mathematica:  DigitSum[n_, b_:10] := Total[IntegerDigits[n, b]]}
        {Maple: digitsum := (n,b) -> convert(convert(n, base, b), `+`);}
  ri: array[1..8] of byte = (3,5,10,36,3,5,10,36);

begin
  f := 0;
  cnt := 0;
  write('mp_digit_sum test: ');

  for i:= 1 to 8 do begin
    write('.');
    if i=1 then mp_set_pow(a,2,456)
    else if i=5 then mp_set_pow(a,2,100000);
    n := mp_digitsum(a,ri[i]);
    if n<>ds[i] then begin
      inc(f);
      writeln('Failed for i=',i);
    end;
    inc(cnt);
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure reverse_check;
var
  i,cnt,f: integer;
  r,q: mp_digit;
  bs,sa,sb: longint;
{$ifdef BIT16}
  const BDIV = 22;
{$else}
  const BDIV = 15;
{$endif}
begin
  f := 0;
  cnt := 1;
  writeln('mp_reverse test: ');
  mp_read_decimal(a, '7544682451002378345969236876234868272345');
  mp_read_decimal(b, '5432728684326786329695438732001542864457');
  mp_reverse(a,10,c);
  if not mp_is_eq(b,c) then inc(f);
  bs := 100;
  for i:=1 to 40 do begin
    write('.');
    for r:=2 to MaxRadix do begin
      mp_rand_bits(a, bs);
      repeat
        mp_div_d(a,r,@b,q);
        if q<>0 then break;
        mp_exch(a,b);
      until a.used=0;
      mp_reverse(a,r,b);
      sa := mp_digitsum(a,r);
      sb := mp_digitsum(b,r);
      mp_reverse(b,r,c);
      inc(cnt);
      if (not mp_is_eq(a,c)) or (sa<>sb) then begin
        inc(f);
        mp_write_radix(output,a,r); writeln;
        mp_write_radix(output,b,r); writeln;
        mp_write_radix(output,c,r); writeln;
      end;
    end;
    bs := bs + bs div BDIV;
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;

end.
