{Test program for MPArith/mp_real, (c) W.Ehrhardt 2007-2017}

program t_chkfp;

{$i STD.INC}
{$i mp_conf.inc}

{$x+}  {pchar I/O}
{$i+}  {RTE on I/O error}

{$ifdef VER70}
{$N+}
{$endif}


{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
    WinCRT,
  {$endif}
  {$ifdef FPC}
    {cmem,}
  {$endif}
  BTypes, mp_types, mp_base, mp_supp, mp_real;



{---------------------------------------------------------------------------}
procedure mpf_randomx(var a: mp_float);
  {-random numbers u*2^x for test, u and x are uniform random}
  { u from -0.5 .. 0.5,  x from -2.5..2.5*2^default_prec}
begin
  mpf_random(a);
  mpf_sub_dbl(a,0.5,a);
  mpf_mul_2k(a,round((random-0.5)*5*mpf_get_default_prec),a);
end;


var
  a,b,c,d,e,s,t,z: mp_float;

const
  LMAX=16000;

var
  tf: text;
  totalfailed: longint;
  line: array[0..LMAX] of char8;
  lin2: array[0..LMAX] of char8;


{---------------------------------------------------------------------------}
{------------------- Recycled mp_int tests ---------------------------------}
{---------------------------------------------------------------------------}

const
  Prefix = 'm#';

const
  fn_add  = Prefix+'add.dat';
  fn_div  = Prefix+'div.dat';
  fn_mul1 = Prefix+'mul1.dat';
  fn_mul2 = Prefix+'mul2.dat';
  fn_expt = Prefix+'expt.dat';
  fn_sqr  = Prefix+'sqr.dat';
  fn_sub  = Prefix+'sub.dat';


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
  writeln('mpf_add test: ', fn_add);
  if open_file(fn_add) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpf_read_decimal(a,line);
    readline;
    mpf_read_decimal(b,line);
    readline;
    mpf_read_decimal(c,line);
    case random(3) of
      0: begin
           mpf_copy(a,t);
           mpf_add(t,b,t);
         end;
      1: begin
           mpf_copy(b,t);
           mpf_add(a,t,t);
         end;
      2: mpf_add(a,b,t);
    end;

    if not mpf_is_eq(t,c) then begin
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
  writeln('mpf_sub test: ', fn_sub);
  if open_file(fn_sub) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpf_read_decimal(a,line);
    readline;
    mpf_read_decimal(b,line);
    readline;
    mpf_read_decimal(c,line);
    case random(3) of
      0: begin
           mpf_copy(a,t);
           mpf_sub(t,b,t);
         end;
      1: begin
           mpf_copy(b,t);
           mpf_sub(a,t,t);
         end;
      2: mpf_sub(a,b,t);
    end;
    if not mpf_is_eq(t,c) then begin
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
  writeln('mpf_mul test 1: ', fn_mul1);
  if open_file(fn_mul1) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpf_read_decimal(a,line);
    readline;
    mpf_read_decimal(b,line);
    readline;
    mpf_read_decimal(c,line);
    case random(3) of
      0: begin
           mpf_copy(a,t);
           mpf_mul(t,b,t);
         end;
      1: begin
           mpf_copy(b,t);
           mpf_mul(a,t,t);
         end;
      2: mpf_mul(a,b,t);
    end;
    if not mpf_is_eq_rel(t,c) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' b: ', mpf_decimal(b,60));
      writeln(' c: ', mpf_decimal(c,60));
      writeln(' t: ', mpf_decimal(t,60));
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
  writeln('mpf_mul test 2: ', fn_mul2);
  if open_file(fn_mul2) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpf_read_decimal(a,line);
    readline;
    mpf_read_decimal(b,line);
    readline;
    mpf_read_decimal(c,line);
    case random(3) of
      0: begin
           mpf_copy(a,t);
           mpf_mul(t,b,t);
         end;
      1: begin
           mpf_copy(b,t);
           mpf_mul(a,t,t);
         end;
      2: mpf_mul(a,b,t);
    end;
    if not mpf_is_eq_rel(t,c) then begin
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
  writeln('mpf_sqr test: ', fn_sqr);
  if open_file(fn_sqr) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpf_read_decimal(a,line);
    readline;
    mpf_read_decimal(c,line);
    case random(2) of
      0: begin
           mpf_copy(a,t);
           mpf_sqr(t,t);
         end;
      1: mpf_sqr(a,t);
    end;
    if not mpf_is_eq_rel(t,c) then begin
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
procedure div_test;
var
  i,f: integer;
begin
  writeln('mpf_div/mul/int test: ', fn_div);
  if open_file(fn_div) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    readln(tf);
    readline;
    mpf_read_decimal(a,line);
    readline;
    mpf_read_decimal(b,line);
    readline;
    mpf_read_decimal(c,line);
    readline;
    mpf_read_decimal(d,line);
    case random(3) of
      0: begin
           mpf_copy(a,s);
           mpf_div(s,b,s);
         end;
      1: begin
           mpf_copy(b,s);
           mpf_div(a,s,s);
         end;
      2: mpf_div(a,b,s);
    end;

    mpf_int(s,s);
    if not mpf_is_eq_rel(c,s)  then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' b: ', mpf_decimal(b,60));
      writeln(' c: ', mpf_decimal(c,60));
      writeln(' s: ', mpf_decimal(s,60));
    end;
    mpf_mul(s,b,t);
    mpf_sub(a,t,t);
    if not mpf_is_eq_rel(d,t) then begin
      writeln('Failed! #',i);
      inc(f);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' b: ', mpf_decimal(b,60));
      writeln(' d: ', mpf_decimal(t,60));
      writeln(' t: ', mpf_decimal(t,60));
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
  i,f: integer;
begin
  writeln('mpf_expt test: ', fn_expt);
  if open_file(fn_expt) then exit;
  i:=1;
  f:=0;
  while not eof(tf) do begin
    {write(i,#13);}
    readln(tf);
    readline;
    mpf_read_decimal(a,line);
    readline;
    mpf_read_decimal(b,line);
    readline;
    mpf_read_decimal(c,line);
    if not s_mpf_is_neg(a) then begin
      case random(3) of
        0: begin
             mpf_copy(a,t);
             mpf_expt(t,b,t);
           end;
        1: begin
             mpf_copy(b,t);
             mpf_expt(a,t,t);
           end;
        2: mpf_expt(a,b,t);
      end;
      if not mpf_is_eq_rel(t,c) then begin
        writeln('Failed! #',i);
        writeln(' a: ', mpf_decimal(a,60));
        writeln(' b: ', mpf_decimal(b,60));
        writeln(' c: ', mpf_decimal(c,60));
        writeln(' t: ', mpf_decimal(t,60));
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
{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure mpf_pi_agm(var a: mp_float);
  {-Calculate a = pi with Brent/Salamin AGM algorithm / s_mpf_agm}
var
  c,s: mp_float;
  aprec: longint;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_pi_agm');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  aprec := a.bitprec;

  mpf_initp2(c,s,a.bitprec+32);
  if mp_error<>MP_OKAY then exit;

  mpf_set_dbl(s,0.5);
  mpf_sqrt(s,s);
  mpf_set1(a);
  s_mpf_agm(a,s,c,@s);

  {s = 0.5  - s}
  mpf_set_dbl(a,0.5);
  mpf_sub(a,s,s);

  {c = 2*AGM^2}
  mpf_sqr(c,c);
  s_mpf_incexp(c,1);

  {pi = 2*AGM^2/(0.5-s)}
  mpf_div(c,s,a);

  mpf_chg_prec(a,aprec);

  mpf_clear2(c,s);
end;


{---------------------------------------------------------------------------}
procedure mpf_pi_chud(var x: mp_float; pp: longint);
  {-Calculate x = pi with Chudnovsky formula}
const
  k1: double = 545140134.0;
  k2: double = 13591409.0;
  k3: double = 1823176476672000.0;
var
  h: mp_float;
  a,b: mp_int;
  alpha,alpha2: double;
  n,n1,prec,lp: longint;
begin
  {Ref: M.Martin[8], LiDIA[22]}
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(x) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_pi_chud');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {Uses increasing precision lp in loop}
  alpha2 := 47.110413;
  prec   := pp + 8;
  if prec >= 128 then lp := 128 else lp := prec;
  alpha := lp;

  mpf_initp(h,lp);
  if mp_error<>MP_OKAY then exit;

  mp_init2(a,b);
  if mp_error<>MP_OKAY then begin
    mpf_clear2(x,h);
    exit;
  end;

  mp_set_int(b,100100025);
  mp_mul_int(b,327843840,b);

  n := trunc(1.0 + prec/alpha2);
  n1:= 6*n-1;

  mpf_set_int(x,n);
  mpf_mul_dbl(x,k1,x);
  mpf_add_dbl(x,k2,x);
  mpf_copy(x,h);

  while (n>0) and (mp_error=MP_OKAY) do begin
    if n1 <= 1290 then begin
      mp_set_int(a,n1*(n1-2)*(n1-4));
      mpf_mul_mpi(x,a,x);
      mp_set_int(a,n*n*n);
    end
    else if n1 <= 46340 then begin
      mp_set_int(a,(n1-4)*(n1-2));
      mp_mul_w(a,n1,a);
      mpf_mul_mpi(x,a,x);
      mp_set_int(a,n*n);
      mp_mul_w(a,n,a);
    end
    else begin
      mp_set_int(a,n1);
      mp_mul_int(a,n1-2,a);
      mp_mul_int(a,n1-4,a);
      mpf_mul_mpi(x,a,x);
      mp_set_int(a,n);
      mp_mul_int(a,n,a);
      mp_mul_int(a,n,a);
    end;
    mp_mul(a,b,a);
    mpf_div_mpi(x,a,x);
    mpf_sub_dbl(h,k1,h);
    mpf_sub(h,x,x);
    {increase precision}
    alpha := alpha + alpha2;
    lp := trunc(1.0 + alpha);
    if lp > prec then lp := prec;
    mpf_chg_prec(x,lp);
    mpf_chg_prec(h,lp);
    dec(n);
    dec(n1,6);
  end;

  mpf_set_dbl(h,k3);
  mpf_sqrt(h,h);
  mpf_div(h,x,x);
  mpf_chg_prec(x,pp);

  mpf_clear(h);
  mp_clear2(a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_pi_machin(var a: mp_float; pp: longint);
  {-Calculate a = pi with Machin/Gauss formula}
var
  sump: mp_int;
  ndigits: word;
const
  Radix = 16;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_pi_machin');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  pp := pp and (not 7);
  ndigits := 16 + pp div 4;
  mp_init(sump);

  {calculate mantissa with mp_arctanw using Machin/Gauss formula}
  {pi/4 = 12*arctan(1/18) + 8*arctan(1/57) - 5*arctan(1/239)}
  mp_arctanw(12, 18, 16, ndigits, a.mantissa);
  mp_arctanw( 8, 57, Radix, ndigits, sump);
  mp_add(a.mantissa, sump, a.mantissa);
  mp_arctanw( 5,239, Radix, ndigits, sump);
  mp_sub(a.mantissa, sump, a.mantissa);

  mp_clear(sump);

  {convert to mp_float}
  a.exponent := 2-mp_bitsize(a.mantissa);
  a.bitprec  := pp;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure test_pi;
  {-Pi calculation: mpf_set_pip. Machin, AGM, Chudnovsky, 6atan(1/û3)}
var
  pic, pim, pis, pia, x: mp_float;
{$ifdef BIT16}
const
  piprec = 4000;
{$else}
const
  piprec = 20000;
{$endif}
begin

  mpf_initp5(pia, pic, pim, pis, x, piprec);
  writeln('Pi test, decimal precision = ',piprec*ln(2)/ln(10):1:1, '  (default for a1+a2+a3)');

  mpf_set_pip(pis, piprec);
  writeln(' mpf_set_pip=',mpf_decimal(pis,60));

  mpf_pi_machin(pim, piprec);
  writeln('      Machin=',mpf_decimal(pim,60));

  mpf_pi_agm(pia);
  writeln('         AGM=',mpf_decimal(pia,60));

  mpf_pi_chud(pic, piprec);
  writeln('  Chudnovsky=',mpf_decimal(pic,60));

  mpf_set_dbl(x,3);
  mpf_sqrt(x,x);
  mpf_div_d(x,3,x);
  mpf_arctan(x,x);
  mpf_mul_int(x,6,x);
  writeln(' 6atan(1/û3)=',mpf_decimal(x,60));

  {pi=arctan(1)+arctan(2)+arctan(3) with def.prec.}
  mpf_set1(a);
  mpf_arctan(a,z);

  mpf_set_dbl(a,2.0);
  mpf_arctan(a,b);
  mpf_add(z,b,z);

  mpf_set_dbl(a,3.0);
  mpf_arctan(a,b);
  mpf_add(z,b,z);
  writeln('    a1+a2+a3=',mpf_decimal(z,60));


  if not mpf_is_eq_rel(x,pis) then begin
    inc(totalfailed);
    mpf_sub(x,pis,x);
    writeln('Failed! Diff atan=',mpf_decimal(x,30));
  end;
  if not mpf_is_eq_rel(pic,pis) then begin
    inc(totalfailed);
    mpf_sub(pic,pis,x);
    writeln('Failed! Diff Chud=',mpf_decimal(x,30));
  end;
  if not mpf_is_eq_rel(pim,pis) then begin
    inc(totalfailed);
    mpf_sub(pim, pis, x);
    writeln('Failed! Diff Mach=',mpf_decimal(x,30));
  end;
  if not mpf_is_eq_rel(pia,pis) then begin
    inc(totalfailed);
    mpf_sub(pia, pis, x);
    writeln('Failed! Diff  AGM=',mpf_decimal(x,30));
  end;

  if not mpf_is_eq_rel(pis,z) then begin
    inc(totalfailed);
    mpf_sub(z, pis, x);
    writeln('Failed! Diff a123=',mpf_decimal(x,30));
  end;

  writeln;

  mpf_clear5(pia, pic, pim, pis, x);
end;


{---------------------------------------------------------------------------}
procedure testradixconv;
  {-Test input/output radix conversion}
var
  cnt,f, radix, ndd, prd: word;
  OK: boolean;
  i: integer;
{$ifndef BIT16}
 const dmax=800;
{$else}
 const dmax=200;
{$endif}
var
  u,v,w,x,y: mp_float;
begin
  mpf_initp5(u,v,w,x,y,mpf_get_default_prec);
  writeln('Test radix conversion (MaxRadix=',MaxRadix,'):  ');
  f := 0;
  cnt := 0;
  ndd := 2;
  while ndd<=dmax do begin
    for radix:=2 to MaxRadix do begin
      mpf_set_int(x,radix);
      mpf_expt_int(x, longint(1)-ndd, x);
      prd := ndd;
      if Radix>2 then prd := round(ndd/lograd[radix]);
      inc(prd,32);

      mpf_chg_prec(u,prd);
      mpf_chg_prec(v,prd);
      mpf_chg_prec(w,prd);
      mpf_chg_prec(x,prd);
      mpf_chg_prec(y,prd);

      mpf_randomx(u);
      mp_show_plus := random(2)=0;
      mp_uppercase := random(2)=0;
      mpf_toradix_n(u, radix, ndd, line, LMAX);
      mpf_read_radix(v,line,radix);
      mpf_toradix_n(v, radix, ndd, lin2, LMAX);
      inc(cnt);


      if mpf_is0(u) then begin
        OK := mpf_is0(v);
      end
      else begin
        mpf_sub(u,v,w);
        mpf_div(w,u,w);
        mpf_abs(w,w);
        OK := mpf_is_le(w,x);
      end;

      if OK then begin
        {In/Out conversion sequence: X1 -> A1 -> X2 -> A2 -> X3 -> A3}
        {Here X1 near X2, A1 <> A2 is possible, but A2 = A3 ... should be stable}
        mpf_read_radix(y, lin2, radix);
        mpf_toradix_n(y, radix, ndd, line, LMAX);
        i := 0;
        while line[i]<>#0 do begin
          if line[i]<>lin2[i] then begin
            OK := false;
            break;
          end;
          inc(i);
        end;
      end;

      if not OK then begin
        writeln('Diff radix=',radix, '  digits=',ndd);
        writeln('   u:', mpf_decimal(u,60));
        writeln('   v:', mpf_decimal(v,60));
        writeln('   w:', mpf_decimal(w,60));
        writeln('out1:', line);
        writeln('out2:', lin2);
        inc(f);
      end;
    end;
    write('.');
    ndd := 2 + 11*ndd div 10;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
  mpf_clear5(u,v,w,x,y);

end;


{---------------------------------------------------------------------------}
procedure testradixconv_alt;
  {-Test input/output radix conversion}
var
  cnt,f, radix, ndd, prd: word;
  OK: boolean;
  i: integer;
{$ifndef BIT16}
 const dmax=800;
{$else}
 const dmax=200;
{$endif}
var
  u,v,w,x,y: mp_float;
begin
  mpf_initp5(u,v,w,x,y,mpf_get_default_prec);
  writeln('Test alternative radix conversion (MaxRadix=',MaxRadix,'):  ');
  f := 0;
  cnt := 0;
  ndd := 2;
  while ndd<=dmax do begin
    for radix:=2 to MaxRadix do begin
      mpf_set_int(x,radix);
      mpf_expt_int(x, longint(1)-ndd, x);
      prd := ndd;
      if Radix>2 then prd := round(ndd/lograd[radix]);
      inc(prd,32);

      mpf_chg_prec(u,prd);
      mpf_chg_prec(v,prd);
      mpf_chg_prec(w,prd);
      mpf_chg_prec(x,prd);
      mpf_chg_prec(y,prd);

      mpf_randomx(u);
      mp_show_plus := random(2)=0;
      mp_uppercase := random(2)=0;
      mpf_toradix_alt(u, radix, ndd, line, LMAX);
      mpf_read_radix(v,line,radix);
      mpf_toradix_alt(v, radix, ndd, lin2, LMAX);
      inc(cnt);


      if mpf_is0(u) then begin
        OK := mpf_is0(v);
      end
      else begin
        mpf_sub(u,v,w);
        {relative error if |u| >= 1, absolute otherwise}
        if s_mpf_ldx(u)>0 then mpf_div(w,u,w);
        mpf_abs(w,w);
        OK := mpf_is_le(w,x);
        {$ifdef debug}
          if not OK then
            OK := OK;
        {$endif}
      end;

      if OK then begin
        {In/Out conversion sequence: X1 -> A1 -> X2 -> A2 -> X3 -> A3}
        {Here X1 near X2, A1 <> A2 is possible, but A2 = A3 ... should be stable}
        mpf_read_radix(y, lin2, radix);
        mpf_toradix_alt(y, radix, ndd, line, LMAX);
        i := 0;
        while line[i]<>#0 do begin
          if line[i]<>lin2[i] then begin
            OK := false;
            break;
          end;
          inc(i);
        end;
      end;

      if not OK then begin
        writeln('Diff radix=',radix, '  digits=',ndd);
        writeln('   u:', mpf_decimal(u,60));
        writeln('   v:', mpf_decimal(v,60));
        writeln('   w:', mpf_decimal(w,60));
        writeln('   x:', mpf_decimal(x,60));
        writeln('out1:', line);
        writeln('out2:', lin2);
        inc(f);
      end;
    end;
    write('.');
    ndd := 2 + 11*ndd div 10;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
  mpf_clear5(u,v,w,x,y);

end;


{---------------------------------------------------------------------------}
procedure transz_test1;
  {-Function test 1 (inverse function of function)}
const
  NTT=60;
{$ifndef BIT16}
const
  FAC = 32;
{$else}
const
  FAC = 4;
{$endif}
var
  f,i,cnt: integer;
  r,r_at,r_le,r_ss,r_11,r_xx: double;
begin
  f := 0;
  cnt := 0;
  writeln('Function test 1 (inverse function of function): ');
  mpf_set_pi2k(e,-1);
  r_at := 0.0;
  r_le := 0.0;
  r_ss := 0.0;
  r_11 := 0.0;
  r_xx := 0.0;
  for i:=1 to NTT*FAC do begin
    if i mod FAC = 0 then write('.');
    mpf_random(a);
    mpf_mul(a,e,a);
    mpf_tan(a,b);
    {a in [0,pi/2), b in [0, inf)}
    inc(cnt);
    mpf_arctan(b,t);
    r := mpf_reldev(t,a);
    if r>r_at then r_at := r;
    if r>8 then begin
      mpf_sub(a,t,d);
      inc(f);
      writeln('Failed arctan(tan(a)) <> a: ',r:1:1);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' t: ', mpf_decimal(t,60));
      writeln(' d: ', mpf_decimal(d,60));
    end;
    if not s_mpf_is_le0(b) then begin
      inc(cnt);
      mpf_ln(b,c);
      mpf_exp(c,t);
      r := mpf_reldev(t,b);
      if r>r_le then r_le := r;
      if r>8 then begin
        mpf_sub(b,t,d);
        inc(f);
        writeln('Failed exp(ln(b)) <> b: ',r:1:1);
        writeln(' b: ', mpf_decimal(b,60));
        writeln(' t: ', mpf_decimal(t,60));
        writeln(' d: ', mpf_decimal(d,60));
      end;
    end;
    inc(cnt);
    mpf_sqrt(b,c);
    mpf_sqr(c,t);
    r := mpf_reldev(t,b);
    if r>r_ss then r_ss := r;
    if r>4 then begin
      mpf_sub(b,t,d);
      inc(f);
      writeln('Failed sqr(sqrt(b)) <> b: ',r:1:1);
      writeln(' b: ', mpf_decimal(b,60));
      writeln(' t: ', mpf_decimal(t,60));
      writeln(' d: ', mpf_decimal(d,60));
    end;

    inc(cnt);
    mpf_random(a);
    mpf_sub_dbl(a,0.5,a);

    mpf_ln1p(a,c);
    mpf_expm1(c,t);
    r := mpf_reldev(t,a);
    if r>r_11 then r_11 := r;
    if r>8 then begin
      mpf_sub(a,t,d);
      inc(f);
      writeln('Failed expm1(ln1p(a)) <> a: ',r:1:1);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' t: ', mpf_decimal(t,60));
      writeln(' d: ', mpf_decimal(d,60));
    end;

    inc(cnt);
    mpf_expm1(a,c);
    mpf_ln1p(c,t);
    r := mpf_reldev(t,a);
    if r>r_11 then r_11 := r;
    if r>8 then begin
      mpf_sub(a,t,d);
      inc(f);
      writeln('Failed ln1p(expm1(a)) <> a: ',r:1:1);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' t: ', mpf_decimal(t,60));
      writeln(' d: ', mpf_decimal(d,60));
    end;

    mpf_random(a);
    mpf_sub_dbl(a,0.5,a);
    mpf_mul_d(a,2,a);

    inc(cnt);
    mpf_arcsin(a,c);
    mpf_sin(c,t);
    r := mpf_reldev(t,a);
    if r>r_xx then r_xx := r;
    if r>8 then begin
      mpf_sub(a,t,d);
      inc(f);
      writeln('Failed sin(arcsin(a)) <> a: ',r:1:1);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' t: ', mpf_decimal(t,60));
      writeln(' d: ', mpf_decimal(d,60));
    end;

    inc(cnt);
    mpf_arccos(a,c);
    mpf_cos(c,t);
    {calculate absolute error, because for small a arccos(a) is near pi/2}
    {and only few low bits change (if any), e.g. for very small arguments}
    {arccos is computationally constant and relative error is meaningless.}
    {Same argument applies to cos because for small a is constant ~ 1.}
    mpf_sub(a,t,d);
    mpf_mul_2k(d,a.bitprec,b);
    r := b.exponent + mp_bitsize(b.mantissa);
    if r > 3 then begin
      inc(f);
      writeln('Failed cos(arccos(a)) <> a: ',r:1:1);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' t: ', mpf_decimal(t,60));
      writeln(' d: ', mpf_decimal(d,60));
    end;

    inc(cnt);
    mpf_arctanh(a,c);
    mpf_tanh(c,t);
    r := mpf_reldev(t,a);
    if r>r_xx then r_xx := r;
    if r>8 then begin
      mpf_sub(a,t,d);
      inc(f);
      writeln('Failed tanh(atanh(a)) <> a: ',r:1:1);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' t: ', mpf_decimal(t,60));
      writeln(' d: ', mpf_decimal(d,60));
    end;

    inc(cnt);
    mpf_arcsinh(a,c);
    mpf_sinh(c,t);
    r := mpf_reldev(t,a);
    if r>r_xx then r_xx := r;
    if r>8 then begin
      mpf_sub(a,t,d);
      inc(f);
      writeln('Failed sinh(arcsinh(a)) <> a: ',r:1:1);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' t: ', mpf_decimal(t,60));
      writeln(' d: ', mpf_decimal(d,60));
    end;

    if not s_mpf_is_le0(b) then begin
      inc(cnt);
      mpf_log2(b,c);
      mpf_exp2(c,t);
      r := mpf_reldev(t,b);
      if r>r_xx then r_xx := r;
      if r>8 then begin
        mpf_sub(b,t,d);
        inc(f);
        writeln('Failed 2^(log2(b)) <> b: ',r:1:1);
        writeln(' b: ', mpf_decimal(b,60));
        writeln(' t: ', mpf_decimal(t,60));
        writeln(' d: ', mpf_decimal(d,60));
      end;
      inc(cnt);
      mpf_log10(b,c);
      mpf_exp10(c,t);
      r := mpf_reldev(t,b);
      if r>r_xx then r_xx := r;
      if r>8 then begin
        mpf_sub(b,t,d);
        inc(f);
        writeln('Failed 10^(log10(b)) <> b: ',r:1:1);
        writeln(' b: ', mpf_decimal(b,60));
        writeln(' t: ', mpf_decimal(t,60));
        writeln(' d: ', mpf_decimal(d,60));
      end;
    end;

  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  writeln('Max Rel Dev (AT/LE/SS/L1E1/other): ', r_at:1:1, ' / ',r_le:1:1, ' / ',r_ss:1:1, ' / ',r_11:1:1, ' / ',r_xx:1:1);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure func_test;
  {-Tests for other elementary functions}
type
  t_testfunc = procedure(const a: mp_float; var b: mp_float);
const
  NUMFUN = 58;
var
  u,v,w: mp_float;
  f,i,k,cnt: integer;
  tf: t_testfunc;
  ts: string[20];
  skip: boolean;
  rd: double;
  fcnt: array[0..NUMFUN-1] of integer;
{$ifndef BIT16}
const
  imax = 32000;
  imsk = $1FF;
{$else}
const
  imax = 4800;
  imsk = $3f;
{$endif}
begin
  writeln('Function test 2 (compare with double precision):  ');
  mpf_initp3(u,v,w,mpf_get_default_prec*2);
  fillchar(fcnt,sizeof(fcnt),0);
  f := 0;
  cnt := 0;
  tf := nil;
  for i:=1 to imax do begin
    skip := false;
    if i and imsk = 0 then write('.');
    mpf_randomx(a);
    k := random(NUMFUN);
    inc(fcnt[k]);
    case k of
      0:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sqrt;
            ts := 'mpf_sqrt';
            s_mpf_abs(a);
          end;
      1:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sqr;
            ts := 'mpf_sqr';
          end;
      2:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_ln;
            ts := 'mpf_ln';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
            s_mpf_abs(a);
          end;
      3:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_exp;
            ts := 'mpf_exp';
            skip := s_mpf_ldx(a) > 22;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 22;
            end;
          end;
      4:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arctan;
            ts := 'mpf_arctan';
          end;
      5:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_cos;
            ts := 'mpf_cos';
          end;
      6:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sin;
            ts := 'mpf_sin';
          end;
      7:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_tan;
            ts := 'mpf_tan';
          end;
      8:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_frac;
            ts := 'mpf_frac';
          end;
      9:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_int;
            ts := 'mpf_int';
          end;
     10:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_expm1;
            ts := 'mpf_expm1';
            skip := s_mpf_ldx(a) > 0;
            while skip do begin
              mpf_random(a);
              skip := s_mpf_ldx(a) > 0;
            end;
          end;
     11:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_ln1p;
            ts := 'mpf_ln1p';
            skip := s_mpf_ldx(a) > 0;
            while skip do begin
              mpf_random(a);
              skip := s_mpf_ldx(a) > 0;
            end;
          end;
     12:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_cosh;
            ts := 'mpf_cosh';
            skip := s_mpf_ldx(a) > 22;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 22;
            end;
          end;
     13:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sinh;
            ts := 'mpf_sinh';
            skip := s_mpf_ldx(a) > 22;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 22;
            end;
          end;
     14:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_tanh;
            ts := 'mpf_tanh';
          end;

     15:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arctanh;
            ts := 'mpf_arctanh';
            skip := s_mpf_ldx(a) > 0;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 0;
            end;
          end;
     16:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arccosh;
            ts := 'mpf_arccosh';
            skip := s_mpf_ldx(a) <= 0;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) <= 0;
            end;
            s_mpf_abs(a);
          end;
     17:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arccosh1p;
            ts := 'mpf_arccosh1p';
            s_mpf_abs(a);
          end;
     18:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arcsinh;
            ts := 'mpf_arcsinh';
          end;

     19:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arcsin;
            ts := 'mpf_arcsin';
            skip := mpf_cmp_mag_dbl(a,1.0)=MP_GT;
            while skip do begin
              mpf_randomx(a);
              skip := mpf_cmp_mag_dbl(a,1.0)=MP_GT;
            end;
          end;
     20:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arccos;
            ts := 'mpf_arccos';
            skip := mpf_cmp_mag_dbl(a,1.0)=MP_GT;
            while skip do begin
              mpf_randomx(a);
              skip := mpf_cmp_mag_dbl(a,1.0)=MP_GT;
            end;
          end;
     21:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_ccell1;
            ts := 'mpf_ccell1';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
            s_mpf_abs(a);
          end;
     22:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_ccell2;
            ts := 'mpf_ccell2';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
            s_mpf_abs(a);
          end;
     23:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_log2;
            ts := 'mpf_log2';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
            s_mpf_abs(a);
          end;
     24:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_log10;
            ts := 'mpf_log10';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
            s_mpf_abs(a);
          end;
     25:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_exp2;
            ts := 'mpf_exp2';
            skip := s_mpf_ldx(a) > 22;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 22;
            end;
          end;
     26:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_exp10;
            ts := 'mpf_exp10';
            skip := s_mpf_ldx(a) > 20;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 20;
            end;
          end;
     27:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_cot;
            ts := 'mpf_cot';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
          end;
     28:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_csc;
            ts := 'mpf_csc';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
          end;
     29:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sec;
            ts := 'mpf_sec';
          end;
     30:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_coth;
            ts := 'mpf_coth';
            skip := mpf_is0(a) or (s_mpf_ldx(a) > 22);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a) or (s_mpf_ldx(a) > 22);
            end;
          end;
     31:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_csch;
            ts := 'mpf_csch';
            skip := mpf_is1(a) or (s_mpf_ldx(a) > 22);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is1(a) or (s_mpf_ldx(a) > 22);
            end;
          end;
     32:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sech;
            ts := 'mpf_sech';
            skip := mpf_is0(a) or (s_mpf_ldx(a) > 22);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a) or (s_mpf_ldx(a) > 22);
            end;
          end;
     33:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arccot;
            ts := 'mpf_arccot';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
          end;
     34:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arccotc;
            ts := 'mpf_arccotc';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
          end;
     35:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arccoth;
            ts := 'mpf_arccoth';
            skip := mpf_is1(a) or (s_mpf_ldx(a) < 1);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is1(a) or (s_mpf_ldx(a) < 1);
            end;
          end;
     36:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arccsc;
            ts := 'mpf_arccsc';
            skip := s_mpf_ldx(a) < 1;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) < 1;
            end;
          end;
     37:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arccsch;
            ts := 'mpf_arccsch';
            skip := mpf_is0(a);
            while skip do begin
              mpf_randomx(a);
              skip := mpf_is0(a);
            end;
          end;
     38:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arcsec;
            ts := 'mpf_arcsec';
            skip := s_mpf_ldx(a) < 1;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) <1;
            end;
          end;
     39:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arcsech;
            ts := 'mpf_arcsech';
            skip := s_mpf_is_le0(a) or (s_mpf_ldx(a) > 0);
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_is_le0(a) or (s_mpf_ldx(a) > 0);
            end;
          end;
     40:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sqrt1pm1;
            ts := 'mpf_sqrt1pm1';
            skip := mpf_cmp_dbl(a,-1.0)=MP_LT;
            while skip do begin
              mpf_randomx(a);
              skip := mpf_cmp_dbl(a,-1.0)=MP_LT;
            end;
          end;
     41:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sinc;
            ts := 'mpf_sinc';
          end;
     42:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_hav;
            ts := 'mpf_hav';
          end;
     43:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_archav;
            ts := 'mpf_archav';
            skip := s_mpf_is_neg(a) or (mpf_cmp_dbl(a,1)=MP_GT);
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_is_neg(a) or (mpf_cmp_dbl(a,1)=MP_GT);
            end;
          end;
     44:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_gd;
            ts := 'mpf_gd';
          end;

     45:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_arcgd;
            ts := 'mpf_arcgd';
            skip := s_mpf_ldx(a) > 1;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 1;
            end;
          end;
     46:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_lambertw;
            ts := 'mpf_lambertw';
            skip := mpf_cmp_dbl(a,-0.367879441171)=MP_LT;
            while skip do begin
              mpf_randomx(a);
              skip := mpf_cmp_dbl(a,-0.367879441171)=MP_LT;
            end;
          end;
     47:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_coshm1;
            ts := 'mpf_coshm1';
            skip := s_mpf_ldx(a) > 4;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 4;
            end;
          end;
     48:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_exp2m1;
            ts := 'mpf_exp2m1';
            skip := s_mpf_ldx(a) > 4;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 4;
            end;
          end;
     49:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_exp10m1;
            ts := 'mpf_exp10m1';
            skip := s_mpf_ldx(a) > 4;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) > 4;
            end;
          end;
     50:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_log10p1;
            ts := 'mpf_log10p1';
            {range -1 < a < 1}
            skip := s_mpf_ldx(a) >= 1;
            while skip do begin
              mpf_random(a);
              skip := s_mpf_ldx(a) >= 1;
            end;
          end;
     51:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_log2p1;
            ts := 'mpf_log2p1';
            skip := s_mpf_ldx(a) >= 1;
            while skip do begin
              mpf_random(a);
              skip := s_mpf_ldx(a) >= 1;
            end;
          end;
     52:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_cospi;
            ts := 'mpf_cospi';
          end;
     53:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_sinpi;
            ts := 'mpf_sinpi';
          end;
     54:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_tanpi;
            ts := 'mpf_tanpi';
            skip := mpf_cmp_mag_dbl(a,0.5)=MP_GT;
            while skip do begin
              mpf_randomx(a);
              skip := mpf_cmp_mag_dbl(a,0.5)=MP_GT;
            end;
          end;
     55:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_cbrt;
            ts := 'mpf_cbrt';
          end;
     56:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_cell1;
            ts := 'mpf_cell1';
            skip := s_mpf_ldx(a) >= 1;
            while skip do begin
              mpf_randomx(a);
              skip := s_mpf_ldx(a) >= 1;
            end;
            s_mpf_abs(a);
          end;
     57:  begin
            tf := {$ifdef FPC_ProcVar}@{$endif}mpf_cell2;
            ts := 'mpf_cell2';
            skip := mpf_cmp_mag_dbl(a,1.0)=MP_GT;
            while skip do begin
              mpf_randomx(a);
              skip := mpf_cmp_mag_dbl(a,1.0)=MP_GT;
            end;
            s_mpf_abs(a);
          end;
     else begin
            ts := '';
            skip := true;
          end;
    end;

    if not skip then begin
      inc(cnt);
      mpf_copyp(a,u);
      tf(a,b);
      if mp_error<>MP_OKAY then begin
        set_mp_error(0);        {clear error variable}
        writeln(i, ', test of ',ts);
        writeln(' a: ', mpf_decimal(a,60));
      end;
      tf(u,v);
      rd := mpf_reldev(v,b);
      if rd>4 then begin
        if rd>8 then begin
          write('Failed! #');
          inc(f);
        end
        else write('Bad precision: #');
        writeln(i, ', test of ',ts, ':  rd=',rd:1:1);
        writeln(' a=', mpf_decimal(a,70));
        writeln(' b=', mpf_decimal(b,70));
        writeln(' v=', mpf_decimal(v,70));
        writeln(' x: ', s_mpf_ldx(a));
      end;
    end;
  end;
  k := imax;
  for i:=0 to NUMFUN-1 do begin
    if fcnt[i]<k then k := fcnt[i];
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f, ',  min fct cnt: ',k);
  inc(totalfailed,f);
  writeln;
  mpf_clear3(u,v,w);
end;


{---------------------------------------------------------------------------}
procedure lnagm_test;
var
  u,v,w: mp_float;
  f,i,k,cnt: integer;
  skip: boolean;
  rd: double;
  maxcp,cp,saveprec,savecut: longint;

  procedure mpf_randomp(var a: mp_float);
    {-random numbers u*2^x for test, u and x are uniform random}
    { u from -0.5 .. 0.5,  x from -2.5..2.5*2^default_prec}
  begin
    mpf_random(a);
    mpf_sub_dbl(a,0.5,a);
    mpf_mul_2k(a,round((random-0.5)*5*cp),a);
  end;

begin
  writeln('s_mpf_lnagm test:  ');
  saveprec := mpf_get_default_prec;
  savecut  := mpf_lna_cutoff;
  mpf_init3(u,v,w);
  f := 0;
  cnt := 0;
  cp := 128;
  {$ifndef BIT16}
    maxcp := MPF_MAX_PREC div 6;
  {$else}
    maxcp := MPF_MAX_PREC div 3;
  {$endif}
  repeat
    {write(cp,#13);}
    write(' ',cp);
    k := maxcp div cp;
    if k=0 then k:=1;
    {$ifndef BIT16}
      if cp<=1000 then k := 2*k;
    {$endif}
    for i:=1 to k do begin
      mpf_chg_prec(u,cp);
      mpf_chg_prec(v,cp);
      mpf_chg_prec(w,cp);
      mpf_randomp(u);
      skip := mpf_is0(u);
      while skip do begin
        mpf_randomp(u);
        skip := mpf_is0(u);
      end;
      s_mpf_abs(u);

      if not skip then begin
        inc(cnt);
        {force lnagm}
        mpf_lna_cutoff := 0;
        mpf_ln(u,v);
        if mp_error<>MP_OKAY then begin
          set_mp_error(0);        {clear error variable}
          writeln(' u: ', mpf_decimal(u,60));
        end;
        {force std ln}
        mpf_lna_cutoff := MaxLongint;
        mpf_ln(u,w);
        rd := mpf_reldev(v,w);
        if rd>2 then begin
          if rd>4 then begin
            write('Failed! #');
            inc(f);
          end
          else write('Bad precision: #');
          writeln(i, ':  rd=',rd:1:1);
          writeln(' u=', mpf_decimal(u,70));
          writeln(' v=', mpf_decimal(v,70));
          writeln(' w=', mpf_decimal(w,70));
          writeln(' x: ', s_mpf_ldx(u));
        end;
      end;
    end;
    cp := round(1.5*cp) and $FFFFFF8;
  until cp>=maxcp;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
  mpf_set_default_prec(saveprec);
  mpf_lna_cutoff := savecut;
  mpf_clear3(u,v,w);
end;


{---------------------------------------------------------------------------}
procedure func3_test;
  {-Function test 3 (special trig values)}
var
  f,cnt: longint;
  procedure singletest(const x: mp_float; n: longint; const msg: mp_string);
    {-Check if x = sign(n)*sqrt(|n|/12)}
  var
    r: double;
  begin
    mpf_set_int(d,abs(n));
    mpf_div_d(d,12,d);
    mpf_sqrt(d,d);
    if n<0 then s_mpf_chs(d);
    inc(cnt);
    r := mpf_reldev(x,d);
    if r>6 then begin
      mpf_sub(x,d,e);
      if r>8 then begin
        inc(f);
        write('Failed ');
      end
      else write('Bad precision ');
      writeln(msg,': diff=', mpf_decimal(e,30), ', r=',r:1:1);
    end;
  end;
begin
  writeln('Function test 3 (special trig values): ');
  cnt := 0;
  f := 0;

  mpf_set0(a);
  mpf_trig(a,@c,@s,@t);
  singletest(c, 12,'cos(0)');
  singletest(s,  0,'sin(0)');
  singletest(t,  0,'tan(0)');

  mpf_set_pi(a);
  mpf_div_d(a,6,a);
  mpf_trig(a,@c,@s,@t);
  singletest(c,  9,'cos(pi/6)');
  singletest(s,  3,'sin(pi/6)');
  singletest(t,  4,'tan(pi/6)');

  mpf_set_pi(a);
  mpf_div_d(a,4,a);
  mpf_trig(a,@c,@s,@t);
  singletest(c,  6,'cos(pi/4)');
  singletest(s,  6,'sin(pi/4)');
  singletest(t, 12,'tan(pi/4)');

  mpf_set_pi(a);
  mpf_div_d(a,3,a);
  mpf_trig(a,@c,@s,@t);
  singletest(c,  3,'cos(pi/3)');
  singletest(s,  9,'sin(pi/3)');
  singletest(t, 36,'tan(pi/3)');

  mpf_set_pi(a);
  mpf_div_d(a,2,a);
  mpf_trig(a,@c,@s,nil);
  singletest(c,  0,'cos(pi/2)');
  singletest(s, 12,'sin(pi/2)');

  mpf_set_pi2k(a,1);
  mpf_div_d(a,3,a);
  mpf_trig(a,@c,@s,@t);
  singletest(c ,-3,'cos(2*pi/3)');
  singletest(s  ,9,'sin(2*pi/3)');
  singletest(t,-36,'tan(2*pi/3)');

  mpf_set_pi2k(a,-2);
  mpf_mul_d(a,3,a);
  mpf_trig(a,@c,@s,@t);
  singletest(c, -6,'cos(3*pi/4)');
  singletest(s,  6,'sin(3*pi/4)');
  singletest(t,-12,'tan(3*pi/4)');

  mpf_set_pi(a);
  mpf_mul_d(a,5,a);
  mpf_div_d(a,6,a);
  mpf_trig(a,@c,@s,@t);
  singletest(c, -9,'cos(5*pi/6)');
  singletest(s,  3,'sin(5*pi/6)');
  singletest(t, -4,'tan(5*pi/6)');

  mpf_set_pi(a);
  mpf_trig(a,@c,@s,@t);
  singletest(c,-12,'cos(pi)');
  singletest(s,  0,'sin(pi)');
  singletest(t,  0,'tan(pi)');

  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;

end;


{---------------------------------------------------------------------------}
procedure arith_test;
  {-Basic arith test (compare with double precision)}
type
  t_testproc = procedure(const a,b: mp_float; var c: mp_float);
var
  u,v,w: mp_float;
  f,i,cnt: longint;
  tp: t_testproc;
  ts: string[10];
  skip: boolean;
{$ifndef BIT16}
const
  imax = 128000;
  imsk = $7FF;
{$else}
const
  imax = 16000;
  imsk = $ff;
{$endif}
begin
  writeln('Basic arith test (compare with double precision):  ');
  mpf_initp3(u,v,w,mpf_get_default_prec*2);
  f := 0;
  cnt := 0;
  tp := nil;
  for i:=1 to imax do begin
    skip := false;
    if i and imsk = 0 then write('.');
    mpf_randomx(a);
    mpf_randomx(b);
    case random(4) of
      0:  begin
           tp := {$ifdef FPC_ProcVar}@{$endif}mpf_add;
           ts := 'mpf_add';
          end;
      1:  begin
           tp := {$ifdef FPC_ProcVar}@{$endif}mpf_sub;
           ts := 'mpf_sub';
          end;
      2:  begin
           tp := {$ifdef FPC_ProcVar}@{$endif}mpf_mul;
           ts := 'mpf_mul';
          end;
      3:  begin
           skip := mpf_is0(b);
           tp := {$ifdef FPC_ProcVar}@{$endif}mpf_div;
           ts := 'mpf_div';
          end;
     else begin
           skip := true;
           ts := '';
          end;
    end;
    if not skip then begin
      inc(cnt);
      mpf_copyp(a,u);
      mpf_copyp(b,v);
      tp(a,b,c);
      tp(u,v,w);
      if not mpf_is_eq_rel(w,c) then begin
        writeln('Failed! #',i, ', test of ',ts);
        writeln(' a: ', mpf_decimal(a,60));
        writeln(' b: ', mpf_decimal(b,60));
        writeln(' c: ', mpf_decimal(c,60));
        writeln(' w: ', mpf_decimal(w,60));
        inc(f);
      end;
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
  mpf_clear3(u,v,w);
end;


{---------------------------------------------------------------------------}
procedure cmp_mag_test;
var
  f,bp: longint;

  procedure singletest(tst: integer);
  begin
    if mpf_cmp_mag(a,b)<>tst then begin
      inc(f);
      writeln(' a: ', mpf_decimal(a,60));
      writeln(' b: ', mpf_decimal(b,60));
      writeln('|a| ? |b| = ', mpf_cmp_mag(b,a), ',  should be =', tst);
    end;
  end;

begin
  f := 0;
  writeln('Test mpf_cmp_mag:  ');
  mpf_set_int(a, 1000000001);
  mpf_set_int(b, 1000010001);
  singletest(MP_LT);

  mpf_set_int(a, 1000000001);
  mpf_set_int(b,-1000010001);
  singletest(MP_LT);

  mpf_set_int(a,-1000000001);
  mpf_set_int(b, 1000010001);
  singletest(MP_LT);

  mpf_set_int(a,-1000000001);
  mpf_set_int(b,-1000010001);
  singletest(MP_LT);


  mpf_set_int(a, 1000010001);
  mpf_set_int(b, 1000000001);
  singletest(MP_GT);

  mpf_set_int(a,-1000010001);
  mpf_set_int(b, 1000000001);
  singletest(MP_GT);

  mpf_set_int(a, 1000010001);
  mpf_set_int(b,-1000000001);
  singletest(MP_GT);

  mpf_set_int(a,-1000010001);
  mpf_set_int(b,-1000000001);
  singletest(MP_GT);


  bp := b.bitprec;
  mpf_chg_prec(b, b.bitprec+1);
  mpf_set_int(a, 1000000001);
  mpf_set_int(b, 1000000001);
  singletest(MP_EQ);

  mpf_set_int(a, 1000000001);
  mpf_set_int(b,-1000000001);
  singletest(MP_EQ);

  mpf_set_int(a,-1000000001);
  mpf_set_int(b, 1000000001);
  singletest(MP_EQ);

  mpf_set_int(a,-1000000001);
  mpf_set_int(b,-1000000001);
  singletest(MP_EQ);

  mpf_chg_prec(b, bp);
  writeln('No. of checks: ',12,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure agm_elli_test;
const
  tva: array[1..20] of pchar8 = ({round(agm(1,n)*2^240), n=1..2. Caclulated with Pari/GP}
        '1766847064778384329583297500742918515827483896875618958121606201292619776',
        '2573926957200703558147887463417964957938985756017045420613346394230663419',
        '3292725843347979949271161057967495955136312416725880990796772911557788311',
        '3963088463295177174332502527717098149043913537447416237536704499296782586',
        '4600884228098463619858436804875659458428469853841701678260799680444852709',
        '5214546525688863867197396831180174227119381327646069826878459014448985188',
        '5809255302810335136162865418664341255007672308141297668321477362879142810',
        '6388488189342210894978001183229490183994993708569758472553850254803734581',
        '6954726145899682414762414990545474705074267835822378918728774983334671624',
        '7509819299794309395330362996330944286026095999027732273228810841904073216',
        '8055194849386119280151234093662571549262150720203836064694201948039599537',
        '8591983588231673480030807893363511924180596724488468582932217361625978792',
        '9121101111043922051272758019410867212013179339136245446917613971403667052',
        '9643302203321723994120661641781723734760378658547639876507028659217826232',
       '10159218566927754261734266992744678707169780923025395402573723350816097846',
       '10669385763965628162075344586229561288885425163449686966094336426213100967',
       '11174262944411611380162533132101796380671697355437999624354432620868460965',
       '11674247602317875538272741952217973482035364311655476490053837539506197100',
       '12169686820612293989049687279088381020518235386301948917116127104356713870',
       '12660885981187313277839321755922689466236688731997910507369599206106733322');
const
  tve: array[1..22] of pchar8 = ({(K'(k) E'(k))*2^240), k=0.1 .. 1.1. Calculated with MapleV, Rel4}

       '6529626027283988547390322592239457152074437505292804297768750294153048002',
       '1795105212861602171032500864587168232886376930667646635133196423068080917',

       '5329009504375548255615714280879284326945263678585632821548179005462985676',
       '1856076776290631611288799432795780765492309447475595172788826407376913996',

       '4642873598696137557645706926138437291825706593099940758511437502510640018',
       '1937308083199946532935467224203654472289624571736813509053891153307098606',

       '4168457886739832137920272456778622429538532861010546817035627679852581592',
       '2033032522053244646716813413301614438535637256400181883182699263429588525',

       '3810233341933401578097184188736413112865811775684941997887526231944615637',
       '2139750787591502804508697108725649529391374175244862394958291646912977339',

       '3525394856061084309896151953091650877452113131713008648419882190371053779',
       '2255115150719806812006117089372758833069336520714870417857543078989588467',

       '3291001434338053873912577671011431207650513657710743132604475957118680902',
       '2377455709116975050970501267228504382196477120855740096568995333130992810',

       '3093314217831291320318439661586246015095366463933904964346067348446792061',
       '2505536483092696199372388560344369728375248087752380979780250941594359900',

       '2923454602345568553406399345806686337519679716205002868897082161842587972',
       '2638415244414219563413106606730788636934293448096451575107905096431876360',

       '2775356879362230867616538564893568898417613047244204940481354912704001575',
       '2775356879362230867616538564893568898417613047244204940481354912704001575',

       '2644697353231021745956034972664051740369391126155904589279830224162767342',
       '2915776955732383395145087586791551003448967138978924941774017467621154981');
const
  tv2: array[0..7] of pchar8 = ({K(k), E(k), k=0, 0.25, 0.5, 0.75}
       '1.5707963267948966192313216916397514420985846996875529104874722961539082031431045',
       '1.5707963267948966192313216916397514420985846996875529104874722961539082031431045',
       '1.5962422221317835101489690714979498795055744578951225772244932732837986996771554',
       '1.5459572561054650349504124399206106120169723661630945379891425808428954404634918',
       '1.6857503548125960428712036577990769895008008941410890441199482978934337028823468',
       '1.4674622093394271554597952669909161360253617523272319605007906364908242272712906',
       '1.9109897807518291965531482187613425592531451316788338626619517081939037995321215',
       '1.3184721079946209973718427944979309306026706470626047161768627051844685382391930');

var
  f,i,cnt: longint;
  r: double;
begin
  {Test values cannot be used with precisions > 240}
  if mpf_get_default_prec>240 then exit;
  write('AGM/CK/CE test: ');
  f := 0;
  cnt := 0;
  for i:=1 to 20 do begin
    write('.');
    mp_read_decimal(a.mantissa,tva[i]);
    mpf_set_mpi2k(a,a.mantissa,-240);
    mpf_set_dbl(b,1);
    mpf_set_int(c,i);
    mpf_agm(b,c,d);
    inc(cnt);
    r := mpf_reldev(d,a);
    if r>4 then begin
      writeln;
      writeln('Failure AGM( 1,',i:2,'); r = ',r:1:1);
      inc(f);
    end;
  end;
  for i:=1 to 11 do begin
    write('.');
    mp_read_decimal(d.mantissa,tve[2*i-1]);
    mpf_set_mpi2k(d,d.mantissa,-240);
    mp_read_decimal(e.mantissa,tve[2*i]);
    mpf_set_mpi2k(e,e.mantissa,-240);
    mpf_set_int(a,i); mpf_div_d(a,10,a);
    inc(cnt,4);
    mpf_ccell1(a,c);
    mpf_ccell2(a,b);
    r := mpf_reldev(c,d);
    if r>4 then begin
      writeln;
      writeln('Failure mp_cell1(',0.1*i:1:1,'); r = ',r:1:1);
      inc(f);
    end;
    r := mpf_reldev(b,e);
    if r>4 then begin
      writeln;
      writeln('Failure mp_cell2(',0.1*i:1:1,'); r = ',r:1:1);
      inc(f);
    end;
    mpf_ccell12(a,b,c);
    r := mpf_reldev(b,d);
    if r>4 then begin
      writeln;
      writeln('Failure mp_cell112(',0.1*i:1:1,'); r(CK) = ',r:1:1);
      inc(f);
    end;
    r := mpf_reldev(c,e);
    if r>4 then begin
      writeln;
      writeln('Failure mp_cell12(',0.1*i:1:1,'); r(CE) = ',r:1:1);
      inc(f);
    end;
  end;
  for i:=0 to 3 do begin
    write('.');
    mpf_read_decimal(d,tv2[2*i]);
    mpf_read_decimal(e,tv2[2*i+1]);
    mpf_set_dbl(a,i*0.25);
    inc(cnt,2);
    mpf_cell1(a,c);
    mpf_cell2(a,b);
    r := mpf_reldev(c,d);
    if r>4 then begin
      writeln;
      writeln('Failure mp_ell1(',0.25*i:1:1,'); r = ',r:1:1);
      inc(f);
    end;
    r := mpf_reldev(b,e);
    if r>4 then begin
      writeln;
      writeln('Failure mp_ell2(',0.25*i:1:1,'); r = ',r:1:1);
      inc(f);
    end;
  end;

  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure atan2_test;
  {-Check appr. values of arctan2}
var
  f,i,j,cnt: longint;
  x,y: double;
  {$ifdef MPC_PurePascal}
    function arctan2(y, x: double): double;
    var
      z: double;
    begin
      if x=0.0 then begin
        if y=0.0 then arctan2 := 0.0
        else if y>0.0 then arctan2 := 0.5*Pi
        else arctan2 := -0.5*Pi;
      end
      else begin
        z := arctan(abs(y/x));
        if x>0 then begin
          if y<0.0 then arctan2 := -z else arctan2 := z;
        end
        else begin
          if y<0.0 then arctan2 := z-Pi else arctan2 := Pi-z;
        end;
      end;
    end;
  {$else}
    function arctan2(const y, x: double): double; assembler;
    asm
      fld     y
      fld     x
      fpatan
      fwait
    end;
  {$endif}
begin
  write('arctan2 test: ');
  f := 0;
  cnt := 0;
  for i:=-5 to 5 do begin
    for j:=-5 to 5 do begin
      if cnt and 3 = 0 then write('.');
      mpf_set_int(a,i);
      mpf_set_int(b,j);
      mpf_arctan2(b,a,c);
      y := mpf_todouble(c);
      x := arctan2(j,i);
      inc(cnt);
      if abs(x-y) > 1e-15*abs(x) then begin
        writeln;
        writeln('Failure arctan2(',j:2, ',',i:2,') = ', arctan2(j,i):9:6,'   ', mpf_decimal(c,40));
        inc(f);
      end;
    end;
  end;
  writeln;
  writeln('No. of checks: ',cnt,', failed: ',f);
  inc(totalfailed,f);
  writeln;
end;

{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

label
  done;

begin
  {mpf_set_default_prec(960);}

  writeln('Test of MP V', MP_VERSION, '   [mp_real]   (c) 2007-2018 W.Ehrhardt ');
  writeln('MAXDigits = ',MAXDigits, ',  MP_MAXBIT = ',MP_MAXBIT, ',  MAXRadix = ', MAXRadix);
  writeln('Karatsuba  cutoffs: mul/sqr = ',mp_mul_cutoff,'/',mp_sqr_cutoff);
  writeln('Toom-3, BZ cutoffs: mul/sqr = ',mp_t3m_cutoff,'/',mp_t3s_cutoff,  ',  div = ',mp_bz_cutoff);
  writeln('Current mp_float default bit precision = ', mpf_get_default_prec,
          ',  decimal precision = ', mpf_get_default_prec*ln(2)/ln(10):1:1);
  writeln;


  mpf_initp5(a,b,c,d,e,mpf_get_default_prec);
  mpf_initp3(s,t,z,mpf_get_default_prec);

  totalfailed := 0;
{
  agm_elli_test;
  func_test;
  func3_test;
  goto done;
}

  testradixconv;
  if totalfailed<>0 then begin
    writeln('Failures in testradixconv! Program aborted!');
    goto done;
  end;

  testradixconv_alt;

  mp_show_plus := true;
  mp_uppercase := true;

  add_test;
  sub_test;
  mul_test1;
  mul_test2;
  div_test;
  sqr_test;
  expt_test;

  cmp_mag_test;
  arith_test;
  atan2_test;
  agm_elli_test;
  lnagm_test;
  transz_test1;
  func_test;
  func3_test;
  test_pi;

{---------------------------------------------------------------------------}
done:
{---------------------------------------------------------------------------}

  mpf_clear5(a,b,c,d,e);
  mpf_clear3(s,t,z);

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
