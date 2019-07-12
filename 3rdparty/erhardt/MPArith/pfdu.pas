unit pfdu;

{simple unit for prime factor decomposition}

interface

{$i STD.INC}

uses
  mp_types;

{$i mp_conf.inc}

{$define UsePollardBrent}      {Use Pollard-Brent, not Pollard rho}


(*************************************************************************

 DESCRIPTION     :  simple driver unit for prime factor decomposition
                    usage:
                       pfd_initialize(ctx);
                       pfd_factor(ctx, n1);
                       pfd_factor(ctx, n2);
                       ...
                       pfd_finalize(ctx);

 REQUIREMENTS    :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA   :  (mp_types)

 MEMORY USAGE    :  heap

 DISPLAY MODE    :  textmode

 Vers  Date      Author      Modification
 ----  --------  -------     ------------------------------------------
 0.10  09.08.04  W.Ehrhardt  Initial version, incl Pollard rho and (-1)
 0.20  15.08.05  we          William's (p+1)
 0.30  17.09.05  we          Brent's ECM
 0.31  17.09.05  we          separate unit forked from t_pfd
 0.32  17.09.05  we          pfd_ctx record
 0.33  18.09.05  we          typed const rando, faster abort handling
 0.34  18.09.05  we          pfd_banner
 0.35  27.09.05  we          WIN/CRT moved to pdfu_crt
 0.36  27.09.05  we          SmallDigitBitSpecial
 0.37  01.10.05  we          SmallDigitBitSpecial up to sqrt(2^31)
 0.38  09.10.05  we          mp_small_factor with parameter fmax
 0.39  30.12.05  we          mp_count_bits/CountBits32 renamed to mp_bitsize/bitsize32
 0.40  05.08.06  we          mp_show_progress set to initial value of trace
 0.41  20.08.06  we          uses mp_max_small, check num=0 or 1, {.$define Use_Brent_ECM}
 0.42  28.08.06  we          uses mp_is_power after small primes test
 0.43  07.09.06  we          ECM1 and ECM2
 0.44  08.09.06  we          CheckWordFactor
 0.45  10.09.06  we          pfd_reset
 0.46  22.10.06  we          New pfd_ctx with exponents, sort factors, boolean ECM2
 0.47  27.10.06  we          Bugfix: mp_clear(t) in smallfactor
 0.48  03.11.06  we          CombSort included
 0.49  06.11.06  we          FPC again: @ operators needed in sortfactors
 0.50  01.05.07  we          Removed MaxLong
 0.51  11.05.07  we          Removed Use_Brent_ECM
 0.52  17.09.08  we          CheckWordFactor with mp_is_longint
 0.53  17.09.08  we          Update PFDUVers, CheckWordFactor: remove get_int
 0.54  29.12.08  we          uses mp_prime
 0.55  06.01.09  we          uses mp_prime moved to implementation
 0.56  21.01.09  we          changes related to (s)mp_divrem
 0.57  29.07.12  we          Use Shanks' SQUFOF
 0.58  15.08.12  we          mp_pollard_brent
 0.59  19.08.12  we          Split ECM2 into two phases, some parameter tuning
 0.60  14.07.14  we          Fermat factorization
 0.61  03.08.14  we          mp_holf: Hart's OneLineFactor
 0.62  13.09.16  we          Print factors/check only if trace is true
 0.63  17.07.18  we          Lim240

**************************************************************************)

(*-------------------------------------------------------------------------
 (C) Copyright 2005-2018 Wolfgang Ehrhardt

 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from
 the use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software in
    a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution.
----------------------------------------------------------------------------*)

const
  MAX_SMALL = 256;             {max number of small (mp_digit) factors}
  MAX_BIG   =  50;             {max number of large (mp_int) factors  }
  PFDUVers  = 'V0.63';         {unit version}

{$ifdef J_OPT}
var
{$else}
const
{$endif}
  trace  : boolean = true;      {trace flow of algorithms}
  rando  : boolean = true;      {call randomize in pfd_initialize}
  UseECM2: boolean = true;      {Use ECM2 phase}
  CMax   : word =    16;        {number of curves for mp_ecm_factor }
  RhoCnt : word =  1000;        {number of iteration for Pollard rho}
  PBrmax : word =  8192;        {Rmax for Pollord-Brent, 0=skip}
  FermCnt: word =   512;        {Fermat steps}
{$ifdef BIT16}
  HolfCnt: word =  5000;        {HOLF steps}
  PP1Bnd : word =  8000;        {Prime bound for Williams p+1}
  PM1Bnd : word = 10000;        {Prime bound for Pollard  p-1}
{$else}
  HolfCnt: word = 10000;        {HOLF steps}
  PP1Bnd : word = 16000;        {Prime bound for Williams p+1}
  PM1Bnd : word = 20000;        {Prime bound for Pollard  p-1}
{$endif}
  lim240 : boolean = true;      {exit is more than 240 digits}

var
  abort: boolean;               {global abort flag}

type
  pfd_ctx = record             {factorization context record}
              bfac : array[1..MAX_BIG] of mp_int;     {big factor}
              sfac : array[1..MAX_SMALL] of mp_digit; {small factors}
              sexp : array[1..MAX_SMALL] of word;     {small exponents}
              smlim: longint;  {mp_max_small_sqr after small factors}
              nextw: word;     {next word to test in CheckWordFactor}
              nb,ns: integer;  {number of big/small factors}
            end;
  PCTX    = ^pfd_ctx;          {Pointer to context, needed for sorting etc}


procedure pfd_banner;
  {-Write lib/unit banner}

procedure pfd_initialize(var ctx: pfd_ctx);
  {-Initialize factorization context}

procedure pfd_factor(var ctx: pfd_ctx; const num: mp_int);
  {-Factorize num}

procedure pfd_finalize(var ctx: pfd_ctx);
  {-Cleanup up factorization context}

procedure pfd_reset(var ctx: pfd_ctx);
  {-Reset context, but do npt clear big array}


implementation


uses
  mp_base, mp_prime, mp_numth, mp_pfu, mp_prng;


{---------------------------------------------------------------------------}
{---------------  CombSort is imported from sort.pas unit ------------------}
{---------------------------------------------------------------------------}

type
  less_funcP = function(i,j: integer; P: pointer): boolean;
                {-Compare function (pointer version), return true if item(i) < item(j)}

  swap_procP = procedure(i,j: integer; P: pointer);
                {-Swap procedure (pointer version), swaps item(i) and item(j)}


{---------------------------------------------------------------------------}
procedure CombSortP(L,R: integer; less: less_funcP; swap: swap_procP; P: pointer);
  {-General CombSort routine (pointer version), sorts items L..R}
var
  i,j,gap: integer;
  swapped: boolean;
begin
  gap := R-L;
  if gap<1 then exit;
  repeat
    gap := longint(gap)*10 div 13;
    if (gap=9) or (gap=10) then gap := 11
    else if gap<1 then gap:=1;
    swapped := false;
    for i:=L to R-gap do begin
      j := i + gap;
      if less(j,i,P) then begin
        swap(i,j,P);
        swapped := true;
      end
    end
  until (gap=1) and not swapped;
end;


{---------------------------------------------------------------------------}
procedure pfd_banner;
  {-Write lib/unit banner}
begin
  writeln('MPArith Version ', MP_VERSION, '   (c) 2004-2010 W.Ehrhardt');
  writeln('PFDU - prime factor decomposition unit ',PFDUVers);
end;


{---------------------------------------------------------------------------}
procedure pfd_initialize(var ctx: pfd_ctx);
  {-Initialize factorization context}
begin
  fillchar(ctx, sizeof(ctx), 0);
  if rando then mp_random_randomize;
  mp_init_multi(ctx.bfac);
end;


{---------------------------------------------------------------------------}
procedure pfd_finalize(var ctx: pfd_ctx);
  {-Cleanup up factorization context}
begin
  mp_clear_multi(ctx.bfac);
end;


{---------------------------------------------------------------------------}
procedure pfd_reset(var ctx: pfd_ctx);
  {-Reset context, but do npt clear big array}
begin
  with ctx do begin
    smlim := 0;
    nextw := 0;
    nb    := 0;
    ns    := 0;
  end;
end;


{---------------------------------------------------------------------------}
procedure smallfactor(var ctx: pfd_ctx; var n: mp_int);
  {-Get small factors (and smlim)}
var
  f,f0,r: mp_digit;
  t: mp_int;
begin
  if MP_Error <> MP_OKAY then exit;
  f0 := 0;
  mp_init(t);
  if MP_Error=MP_OKAY then with ctx do begin
    repeat
      mp_small_factor(n,f0,mp_max_small,f);
      if f=0 then break;
      if ns<MAX_SMALL then begin
        inc(ns);
        sfac[ns] := f;
        sexp[ns] := 0;
      end;
      mp_mod_d(n,f,r);
      if r<>0 then begin
        writeln('Internal error: non zero r in smallfactor');
        halt;
      end;
      while r=0 do begin
        mp_div_d(n,f,@t, r);
        if r=0 then begin
          inc(sexp[ns]);
          mp_exch(t,n);
        end;
        if mp_is1(n) then break;
      end;
      f0 := f;

      if trace then writeln(f:5, ' ', mp_decimal(n));
    until false;
    smlim := mp_max_small_sqr;
    if trace then writeln('SLim: ', smlim);
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
function prime(const n: mp_int): boolean;
  {-Check if n is prime}
var
  p: boolean;
begin
  if trace then write('PRIM: ',mp_decimal(n));
  p := mp_is_pprime(n);
  prime := p;
  if trace then writeln(' - ',p);
end;


{---------------------------------------------------------------------------}
procedure push(var ctx: pfd_ctx; const n: mp_int);
  {-Insert n into factor list}
begin
  with ctx do begin
    if nb<MAX_BIG then begin
      inc(nb);
      mp_copy(n,bfac[nb]);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure CheckWordFactor(var ctx: pfd_ctx; const N: mp_int; var f: mp_int);
  {-Check for word size factors; factor in f if <>0}
var
  nl: longint;
  i,i0,r: word;
const
  MaxW = $FFF2;
begin
  mp_zero(f);
  i0 := ctx.nextw;
  if mp_is_longint(N,nl) then begin
    {use mp_max_small_sqr to avoid irrelevant hints and warnings}
    if (mp_max_small_sqr<sqr(longint($7FFF))) and (i0<=46341) then begin
      {uses 15 bit integer/longint code}
      if trace then writeln('Word: ',nl);
      {Check the range MP_DIGIT_MAX+1 up to sqrt(2^31)}
      for i:=i0 to 46341 do begin
        if IsPrime16(i) then begin
          if mp_error <> MP_OKAY then exit;
          if nl mod i = 0 then begin
            mp_set_w(f,i);
            exit;
          end;
        end;
        ctx.nextw := i;
      end;
    end;
  end
  else begin
    if i0<MaxW then begin
      if trace then write('Word: ');
      for i:=i0 to $FFF2 do begin
        if trace and (i and $7FF = 0) then write('.');
        if IsPrime16(i) then begin
          if mp_error <> MP_OKAY then exit;
          mp_mod_w(N,i,r);
          if r=0 then begin
            mp_set_w(f,i);
            break;
          end;
        end;
        ctx.nextw := i;
      end;
      if trace then writeln;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function lesss(i,j: integer; P: pointer): boolean; {$ifdef BIT16} far; {$endif}
  {-Compare function for small factors}
begin
  with PCTX(P)^ do lesss := sfac[i]<sfac[j];
end;


{---------------------------------------------------------------------------}
procedure swaps(i,j: integer; P: pointer);  {$ifdef BIT16} far; {$endif}
  {-Swap procedure for small factors}
var
  tf: mp_digit;
  tx: word;
begin
  with PCTX(P)^ do begin
    tf := sfac[i];
    sfac[i] := sfac[j];
    sfac[j] := tf;
    tx := sexp[i];
    sexp[i] := sexp[j];
    sexp[j] := tx;
  end;
end;


{---------------------------------------------------------------------------}
function lessb(i,j: integer; P: pointer): boolean; {$ifdef BIT16} far; {$endif}
  {-Compare function for big factors}
begin
  with PCTX(P)^ do begin
    if (bfac[i].sign=MP_ZPOS) and (bfac[j].sign=MP_ZPOS) then lessb := mp_is_lt(bfac[i], bfac[j])
    else lessb := mp_is_gt(bfac[i], bfac[j]);
  end;
end;


{---------------------------------------------------------------------------}
procedure swapb(i,j: integer; P: pointer);  {$ifdef BIT16} far; {$endif}
  {-Swap procedure for big factors}
begin
  with PCTX(P)^ do mp_exch(bfac[i], bfac[j]);
end;


{---------------------------------------------------------------------------}
procedure sortfactors(var ctx: pfd_ctx);
  {-Sort small and big factors}
begin
  {$ifdef FPC}
    combsortp(1, ctx.ns, @lesss, @swaps, @ctx);
    combsortp(1, ctx.nb, @lessb, @swapb, @ctx);
  {$else}
    combsortp(1, ctx.ns, lesss, swaps, @ctx);
    combsortp(1, ctx.nb, lessb, swapb, @ctx);
  {$endif}
end;


{---------------------------------------------------------------------------}
procedure pfd1(var ctx: pfd_ctx; var n: mp_int);
  {-Prime factorization core (recursive)}
var
  n1,n2: mp_int;
  k,res,seed: longint;
begin
  if trace then writeln('PFD1: ',mp_decimal(n));
  if abort then exit;

  with ctx do begin
    {Test if n=1 or less than smlim}
    if mp_cmp_d(n,1)=MP_EQ then exit;
    if (mp_cmp_int(n,smlim)=MP_LT) or prime(n) then begin
      push(ctx, n);
      exit;
    end;

    mp_init2(n1,n2);
    {n=n1, d=n2}

    mp_copy(n,n1);
    if nextw=0 then nextw := succ(mp_max_small);
    CheckWordFactor(ctx,n1,n2);

    if mp_iszero(n2) and (not abort) then begin
      {test perfekt power}
      if trace then write('PPow: ');
      mp_is_power(n1,n2,k);
      if trace then begin
        if k>1 then write(k);
        writeln;
      end;
      if k>1 then begin
        if prime(n2) then begin
          while k>0 do begin
            push(ctx,n2);
            dec(k);
          end;
          exit;
        end
      end;
    end;

    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (not abort) then begin
      if mp_bitsize(n1) <= 60 then begin
        {try Shanks' SQUFOF}
        if trace then write('SQUF: ');
        mp_squfof(n1,k,res);
        if (res=0) and (k>1) then mp_set_int(n2,k);
        if trace then writeln;
      end;
    end;

    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (FermCnt>0) and (not abort) then begin
      {try Fermat factorization}
      if trace then write('Ferm: ');
      mp_fermat_factor(n1,n2, FermCnt);
      if trace then writeln;
    end;

    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (HolfCnt>0) and (not abort) then begin
      {try Fermat factorization}
      if trace then write('HOLF: ');
      mp_holf(n1,n2, HolfCnt);
      if trace then writeln;
    end;

{$ifdef UsePollardBrent}
    {Try Pollord-Brent before p-1 and p+1}
    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (PBrmax>0) and (not abort) then begin
      {try Pollard-Brent}
      if trace then write('PR-B: ');
      mp_pollard_brent_ex(n1,n2,1,PBrmax);
      if trace then writeln;
    end;
    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (PM1Bnd>0) and (not abort) then begin
      {try Pollard (p-1)}
      if trace then write('PP-1: ');
      mp_pollard_pm1(n1,n2,PM1Bnd);
      if trace then writeln;
    end;
    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (PP1Bnd>0) and (not abort) then begin
      {try Williams (p+1)}
      if trace then write('WP+1: ');
      mp_williams_pp1(n1,n2,PP1Bnd,3);
      if trace then writeln;
    end;
{$else}
    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (PM1Bnd>0) and (not abort) then begin
      {try Pollard (p-1)}
      if trace then write('PP-1: ');
      mp_pollard_pm1(n1,n2,PM1Bnd);
      if trace then writeln;
    end;
    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (PP1Bnd>0) and (not abort) then begin
      {try Williams (p+1)}
      if trace then write('WP+1: ');
      mp_williams_pp1(n1,n2,PP1Bnd,3);
      if trace then writeln;
    end;
    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (RhoCnt>0) and (not abort) then begin
      {try Pollard rho}
      if trace then write('PRho: ');
      mp_pollard_rho(n1,n2,RhoCnt);
      if trace then writeln;
    end;
{$endif}

    if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (CMax>3) and (not abort) then begin
      {try ECM}
      if trace then write('ECM1: ');
      seed := 0;
      mp_ecm_factor(n1,n2, 3,400,seed,k);
      if trace then writeln;
    end;
    if UseECM2 then begin
      if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (CMax>0) and (not abort) then begin
        {try ECM}
        if trace then write('ECM2: ');
        seed := 0;
        mp_ecm_factor(n1,n2,6,3000,seed,k);
        if trace then writeln;
      end;
      if ((mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2)) and (CMax>0) and (not abort) then begin
        {try ECM}
        if trace then write('ECM3: ');
        seed := 0;
        mp_ecm_factor(n1,n2,CMax,ECM_C1Max,seed,k);
        if trace then writeln;
      end;
    end;

    if (mp_cmp(n1, n2)=MP_EQ) or mp_iszero(n2) or abort then begin
      {no factor found, insert negative n in to list}
      mp_chs(n1, n1);
      push(ctx, n1);
    end
    else begin
      {found factor d=n2, factorize d and n/d}
      mp_div(n1,n2,n1);
      pfd1(ctx, n2);
      pfd1(ctx, n1);
    end;
    mp_clear2(n1,n2);
  end;
end;


{---------------------------------------------------------------------------}
procedure pfd_factor(var ctx: pfd_ctx; const num: mp_int);
  {-Factorize num}
var
  i: integer;
  n,x: mp_int;
begin
  if mp_not_init(num) then begin
    writeln('Input number not initialized');
    exit;
  end;
  if lim240 and (mp_radix_size(num,10)>240)  then begin
    if trace then writeln('Input number with more that 240 decimal digits');
    exit;
  end;
  abort := false;
  with ctx do begin
    nb := 0;
    ns := 0;
    mp_init2(n,x);
    if mp_error<>0 then begin
      writeln('Error initializing copy of input');
      exit;
    end;
    mp_abs(num,n);
    if mp_cmp_d(n,1)=MP_GT then begin
      {Get small (mp_digit) factors}
      smallfactor(ctx, n);
      if mp_cmp_d(n,1)=MP_GT then begin
        {find large factor}
        pfd1(ctx, n);
      end;
      sortfactors(ctx);
      if trace then begin
        {Print and multiply factors}
        mp_set(n,1);
        for i:=1 to ns do begin
          if sexp[i]=0 then sexp[i]:=1;
          mp_set_pow(x,sfac[i],sexp[i]);
          mp_mul(n, x, n);
          if sexp[i]=1 then write(sfac[i])
          else write(sfac[i],'^',sexp[i]);
          if i<ns then write(' * ');
        end;
        for i:=1 to nb do begin
          if (i>1) or (ns>0) then write(' * ');
          write(mp_decimal(bfac[i]));
          mp_mul(n, bfac[i], n);
        end;
      end;
    end;
    if trace then begin
      writeln;
      writeln('Start: ',mp_decimal(num));
      writeln('Check: ',mp_decimal(n));
    end;
  end;
  mp_clear2(n,x);
end;


begin
  mp_show_progress := trace;
end.
