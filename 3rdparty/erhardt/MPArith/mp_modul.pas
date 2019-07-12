unit mp_modul;

{MP basic modular arithmetic, GCD and LCM, Jacobi/Kronecker}

interface

{$i STD.INC}

uses
  mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION   :  MP basic modular arithmetic, GCD and LCM, Jacobi/Kronecker

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  heap

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] LibTomMath 0.30+ by Tom St Denis
                  [2] MPI by M.J. Fromberger
                  [3] Knuth, D.E.: The Art of computer programming. Vol 2
                      Seminumerical Algorithms, 3rd ed., 1998
                  [5] (HAC) Menezes,A., von Oorschot,P., Vanstone, S: Handbook of
                      Applied Cryptography, 1996, www.cacr.math.uwaterloo.ca/hac
                 [12] PARI/GP at http://pari.math.u-bordeaux.fr/
                 [23] J. Shallit, J. Sorenson, A Binary Algorithm for the Jacobi Symbol,
                      SIGSAM Bulletin, 27(1), 4-11, 1993; available online at
                      http://blue.butler.edu/~jsorenso/papers/binjac.ps or
                      http://citeseer.ist.psu.edu/article/shallit93binary.html
                 [24] H. Cohen, A Course in Computational Algebraic Number Theory
                      4th printing, 2000
                 [28] J.P. Sorenson, An analysis of Lehmer's Euclidean GCD algorithm,
                      ACM International Symposium on Symbolic and Algebraic Computation, 1995
                      http://blue.butler.edu/~jsorenso/papers/lehmer.pdf


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.23.00  24.09.12  W.Ehrhardt  Split from old mp_numth
 1.23.01  24.10.12  we          mp_initial_mod initialized (removed from mp_numth)

 1.24.00  15.12.12  we          Some word/integer types changed to TNInt

**************************************************************************)

(*-------------------------------------------------------------------------
  This code uses material/ideas from the following 3rd party libraries:
   - LibTomMath 0.30+ by Tom St Denis
   - MPI 1.8.6 by Michael J. Fromberger
   - NX V0.18 and V0.9+ by Marcel Martin
  See the file '3rdparty.mpa' for the licenses.
----------------------------------------------------------------------------*)


(*-------------------------------------------------------------------------
 (C) Copyright 2004-2012 Wolfgang Ehrhardt

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


var
  mp_initial_mod  : boolean;  {Single initial modular reduction in mp_gcd}
                              {and in mp_kronecker/mp_jacobi}


procedure mp_addmod(const a,b,c: mp_int; var d: mp_int);
  {-Calculate d = a + b (mod c)}

procedure mp_exptmod(const a,b,c: mp_int; var d: mp_int);
  {-Compute d = a^b mod c, c>0. If b<0, a must have an inverse mod c}

procedure mp_exptmod_d(const a: mp_int; b: longint; const c: mp_int; var d: mp_int);
  {-Compute d = a^b mod c, c>0, b longint. If b<0, a must have an inverse mod c}

procedure mp_gcd(const a,b: mp_int; var c: mp_int);
  {-Calculate c = gcd(a,b) using the binary method}

function  mp_gcd1(const a,b: mp_int; var c: mp_int): boolean;
  {-Calculate c = gcd(a,b) using the binary method, return true if c=1 and no error}

procedure mp_gcd_euclid(const a,b: mp_int; var c: mp_int);
  {-Calculate c = gcd(a,b), non optimized Euclid}

procedure mp_gcd_ml(const a,b: mp_int; var u: mp_int);
  {-Calculate u = gcd(a,b) using the Sorenson's Modified Lehmer method}

procedure mp_invmod(const a, b: mp_int; var c: mp_int);
  {-Compute c = a^-1 (mod b), b>0, via mp_xgcd, MP_UNDEF error if there is no inverse}

function  mp_invmodf(const a, b: mp_int; var c: mp_int): boolean;
  {-Compute c = a^-1 (mod b), b>0, via mp_xgcd, return true if inverse exists}

function  mp_jacobi(const a,n: mp_int): integer;
  {-Compute the Jacobi/Legendre symbol (a|n), n: odd and > 2}

function  mp_jacobi_lm(a: longint; const n: mp_int): longint;
  {-Compute the Jacobi/Legendre symbol (a|n), n: odd and > 2}

function  mp_jacobi_ml(const a: mp_int; n: longint): longint;
  {-Compute the Jacobi/Legendre symbol (a|n), n: odd and > 2}

function  mp_kronecker(const a,n: mp_int): integer;
  {-Compute the Kronecker symbol (a|n)}

procedure mp_lcm(const a,b: mp_int; var c: mp_int);
  {-Compute least common multiple as |a*b|/(a, b)}

procedure mp_mulmod(const a,b,c: mp_int; var d: mp_int);
  {-Calculate d = a * b (mod c)}

procedure mp_sqrmod(const a,c: mp_int; var d: mp_int);
  {-Calculate d = a * a (mod c)}

procedure mp_submod(const a,b,c: mp_int; var d: mp_int);
  {-Calculate d = a - b (mod c)}

procedure mp_xgcd(const a,b: mp_int; p1,p2,p3: pmp_int);
  {-Extended gcd algorithm, calculate  a*p1^ + b*p2^ = p3^ = gcd(a,b)}
  { p1,p2,p3 may be nil if the values are not required}

procedure mp_xgcd_bin(const a,b: mp_int; p1,p2,p3: pmp_int);
  {-Extended binary gcd, calculate a*p1^ + b*p2^ = p3^ = gcd(a,b)  }
  { p1,p2,p3 may be nil if the values are not needed. Note that p1^}
  { and p2^ are NOT unique and may differ from those of mp_xgcd!!  }

procedure mp_xlcm(const a,b: mp_int; var c,x,y: mp_int);
  {-Calculate c,x,y with lcm(a,b)=c=x*y and x|a, y|b, gcd(x,y)=1.}
  { c,x,y should be different variables, otherwise result may be inconsistent}


implementation

uses
  mp_base;

type
  TRedType = (MR_Barret, MR_Montgomery
             {$ifdef MPC_Reduce_2k} ,MR_Reduce2k {$endif}); {supported reduction types for exptmod}


{---------------------------------------------------------------------------}
procedure mp_addmod(const a,b,c: mp_int; var d: mp_int);
  {-Calculate d = a + b (mod c)}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) or mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_addmod');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init(t);
  if mp_error<>MP_OKAY then exit;

  mp_add(a, b, t);
  mp_mod(t, c, d);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_exptmod_win(const g,e,p: mp_int; var b: mp_int; redmode: TRedType);
  {-Internal: Compute y=g^|e| mod p, p>0, internal sliding windows}
label
  __M, __MU, __RES;
var
  bitbuf, bitcpy, bitcnt, mode, x, y, winsize, wmax1,wmax2: integer;
  digidx: TNint;
  ps2: longint;
  buf: mp_digit;
  mp : mp_digit;
  {$ifdef MPC_Reduce_2k}
  d2k: mp_digit;
  {$endif}
  bc: longint;
  res, mu: mp_int;
  M: array[1..WEXP_TABSIZE] of mp_int;

  {---------------------------------------------}
  procedure Gen_Redux(var mpi: mp_int);
    {-General modular reduction of mpi driven by redmode}
  begin
    case redmode of
      MR_Montgomery: mp_montgomery_reduce(mpi, p, mp);
          MR_Barret: mp_reduce(mpi, p, mu);
    {$ifdef MPC_Reduce_2k}
        MR_Reduce2k: begin
                       mp_reduce_2k(mpi, p, d2k);
                     end;
    {$endif}
               else  if mp_error=MP_OKAY then set_mp_error(MP_BADARG);
    end;
  end;

begin
  {Uses a left-to-right k-ary sliding window to compute the modular exponentiation}
  {The value of k changes based on the size of the exponent.}
  if mp_error<>MP_OKAY then exit;
  {No checks}
  {find window size}
  bc := mp_bitsize(e);
  if bc<=7 then winsize := 2
  else if bc<=36 then winsize := 3
  else if bc<=140 then winsize := 4
  else if bc<=450 then winsize := 5
  else if bc<=1303 then winsize := 6
  else if bc<=3529 then winsize := 7
  else winsize := 8;

  if winsize>WEXP_MAX then winsize := WEXP_MAX;
  wmax1 := 1 shl (winsize-1);
  wmax2 := 1 shl winsize;

  ps2 := 2*(p.used+1);
  {initialize M array}
  {initialize first cell}
  mp_init_size(M[1],ps2);
  if mp_error<>MP_OKAY then exit;

  {now initialize the second half of the array}
  for x:=wmax1 to wmax2-1 do begin
    mp_init_size(M[x],ps2);
    if mp_error<>MP_OKAY then begin
      for y:=wmax2 to x-1 do mp_clear(M[y]);
      mp_clear(M[1]);
      exit;
    end;
  end;

  {create mu, used for Barrett reduction}
  mp_init_size(mu,ps2);
  if mp_error<>MP_OKAY then goto __M;

  {setup result}
  mp_init(res);
  if mp_error<>MP_OKAY then goto __MU;

  {Do initial setup depending on reduction type}
  if Redmode=MR_Montgomery then begin
    mp_montgomery_setup(p, mp);
    mp_montgomery_calcnorm(res, p);
    {now set M[1] to G * R mod m}
    mp_mulmod(g, res, p, M[1]);
  end
  else begin
    {The M table contains powers of the base, }
    {e.g. M[x] = g^x mod p                   }
    mp_set(res, 1);
    mp_mod(g, p, M[1]);
    if Redmode=MR_Barret then mp_reduce_setup(mu, p)
  {$ifdef MPC_Reduce_2k}
    else if Redmode=MR_Reduce2k then begin
      mp_reduce_2k_setup(p, d2k)
    end
  {$endif}
    else begin
      {unsupported reduction}
      set_mp_error(MP_BADARG);
    end;
  end;
  if mp_error<>MP_OKAY then goto __RES;

  {Create M table                           }
  {The first half of the table is not       }
  {computed though accept for M[0] and M[1] }

  {compute the value at M[wmax1] by squaring M[1] (winsize-1) times}
  mp_copy(M[1], M[wmax1]);
  for x:=0 to winsize-2 do begin
    mp_sqr(M[wmax1], M[wmax1]);
    Gen_Redux(M[wmax1]);
    if mp_error<>MP_OKAY then goto __RES;
  end;

  {create upper table, that is M[x] = M[x-1] * M[1] (mod p)}
  for x:=wmax1+1 to wmax2-1 do begin
    mp_mul(M[x-1], M[1], M[x]);
    Gen_Redux(M[x]);
    if mp_error<>MP_OKAY then goto __RES;
  end;

  {set initial mode and bit cnt}
  mode   := 0;
  bitcnt := 1;
  buf    := 0;
  digidx := e.used - 1;
  bitcpy := 0;
  bitbuf := 0;

  repeat
    if mp_error<>MP_OKAY then goto __RES;
    {grab next digit as required }
    dec(bitcnt);
    if bitcnt=0 then begin
      {if digidx = -1 we are out of digits}
      if digidx = -1 then break;
      {read next digit and reset the bitcnt}
      buf    := e.pdigits^[digidx]; dec(digidx);
      bitcnt := DIGIT_BIT;
    end;

    {grab the next msb from the exponent}
    y := (buf shr mp_digit(DIGIT_BIT - 1)) and 1;

    {temporarly turn range checks off}
    {$R-}
    buf := buf shl 1;
    {$ifdef RangeChecks_on} {$R+} {$endif}

    {if the bit is zero and mode = 0 then we ignore it             }
    {These represent the leading zero bits before the first 1 bit  }
    {in the exponent.  Technically this opt is not required but it }
    {does lower the # of trivial squaring/reductions used          }

    if (mode=0) and (y=0) then continue;

    {if the bit is zero and mode == 1 then we square}
    if (mode=1) and (y=0) then begin
      mp_sqr(res, res);
      Gen_Redux(res);
      if mp_error<>MP_OKAY then goto __RES;
      continue;
    end;

    {else we add it to the window}
    inc(bitcpy);
    bitbuf := bitbuf or (y shl (winsize - bitcpy));
    mode   := 2;

    if bitcpy=winsize then begin
      {ok window is filled so square as required and multiply}
      {square first}
      for x:=0 to winsize-1 do begin
        mp_sqr(res, res);
        Gen_Redux(res);
        if mp_error<>MP_OKAY then goto __RES;
      end;

      {then multiply}
      mp_mul(res, M[bitbuf], res);
      {and reduce}
      Gen_Redux(res);
      if mp_error<>MP_OKAY then goto __RES;

      {empty window and reset}
      bitcpy := 0;
      bitbuf := 0;
      mode   := 1;
    end;
  until false;

  {if bits remain then square/multiply}
  if (mode=2) and (bitcpy > 0) then begin
    {square then multiply if the bit is set}
    for x:=0 to bitcpy-1 do begin
      mp_sqr(res, res);
      Gen_Redux(res);
      if mp_error<>MP_OKAY then goto __RES;
      bitbuf := bitbuf shl 1;
      if (bitbuf and wmax2) <> 0 then begin
        {then multiply}
        mp_mul(res, M[1], res);
        {and reduce}
        Gen_Redux(res);
        if mp_error<>MP_OKAY then goto __RES;
      end;
    end;
  end;

  if Redmode=MR_Montgomery then begin
    {fix result if Montgomery reduction is used recall that any value}
    {in a Montgomery system is actually multiplied by R mod n.  So we}
    {have to reduce one more time to cancel out the factor of R.     }
    mp_montgomery_reduce(res, P, mp);
    if mp_error<>MP_OKAY then goto __RES;
  end;

  mp_exch(res, b);

__RES: mp_clear(res);
 __MU: mp_clear(mu);
  __M: mp_clear(M[1]);

  for x:=wmax1 to wmax2-1 do mp_clear(M[x]);

end;


{---------------------------------------------------------------------------}
procedure mp_exptmod(const a,b,c: mp_int; var d: mp_int);
  {-Compute d = a^b mod c, c>0. If b<0, a must have an inverse mod c}
var
  rt: TRedType;
  t: array[0..1] of mp_int; {a:t[0], b:t[1]}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) or mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_exptmod');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {modulus c must be positive}
  if s_mp_is_le0(c) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_exptmod: c<=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {if exponent b is negative we have to recourse}
  if b.sign=MP_NEG then begin
     {first compute 1/a mod c}
     mp_init_multi(t);
     if mp_error=MP_OKAY then begin
       mp_invmod(a, c, t[0]);
       {now get |b|}
       mp_abs(b, t[1]);
       {and now compute (1/a)^|b| instead of a^b [b < 0]}
       mp_exptmod(t[0], t[1], c, d);
       mp_clear_multi(t);
     end;
     exit;
  end;

  {easy outs}
  if mp_is1(c) then begin
    mp_zero(d);
    exit;
  end;
  if mp_is1(b) then begin
    mp_mod(a,c,d);
    exit;
  end;

  {Default: Barrett reduction}
  rt := MR_Barret;
  if mp_isodd(c) then rt := MR_Montgomery;

{$ifdef MPC_Reduce_2k}
  if mp_reduce_is_2k(c) then rt := MR_Reduce2k;
  {*tbd: DR module variants}
{$endif}

  {Use sliding window routine to compute the modular exponentiation}
  mp_exptmod_win(a, b, c, d, rt)

end;


{---------------------------------------------------------------------------}
procedure mp_exptmod_d(const a: mp_int; b: longint; const c: mp_int; var d: mp_int);
  {-Compute d = a^b mod c, c>0, b longint. If b<0, a must have an inverse mod c}
var
  tmp: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {arg check in mp_exptmod}
  mp_init_set_int(tmp, b);
  if mp_error=MP_OKAY then begin
    mp_exptmod(a,tmp,c,d);
    mp_clear(tmp);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_gcd(const a,b: mp_int; var c: mp_int);
  {-Calculate c = gcd(a,b) using the binary method}
var
  u,v: mp_int;
  k,u_lsb,v_lsb: longint;
{$ifdef MPC_UseGCD32}
  ui,vi: longint;
{$endif}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_gcd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if one arg is zero, gcd is the abs of the other}
  if mp_iszero(a) then begin
    mp_abs(b, c);
    exit;
  end
  else if mp_iszero(b) then begin
    mp_abs(a, c);
    exit;
  end;

  {get copies of a and b we can modify}
  mp_init2(u,v);
  if mp_error<>MP_OKAY then exit;

  {u,v are positive for the remainder of the algorithm}
  mp_abs(a,u);
  mp_abs(b,v);

  if mp_error=MP_OKAY then begin
    {Do a single initial modular reduction if the sizes of u and v}
    {are very different, see H.Cohen's Alg. 1.3.5 in [24]}
    if mp_initial_mod then begin
      {u>0, v>0. Make u >= v}
      if (u.used < v.used) or (mp_cmp_mag(u,v)=MP_LT) then mp_exch(u,v);
      if u.used > v.used+1 then begin
        {do the initial reduction}
        mp_mod(u,v,u);
        {done in u mod v = 0}
        if u.used=0 then begin
          {store result and clear local vars}
          mp_exch(c,v);
          mp_clear2(u,v);
          exit;
        end;
        {$ifdef MPC_UseGCD32}
          if mp_is_longint(u,ui) then begin
            {u <>0 but fits into longint}
            mp_mod_int(v,ui,vi);
            mp_set_int(c,gcd32(ui,vi));
            mp_clear2(u,v);
            exit;
          end;
        {$endif}
      end;
    end;

    {Here both u and v are > 0; find the common power of two for u and v}
    u_lsb := mp_cnt_lsb(u);
    v_lsb := mp_cnt_lsb(v);
    if u_lsb<v_lsb then k:=u_lsb else k:=v_lsb;

    {Make both u and v odd}
    if u_lsb>0 then mp_shr(u, u_lsb, u);
    if v_lsb>0 then mp_shr(v, v_lsb, v);

    while (mp_error=MP_OKAY) and (v.used>0) do begin
      {$ifdef MPC_UseGCD32}
        if mp_is_longint(v,vi) then begin
          mp_mod_int(u,vi,ui);
          mp_set_int(u,gcd32(ui,vi));
          break;
        end;
      {$endif}
      {done if u=v, if u>v swap u and v}
      if u.used>=v.used then begin
         case mp_cmp_mag(u,v) of
           MP_GT: mp_exch(u, v);
           MP_EQ: break;
         end;
      end;
      {subtract smaller from larger, resulting v is > 0}
      s_mp_sub(v, u, v);
      if v.pdigits^[0] and 2 <>0 then begin
        {Only one factor two can be removed, so reduce overhead. This}
        {occurs about one-half of the time, see Knuth [3] 4.5.2, p348}
        mp_shr1(v);
      end
      else begin
        {Divide out all factors of two}
        mp_shr(v, mp_cnt_lsb(v), v);
      end;
    end;

    {multiply by 2^k which we divided out at the beginning}
    mp_shl(u, k, c);
  end;
  mp_clear2(u,v);
end;


{---------------------------------------------------------------------------}
function mp_gcd1(const a,b: mp_int; var c: mp_int): boolean;
  {-Calculate c = gcd(a,b) using the binary method, return true if c=1 and no error}
begin
  mp_gcd(a,b,c);
  mp_gcd1 := mp_is1(c);
end;


{---------------------------------------------------------------------------}
procedure mp_gcd_euclid(const a,b: mp_int; var c: mp_int);
  {-Calculate c = gcd(a,b), non optimized Euclid}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_gcd_euclid');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Simple Euclid algorithm, Knuth A}
  mp_init(t);
  if mp_error=MP_OKAY then begin
    {make both mp_ints positive}
    {handle b first (may be @b=@c)}
    mp_abs(b,t);
    mp_abs(a,c);
    {ggc(c,t) = gcd(t, c mod t)}
    while (mp_error=MP_OKAY) and (not mp_iszero(t)) do begin
      mp_mod(c,t,c);
      mp_exch(t,c);
    end;
    mp_clear(t);
  end;
end;


{$ifdef MP_32BIT}
type
  tml_int = longint; {FPC has 16 bit integer/maxint with 32 bit code!!!!!}
  tml_dbl = int64;
const
  Lehmer_MaxSingle = MP_DIGIT_MAX;        {use with mul_d()}
  Lehmer_MaxK      = 2*DIGIT_BIT-1;
{$else}
type
  tml_int = integer;
  tml_dbl = longint;
const
  Lehmer_MaxSingle = $7FFF;               {use with mul_w()}
  Lehmer_MaxK      = 29;
{$endif}


{---------------------------------------------------------------------------}
function Lehmer(const u,v: mp_int; k: integer; var x0,x1,y0,y1: tml_int): boolean;
  {-Single precision Lehmer steps to reduce u,v based on highest k bits of}
  { u,v. Return true if successful and reduction using xi,yi should be done}
var
  hu,hv,hb: mp_word;
  i: TNInt;
  uh,vh,q,x2,y2: tml_dbl;

begin

  {This is an implementation of Sorenson's Modified Lehmer procedure [28]. }
  {It is based on results of Lehmer, Collins, Jebelean, and Sorenson. The  }
  {procedure performs single digit calculations that allow to combine some }
  {Euclidean GCD steps into a single multiprecision calculation. It returns}
  {Sorenson's x[i], y[i], x_[i-1], and y_[i-1] in x1, y1, x0, and y0. These}
  {values are used in both GCD procedures (standard and extended).         }

  Lehmer := false;

  {Check the technical conditions and exit if they are not fulfilled}
  if (k>Lehmer_MaxK) or (u.sign=MP_NEG) or (v.sign=MP_NEG) or (u.used<2) or (v.used>u.used) or (u.used>v.used+1) then exit;

  {Get the leading two digits of u and v; here v[u.used-1] may be zero.}
  i  := pred(u.used);
  hu := (mp_word(u.pdigits^[i]) shl DIGIT_BIT) + u.pdigits^[pred(i)];
  if v.used<u.used then hv := v.pdigits^[pred(i)]
  else hv := (mp_word(v.pdigits^[i]) shl DIGIT_BIT) + v.pdigits^[pred(i)];

  {Get the highest k bits of u and the corresponding bits of v}
  hb := mp_word(1) shl k;
  i  :=0;
  while hu>=hb do begin
    hu := hu shr 1;
    inc(i);
  end;
  hv := hv shr i;
  if hv=0 then exit;

  i  := 0;
  uh := hu;
  vh := hv;
  x0 := 1;
  y0 := 0;
  x1 := 0;
  y1 := 1;
  repeat
    q  := uh div vh;
    x2 := vh;
    vh := uh - q*vh;
    uh := x2;
    x2 := x0 - q*x1;
    y2 := y0 - q*y1;
    {In addition to the standard break conditions, exit if the next  }
    {new x2/y2 values are greater than the maximum allowed values. If}
    {this happens during the first iteration, Lehmer is still false  }
    if (abs(x2) > Lehmer_MaxSingle) or (abs(y2) > Lehmer_MaxSingle) then break;
    inc(i);
    if odd(i) then begin
      if (vh < -y2) or (uh-vh < x2-x1) then break;
    end
    else begin
      if (vh < -x2) or (uh-vh < y2-y1) then break;
    end;
    x0 := x1;
    x1 := x2;
    y0 := y1;
    y1 := y2;
    Lehmer := true;
  until false;
end;


{---------------------------------------------------------------------------}
procedure mp_gcd_ml(const a,b: mp_int; var u: mp_int);
  {-Calculate u = gcd(a,b) using the Sorenson's Modified Lehmer method}
var
  s,t,v: mp_int;
  tu,tv: longint;
  bsu,bsv: longint;
  x0,x1,y0,y1: tml_int;
const
  k=Lehmer_MaxK; k2=(k+1) div 2;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(u) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_gcd_ml');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init3(s,t,v);
  if mp_error=MP_OKAY then begin
    {make both mp_ints positive}
    {handle b first (may be @b=@u)}
    mp_abs(b,v);
    mp_abs(a,u);
    {make u>=v, this will be a loop invariant}
    if (u.used<v.used) or (mp_cmp_mag(u,v)=MP_LT) then mp_exch(u,v);

    bsu := mp_bitsize(u);
    bsv := mp_bitsize(v);
    {Iterate until v is in longint range, then use gcd32}
    while (mp_error=MP_OKAY) and (bsv>31) do begin
      if (bsu-bsv <= k2) and Lehmer(u,v,k,x0,x1,y0,y1) then begin
        {calculate new (u, v):  v := x1*u + y1*v;  u := x0*u + y0*v;}
        {$ifdef MP_16BIT}
           {v := x1*u + y1*v}
           mp_mul_w(u,abs(x1),s);  if x1<0 then s_mp_chs(s);
           mp_mul_w(v,abs(y1),t);  if y1<0 then s_mp_chs(t);
           mp_add(s,t,t);
           mp_exch(t,v);
           {u := x0*u + y0*v}
           mp_mul_w(u,abs(x0),s);  if x0<0 then s_mp_chs(s);
           mp_mul_w(t,abs(y0),t);  if y0<0 then s_mp_chs(t);
           mp_add(s,t,u);
        {$else}
           {v := x1*u + y1*v}
           mp_mul_d(u,abs(x1),s);  if x1<0 then s_mp_chs(s);
           mp_mul_d(v,abs(y1),t);  if y1<0 then s_mp_chs(t);
           mp_add(s,t,t);
           mp_exch(t,v);
           {u := x0*u + y0*v}
           mp_mul_d(u,abs(x0),s);  if x0<0 then s_mp_chs(s);
           mp_mul_d(t,abs(y0),t);  if y0<0 then s_mp_chs(t);
           mp_add(s,t,u);
        {$endif}
      end;
      mp_mod(u,v,u);
      mp_exch(v,u);
      bsu := bsv;
      bsv := mp_bitsize(v);
    end;
    if (mp_error=MP_OKAY) and (v.used>0) then begin
      {v<>0 and fits into longint}
      tv := mp_get_int(v);
      mp_mod_int(u,tv,tu);
      mp_set_int(u,gcd32(tv,tu));
    end;
    mp_clear3(s,t,v);
  end;
end;


{---------------------------------------------------------------------------}
function mp_invmodf(const a, b: mp_int; var c: mp_int): boolean;
  {-Compute c = a^-1 (mod b), b>0, via mp_xgcd, return true if inverse exists}
var
  u,v: mp_int;
  la, lb, lc: longint;
begin

  {Note: this function is normally faster than the old fast_invmod}
  {for odd b due to the added Lehmer steps in mp_xgcd. It is much }
  {faster than the old mp_invmodf for even b.}

  mp_invmodf := false;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_invmodf');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if a and b are even or b<=0 then return no success}
  if (mp_iseven(b) and mp_iseven(a)) or (b.sign=MP_NEG) or mp_iszero(b) then exit;

  {use invmod32 if b has less than 32 bits}
  if mp_is_longint(b,lb) then begin
    if not mp_is_longint(a,la) then mp_mod_int(a,lb,la);
    lc := invmod32(la,lb);
    if lc<>0 then begin
      mp_set_int(c,lc);
      mp_invmodf := true;
    end;
    exit;
  end;

  {initialize temps}
  mp_init2(u,v);
  if mp_error<>MP_OKAY then exit;

  {calculate a*u + b*? = v = gcd(a,b)}
  mp_xgcd(a,b,@u,nil,@v);

  {if gcd(a,b><>1 then there is no inverse}
  if mp_is1(v) then begin
    {Make 0 <= a^-1 < b}
    while (mp_error=MP_OKAY) and (u.sign=MP_NEG) do mp_add(u,b,u);
    while (mp_error=MP_OKAY) and (mp_cmp_mag(u, b)<>MP_LT) do mp_sub(u,b,u);
    mp_exch(u, c);
    mp_invmodf := mp_error=MP_OKAY;
  end;
  mp_clear2(u,v);
end;


{---------------------------------------------------------------------------}
procedure mp_invmod(const a, b: mp_int; var c: mp_int);
  {-Compute c = a^-1 (mod b), b>0, via mp_xgcd, MP_UNDEF error if there is no inverse}
begin
  if not mp_invmodf(a,b,c) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXUndef.Create('mp_invmod: no inverse');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$endif}
    set_mp_error(MP_UNDEF);
  end;
end;


{---------------------------------------------------------------------------}
function kron_intern(var x,y: mp_int): integer;
  {-Internal kronecker: no mp_init, initial range/sanity checks etc}
var
  res: integer;
  m8: mp_digit;
  k: longint;
begin
  (*
  Here the PARI/GP definition is used: The Kronecker symbol is the Legendre
  symbol (x|y) extended to all integers by complete multiplicativity,
  plus special rules for y = 0, -1, or 2:

  y = 0:  (x|0)  =  1 if |x| = 1
                    0 otherwise

  y =-1: (x|-1)  =  1 if x >= 0,
                   -1 if x <  0.

  y = 2:  (x|2)  =  0 if x is even
                 =  1 if x = 1, 7 mod 8
                 = -1 if x = 3, 5 mod 8
  *)

  {initialize return val for zero/error}
  kron_intern := 0;

  {easy case (x|0)}
  if mp_iszero(y) then begin
    if mp_is1a(x) then kron_intern := 1;
    exit;
  end;

  {initialize accumulated result}
  res := 1;

  {here y<>0, make y positive}
  if y.sign=MP_NEG then begin
    {(x|y) = (x|-1)*(x|-y)}
    if x.sign=MP_NEG then res := -res;
    mp_abs(y,y);
  end;

  {if y even, reduce to odd case}
  if y.pdigits^[0] and 1 = 0 then begin
    {(x|2)=0 if x is even}
    if mp_iseven(x) then exit;
    {y is even; divide out highest power of two: y = 2^k*y'}
    {(x|y) = (x|2)^k * (x|y')}
    mp_makeodd(y,y,k);
    if odd(k) then begin
      {note: x is odd and <> 0, so pdigits^[0] is valid}
      m8 := x.pdigits^[0] and 7;
      if (m8=3) or (m8=5) then res := -res;
    end;
  end;

  {Here y is positive and odd, and we actually calculate a Jacobi symbol.}
  {x>=0 is needed for binary algorithm, so take abs and adjust sign.}
  if x.sign=MP_NEG then begin
    {(-x|y) = (x|y)*(-1|y); change sign if y=3 (mod 4)}
    if y.pdigits^[0] and 3 = 3 then res := -res;
    x.sign := MP_ZPOS;
  end;

  if mp_initial_mod and (y.used<=x.used) then begin
    {Here x,y>0, y is odd. Reduce size once, see Cohen [24], Alg.1.4.12}
    mp_mod(x,y,x);
  end;

  {This is the Jacobi loop from Shallit/Sorenson [23]}
  while (x.used>0) and (mp_error=MP_OKAY) do begin
    {y odd, x odd or even. divide out highest power of two: 2^k}
    mp_makeodd(x,x,k);
    {only if odd exponent check y mod 8 and change sign}
    if odd(k) then begin
      { (2|y) = -1 if y = +-3 (mod 8)}
      m8 := y.pdigits^[0] and 7;
      if (m8=3) or (m8=5) then res := -res;
    end;
    {x,y odd}
    if x.used<=y.used then begin
      case mp_cmp_mag(x,y) of
        MP_LT: begin
                 {if x<y then interchange(x,y) and adjust res}
                 mp_exch(x,y);
                 if (y.pdigits^[0] and 3 = 3) and (x.pdigits^[0] and 3 = 3) then res := -res;
               end;
        MP_EQ: begin
                 {done if x=y}
                 break;
               end;
      end;
    end;
    {Here x,y odd and x>y so we can use s_mp_sub}
    s_mp_sub(x,y,x);
    {y odd, x even, x > 0.  Division by 2 and sign adjustment from the}
    {Shallit/Sorenson code is not needed here because these operations}
    {will be done in the next pass of the loop.}
  end;
  if (mp_error=MP_OKAY) and mp_is1(y) then kron_intern := res;
end;


{---------------------------------------------------------------------------}
function mp_kronecker(const a,n: mp_int): integer;
  {-Compute the Kronecker symbol (a|n)}
var
  a2, n2: mp_int;
begin
  mp_kronecker := 0;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_kronecker');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Initialize local mp_int copies}
  mp_init_copy(a2,a);
  if mp_error<>MP_OKAY then exit;

  mp_init_copy(n2, n);
  if mp_error=MP_OKAY then begin
    mp_kronecker := kron_intern(a2,n2);
    mp_clear(n2);
  end;
  mp_clear(a2);
end;


{---------------------------------------------------------------------------}
function mp_jacobi(const a,n: mp_int): integer;
  {-Compute the Jacobi/Legendre symbol (a|n), n: odd and > 2}
var
  a2, n2: mp_int;
begin
  mp_jacobi := 0;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_jacobi');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Error if n<3 or even}
  if (n.sign=MP_NEG) or mp_is1(n) or mp_iseven(n) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_jacobi: n<3 or even');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {Initialize local mp_int copies}
  mp_init_copy(a2,a);
  if mp_error<>MP_OKAY then exit;

  mp_init_copy(n2, n);
  if mp_error=MP_OKAY then begin
    {Invalid n values are excluded above, use Kronecker symbol}
    mp_jacobi := kron_intern(a2,n2);;
    mp_clear(n2);
  end;
  mp_clear(a2);
end;


{---------------------------------------------------------------------------}
function mp_jacobi_lm(a: longint; const n: mp_int): longint;
  {-Compute the Jacobi/Legendre symbol (a|n), n: odd and > 2}
var
  j,k: integer;
  t: longint;
  m8: mp_digit;
begin
  mp_jacobi_lm := 0;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_jacobi_lm');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Error if n<3 or even}
  if (n.sign=MP_NEG) or mp_is1(n) or mp_iseven(n) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_jacobi_lm: n<3 or even');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if a=0 then exit;

  {initialize accumulated result}
  j := 1;
  m8 := n.pdigits^[0] and 7;

  if a<0 then begin
    {(-a|b) = (a|b)*(-1|b); change sign if b=3 (mod 4)}
    if m8 and 3 = 3 then j := -j;
    a := -a;
  end;

  {if a is even divide out highest power of two: 2^k}
  if a and 1 = 0 then begin
    k:=0;
    repeat
      inc(k);
      a := a shr 1;
    until odd(a);
    {if odd exponent use (2|b) = -1 if b = +-3 (mod 8)}
    if odd(k) and ((m8=3) or (m8=5)) then j:=-j;
  end;

  {done if a=1}
  if a=1 then mp_jacobi_lm := j
  else begin
    {adjust result using quadratic reciprocity}
    if (a and 3 = 3) and (m8 and 3 = 3) then j := -j;
    {swap variables, reduce mp_int mod longint, and use jacobi32}
    mp_mod_int(n,a,t);
    if t<>0 then mp_jacobi_lm := j*jacobi32(t,a);
  end;
end;


{---------------------------------------------------------------------------}
function mp_jacobi_ml(const a: mp_int; n: longint): longint;
  {-Compute the Jacobi/Legendre symbol (a|n), n: odd and > 2}
var
  am: longint;
begin
  mp_jacobi_ml := 0;
  if mp_error<>MP_OKAY then exit;

  {Error if n<3 or even}
  if (n<3) or (n and 1 = 0) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_jacobi_ml: n<3 or even');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {mp_mod_int does init check on a and returns am = a mod n}
  mp_mod_int(a,n,am);
  mp_jacobi_ml := jacobi32(am,n);
end;


{---------------------------------------------------------------------------}
procedure mp_lcm(const a,b: mp_int; var c: mp_int);
  {-Compute least common multiple as |a*b|/(a, b)}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_lcm');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if a=0 or b=0 then c:=0}
  if mp_is0(a) or mp_is0(b) then begin
    mp_zero(c);
    exit;
  end;

  mp_init(t);
  if mp_error<>MP_OKAY then exit;

  {get the gcd of the two inputs}
  if mp_gcd1(a,b,t) then mp_mul(a,b,c)
  else begin
    {divide the largest by the GCD}
    if mp_cmp_mag(a,b)=MP_GT then begin
      {lcm  = (a div t)*b}
      mp_div(a, t, t);
      mp_mul(b, t, c);
    end
    else begin
      {lcm  = a * (b div t)}
      mp_div(b, t, t);
      mp_mul(a, t, c);
    end;
  end;
  {fix the sign to positive}
  c.sign := MP_ZPOS;
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_mulmod(const a,b,c: mp_int; var d: mp_int);
  {-Calculate d = a * b (mod c)}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) or mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mulmod');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init(t);
  if mp_error<>MP_OKAY then exit;

  mp_mul(a, b, t);
  mp_mod(t, c, d);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_sqrmod(const a,c: mp_int; var d: mp_int);
  {-Calculate d = a * a (mod c)}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(c) or mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sqrmod');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_init(t);
  if mp_error=MP_OKAY then begin
    mp_sqr(a, t);
    mp_mod(t, c, d);
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_submod(const a,b,c: mp_int; var d: mp_int);
  {-Calculate d = a - b (mod c)}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) or mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_submod');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init(t);
  if mp_error<>MP_OKAY then exit;

  mp_sub(a, b, t);
  mp_mod(t, c, d);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_xgcd(const a,b: mp_int; p1,p2,p3: pmp_int);
  {-Extended gcd algorithm, calculate  a*p1^ + b*p2^ = p3^ = gcd(a,b)}
  { p1,p2,p3 may be nil if the values are not required.}
var
  u1,u3,v1,v3,t1,t3,q: mp_int;
  x0,x1,y0,y1: tml_int;
  TryLehmer: boolean;
const
  k=Lehmer_MaxK; k2=(k+1) div 2;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or ((p1<>nil) and mp_not_init(p1^)) or
       ((p2<>nil) and mp_not_init(p2^)) or ((p3<>nil) and mp_not_init(p3^))
    then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_xgcd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if ((p1=p2) and (p1<>nil)) or ((p1=p3) and (p1<>nil)) or ((p2=p3) and (p2<>nil)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_xgcd: identical pointers');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}

  if (p1=nil) and (p2=nil) then begin
    {easy out for trivial cases}
    if p3<>nil then mp_gcd(a,b,p3^);
    exit;
  end;

  mp_init7(u1,u3,v1,v3,t1,t3,q);
  if mp_error<>MP_OKAY then exit;

  {Note: the following setup returns 1*0 + 0*0 = 0, for a=b=0! This}
  {is somewhat strange but OK. Knuth's algorithm X gives the same. }
  {Usage of u2,v2,t2 suppressed, p2^ = u2 will be = (u3 - a*u1)/b. }
  {We need positive working variables u3,v3. So init via abs, the  }
  {sign of u3 is corrected after loop and before calculation of u2.}

  {initialize, (u1,u3) = (1,a)}
  mp_set1(u1);
  mp_abs(a, u3);

  {initialize, (v1,v3) = (0,b), v1=0 after init}
  mp_abs(b, v3);

  {The Lehmer step needs u3 >= v3 > 0. This condition may no hold}
  {initially. For extended GCD we cannot simply swap u3 and v3.  }
  {Therefore we skip the step in the first iteration, after that }
  {the condition becomes true via the Euclidean mod operation.   }
  TryLehmer := mp_cmp_mag(u3,v3)<>MP_LT;

  {loop while v3 != 0}
  while (mp_error=MP_OKAY) and (not mp_iszero(v3)) do begin
    if TryLehmer and (mp_bitsize(u3)-mp_bitsize(v3) <= k2) then begin
      if Lehmer(u3,v3,k,x0,x1,y0,y1) then begin
        {calculate new (ui, vi):  vi := x1*ui + y1*vi;  ui := x0*ui + y0*vi;}
        {$ifdef MP_16BIT}
          {v3 := x1*u3 + y1*v3}
          mp_mul_w(u3,abs(x1), q);  if x1<0 then s_mp_chs(q);
          mp_mul_w(v3,abs(y1),t1);  if y1<0 then s_mp_chs(t1);
          mp_add(q,t1,t1);
          {u3 := x0*u3 + y0*v3}
          mp_mul_w(u3,abs(x0), q);  if x0<0 then s_mp_chs(q);
          mp_mul_w(v3,abs(y0),t3);  if y0<0 then s_mp_chs(t3);
          mp_add(q,t3,t3);
          {the new u3,v3 should be >= 0, skip the Lehmer step if not}
          if (t1.sign=MP_ZPOS) and (t3.sign=MP_ZPOS) then begin
            mp_exch(t3,u3);
            mp_exch(t1,v3);
            {v1 := x1*u1 + y1*v1}
            mp_mul_w(u1,abs(x1), q);  if x1<0 then s_mp_chs(q);
            mp_mul_w(v1,abs(y1),t1);  if y1<0 then s_mp_chs(t1);
            mp_add(q,t1,t1);
            mp_exch(t1,v1);
            {u1 := x0*u1 + y0*v1}
            mp_mul_w(u1,abs(x0), q);  if x0<0 then s_mp_chs(q);
            mp_mul_w(t1,abs(y0),t1);  if y0<0 then s_mp_chs(t1);
            mp_add(q,t1,u1);
          end;
        {$else}
          {v3 := x1*u3 + y1*v3}
          mp_mul_d(u3,abs(x1), q);  if x1<0 then s_mp_chs(q);
          mp_mul_d(v3,abs(y1),t1);  if y1<0 then s_mp_chs(t1);
          mp_add(q,t1,t1);
          {u3 := x0*u3 + y0*v3}
          mp_mul_d(u3,abs(x0), q);  if x0<0 then s_mp_chs(q);
          mp_mul_d(v3,abs(y0),t3);  if y0<0 then s_mp_chs(t3);
          mp_add(q,t3,t3);
          {the new u3,v3 should be >= 0, skip the Lehmer step if not}
          if (t1.sign=MP_ZPOS) and (t3.sign=MP_ZPOS) then begin
            mp_exch(t3,u3);
            mp_exch(t1,v3);
            {v1 := x1*u1 + y1*v1}
            mp_mul_d(u1,abs(x1), q);  if x1<0 then s_mp_chs(q);
            mp_mul_d(v1,abs(y1),t1);  if y1<0 then s_mp_chs(t1);
            mp_add(q,t1,t1);
            mp_exch(t1,v1);
            {u1 := x0*u1 + y0*v1}
            mp_mul_d(u1,abs(x0), q);  if x0<0 then s_mp_chs(q);
            mp_mul_d(t1,abs(y0),t1);  if y0<0 then s_mp_chs(t1);
            mp_add(q,t1,u1);
          end;
        {$endif}
      end;
    end;

    {q = u3/v3, t3=u3-q*v3}
    mp_divrem(u3, v3, @q, @t3);

    {(t1,t3) = (u1,u3) - q*(v1,v3)}
    mp_mul(v1, q,  t1);
    mp_sub(u1, t1, t1);

    {(u1,u3) = (v1,v3)}
    mp_exch(v1, u1);
    mp_exch(v3, u3);

    {(v1,v3) = (t1,t3)}
    mp_exch(t1, v1);
    mp_exch(t3, v3);

    TryLehmer := true;
  end;

  if mp_error=MP_OKAY then begin
    {adjust sign if a<0}
    if a.sign=MP_NEG then s_mp_chs(u1);
    {copy results if requested}
    if p2<>nil then begin
      {here v3=0}
      if mp_iszero(b) then mp_exch(p2^, v3)
      else begin
        {use v3 to calculate p2^ = (u3 - a*u1)/b}
        mp_mul(a,u1,v1);
        mp_sub(u3,v1,v3);
        mp_div(v3,b,p2^);
      end;
    end;
    if p1<>nil then mp_exch(p1^, u1);
    if p3<>nil then mp_exch(p3^, u3);
  end;
  mp_clear7(u1,u3,v1,v3,t1,t3,q);
end;


{---------------------------------------------------------------------------}
procedure mp_xgcd_bin(const a,b: mp_int; p1,p2,p3: pmp_int);
  {-Extended binary gcd, calculate a*p1^ + b*p2^ = p3^ = gcd(a,b)  }
  { p1,p2,p3 may be nil if the values are not needed. Note that p1^}
  { and p2^ are NOT unique and may differ from those of mp_xgcd!!  }
var
  x,y,u,v,c,d: mp_int;
  g: longint;
  pt: pmp_int;
  swapped: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or ((p1<>nil) and mp_not_init(p1^)) or
       ((p2<>nil) and mp_not_init(p2^)) or ((p3<>nil) and mp_not_init(p3^))
    then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_xgcd_bin');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if ((p1=p2) and (p1<>nil)) or ((p1=p3) and (p1<>nil)) or ((p2=p3) and (p2<>nil)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_xgcd_bin: identical pointers');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}

  {Handle trivial cases}
  if mp_iszero(a) then begin
    {0*a + sign(b)*b = |b|}
    if p3<>nil then mp_abs(b,p3^);
    if p1<>nil then mp_zero(p1^);
    if p2<>nil then mp_set_short(p2^,mp_sign(b));
    exit;
  end;
  if mp_iszero(b) then begin
    {sign(a)*a + 0*b = |a|}
    if p3<>nil then mp_abs(a,p3^);
    if p2<>nil then mp_zero(p2^);
    if p1<>nil then mp_set_short(p1^,mp_sign(a));
    exit;
  end;
  if (p1=nil) and (p2=nil) then begin
    if p3<>nil then mp_gcd(a,b,p3^);
    exit;
  end;

  {Here a<>0, b<>0}
  mp_init6(x,y,u,v,c,d);

  mp_copy(a,x);
  mp_copy(b,y);
  g := 0;

  {Note mp_iseven() returns false if mp_error<>MP_OKAY}

  {Get common power of 2}
  while mp_iseven(x) and mp_iseven(y) do begin
    mp_shr1(x);
    mp_shr1(y);
    inc(g);
  end;

  {Core routine needs odd x, use symmetry of x*p1^ + y*p2^ = gcd(x,y)}
  {and swap x,y and the pointers to the corresponding factors}
  if mp_iseven(x) then begin
    mp_exch(x,y);
    pt := p1;
    p1 := p2;
    p2 := pt;
    swapped := true;
  end
  else swapped := false;

  {Here x odd, perform binary extended core routine v=gcd(x,y).}
  {This is essentially the same routine as in the old}
  {fast_mp_invmod from LTM / HAC 14.61, pp608}

  mp_abs(y,v);
  mp_abs(x,u);
  mp_set1(d);

  repeat
    {u>0, v>0}
    while mp_iseven(u) do begin
      mp_shr1(u);
      if mp_isodd(c) then mp_sub(c, x, c);
      mp_shr1(c);
    end;
    while mp_iseven(v) do begin
      mp_shr1(v);
      if mp_isodd(d) then mp_sub(d, x, d);
      mp_shr1(d);
    end;
    if (u.used>v.used) or (mp_cmp(u, v)<>MP_LT) then begin
      {u >= v, v unchanged > 0}
      mp_sub(u, v, u);
      mp_sub(c, d, c);
    end
    else begin
      {u<v, new v remains > 0}
      mp_sub(v, u, v);
      mp_sub(d, c, d);
    end;
  until (mp_error<>MP_OKAY) or mp_iszero(u);

  {Multiply by common power of 2}
  mp_shl(v,g,v);

  {Adjust coefficient of y}
  if y.sign=MP_NEG then s_mp_chs(d);

  {Here (??)*p1^ + d*p2^ = v = gcd(a,b)}
  {Calculate p1^ if requested}
  if p1<>nil then begin
    if swapped then begin
      if mp_iszero(b) then mp_zero(p1^)
      else begin
        {p1^ = (v - d*a)/b}
        mp_mul(d, a, y);
        mp_sub(v, y, y);
        mp_div(y, b, p1^);
      end;
    end
    else begin
      if mp_iszero(a) then mp_zero(p1^)
      else begin
        {p1^ = (v - d*b)/a}
        mp_mul(d, b, y);
        mp_sub(v, y, y);
        mp_div(y, a, p1^);
      end;
    end;
  end;

  {Store p2^ and p3^ if requested}
  if p2<>nil then mp_exch(p2^, d);
  if p3<>nil then mp_exch(p3^, v);

  mp_clear6(x,y,u,v,c,d);
end;


{---------------------------------------------------------------------------}
procedure mp_xlcm(const a,b: mp_int; var c,x,y: mp_int);
  {-Calculate c,x,y with lcm(a,b)=c=x*y and x|a, y|b, gcd(x,y)=1.}
  { c,x,y should be different variables, otherwise result may be inconsistent}
var
  g,t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) or mp_not_init(x) or mp_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_xlcm');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mp_is0(a) or mp_is0(b) then begin
    {if a=0 or b=0 return (0,|a|,|b|)}
    {The following code looks unnecessary complicated but is}
    {needed if the variables 'overlap' @a=@y, @b=@x etc}
    if mp_is0(a) then begin
      mp_abs(b,y);
      mp_zero(x);
    end
    else begin
      mp_abs(a,x);
      mp_zero(y);
    end;
    mp_zero(c);
    exit;
  end;

  mp_init2(g,t);
  if mp_error<>MP_OKAY then exit;

  {t = |a|}
  mp_abs(a,t);
  {first get gcd(a,b)}
  if mp_gcd1(a,b,g) then begin
    {gcd(a,b)=1:  x=|a|, y=|b|, c=x*y}
    mp_abs(b,g);
    mp_mul(t,g,c);
    mp_exch(t,x);
    mp_exch(g,y);
  end
  else begin
    {g = |b|/gcd(a,b)}
    mp_div(b,g,g); g.sign := MP_ZPOS;
    {c = lcm(a,b) = t*g}
    mp_mul(t,g,c);
    {now calculate x alias t}
    while (mp_error=MP_OKAY) and not mp_gcd1(t,g,g) do mp_div(t,g,t);
    {get y=lcm(a,b)/x}
    mp_div(c,t,y);
    mp_exch(t,x);
  end;
  mp_clear2(g,t);
end;

begin
  mp_initial_mod := true;
end.
