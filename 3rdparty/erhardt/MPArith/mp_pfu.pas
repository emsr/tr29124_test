unit mp_pfu;

{Multi precision prime factorization routines}

interface

{$i STD.INC}

uses
  mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION   :  Multi precision prime factorization routines

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  heap

 DISPLAY MODE  :  (text if mp_show_progress)

 REFERENCES    :  [3] Knuth, D.E.: The Art of computer programming. Vol 2
                      Seminumerical Algorithms, 3rd ed., 1998
                  [4] Forster, O.: Algorithmische Zahlentheorie, 1996
                  [6] R. P. Brent, Factor: an integer factorization program for
                      the IBM PC, Report TR-CS-89-23, October 1989, 7 pp.
                      http://maths-people.anu.edu.au/~brent/pub/pub117.html
                      http://maths-people.anu.edu.au/~brent/ftp/rpb117/rpb117.exe
                  [8] Marcel Martin: NX - Numerics library of multiprecision
                      numbers for Delphi and Free Pascal, 2006-2009
                      www.ellipsa.eu/public/nx/index.html
                 [10] Crandall,R., C.Pomerance: Prime Numbers, A Computational
                      Perspective, 2nd ed., 2005
                 [40] H. Riesel, Prime Numbers and Computer Methods for Factorization,
                      Vol. 126 of Progress in Mathematics, Boston, 2nd ed. 1994.
                      Paperback reprint 2012 in Modern Birkh„user Classics Series.


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.23.00  24.09.12  W.Ehrhardt  Factorization routines from mp_numth
 1.29.00  14.07.14  we          mp_fermat_factor
 1.30.00  03.08.14  we          mp_holf: Hart's OneLineFactor
**************************************************************************)

(*-------------------------------------------------------------------------
  This code uses material/ideas from the following 3rd party libraries:
   - LibTomMath 0.30+ by Tom St Denis
   - MPI 1.8.6 by Michael J. Fromberger
   - NX V0.18 and V0.9+ by Marcel Martin
  See the file '3rdparty.mpa' for the licenses.
----------------------------------------------------------------------------*)


(*-------------------------------------------------------------------------
 (C) Copyright 2004-2014 Wolfgang Ehrhardt

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
  mp_max_small = mp_digit(MP_DIGIT_MAX and $7FFF); {max. tested small factor}


type
  TProgressProc = procedure(checkonly: boolean; cnt,maxcnt: longint; var abort: boolean);
                    {procedure to show progress in factorization routines.}
                    {If on return abort is true, factorization is aborted.}
                    {checkonly: check only for abort, not recommended to show progress,}
                    {cnt, maxcnt: current and maximum iteration cnt, maxcnt=0; no maximum known,}
                    {abort: set to true to abort calculation, proc is called with abort=false}

var
  mp_show_progress: boolean;  {Factorization routines rho, (p+1), (p-1) call}
                              {progress procedure after one accum cycle}

var
  mp_isconsole    : boolean;  {If true use writeln if mp_show_progress is}
                              {true, but no progress procedure is defined}

var
  mp_max_small_sqr: longint;  {square of max small factor, ie numbers less}
                              {are prime if mp_small_factor returns f=0}

const
  ECM_C1Min =  100;       {Minimum and maximum of C1 in mp_ecm_factor  }
  ECM_C1Max = 6500;       {Phase #1 makes use of primes up to Prime[C1] <= NumPrimes16 }


procedure mp_set_progress(const PP: TProgressProc);
  {-Make PP the new progress proc}


procedure mp_ecm_factor(const n: mp_int; var f: mp_int; CMax,C1: word; var seed,phase: longint);
  {-Find a factor f of n by inversionless ECM method, f=0 if no factor found.}
  { n must be > 1 and @f<>@n otherwise a runtime error / exception is raised.}
  { Two-phase elliptic curve algorithm, C1 primes are used in phase 1, C1 is}
  { clipped to [ECM_C1Min,ECM_C1Max]. Up to CMax curves (32 if CMax=0), if}
  { input seed=0, mp_random_int is used as initial value. The code is based}
  { on Marcel Martin's NX v0.18/FPC implementation of Algorithm 7.4.4 from}
  { Crandall/Pomerance, incorporating various enhancements of Brent,}
  { Crandall, Montgomery, Woltman, and Zimmermann. WE changes: with Barrett}
  { reduction and intermediate GCDs, 32 bit prime table via config.}

procedure mp_ecm_simple(const n: mp_int; var f: mp_int; seed: longint);
  {-Simplified version of mp_ecm_factor with CMax=0, C1=ECM_C1Min}

procedure mp_fermat_factor(const N: mp_int; var f: mp_int; cnt: longint);
  {-Find a factor f of N with Fermat's method, f=0 if no factor found; cnt: number of tests}

procedure mp_holf(const n: mp_int; var f: mp_int; cnt: longint);
  {-Find a factor f of n with Hart's OneLineFactor, f=0 if no factor found; cnt: number of tests}

procedure mp_pollard_brent(const N: mp_int; var f: mp_int);
  {-Find a factor f of N with Brent's version of Pollard rho, f=0 if no factor found; c=1, rmax=8192}

procedure mp_pollard_brent_ex(const N: mp_int; var f: mp_int; c: mp_digit; rmax: word);
  {-Find a factor f of N with Brent's version of Pollard rho, f=0 if no factor found;}
  { c = constant in iteration formula, rmax <= 16384 is the maximal r}

procedure mp_pollard_pm1(const N: mp_int; var f: mp_int; bound: word);
  {-Find a factor f of N with p-1 method, f=0 if no factor found}
  { primes <= bound <= 65000 will be included in pp_expo product}

procedure mp_pollard_rho(const N: mp_int; var f: mp_int; cnt: longint);
  {-Find a factor f of N with Pollard's rho method, f=0 if no factor found; cnt: number accumulation cycles}

procedure mp_squfof(const n: mp_int; var f, res: longint);
  {-Find a factor f of n < 2^60 with Shanks' SQUFOF method; n must be}
  { composite. If res=0 then f is a factor; if res<>0, no factor is  }
  { found. mp_squfof is based on Riesel's SQUFOF code. Error values: }
  {#F}
  { res = -1:  general error                                         }
  { res = -2:  n is < 2 or >= 2^60                                   }
  { res = -3:  forward  cycle: no square form found                  }
  { res = -4:  backward cycle: no factor found                       }
  { res = -5:  queue overflow                                        }
  { res >  0:  search period (=res) for CF of sqrt(n) is reached     }
  {#F}

procedure mp_williams_pp1(const N: mp_int; var f: mp_int; bound,numtest: word);
  {-Find a factor f of N with William's p+1 method, f=0 if no success. numtest}
  { random seeds will be tried, should be about 3 or 4, if 0 then 3 is used}
  { primes <= bound <= 65000 will be included in pp_expo product}


procedure s_mp_coshmult(const a,b,c,mc: mp_int; var d: mp_int);
  {-Internal coshmult, d=mod_coshmult(a,b,c), mc: Barrett parameter for c}


implementation


uses
  mp_base, mp_prime, mp_modul, mp_numth, mp_prng;

var
  mp_progress: TProgressProc;  {Progress function for factoring procedures}


{---------------------------------------------------------------------------}
function ProgressAssigned: boolean;
  {-Check if progress is <> nil}
begin
  {$ifdef FPC}
    {$ifdef FPC_DELPHI}
      ProgressAssigned := @mp_progress<>nil;
    {$else}
      ProgressAssigned := mp_progress<>nil;
    {$endif}
  {$else}
    ProgressAssigned := @mp_progress<>nil;
  {$endif}
end;


{---------------------------------------------------------------------------}
procedure mp_set_progress(const PP: TProgressProc);
  {-Make PP the new progress proc}
begin
  mp_progress := PP;
end;


{---------------------------------------------------------------------------}
procedure mp_ecm_factor(const n: mp_int; var f: mp_int; CMax,C1: word; var seed,phase: longint);
  {-Find a factor f of n by inversionless ECM method, f=0 if no factor found.}
  { n must be > 1 and @f<>@n otherwise a runtime error / exception is raised.}
  { Two-phase elliptic curve algorithm, C1 primes are used in phase 1, C1 is}
  { clipped to [ECM_C1Min,ECM_C1Max]. Up to CMax curves (32 if CMax=0), if}
  { input seed=0, mp_random_int is used as initial value. The code is based}
  { on Marcel Martin's NX v0.18/FPC implementation of Algorithm 7.4.4 from}
  { Crandall/Pomerance, incorporating various enhancements of Brent,}
  { Crandall, Montgomery, Woltman, and Zimmermann. WE changes: with Barrett}
  { reduction and intermediate GCDs, 32 bit prime table via config.}
label
  nextcurve, gcd2, final, err1;
const
  D = 200;            {3*(D+2) mp_ints of order n in vectors xx,zz,xz}
type
  TECMVec = array[0..D+1] of mp_int;
  PECMVec = ^TECMVec;
var
  i: integer;
  curve: word;                 {curve counter}
  j,k,p1,p,q,max,kD1,kDBound: longint;
  C2,CC,CP: longint;           {C2 and progress bound calculated from C1}
  an,ad: mp_int;               {constants for curve, both reduced mod n}
  xr,zr,xs,zs: mp_int;         {Q=[xr:zr] in phase 1}
  x0,z0: mp_int;               {used only in ecm_mul}
  mu: mp_int;                  {Barrett parameter for n}
  t1,t2,t3: mp_int;            {temporary}
  xx,zz,xz: PECMVec;           {Phase 2 vectors}
  fs: mp_digit;                {small factor}
  allocated: boolean;          {Phase 2 vectors allocated and initialized}
  cancel: boolean;             {Progress cancel flag}

{$ifdef MPC_ECM_Primetable}
const
  C2Max = ECM_C1Max*100;  {Prime[C2Max] should be smaller than 2^31}
type
  TPrimeTable = array[1..C2MAX] of longint;
var
  ppt: ^TPrimeTable;
{$endif}

  {----------------------------------------------}
  procedure ecm_double(var x,z: mp_int);
    {-Double the point: (x,z) := 2*(x,z)}
    { assume input in [0,n), return output in [0,n)}
  begin
    {t1 := (x + z)^2}
    mp_add(x,z,t1);        if mp_cmp(n,t1)<>MP_GT then mp_sub(t1,n,t1);
    mp_sqr(t1,t1);         mp_reduce(t1,n,mu);
    {t2 := (x - z)^2}
    mp_sub(x,z,t2);        {no need to make positive (x-z)^ < n^2}
    mp_sqr(t2,t2);         mp_reduce(t2,n,mu);
    {t3 := t1 - t2}
    mp_sub(t1,t2,t3);      if t3.sign=MP_NEG then mp_add(t3,n,t3);
    {t2 := 4*t2}
    mp_mul_2(t2,t2);       if mp_cmp(n,t2)<>MP_GT then mp_sub(t2,n,t2);
    mp_mul_2(t2,t2);       if mp_cmp(n,t2)<>MP_GT then mp_sub(t2,n,t2);
    mp_mul(t1,t2,x);       mp_reduce(x,n,mu);
    {x := t1*t2*ad}
    mp_mul(x,ad,x);        mp_reduce(x,n,mu);
    mp_mul(t2,ad,t2);      mp_reduce(t2,n,mu);
    {t1 := t3*an}
    mp_mul(t3,an,t1);      mp_reduce(t1,n,mu);
    {z := (t2+t1)*t3}
    mp_add(t2,t1,z);       if mp_cmp(n,z)<>MP_GT then mp_sub(z,n,z);
    mp_mul(z,t3,z);        mp_reduce(z,n,mu);
  end;


  {----------------------------------------------}
  procedure ecm_add(var x2,z2,x1,z1,xd,zd: mp_int);
    {-Add (x2,z2) = (x2,z2) + (x1,z1),  with (xd,zd) = (x2,z2) - (x1,z1)}
    { assume input in [0,n), return output in [0,n)}
  begin
    {t1 := x1-z1}
    mp_sub(x1,z1,t1);  if t1.sign=MP_NEG      then mp_add(t1,n,t1);
    {t2 := x2+z2}
    mp_add(x2,z2,t2);  if mp_cmp(n,t2)<>MP_GT then mp_sub(t2,n,t2);
    {t2 := (x2+z2)(x1-z1)}
    mp_mul(t2,t1,t2);  mp_reduce(t2,n,mu);
    {t1 := x1+z1}
    mp_add(x1,z1,t1);  if mp_cmp(n,t1)<>MP_GT then mp_sub(t1,n,t1);
    {t3 := x2-z2}
    mp_sub(x2,z2,t3);  if t3.sign=MP_NEG      then mp_add(t3,n,t3);
    {t1 := (x1+z1)(x2-z2)}
    mp_mul(t1,t3,t1);  mp_reduce(t1,n,mu);
    {x+ := (x2+z2)(x1-z1)+(x1+z1)(x2-z2) = 2(x1x2-z1z2)}
    mp_add(t2,t1,x2);  if mp_cmp(n,x2)<>MP_GT then mp_sub(x2,n,x2);
    mp_sqr(x2,x2);     mp_reduce(x2,n,mu);
    {x+ := 4zd*(x1x2-z1z2)^2}
    mp_mul(x2,zd,x2);  mp_reduce(x2,n,mu);
    {z+ := (x2+z2)(x1-z1)-(x1+z1)(x2-z2) = 2(x1z2-x2z1)}
    mp_sub(t2,t1,z2);  {no need to make positive (t2-t1)^2 < n^2}
    mp_sqr(z2,z2);     mp_reduce(z2,n,mu);
    {z+ := 4xd*(x1z2-x2z1)^2}
    mp_mul(z2,xd,z2);  mp_reduce(z2,n,mu);
    {cf C/P [10] (7.6), with B=0, A=1. Factors 4 are irrelevant for x/z}
  end;

  {----------------------------------------------}
  procedure ecm_mul_int(var x,z: mp_int; e: longint);
    {-Elliptic multiplication (Montgomery method) (x,y) := e*(x,y), e>0}
    { assume input x,z in (0,n), return output x,z in [0,n)}
  var
    c: longint;
  begin
    if e > 2 then begin
      mp_copy(x,x0);
      mp_copy(z,z0);
      mp_copy(x,xs);
      mp_copy(z,zs);
      ecm_double(xs,zs);
      c := longint($40000000);
      while (c and e) = 0 do c := c shr 1;
      c := c shr 1;
      repeat
        if (e and c) = 0 then begin
          ecm_add(xs,zs,x,z,x0,z0);
          ecm_double(x,z);
        end
        else begin
          ecm_add(x,z,xs,zs,x0,z0);
          ecm_double(xs,zs);
        end;
        c := c shr 1;
      until c = 0;
    end
    else if e=2 then ecm_double(x,z);
    {else if e=1 return x,z unchanged}
  end;

  {----------------------------------------------}
  procedure ecm_setup(var seed: longint);
    {-Calculate ad,an from random seed in [6,n-1], and point Q on curve}
    { C = an/ad - 2  mod n, from curve y^2 = x^3 + Cx^2 + x)}
    { Q : a point on the curve, with coordinates [xr:zr]}
  begin
    inc(seed);
    if (seed<2) or (seed=5) then seed := 6;

    {initial point}
    mp_set_int(xs, seed);
    mp_sqr(xs,xs);
    mp_sub_d(xs,5,xs);        {u  = s^2-5}

    mp_set_int(zs, seed);
    mp_shl(zs,2,zs);          {v  = 4s}

    mp_exptmod_d(xs,3,n,xr);  {xr = u^3}
    mp_exptmod_d(zs,3,n,zr);  {zr = v^3}

    {coef C = an/ad - 2}
    mp_sub(zs,xs,an);
    mp_exptmod_d(an,3,n,an);  { (v-u)^3}
    mp_mul_d(xs,3,xs);
    mp_add(xs,zs,ad);
    mp_mulmod(an,ad,n,an);    { (v-u)^3*(3u+v)}

    mp_shl(zs,2,zs);
    mp_mulmod(xr,zs,n,ad);    { 4*u^3*v}
  end;

begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) or mp_not_init(f)  then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_ecm_factor');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Check n>1 and @n<>@f}
  if mp_cmp_d(n,2)=MP_LT then begin
    {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_ecm_factor: n < 2');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
  if @n=@f then begin
    {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_ecm_factor: @n=@f');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {check small factor up to 127}
  mp_small_factor(n,2,127,fs);
  if fs<>0 then begin
    mp_set(f,fs);
    exit;
  end;

  if C1<ECM_C1Min then C1 := ECM_C1Min;
  if C1>ECM_C1Max then C1 := ECM_C1Max;
  C2 := 100*C1;

  if CMax=0 then CMax := 32;
  curve := 0;

  {Not yet allocated, set pointers to nil. Allocation and initialization}
  {of dynamic arrays will be done if phase 2 is entered the first time}
  allocated := false;
  xx := nil;
  zz := nil;
  xz := nil;
  {$ifdef MPC_ECM_Primetable}
    ppt := nil;
  {$endif}

  {initialize local mp_ints and calculate Barrett parameter for n}
  mp_init6(xr,zr,xs,zs,an,ad); if mp_error<>0 then exit;
  mp_init6(x0,z0,t1,t2,t3,mu); if mp_error<>0 then goto err1;
  mp_reduce_setup(mu,n);

  {generate seed for first curve, use random if not in (5,n)}
  while (seed<2) or (seed=5) or not (mp_cmp_int(n,seed)=MP_GT) do begin
    if mp_bitsize(n)<31 then seed := 6 + mp_random_digit mod (MP_DIGIT_MAX-6)
    else seed := 6+mp_random_int;
    if mp_error<>MP_OKAY then goto final;
  end;

  {Marcel Martin uses a prime table of C2 longints. Here the Primes[j]}
  {can be generated via next/prevprime32 as a space/time tradeoff:}
  {Always for 16 bit, used for 32 bit if MPC_ECM_Primetable is not defined.}

  p1 := prime32(C1);

nextcurve:

  mp_set1(f);

  inc(curve);
  ecm_setup(seed);

  {---------------------------}
  {---- Perform phase 1 ------}
  {---------------------------}

  phase := 1;
  CP:= C1 div 4; {Progress call}
  CC:= CP;

  {loop over primes}
  for i:=1 to C1 do begin
    if mp_error<>MP_OKAY then goto final;
    p := Primes16[i];
    q := p;
    {find max exponent max with p^max<B1=Prime[C1], p=Pime[i]}
    max := p1 div p;
    while p <= max do p := p*q;
    {Q := [p^max]*Q}
    ecm_mul_int(xr,zr,p);
    {Check if cancel, no progress indicator}
    if odd(i) and (ProgressAssigned) then begin
      cancel := false;
      mp_progress(true, 0, C1, cancel);
      if cancel then begin
        mp_zero(f);
        goto final;
      end;
    end;
    if i>CC then begin
      inc(CC,CP);
      if mp_show_progress then begin
        {show progress and check for cancel}
        if ProgressAssigned then begin
          cancel := false;
          mp_progress(false, 1, C1, cancel);
          if cancel then begin
            mp_zero(f);
            goto final;
          end;
        end
        else if mp_isconsole then write('.');
      end;
      {WE: check non-trivial GCD after CP primes, this is not done in C/P/M}
      mp_gcd(zr,n,f);
      if (mp_cmp_d(f,1)=MP_GT) and mp_is_lt(f,n) then goto final;
    end;
  end;

  {gcd(z,n)}
  if not mp_gcd1(zr,n,f) then begin
    {GCD not 1, test other trivial case and exit if non-trivial}
    if mp_is_eq(f,n) then begin
      if curve<CMAX then goto nextcurve  {try next curve}
      else mp_zero(f); {indicate failure}
    end;
    goto final;
  end;

  {---------------------------}
  {---- Perform phase 2 ------}
  {---------------------------}

  phase := 2;
  if not allocated then begin
    allocated := true;
    {$ifdef MPC_ECM_Primetable}
      ppt := IAlloc(C2*sizeof(longint));
      nextprime32_array(1,C2,ppt^);
    {$endif}
    xx := mp_getmem(sizeof(TECMVec));
    zz := mp_getmem(sizeof(TECMVec));
    xz := mp_getmem(sizeof(TECMVec));
    if xx<>nil then mp_init_multi(xx^);
    if zz<>nil then mp_init_multi(zz^);
    if xz<>nil then mp_init_multi(xz^);
  end;

  k := (p1+1) div D;
  mp_copy(xr,xx^[0]);
  mp_copy(zr,zz^[0]);
  kD1 := k*D+1;
  ecm_mul_int(xx^[0],zz^[0],kD1);

  i := D+1;
  mp_copy(xr,xx^[i]);
  mp_copy(zr,zz^[i]);
  ecm_mul_int(xx^[i],zz^[i],kD1+D+D);

  {i=1}
  mp_copy(xr,xx^[1]);
  mp_copy(zr,zz^[1]);
  ecm_double(xx^[1],zz^[1]);
  mp_mul(xx^[1],zz^[1],xz^[1]);
  mp_reduce(xz^[1],n, mu);

  {i=2}
  mp_copy(xx^[1],xx^[2]);
  mp_copy(zz^[1],zz^[2]);
  ecm_double(xx^[2],zz^[2]);
  mp_mul(xx^[2],zz^[2],xz^[2]);
  mp_reduce(xz^[2],n,mu);

  for i:=3 to D do begin
    if mp_error<>MP_OKAY then goto final;
    mp_copy(xx^[i-1],xx^[i]);
    mp_copy(zz^[i-1],zz^[i]);
    ecm_add(xx^[i],zz^[i],xx^[1],zz^[1],xx^[i-2],zz^[i-2]);
    mp_mul(xx^[i],zz^[i],xz^[i]);
    mp_reduce(xz^[i],n,mu);
  end;

  CP:= C2 div 5; {Progress call and check GCD five times in phase 2}
  CC:= CP;
  j := C1+1;
  {$ifdef MPC_ECM_Primetable}
    if (ppt<>nil) and (j<=C2) then p:=ppt^[j]
    else p := nextprime32(p1+1);
  {$else}
    p := nextprime32(p1+1);
  {$endif}

  repeat
    if mp_error<>MP_OKAY then goto final;
    mp_mul(xx^[0],zz^[0],xz^[0]);
    mp_reduce(xz^[0],n,mu);
    kD1 := k*D + 1;
    if p=kD1 then begin
      if j=C2 then goto gcd2;
      inc(j);
      {$ifdef MPC_ECM_Primetable}
        if (ppt<>nil) and (j<=C2) then p:=ppt^[j]
        else p := nextprime32(p+1);
      {$else}
        p := nextprime32(p+1);
      {$endif}
    end;
    kDBound := kD1 + D + D - 2;
    repeat
      if mp_error<>MP_OKAY then goto final;
      i := (p-kD1) shr 1;
      {accumulate (xx0 - xxi)(zz0 + zzi) - xx0*zz0 + xxi*zzi}
      mp_sub(xx^[0],xx^[i],t1);   if t1.sign=MP_NEG      then mp_add(t1,n,t1);
      mp_add(zz^[0],zz^[i],t2);   if mp_cmp(n,t2)<>MP_GT then mp_sub(t2,n,t2);
      mp_mul(t2,t1,t2);           mp_reduce(t2,n,mu);
      mp_sub(t2,xz^[0],t2);       if t2.sign=MP_NEG then mp_add(t2,n,t2);
      mp_add(t2,xz^[i],t2);       if mp_cmp(n,t2)<>MP_GT then mp_sub(t2,n,t2);
      mp_mul(f,t2,f);             mp_reduce(f,n,mu);
      if j=C2 then goto gcd2;
      inc(j);
      {$ifdef MPC_ECM_Primetable}
        if (ppt<>nil) and (j<=C2) then p:=ppt^[j]
        else p := nextprime32(p+1);
      {$else}
        p := nextprime32(p+1);
      {$endif}
    until p > kDBound;

    if ProgressAssigned then begin
      cancel := false;
      if j>CC then begin
        mp_progress(false, j, C2, cancel);
        inc(CC,CP);
        if not mp_gcd1(f,n,t1) and not mp_is_eq(f,t1) then begin
          {non-trivial intermediate GCD found, set factor and leave}
          mp_exch(f,t1);
          goto final;
        end;
      end
      else mp_progress(true, j, C2, cancel);
      if cancel then begin
        mp_zero(f);
        goto final;
      end;
    end;

    inc(k,2);
    i := D+1;
    mp_copy(xx^[i],xs);
    mp_copy(zz^[i],zs);
    ecm_add(xx^[i],zz^[i],xx^[i-1],zz^[i-1],xx^[0],zz^[0]);
    mp_exch(xx^[0],xs);
    mp_exch(zz^[0],zs);
  until false;

gcd2:

  mp_gcd(f,n,f);
  if (mp_cmp_d(f,1)=MP_GT) and mp_is_lt(f,n) then begin
    {done if non-trivial GCD}
    goto final;
  end;

  {Try another curve if maximum curve count not exceeded}
  if curve < CMAX then goto nextcurve;

  {indicate failure}
  mp_zero(f);

final:

  {"finally" part, deallocate and clear resources}
  if allocated then begin
    {$ifdef MPC_ECM_Primetable}
      freemem(pointer(ppt),C2*sizeof(longint));
    {$endif}
    if xx<>nil then begin
      mp_clear_multi(xx^);
      mp_freemem(pointer(xx),sizeof(TECMVec));
    end;
    if zz<>nil then begin
      mp_clear_multi(zz^);
      mp_freemem(pointer(zz),sizeof(TECMVec));
    end;
    if xz<>nil then begin
      mp_clear_multi(xz^);
      mp_freemem(pointer(xz),sizeof(TECMVec));
    end;
  end;
  mp_clear6(x0,z0,t1,t2,t3,mu);

err1:
  mp_clear6(xr,zr,xs,zs,an,ad);

end;


{---------------------------------------------------------------------------}
procedure mp_ecm_simple(const n: mp_int; var f: mp_int; seed: longint);
  {-Simplified version of mp_ecm_factor with CMax=0, C1=ECM_C1Min}
var
  phase: longint;
begin
  mp_ecm_factor(n,f,0,ECM_C1Min,seed,phase);
end;


{---------------------------------------------------------------------------}
procedure s_mp_ppexpo(var x: mp_int; B0,B1: word);
  {-Product of primes B0 < p <= B1 and integers isqrt(B0) < n <= isqrt(B1)}
var
  m0, m1, i: word;
  function isqrt(w: word): word;
  begin
    {No error possible}
    isqrt := word(trunc(sqrt(w)));
  end;
begin
  if mp_error<>MP_OKAY then exit;
  {Ref: Forster[4], 14, "Die (p-1)-Faktorisierungs-Methode"}

  {$ifdef MPC_ArgCheck}
    if mp_not_init(x) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_ppexpo');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_set1(x);
  if mp_error=MP_OKAY then begin
    m0 := isqrt(B0)+1; if m0<2 then m0:=2;
    m1 := isqrt(B1);
    for i := m0 to m1 do begin
      mp_mul_w(x,i,x);
      if mp_error<>MP_OKAY then break;
    end;
    if odd(B0) then inc(B0);
    i := B0+1;
    while (i<=B1) and (mp_error=MP_OKAY) do begin
      if IsPrime16(i) then mp_mul_w(x,i,x);
      inc(i,2);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_williams_pp1(const N: mp_int; var f: mp_int; bound,numtest: word);
  {-Find a factor f of N with William's p+1 method, f=0 if no success. numtest}
  { random seeds will be tried, should be about 3 or 4, if 0 then 3 is used}
  { primes <= bound <= 65000 will be included in pp_expo product}
const
  anz0=128;
label
  leave;
var
  B0,B1: word;
  a,r,x: mp_int;
  cancel: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {Ref: Forster[4], 18, "Die (p+1)-Faktorisierungs-Methode"}

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(N) or mp_not_init(f)  then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_williams_pp1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init3(r,a,x); if mp_error<>MP_OKAY then exit;

  {Setup Barrett reduction}
  mp_reduce_setup(r, N);

  if numtest=0 then numtest := 3;
  cancel := false;

  while (numtest>0) and (not cancel) do begin
    dec(numtest);

    {a random >= 2}
    repeat
      mp_rand(a, N.used);
      mp_mod(a,N,a);
      if mp_error<>MP_OKAY then goto leave;
    until mp_cmp_d(a,1)=MP_GT;

    {x = a*a-1}
    mp_sqr(a,x);
    mp_dec(x);
    mp_gcd(x,N,f);

    {found non trivial factor, goto leave with res=MP_OKAY}
    if (mp_error<>MP_OKAY) or (mp_cmp_d(f,1)=MP_GT) then goto leave;
    if bound>65000 then bound:=65000;

    B0 := 0;
    while B0<bound do begin
      B1 := B0+anz0;
      if B1>bound then B1:=bound;
      s_mp_ppexpo(x, B0, B1);
      s_mp_coshmult(a,x,N,r,a);
      if mp_error<>MP_OKAY then goto leave;
      if mp_is1(a) then break;
      if mp_show_progress and (B0 and $100 = $100) then begin
        if ProgressAssigned then begin
          cancel := false;
          mp_progress(false, B0, 0, cancel);
          if cancel then break;
        end
        else if mp_isconsole then write('.');
      end;
      mp_sub_d(a,1,x);
      mp_gcd(x,n,f);
      {error or non trivial factor?}
      if (mp_error<>MP_OKAY) or (mp_cmp_d(f,1)=MP_GT) then goto leave;
      inc(B0, anz0);
    end;
    mp_zero(f);
  end;

leave:
   mp_clear3(x,a,r);
end;


{---------------------------------------------------------------------------}
procedure s_mp_coshmult(const a,b,c,mc: mp_int; var d: mp_int);
  {-Internal coshmult, d=mod_coshmult(a,b,c), mc: Barrett parameter for c}
var
  v,w,t,q: mp_int;
  k: longint;        {bit loop counter}
  cmp: integer;      {cmp(b,1) result}
  bk,b1: boolean;    {status bit[k], bit[k+1]}
begin

  if mp_error<>MP_OKAY then exit;

  {s_mp_coshmult = v[b] = v[b;a] mod c}
  {v[k+2] = 2|a|*v[k+1] - v[k], v[0]=1, v[1]=|a| }

  {easy outs for b<=1, b<0 is handled as b=0}
  cmp := mp_cmp_d(b,1);
  if cmp=MP_LT then begin
    {v[0] = 1}
    mp_set1(d);
    exit;
  end
  else if cmp=MP_EQ then begin
    {v[1] = |a|}
    mp_abs(a,d);
    exit;
  end;

  {create local mp_int}
  mp_init4(v,q,w,t);  if mp_error<>0 then exit;

  {compatibility with Forster: silently set q = abs(a) mod c}
  mp_abs(a,q);
  mp_mod(q,c,q);

  {initialize}
  mp_set1(v);     {v = 1}
  mp_copy(q,w);   {w = |a|}
  bk := false;    {bit[k+1] = 0}

  {v[2k]   = 2*v[k]*v[k]   - a}
  {v[2k+1] = 2*v[k]*v[k+1] - 1}

  for k:=mp_bitsize(b) downto 0 do begin
    {Here: 0 <= v,w < c}
    {exit if error}
    if mp_error<>MP_OKAY then break;
    {get bit[k], swap v,w if different from b1[k+1]}
    b1 := bk;
    bk := mp_isbit(b,k);
    if b1<>bk then mp_exch(v,w);
    {t = 2*v}
    mp_mul_2(v,t);
    if mp_cmp_mag(c,t)<>MP_GT then mp_sub(t,c,t);
    {w = 2*v*w - a}
    mp_mul(t,w,w);
    mp_reduce(w,c,mc);
    mp_sub(w,q,w);
    if w.sign=MP_NEG then mp_add(w,c,w);
    {v = 2*v*v - 1}
    mp_mul(t,v,v);
    mp_reduce(v,c,mc);
    mp_dec(v);
    if v.sign=MP_NEG then mp_add(v,c,v);
  end;
  if bk then mp_exch(w,d) else mp_exch(v,d);
  mp_clear4(t,w,q,v);
end;


{---------------------------------------------------------------------------}
procedure bigprime_pm1(const N: mp_int; var y,d: mp_int; bound: longint);
  {-Big prime stage for Pollard p-1, find factor d or d=0 if no success}
label
  found;
{$ifndef BIT16}
const
  MAXHDIFF=77;
  MAXBOUND=5000000;
{$else}
const
  MAXHDIFF=66;
  MAXBOUND=2000000;
{$endif}
var
  i,count: integer;
  q,q0: longint;                    {current, prev prime}
  x: array[0..MAXHDIFF] of mp_int;  {x[i] = y^(2i) mod N}
  t,z,mu: mp_int;
  show,cancel: boolean;
begin
  {Ref: Forster[4], 14, "Die (p-1)-Faktorisierungs-Methode"}
  {Note: Forster uses a file with the halved differences of consecutive}
  {primes, here the primes are obtained by calls to nextprime32}
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(N) or mp_not_init(y) or mp_not_init(d)  then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('bigprime_pm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init3(t,z,mu);
  if mp_error<>MP_OKAY then exit;

  mp_init_multi(x);
  if mp_error<>MP_OKAY then begin
    mp_clear3(t,z,mu);
    exit;
  end;

  {Get Barret parameter}
  mp_reduce_setup(mu, N);

  {z=y^2}
  mp_sqrmod(y,N,z);
  {calculate x[i]=y^(2i) mod N}
  mp_set1(x[0]);
  for i:=1 to MAXHDIFF do mp_mulmod(x[i-1],z,N,x[i]);

  {Note: if MAXBOUND is increased the maximum halved difference MAXHDIFF}
  {      should also be adjusted}
  if bound>MAXBOUND then bound:=MAXBOUND;
  q := 3;
  count := 0;
  show  := true;
  mp_exptmod_d(y,q,N,y);
  mp_sub_d(y,1,z);
  {Max. number of GCDs about 150 for 16Bit, 350 for > 16Bit}
  while q<bound do begin
    q0 := q;
    q := nextprime32(q0+1);
    i := (q-q0) shr 1;
    {q0 = prev prime, q= current prime, i = difference/2}
    if i>MAXHDIFF then break;
    mp_mul(y,x[i],y);   mp_reduce(y,N,mu);
    {z = z*(y-1) mod N = (z*y)-z mod N}
    mp_mul(z,y,t);      mp_reduce(t,N,mu);
    mp_sub(t,z,z);      if z.sign=MP_NEG then mp_add(z,N,z);
    inc(count);
    {Note: Forster only calculates the gcd if q>=1000}
    if (count>=1000) or (q>=bound) then begin
      mp_gcd(z,N,d);
      if mp_cmp_d(d,1)=MP_GT then begin
        {only nontrivial factors 1 < d < N}
        if mp_cmp_mag(d,N)=MP_LT then goto found;
      end;
      if mp_show_progress and show then begin
        if ProgressAssigned then begin
          cancel := false;
          mp_progress(false, q, bound, cancel);
          if cancel then break;
        end
        else if mp_isconsole then write('.');
      end;
      count := 0;
      show  := not show;
    end;
  end;

  {loop exceeded or i>MAXHDIFF}
  mp_zero(d);

found:
  mp_clear_multi(x);
  mp_clear3(t,z,mu);
end;


{---------------------------------------------------------------------------}
procedure mp_pollard_pm1(const N: mp_int; var f: mp_int; bound: word);
  {-Find a factor f of N with p-1 method, f=0 if no factor found}
  { primes <= bound <= 65000 will be included in pp_expo product}
const
  anz0=128;
label
  leave;
var
  B0,B1: word;
  base,x: mp_int;
  cancel: boolean;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(N) or mp_not_init(f)  then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_pollard_pm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Easy out}
  if mp_cmp_mag_d(N,4)=MP_LT then begin
    mp_abs(N,f);
    exit;
  end;

  {Ref: Forster[4], 14, "Die (p-1)-Faktorisierungs-Methode"}
  mp_init2(base,x);
  if mp_error<>MP_OKAY then exit;

  mp_rand(base, pred(N.used));
  mp_add_d(base,2,base);
  mp_gcd(base,N,f);

  {exit if non trivial factor}
  if mp_cmp_d(f,1)=MP_GT then goto leave;

  if bound>65000 then bound:=65000;
  B0 := 0;
  while (mp_error=MP_OKAY) and (B0<bound) do begin
    B1 := B0+anz0;
    if B1>bound then B1:=bound;
    s_mp_ppexpo(x, B0, B1);
    mp_exptmod(base,x,N,base);
    if mp_is1(base) then begin
      mp_zero(f);
      goto leave;
    end;
    if mp_show_progress and (B0 and $0100 = $100) then begin
      if ProgressAssigned then begin
        cancel := false;
        mp_progress(false, B0, 0, cancel);
        if cancel then break;
      end
      else if mp_isconsole then write('.');
    end;
    mp_sub_d(base,1,x);
    mp_gcd(x,N,f);
    if mp_cmp_d(f,1)=MP_GT then goto leave;
    inc(B0, anz0);
  end;
  if mp_is1(f) then bigprime_pm1(N, base, f, longint(bound)*100);

leave:

 mp_clear2(x,base);
end;


{---------------------------------------------------------------------------}
procedure mp_pollard_rho(const N: mp_int; var f: mp_int; cnt: longint);
  {-Find a factor f of N with Pollard's rho method, f=0 if no factor found; cnt: number accumulation cycles}
label
  leave;
const
  acc = 256;
var
  p,x,y,d,m: mp_int;
  k: integer;
  i: longint;
  cancel: boolean;
begin
  if mp_error<>MP_OKAY then exit;

  {Ref: Forster[4], 13, "Die Pollard'sche-Rho-Methode"}
  {New WE version using Barrett reduction}

  {$ifdef MPC_ArgCheck}
    if mp_not_init(N) or mp_not_init(f) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_pollard_rho');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Easy out}
  if mp_cmp_mag_d(N,4)=MP_LT then begin
    mp_abs(N,f);
    exit;
  end;

  mp_init5(m,p,x,y,d);
  if mp_error<>MP_OKAY then exit;

  {See Forster's poll_rho}
  {רררררררררררררררררררררר}

  {N.used>0!, x<N}
  mp_rand(x, pred(N.used));
  mp_copy(x,y);
  if mp_error=MP_OKAY then begin
    {Setup Barrett reduction, Forster does not use Barrett}
    mp_reduce_setup(m, N);
    for i:=1 to cnt do begin
      {accumulate differences}
      mp_set1(p);
      {if error burn some cycles}
      for k:=1 to acc do begin
        if mp_error<>MP_OKAY then goto leave;
        {x*x+2 <= (N-1)^2 + 2 = N*N - 2*N + 3 <= N*N - 1 for N>3}
        {same for y, so Barrett can be used}
        {x = x*x + 2 mod N}
        mp_sqr(x,x);
        mp_add_d(x,2,x);
        mp_reduce(x,N,m);
        {y = y*y + 2 mod N}
        mp_sqr(y,y);
        mp_add_d(y,2,y);
        mp_reduce(y,N,m);
        {y = y*y + 2 mod N}
        mp_sqr(y,y);
        mp_add_d(y,2,y);
        mp_reduce(y,N,m);
        {p = p*(y-x) mod N}
        mp_sub(y,x,d); d.sign := MP_ZPOS;
        mp_mul(p,d,p);
        mp_reduce(p,N,m);
      end;

      mp_gcd(p,N,f);
      if mp_show_progress then begin
        if ProgressAssigned then begin
          cancel := false;
          mp_progress(false, i, 0, cancel);
          if cancel then break;
        end
        else if mp_isconsole then write('.');
      end;

      {factor found if 1 < f < N}
      if (mp_cmp_d(f,1)=MP_GT) and (mp_cmp_mag(f,N)=MP_LT) then goto leave;
    end;
    {f=0 indicates no factor found}
    mp_zero(f);
  end;

leave:
  mp_clear5(m,p,x,y,d);
end;


{---------------------------------------------------------------------------}
procedure mp_pollard_brent_ex(const N: mp_int; var f: mp_int; c: mp_digit; rmax: word);
  {-Find a factor f of N with Brent's version of Pollard rho, f=0 if no factor found;}
  { c = constant in iteration formula, rmax <= 16384 is the maximal r}
label
  leave;
const
  m = 64;
var
  q,x,y,ys,d,mu: mp_int;
  i,j,k,r,rm: integer;
  cancel: boolean;
begin
  {Ref: R.P. Brent, An Improved Monte Carlo Factorization Algorithm,}
  {BIT Vol.20, p.176-184, 1980, available from the author's site at }
  {http://wwwmaths-people.anu.edu.au/~brent/pub/pub051.html         }
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(N) or mp_not_init(f) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_pollard_brent_ex');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Easy out}
  if mp_cmp_mag_d(N,4)=MP_LT then begin
    mp_abs(N,f);
    exit;
  end;

  if rmax>16384 then rm := 16384 else rm := rmax;
  if c=0 then c:=1;

  mp_init6(mu,q,x,y,ys,d);
  if mp_error<>MP_OKAY then exit;

  {Setup Barrett reduction}
  mp_reduce_setup(mu, N);

  {The following code is essentially the algorithm P2'' from}
  {Brent's section 7 with m=64, x0=0, and r limited to rmax.}
  r := 1;
  mp_zero(y);  {y := x0 = 0}
  mp_set1(q);  {accumulate differences}

  repeat
    mp_copy(y,x);
    for i:=1 to r do begin
      {y = y*y + c mod N}
      mp_sqr(y,y);
      mp_add_d(y,c,y);
      mp_reduce(y,N,mu);
    end;
    k := 0;
    repeat
      mp_copy(y,ys);
      j := r-k;
      if j<m then j := m;
      for i:=1 to j do begin
        {y = y*y + c mod N}
        mp_sqr(y,y);
        mp_add_d(y,c,y);
        mp_reduce(y,N,mu);
        {q = q*|y-x| mod N}
        mp_sub(y,x,d); d.sign := MP_ZPOS;
        mp_mul(q,d,q);
        mp_reduce(q,N,mu);
      end;
      if mp_show_progress and (k and $ff = 0) then begin
        if ProgressAssigned then begin
          cancel := false;
          mp_progress(false, k, 0, cancel);
          if cancel then begin
            mp_zero(f);
            goto leave;
          end;
        end
        else if mp_isconsole then write('.');
      end;
      mp_gcd(q,N,f);
      k := k + m;
      if mp_error<>MP_OKAY then goto leave;
    until (k>=r) or (mp_cmp_d(f,1)=MP_GT);
    if r>=rm then break;
    r := 2*r;
  until mp_cmp_d(f,1)=MP_GT;

  if mp_is_eq(f,N) then begin
    {use saved y=ys to backtrack}
    repeat
      {ys = ys*ys + c mod N}
      mp_sqr(ys,ys);
      mp_add_d(ys,c,ys);
      mp_reduce(ys,N,mu);
      mp_sub(ys,x,d); d.sign := MP_ZPOS;
      mp_gcd(d,N,f);
    until (mp_cmp_d(f,1)=MP_GT) or (mp_error<>MP_OKAY);
  end;

  if (mp_cmp_d(f,1)=MP_EQ) or (mp_is_eq(f,N)) then begin
    {f=0 indicates no factor found}
    mp_zero(f);
  end;

leave:
  if mp_error<>MP_OKAY then mp_zero(f);
  mp_clear6(mu,q,x,y,ys,d);
end;


{---------------------------------------------------------------------------}
procedure mp_pollard_brent(const N: mp_int; var f: mp_int);
  {-Find a factor f of N with Brent's version of Pollard rho, f=0 if no factor found; c=1, rmax=8192}
begin
  mp_pollard_brent_ex(N,f,1,8192);
end;


{---------------------------------------------------------------------------}
procedure mp_squfof(const n: mp_int; var f, res: longint);
  {-Find a factor f of n < 2^60 with Shanks' SQUFOF method; n must be}
  { composite. If res=0 then f is a factor; if res<>0, no factor is  }
  { found. mp_squfof is based on Riesel's SQUFOF code. Error values: }
  { res = -1:  general error                                         }
  { res = -2:  n is < 2 or >= 2^60                                   }
  { res = -3:  forward  cycle: no square form found                  }
  { res = -4:  backward cycle: no factor found                       }
  { res = -5:  queue overflow                                        }
  { res >  0:  search period (=res) for CF of sqrt(n) is reached     }
const
  MAXLIST = 30;
  MAX2    = 1000000;
  MAX1    = 2*MAX2;
label
  leave;
var
  q,q0,q1,q2,p1,p2,i,u,r,sqn,sq2sqn: longint;
  k: integer;
  x,y: mp_int;
  List: array[0..MAXLIST] of longint;
  found,usesqr,cancel: boolean;
begin
  {This is based on Riesel's program SQUFOF in [40], p.186-193. Main}
  {changes: Readable formating using while/repeat etc, result codes,}
  {no (direct) floating point routines, check queue overflow.       }
  f := 0;
  res := -1;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_squfof');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if (mp_cmp_d(n,2)=MP_LT) or (mp_bitsize(n) > 60) then begin
    {n is < 2 or >= 2^60}
    res := -2;
    exit;
  end;

  if mp_iseven(n) then begin
    f := 2;
    res := 0;
    exit;
  end;

  mp_init2(x,y);
  if mp_error<>MP_OKAY then exit;

  mp_sqrtrem(n,x,y);
  if not mp_is_longint(x,sqn) then goto leave; {Paranoia: should no happen}
  if not mp_is_longint(y,q1)  then goto leave; {Paranoia: should no happen}

  if q1=0 then begin
    {n is a square}
    f := sqn;
    res := 0;
    goto leave;
  end;

  u := 2*sqn;
  if q1 >= sqn then u := u+1;
  sq2sqn := isqrt32(u);

  q0 := 1;
  p1 := sqn;
  r  := 1;
  List[0] := 0;
  found := false;

  {Continued fraction expansion of sqrt(n)}
  for i:=1 to MAX1 do begin
    if mp_show_progress and (i and $7FFF = 0) then begin
      if ProgressAssigned then begin
        cancel := false;
        mp_progress(false, i, MAX1, cancel);
        if cancel then begin
          break;
        end;
      end
      else if mp_isconsole then write('.');
    end;
    q  := (sqn + p1) div q1;
    p2 := q*q1 - p1;
    q2 := q0 + q*(p1-p2);
    q0 := q1;
    q1 := q2;
    p1 := p2;

    u  := q0;
    if not odd(u) then u := u div 2;
    if (u < sq2sqn) and (u > 1) then begin
      k := List[0];
      if List[0] >= MAXLIST then begin
        res := -5;
        goto leave;
      end;
      inc(k);
      List[0] := k;
      List[k] := u;
    end;

    {Here a small denominator is put in the list}
    if odd(i) and (q1 and 3 < 2) and is_square32ex(q1,r) then begin
      {A square is found}
      usesqr := true;
      for k:=1 to List[0] do begin
        if r=List[k] then begin
          {Square of no use, next i}
          usesqr := false;
          break;
        end;
      end;
      if usesqr then begin
        if r > 1 then begin
          found := true;
          break;
        end
        else begin
          {The period has been searched without finding any useful form}
          res := i;
          goto leave;
        end
      end;
    end;
  end;

  if not found then res := -3
  else begin
    {Computation of the inverse square root of the square form}
    q0 := r;
    p1 := sqn - (sqn - p1) mod r;
    {Compute q1 = (N - p1^2) div q0, division must be exact!}
    mp_set_int(x,p1);
    mp_sqr(x,x);
    mp_sub(n,x,x);
    mp_div_int(x,r,@y,u);
    if (u<>0) or not mp_is_longint(y,q1) then goto leave; {Paranoia: should no happen}
    {Reverse cycle to find a factor}
    for i:=1 to MAX2 do begin
      q  := (sqn + p1) div q1;
      p2 := q*q1 - p1;
      q2 := q0 + q*(p1-p2);
      q0 := q1;
      u  := p1;
      p1 := p2;
      q1 := q2;
      if u = p2 then begin
        {factor found}
        if not odd(q0) then q0 := q0 shr 1;
        f := q0;
        res := 0;
        goto leave;
      end;
    end;
    res := -4;
  end;

leave:
  mp_clear2(x,y);

end;


{---------------------------------------------------------------------------}
procedure mp_fermat_factor(const N: mp_int; var f: mp_int; cnt: longint);
  {-Find a factor f of N with Fermat's method, f=0 if no factor found; cnt: number of tests}
var
  a,b,c: mp_int;
  i: longint;
  cancel: boolean;
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(N) or mp_not_init(f) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_fermat_factor');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Easy outs}
  if mp_cmp_mag_d(N,4)=MP_LT then begin
    mp_abs(N,f);
    if mp_is1(f) then mp_zero(f);
    exit;
  end
  else if mp_iseven(N) then begin
    mp_set(f,2);
    exit;
  end;

  mp_init3(a,b,c);
  if mp_error<>MP_OKAY then exit;

  mp_sqrtrem(n,a,b);
  if mp_is0(b) then begin
    {n=a^2}
    mp_copy(a,f);
  end
  else begin
    for i:=1 to cnt do begin
      if mp_show_progress and (i and $FF = 0) then begin
        if ProgressAssigned then begin
          cancel := false;
          mp_progress(false, i, cnt, cancel);
          if cancel then begin
            break;
          end;
        end
        else if mp_isconsole then write('.');
      end;
      mp_inc(a);
      mp_sqr(a,c);
      mp_sub(c,n,c);
      if mp_is_square2(c,@b) then begin
        {n = (a+b)*(a-b)}
        mp_sub(a,b,c);
        if mp_is1(c) then mp_zero(f) {only trival factor 1}
        else mp_add(a,b,f);
        break;
      end;
      if (MP_Error<>MP_OKAY) or (i=cnt) then begin
        mp_zero(f);
        break;
      end;
    end;
  end;
  mp_clear3(a,b,c);
end;


{---------------------------------------------------------------------------}
procedure mp_holf(const n: mp_int; var f: mp_int; cnt: longint);
  {-Find a factor f of n with Hart's OneLineFactor, f=0 if no factor found; cnt: number of tests}
var
  s,m,a,mu: mp_int;
  i: longint;
  r: mp_digit;
  cancel: boolean;
const
  PM = 480; {Pre-multiplier}
begin
  if mp_error<>MP_OKAY then exit;

  {W.B. Hart, A one line factoring algorithm, 2012, Journal of the   }
  {Australian Mathematical Society, Vol.92, pp. 61-69. Available from}
  {http://wrap.warwick.ac.uk/54707/1/WRAP_Hart_S1446788712000146a.pdf}

  {W.B. Hart, A One Line Factoring Algorithm, 2009}
  {http://selmer.warwick.ac.uk/onelinefactor.pdf  }

  {Test cases for HOLF: structure nx=x*nextprime(q*x), x prime, q small}
  {n1=658777688124265765008727722820674755897063681968221644850784763800801}
  {n2=93899836061477776150307035713278480058525235649004392384791771141}
  {n3=10341326324326594086808731594429812699683531021516364916988043}
  {n4=3756150573033065890153057446046108450497850840115161}
  {n5=66951327270230683487974977703880821429755937319}
  {n6=63599406461425156216302869732000579708026183660247, cnt=30000, auch ECM}

  {$ifdef MPC_ArgCheck}
    if mp_not_init(N) or mp_not_init(f) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_holf');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Check factors 2,3,5,7}
  mp_mod_d(n,210,r);
  if r and 1 = 0 then begin
    mp_set(f,2);
    exit;
  end
  else if r mod 3 = 0 then begin
    mp_set(f,3);
    exit;
  end
  else if r mod 5 = 0 then begin
    mp_set(f,5);
    exit;
  end
  else if r mod 7 = 0 then begin
    mp_set(f,7);
    exit;
  end;

  mp_init4(a,s,m,mu);
  if mp_error<>MP_OKAY then exit;

  {Pre-multiply a = PM*|n|}
  mp_mul_d(n,PM,a);
  mp_abs(a,a);

  {Setup Barrett reduction for a}
  mp_reduce_setup(mu, a);

  for i:=1 to cnt do begin
    if mp_error<>MP_OKAY then break;
    if mp_show_progress and (i and $FF = 0) then begin
      if ProgressAssigned then begin
        cancel := false;
        mp_progress(false, i, cnt, cancel);
        if cancel then begin
          mp_zero(f);
          break;
        end;
      end
      else if mp_isconsole then write('.');
    end;
    mp_mul_int(a,i,s);
    mp_sqrtrem(s,s,m);
    if mp_is0(m) then begin
      mp_gcd(n,s,f);
      break;
    end;
    {s = ceil(sqrt(PM*n*i))}
    {m = s^2 mod a}
    mp_inc(s);
    mp_sqr(s,m);
    mp_reduce(m,a,mu);
    if mp_is_square2(m,@f) then begin
      mp_sub(s,f,m);
      mp_gcd(n,m,f);
      break;
    end;
  end;
  if mp_is1(f) or (mp_error<>MP_OKAY) then mp_zero(f);
  mp_clear4(a,s,m,mu);
end;


(*
{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}
*)

begin

  mp_show_progress   := false;  {rho, (p+1), (p-1) after one accum cycle}
  mp_progress        := nil;
  mp_max_small_sqr   := sqr(longint(MP_DIGIT_MAX and $7FFF));

  {$ifndef BIT16}
     mp_isconsole := isconsole;
  {$else}
    {$ifdef WINDOWS}
      mp_isconsole := false;
    {$else}
      mp_isconsole := true;
    {$endif}
  {$endif}
end.
