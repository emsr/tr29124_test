unit mp_cmplx;

{Multi precision complex floating point arithmetic routines}

interface

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
{$endif}

uses
  BTypes, mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION   :  Multi precision complex floating point arithmetic routines

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  heap

 DISPLAY MODE  :  ---

 REFERENCES    : [41] [HMF]: M. Abramowitz, I.A. Stegun. Handbook of Mathematical
                      Functions. New York, 1970, http://www.math.sfu.ca/~cbm/aands/
                 [42] W. Kahan, "Branch Cuts for Complex Elementary Functions, or Much Ado
                      About Nothing's Sign Bit", in The State of Art in Numerical Analysis,
                      ed. by A. Iserles and M.J.D. Powell, 1987, pp. 165-211.
                      Available as http://people.freebsd.org/~das/kahan86branch.pdf

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------

 0.00.10  15.12.13  W.Ehrhardt  Basic type definitions, mpc_init[*]
 0.00.11  16.12.13  we          mpc_clear[*]
 0.00.12  02.01.14  we          mpc_set0/1/i/dbl/ext/mpf
 0.00.13  02.01.14  we          mpc_setp_mpf
 0.00.14  08.02.14  we          mpc_add*, mpc_sub*, mpc_chs
 0.00.15  09.02.14  we          mpc_abs/2, mpc_mul*, mpc_div*
 0.00.16  10.02.14  we          mpc_checksum, mpc_copy/p, mpc_is0/1/i
 0.00.16  10.02.14  we          mpc_checksum, mpc_copy/p, mpc_is0/1/i
 0.00.17  11.02.14  we          mpc_sqrt, mpc_sqr
 0.00.18  11.02.14  we          mpc_arg, mpc_ln
 0.00.19  11.02.14  we          mpc_cis, mpc_exp
 0.00.20  12.02.14  we          mpc_sin/cos/sincos
 0.00.21  12.02.14  we          mpc_sinh/cosh/sincosh
 0.00.21  13.02.14  we          mpc_coth, mpc_tanh
 0.00.22  13.02.14  we          mpc_pow, mpc_cot, mpc_tan
 0.00.23  14.02.14  we          mpc_arctan/h
 0.00.24  15.02.14  we          mpc_arccos/h, mpc_arcsin/h
 0.00.25  16.02.14  we          mpc_nroot/1
 0.00.26  18.02.14  we          mpc_is1a, improved mpc_coth

 1.27.00  19.02.14  we          mp_complex types moved to unit mp_types

 1.28.00  31.03.14  we          mpc_csc/h, mpc_sec/h
 1.28.01  01.04.14  we          mpc_arccot/c
 1.28.02  01.04.14  we          mpc_arccoth/c
 1.28.03  02.04.14  we          mpc_arcsec/h
 1.28.04  02.04.14  we          mpc_arccsc/h
 1.28.05  03.04.14  we          fix mpc_tanh: use s_mpf_ldx(a.re)

 1.31.00  24.11.14  we          mpc_agm1
 1.31.01  25.11.14  we          mpc_agm

 1.33.00  09.10.16  we          mpc_exp2, mpc_exp10
 1.33.01  13.10.16  we          fix arg check in mpc_agm1
 1.33.02  14.10.16  we          mpc_coth without vers()
 1.33.03  14.10.16  we          mpc_ln1p
 1.33.04  15.10.16  we          mpc_expm1

 1.38.00  26.06.18  we          mpc_??_ext changed to mpc_??_dbl

 1.39.00  12.11.18  we          mpc_log10
 1.39.01  13.11.18  we          mpc_is_ia
 1.39.02  13.11.18  we          fix mpc_ln1p  for small imaginary a
 1.39.03  14.11.18  we          fix mpc_expm1 for small imaginary a
 1.39.04  14.11.18  we          fix mpc_arctanh for real |a| > 1
 1.39.05  16.11.18  we          mpc_arctanh: adjust sign on the branch cut
**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2013-2018 Wolfgang Ehrhardt

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


procedure mpc_abs(const a: mp_complex; var b: mp_float);
  {-Calculate the absolute value, b = |a|}

procedure mpc_abs2(const a: mp_complex; var b: mp_float);
  {-Calculate the squared absolute value, b = |a|^2}

procedure mpc_add(const a,b: mp_complex; var c: mp_complex);
  {-Calculate c = a+b}

procedure mpc_add_dbl(const a: mp_complex; b: double; var c: mp_complex);
  {-Calculate c = a+b}

procedure mpc_add_mpf(const a: mp_complex; const b: mp_float; var c: mp_complex);
  {-Calculate c = a+b}

procedure mpc_agm(const x,y: mp_complex; var w: mp_complex);
  {-Calculate the 'optimal' arithmetic-geometric mean w = AGM(x,y)}

procedure mpc_agm1(const z: mp_complex; var w: mp_complex);
  {-Calculate the 'optimal' arithmetic-geometric mean w = AGM(1,z)}

procedure mpc_arccos(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular cosine b = arccos(a)}

procedure mpc_arccosh(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic cosine b = arccosh(a)}

procedure mpc_arccot(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular cotangent b = arccot(a) = arctan(1/a)}

procedure mpc_arccotc(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular cotangent b = arccotc(a) = Pi/2 - arctan(a)}

procedure mpc_arccoth(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic cotangent b = arccoth(a) = arctanh(1/a)}

procedure mpc_arccothc(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic cotangent b = arccothc(a) = arctanh(a) + i*Pi/2}

procedure mpc_arccsc(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular cosecant b = arccsc(a) = arcsin(1/a)}

procedure mpc_arccsch(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic cosecant b = arccsch(a) = arcsinh(1/a)}

procedure mpc_arcsec(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular secant b = arcssec(a) = arccos(1/a)}

procedure mpc_arcsech(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic secant b = arcssech(a) = arccosh(1/a)}

procedure mpc_arcsin(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular sine b = arcsin(a)}

procedure mpc_arcsinh(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic sine b = arcsinh(a)}

procedure mpc_arctan(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular tangent w = arctan(z)}

procedure mpc_arctanh(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic tangent b = arctanh(a)}

procedure mpc_arg(const a: mp_complex; var b: mp_float);
  {-Calculate the principle value of the argument: b = arg(a) = arctan2(a.im, a.re)}

function  mpc_checksum(const a: mp_complex): longint;
  {-Return a checksum for a, -1 if mp_error<>MP_OKAY, -2 if not initialized}

procedure mpc_chs(const a: mp_complex; var b: mp_complex);
  {-Change sign, b = -a}

procedure mpc_cis(const x: mp_float; var a: mp_complex);
  {-Calculate a = exp(i*x) = cos(x) + i*sin(x)}

procedure mpc_clear(var a: mp_complex);
  {-Clear an mp_complex}

procedure mpc_clear2(var a,b: mp_complex);
  {-Clear 2 mp_complex}

procedure mpc_clear3(var a,b,c: mp_complex);
  {-Clear 3 mp_complex}

procedure mpc_clear4(var a,b,c,d: mp_complex);
  {-Clear 4 mp_complex}

procedure mpc_clear5(var a,b,c,d,e: mp_complex);
  {-Clear 5 mp_complex}

procedure mpc_conj(const a: mp_complex; var b: mp_complex);
  {-Return the complex conjugate b = a.re - i*a.im}

procedure mpc_copy(const a: mp_complex; var b: mp_complex);
  {-Copy a to b including bitprecs}

procedure mpc_copyp(const a: mp_complex; var b: mp_complex);
  {-Copy a to b, preserve b's bitprecs}

procedure mpc_cos(const a: mp_complex; var b: mp_complex);
  {-Calculate b = cos(a)}

procedure mpc_cosh(const a: mp_complex; var b: mp_complex);
  {-Calculate b = cosh(a)}

procedure mpc_cot(const a: mp_complex; var b: mp_complex);
  {-Calculate b = cot(a)}

procedure mpc_coth(const a: mp_complex; var b: mp_complex);
  {-Calculate b = coth(a)}

procedure mpc_csc(const a: mp_complex; var b: mp_complex);
  {-Calculate the complex circular cosecant b = csc(a) = 1/sin(a)}

procedure mpc_csch(const a: mp_complex; var b: mp_complex);
  {-Calculate the complex hyperbolic cosecant b = csch(a) = 1/sinh(a)}

procedure mpc_div(const a,b: mp_complex; var c: mp_complex);
  {-Calculate c = a/b}

procedure mpc_div_dbl(const a: mp_complex; b: double; var c: mp_complex);
  {-Calculate c = a/b}

procedure mpc_div_mpf(const a: mp_complex; const b: mp_float; var c: mp_complex);
  {-Calculate c = a/b}

procedure mpc_exch(var a,b: mp_complex);
  {-Exchange two mp_complexes (including bitprec)}

procedure mpc_exp(const a: mp_complex; var b: mp_complex);
  {-Calculate b = exp(a)}

procedure mpc_exp2(const a: mp_complex; var b: mp_complex);
  {-Calculate b = 2^a = exp(a*ln(2))}

procedure mpc_exp10(const a: mp_complex; var b: mp_complex);
  {-Calculate b = 10^a = exp(a*ln(10))}

procedure mpc_expm1(const a: mp_complex; var b: mp_complex);
  {-Calculate b = exp(a)-1}

procedure mpc_init(var a: mp_complex);
  {-Initialize an mp_complex with default precision}

procedure mpc_init2(var a,b: mp_complex);
  {-Initialize two mp_complexes with default precision}

procedure mpc_init3(var a,b,c: mp_complex);
  {-Initialize 3 mp_complexes with default precision}

procedure mpc_init4(var a,b,c,d: mp_complex);
  {-Initialize 4 mp_complexes with default precision}

procedure mpc_init5(var a,b,c,d,e: mp_complex);
  {-Initialize 5 mp_complexes with default precision}

procedure mpc_initp(var a: mp_complex; prec: longint);
  {-Initialize an mp_complex with bit precision prec}

procedure mpc_initp2(var a,b: mp_complex; prec: longint);
  {-Initialize two mp_complexes with bit precision prec}

procedure mpc_initp3(var a,b,c: mp_complex; prec: longint);
  {-Initialize 3 mp_complexes with bit precision prec}

procedure mpc_initp4(var a,b,c,d: mp_complex; prec: longint);
  {-Initialize 4 mp_complexes with bit precision prec}

procedure mpc_initp5(var a,b,c,d,e: mp_complex; prec: longint);
  {-Initialize 5 mp_complexes with bit precision prec}

procedure mpc_initp_multi_p(var pv: array of pmp_complex; prec: longint);
  {-Initialize with bit precision prec a list of mp_complexes given as a pointer}
  { vector; on error the already initialized mp_complexes will be cleared}

procedure mpc_inv(const a: mp_complex; var b: mp_complex);
  {-Calculate b = 1/a}

function  mpc_is0(const a: mp_complex): boolean;
  {-Return true if a=0}

function  mpc_is1(const a: mp_complex): boolean;
  {-Return true if a=1}

function  mpc_is1a(const a: mp_complex): boolean;
  {-Return true if a = +1 or -1}

function  mpc_is_i(const a: mp_complex): boolean;
  {-Return true if a=I}

function mpc_is_ia(const a: mp_complex): boolean;
  {-Return true if a=I or a= -I}

procedure mpc_ln(const a: mp_complex; var b: mp_complex);
  {-Calculate the natural logarithm b = ln(a); principal branch ln(|a|) + i*arg(a), accurate near |a|=1}

procedure mpc_ln1p(const a: mp_complex; var b: mp_complex);
  {-Calculate the natural logarithm b = ln(1+a)}

procedure mpc_log10(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal branch of the base 10 logarithm of a, b=ln(z)/ln(10)}

procedure mpc_mul(const a,b: mp_complex; var c: mp_complex);
  {-Calculate c = a*b}

procedure mpc_mul_2k(const a: mp_complex; k: longint; var b: mp_complex);
  {-Calculate b = a*2^k}

procedure mpc_mul_dbl(const a: mp_complex; b: double; var c: mp_complex);
  {-Calculate c = a*b}

procedure mpc_mul_mpf(const a: mp_complex; const b: mp_float; var c: mp_complex);
  {-Calculate c = a*b}

function  mpc_not_init(const a: mp_complex): boolean;
  {-Sanity check if a is initialized, does not catch all cases!}

procedure mpc_pow(const a,b: mp_complex; var c: mp_complex);
  {-Calculate the principal value of the complex power c = a^b = exp(b*ln(a))}

procedure mpc_nroot(const a: mp_complex; n: longint; var b: mp_complex);
  {-Calculate the nth principal root b = a^(1/n) = exp(ln(a)/n)}

procedure mpc_nroot1(n: longint; var a: mp_complex);
  {-Calculate the principal nth root of unity a = exp(2*Pi*i/n)}

procedure mpc_sec(const a: mp_complex; var b: mp_complex);
  {-Calculate the complex circular secant b = sec(a) = 1/cos(a)}

procedure mpc_sech(const a: mp_complex; var b: mp_complex);
  {-Calculate the complex hyperbolic secant b = sech(a) = 1/cosh(a)}

procedure mpc_set0(var a: mp_complex);
  {-Set a=0}

procedure mpc_set1(var a: mp_complex);
  {-Set a=1}

procedure mpc_seti(var a: mp_complex);
  {-Set a=i}

procedure mpc_set_dbl(var a: mp_complex; x,y: double);
  {-Set a = x + iy with x,y double, error if x,y are NAN or INF}

{$ifndef EXT64}
procedure mpc_set_ext(var a: mp_complex; x,y: extended);
  {-Set a = x + iy with x,y extended, error if x,y are NAN or INF}
{$endif}

procedure mpc_set_mpf(var a: mp_complex; const x,y: mp_float);
  {-Set a = x + iy}

procedure mpc_setp_mpf(var a: mp_complex; const x,y: mp_float);
  {-Set a = x + iy, bitprecs of a are preserved}

procedure mpc_sin(const a: mp_complex; var b: mp_complex);
  {-Calculate b = sin(a)}

procedure mpc_sincos(const a: mp_complex; var s,c: mp_complex);
  {-Calculate s = sin(a), c = cos(a)}

procedure mpc_sinh(const a: mp_complex; var b: mp_complex);
  {-Calculate b = sinh(a)}

procedure mpc_sinhcosh(const a: mp_complex; var s,c: mp_complex);
  {-Calculate s = sinh(a), c = cosh(a)}

procedure mpc_sqr(const a: mp_complex; var b: mp_complex);
  {-Calculate b = a*a}

procedure mpc_sqrt(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal square root b = sqrt(a)}

procedure mpc_sub(const a,b: mp_complex; var c: mp_complex);
  {-Calculate c = a-b}

procedure mpc_sub_dbl(const a: mp_complex; b: double; var c: mp_complex);
  {-Calculate c = a-b}

procedure mpc_sub_mpf(const a: mp_complex; const b: mp_float; var c: mp_complex);
  {-Calculate c = a-b}

procedure mpc_tan(const a: mp_complex; var b: mp_complex);
  {-Calculate b = tan(a)}

procedure mpc_tanh(const a: mp_complex; var b: mp_complex);
  {-Calculate b = tanh(a)}

procedure mpc_xdivc(a: double; const b: mp_complex; var c: mp_complex);
  {-Calculate c = a/b}

procedure s_mpc_mul(const a,b: mp_complex; bc: boolean; var c: mp_complex);
  {-Calculate c = a*b or c=a*conj(b) if bc}


implementation


uses
  mp_base, mp_real;


{---------------------------------------------------------------------------}
procedure mpc_abs(const a: mp_complex; var b: mp_float);
  {-Calculate the absolute value, b = |a|}
begin
  mpf_hypot(a.re, a.im, b);
end;


{---------------------------------------------------------------------------}
procedure mpc_abs2(const a: mp_complex; var b: mp_float);
  {-Calculate the squared absolute value, b = |a|^2}
var
  x: mp_float;
begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b);
  {$endif}
  mpf_initp(x,b.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_sqr(a.re,x);
    mpf_sqr(a.im,b);
    mpf_add(x,b,b);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_add(const a,b: mp_complex; var c: mp_complex);
  {-Calculate c = a+b}
begin
  mpf_add(a.re, b.re, c.re);
  mpf_add(a.im, b.im, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_add_dbl(const a: mp_complex; b: double; var c: mp_complex);
  {-Calculate c = a+b}
begin
  mpf_add_dbl(a.re, b, c.re);
  mpf_copy(a.im, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_add_mpf(const a: mp_complex; const b: mp_float; var c: mp_complex);
  {-Calculate c = a+b}
begin
  mpf_add(a.re, b, c.re);
  mpf_copy(a.im, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_agm1(const z: mp_complex; var w: mp_complex);
  {-Calculate the 'optimal' arithmetic-geometric mean w = AGM(1,z)}
var
  a,b,u: mp_complex;
  lim,cnt,prec: longint;
  r: integer;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(z) or mpc_not_init(w) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpc_agm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mpf_is0(z.im) and (mpf_is0(z.re) or (mpf_is1a(z.re) and s_mpf_is_neg(z.re))) then begin
    mpc_set0(w);
    exit;
  end;

  prec := w.re.bitprec+8;
  lim  := prec div 2;
  mpc_initp3(a,b,u,prec);
  if mp_error<>MP_OKAY then exit;

  {Compute suitable starting values a, b for optimal AGM, see Pari/GP }
  {source code and the discussion in the pari-dev thread 'Complex AGM'}
  {http://pari.math.u-bordeaux.fr/archives/pari-dev-1202/msg00045.html}
  {or http://comments.gmane.org/gmane.comp.mathematics.pari.devel/3543}
  if s_mpf_is_neg(z.re) then begin
    if s_mpf_is_neg(z.im) then begin
      {a := +I*a}
      r := -1;
      mpf_mul_2k(z.im,-1, a.re);
      s_mpf_chs(a.re);
      mpf_add_dbl(z.re, 1.0, a.im);
      mpf_mul_2k(a.im,-1, a.im);
    end
    else begin
      {a := -I*a}
      r := 1;
      mpf_mul_2k(z.im,-1, a.re);
      mpf_add_dbl(z.re, 1.0, a.im);
      mpf_mul_dbl(a.im, -0.5, a.im);
    end;
    mpc_chs(z,b);
    mpc_sqrt(b,b);
  end
  else begin
    {no rotation}
    r := 0;
    mpc_mul_dbl(z, 0.5, a);
    mpf_add_dbl(a.re, 0.5, a.re);
    mpc_sqrt(z,b);
  end;

  {Here a.re >= 0 and b.re >= 0, do standard AGM iteration}
  for cnt:=0 to lim do begin
    {u = a*b}
    mpc_mul(a,b,u);
    {a = (a+b)/2}
    mpc_add(a,b,a);
    mpc_sqrt(u,b);
    {b = sqrt(a*b)}
    mpc_mul_2k(a,-1,a);
    {check difference between arith and geo part}
    mpc_sub(b,a,u);
    if mpc_is0(u) then break
    else begin
      mpc_abs(u,u.re);
      mpc_abs(a,u.im);
      if lim < s_mpf_ldx(u.im)-s_mpf_ldx(u.re) then break;
    end;
  end;

  {iteration difference is < prec/2, do final step to make it < prec}
  mpc_add(a,b,a);
  mpc_mul_2k(a,-1,a);

  {undo rotation}
  case r of
      1: begin
           {w=a*I}
           mpf_copyp(a.im, w.re);
           mpf_copyp(a.re, w.im);
           s_mpf_chs(w.re);
         end;
     -1: begin
           {w=-a*I}
           mpf_copyp(a.im, w.re);
           mpf_copyp(a.re, w.im);
           s_mpf_chs(w.im);
         end;
    else begin
           {w=a}
           mpf_copyp(a.re, w.re);
           mpf_copyp(a.im, w.im);
         end;
  end; {case}

  mpc_clear3(a,b,u);
end;


{---------------------------------------------------------------------------}
procedure mpc_agm(const x,y: mp_complex; var w: mp_complex);
  {-Calculate the 'optimal' arithmetic-geometric mean w = AGM(x,y)}
var
  u,v: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(w.re);
  {$endif}

  if mpc_is0(x) or mpc_is0(y) then begin
    mpc_set0(w);
    exit;
  end;

  mpc_initp2(u,v,w.re.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpc_abs(x,u.re);
    mpc_abs(y,u.im);
    if mpf_is_gt(u.re,u.im) then begin
      {|x| > |y|}
      mpc_div(y,x,v);
      mpc_copyp(x,u);
    end
    else begin
      {|x| <= |y|}
      mpc_div(x,y,v);
      mpc_copyp(y,u);
    end;
    mpc_agm1(v,v);
    mpc_mul(u,v,w);
    mpc_clear2(u,v);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arccos(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular cosine b = arccos(a)}
var
  zm,zp: mp_complex;
  x,y: mp_float;
  prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {Ref: Kahan[42], procedure CACOS}
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  prec := b.re.bitprec+8;
  mpc_initp2(zm,zp, prec);
  if mp_error=MP_OKAY then begin
    mpf_initp2(x,y,prec);
    if mp_error=MP_OKAY then begin
      {zp = sqrt(1+a)}
      mpc_add_dbl(a,1.0,zp);
      mpc_sqrt(zp,zp);
      {zm = sqrt(1-a)}
      mpc_chs(a,zm);
      s_mpf_inc1(zm.re);
      mpc_sqrt(zm,zm);
      {b.re = 2.0*arctan2(re(sqrt(1-a), re(sqrt(1+a)))) }
      mpf_arctan2(zm.re, zp.re, x);
      mpf_mul_2k(x,1,b.re);
      {b.im = arcsinh(im(sqrt(1+conj(a))*sqrt(1-a)))    }
      {use im(sqrt(1+conj(a))) = -im(sqrt(1+a))         }
      mpf_mul(zp.re, zm.im, y);
      mpf_mul(zp.im, zm.re, x);
      mpf_sub(y,x,y);
      mpf_arcsinh(y,b.im);
      mpf_clear2(x,y);
    end;
    mpc_clear2(zm,zp);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arccosh(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic cosine b = arccosh(a)}
var
  zm,zp: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  {Ref: Kahan[42], procedure CACOSH}
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpc_initp2(zm,zp, b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    {zp = sqrt(a+1)}
    mpc_add_dbl(a,1.0,zp);
    mpc_sqrt(zp,zp);
    {zm = sqrt(a-1)}
    mpc_sub_dbl(a,1.0,zm);
    mpc_sqrt(zm,zm);
    {b.im = 2.0*arctan2(im(sqrt(a-1.0), re(sqrt(a+1.0)))) }
    mpf_arctan2(zm.im, zp.re, b.im);
    s_mpf_incexp(b.im,1);
    {b.re = arcsinh(re(sqrt(conj(a)-1.0)*sqrt(a+1.0)))    }
    {use im(sqrt(conj(a)-1)) = -im(sqrt(a-1))             }
    mpf_mul(zp.re, zm.re, zm.re);
    mpf_mul(zp.im, zm.im, zm.im);
    mpf_add(zm.re, zm.im, zm.re);
    mpf_arcsinh(zm.re,b.re);
    mpc_clear2(zm,zp);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arcsin(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular sine b = arcsin(a)}
var
  zm,zp: mp_complex;
  x,y: mp_float;
  prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {Ref: Kahan[42], procedure CASIN}
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  prec := b.re.bitprec+8;
  mpc_initp2(zm,zp, prec);
  if mp_error=MP_OKAY then begin
    mpf_initp2(x,y,prec);
    if mp_error=MP_OKAY then begin
      {zp = sqrt(1+a)}
      mpc_add_dbl(a,1.0,zp);
      mpc_sqrt(zp,zp);
      {zm = sqrt(1-a)}
      mpc_chs(a,zm);
      s_mpf_inc1(zm.re);
      mpc_sqrt(zm,zm);
      {y = im(sqrt(1-conj(a))*sqrt(1+a))}
      mpf_mul(zm.re, zp.im, y);
      mpf_mul(zp.re, zm.im, x);
      mpf_sub(y,x,y);
      {x = re(sqrt(1-a)*sqrt(1+a))}
      mpf_mul(zp.re, zm.re, x);
      mpf_mul(zp.im, zm.im, zm.im);
      mpf_sub(x,zm.im,x);
      {b.re = arctan2(re(a), re(sqrt(1-a)*sqrt(1+a))) }
      {b.im = arcsinh(im(sqrt(1-conj(a))*sqrt(1+a)))  }
      mpf_arctan2(a.re, x, b.re);
      mpf_arcsinh(y, b.im);
      mpf_clear2(x,y);
    end;
    mpc_clear2(zm,zp);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arcsinh(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic sine b = arcsinh(a)}
begin
  {arcsinh(a) = -i*arcsin(i*a)}
  mpc_copyp(a,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
  mpc_arcsin(b,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_arctan(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular tangent w = arctan(z)}
begin
  {Ref HMF[41], 4.4.22: arctan(z) = -i*arctanh(iz)}
  mpc_copyp(a,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
  mpc_arctanh(b,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_arctanh(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic tangent b = arctanh(a)}
var
  x,x1,y,u,v: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if s_mpf_is0(a.re) then begin
    {arctanh(iy)=i*arctan(y);}
    mpf_set0(b.re);
    mpf_arctan(a.im, b.im);
    exit;
  end;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpf_initp5(x,x1,y,u,v, b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    {Ref: Kahan[42], procedure CATANH, normal case with rho=0}
    mpf_abs(a.re, x);
    mpf_abs(a.im, y);
    {v = y^2}
    mpf_sqr(y,v);
    {x1 = 1-x}
    mpf_chs(x,x1);
    s_mpf_inc1(x1);
    {u := ln1p(4.0*x/(sqr(x1) + v));}
    mpf_sqr(x1,u);
    mpf_add(u,v,u);
    mpf_div(x,u,u);
    s_mpf_incexp(u,2);
    mpf_ln1p(u,u);
    {v := arctan2(2.0*y, x1*(1.0+x) - v)}
    s_mpf_incexp(y,1);
    s_mpf_inc1(x);
    mpf_mul(x,x1,x);
    mpf_sub(x,v,x);
    if mpf_is0(y) and (mpf_cmp_dbl(a.re,1) > 0) then begin
      {Adjust sign on the branch cut}
      mpf_set_pi(v);
      mpf_chs(v,v);
    end
    else mpf_arctan2(y,x,v);
    {b.re := copysign(0.25,a.re)*u;}
    {b.im := copysign(0.5 ,a.im)*v;}
    if s_mpf_is_neg(a.re) then s_mpf_chs(u);
    if s_mpf_is_neg(a.im) then s_mpf_chs(v);
    mpf_mul_2k(u,-2,b.re);
    mpf_mul_2k(v,-1,b.im);
    mpf_clear5(x,x1,y,u,v);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arccoth(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic cotangent b = arccoth(a) = arctanh(1/a)}
var
  z: mp_complex;
  lr,prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) or mpc_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpc_arccoth');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  prec := b.re.bitprec+8;
  if s_mpf_ldx(a.im) < 0 then begin
    {|a.im| < 0.5}
    lr := s_mpf_ldx(a.re);
    if (lr >=0 ) and (lr <= 1) then begin
      {0.5 <= |a.re| < 2: add some more precision near poles at +- 1}
      inc(prec,8);
    end;
  end;
  mpc_initp(z, prec);
  if mp_error=MP_OKAY then begin
    mpc_inv(a,z);
    mpc_arctanh(z,z);
    mpc_copyp(z,b);
    mpc_clear(z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arccothc(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic cotangent b = arccothc(a) = arctanh(a) + i*Pi/2}
var
  z: mp_complex;
  prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  if s_mpf_is_neg(a.im) then begin
    {if z.im < 0 then  arccothc(a) = arccoth(a)}
    mpc_arccoth(a,b);
  end
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(b.re);
    {$endif}
    prec := b.re.bitprec+8;
    mpc_initp(z, prec);
    if mp_error=MP_OKAY then begin
      mpc_arctanh(a,z);
      mpf_copyp(z.re,b.re);
      mpf_set_pi2k(z.re,-1);
      mpf_add(z.re, z.im, z.im);
      mpf_copyp(z.im,b.im);
      mpc_clear(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arccot(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular cotangent b = arccot(a) = arctan(1/a)}
begin
  {arccot(z) = i*arccoth(i*z))}
  mpc_copyp(a,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
  mpc_arccoth(b,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
end;


{---------------------------------------------------------------------------}
procedure mpc_arccotc(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular cotangent b = arccotc(a) = Pi/2 - arctan(a)}
var
  z: mp_complex;
  prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  if s_mpf_is_gt0(a.re) then begin
    {if z.re < 0 then  arccotc(a) = arccot(a)}
    mpc_arccot(a,b);
  end
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(b.re);
    {$endif}
    prec := b.re.bitprec+8;
    mpc_initp(z, prec);
    if mp_error=MP_OKAY then begin
      mpc_arctan(a,z);
      mpf_chs(z.im,b.im);
      mpf_set_pi2k(z.im,-1);
      mpf_sub(z.im, z.re, z.re);
      mpf_copyp(z.re,b.re);
      mpc_clear(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arccsc(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular cosecant b = arccsc(a) = arcsin(1/a)}
var
  z: mp_complex;
  lr,prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) or mpc_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpc_arccsc');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  prec := b.re.bitprec+8;
  if s_mpf_ldx(a.im) < 0 then begin
    {|a.im| < 0.5}
    lr := s_mpf_ldx(a.re);
    if (lr >=0 ) and (lr <= 1) then begin
      {0.5 <= |a.re| < 2: add some more precision near z = +- 1}
      inc(prec,8);
    end;
  end;
  mpc_initp(z, prec);
  if mp_error=MP_OKAY then begin
    mpc_inv(a,z);
    mpc_arcsin(z,z);
    mpc_copyp(z,b);
    mpc_clear(z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arccsch(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic cosecant b = arccsch(a) = arcsinh(1/a)}
begin
  if mp_error<>MP_OKAY then exit;
  {Use arccsch(z) = i*arccsc(i*z)}
  mpc_copyp(a,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
  mpc_arccsc(b,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
end;


{---------------------------------------------------------------------------}
procedure mpc_arcsec(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse circular secant b = arcssec(a) = arccos(1/a)}
var
  z: mp_complex;
  lr,prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) or mpc_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpc_arcsec');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  prec := b.re.bitprec+8;
  if s_mpf_ldx(a.im) < 0 then begin
    {|a.im| < 0.5}
    lr := s_mpf_ldx(a.re);
    if (lr >=0 ) and (lr <= 1) then begin
      {0.5 <= |a.re| < 2: add some more precision near z = +- 1}
      inc(prec,8);
    end;
  end;
  mpc_initp(z, prec);
  if mp_error=MP_OKAY then begin
    mpc_inv(a,z);
    mpc_arccos(z,z);
    mpc_copyp(z,b);
    mpc_clear(z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arcsech(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal value of the complex inverse hyperbolic secant b = arcssech(a) = arccosh(1/a)}
var
  z: mp_complex;
  lr,prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) or mpc_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpc_arcsech');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  prec := b.re.bitprec+8;
  if s_mpf_ldx(a.im) < 0 then begin
    {|a.im| < 0.5}
    lr := s_mpf_ldx(a.re);
    if (lr >=0 ) and (lr <= 1) then begin
      {0.5 <= |a.re| < 2: add some more precision near z = +- 1}
      inc(prec,8);
    end;
  end;
  mpc_initp(z, prec);
  if mp_error=MP_OKAY then begin
    mpc_inv(a,z);
    mpc_arccosh(z,z);
    mpc_copyp(z,b);
    mpc_clear(z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_arg(const a: mp_complex; var b: mp_float);
  {-Calculate the principle value of the argument: b = arg(a) = arctan2(a.im, a.re)}
begin
  mpf_arctan2(a.im, a.re, b);
end;


{---------------------------------------------------------------------------}
function mpc_checksum(const a: mp_complex): longint;
  {-Return a checksum for a, -1 if mp_error<>MP_OKAY, -2 if not initialized}
var
  adler: longint;
begin
  mpc_checksum := -1;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) then begin
      mpc_checksum := -2;
      exit;
    end;
  {$endif}
  with a.re do begin
    adler := mp_checksum(mantissa);
    s_mp_checksum(adler,@exponent,sizeof(exponent));
    s_mp_checksum(adler,@bitprec, sizeof(bitprec));
  end;
  with a.im do begin
    s_mp_checksum(adler,@mantissa,sizeof(mantissa));
    s_mp_checksum(adler,@exponent,sizeof(exponent));
    s_mp_checksum(adler,@bitprec, sizeof(bitprec));
  end;
  mpc_checksum := adler;
end;


{---------------------------------------------------------------------------}
procedure mpc_chs(const a: mp_complex; var b: mp_complex);
  {-Change sign, b = -a}
begin
  mpf_chs(a.re, b.re);
  mpf_chs(a.im, b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_cis(const x: mp_float; var a: mp_complex);
  {-Calculate a = exp(i*x) = cos(x) + i*sin(x)}
begin
  mpf_sincos(x, a.im, a.re);
end;


{---------------------------------------------------------------------------}
procedure mpc_clear(var a: mp_complex);
  {-Clear an mp_complex}
begin
  mpf_clear(a.re);
  mpf_clear(a.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_clear2(var a,b: mp_complex);
  {-Clear 2 mp_complex}
begin
  mpf_clear(a.re);
  mpf_clear(a.im);
  mpf_clear(b.re);
  mpf_clear(b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_clear3(var a,b,c: mp_complex);
  {-Clear 3 mp_complex}
begin
  mpc_clear2(a,b);
  mpc_clear(c);
end;


{---------------------------------------------------------------------------}
procedure mpc_clear4(var a,b,c,d: mp_complex);
  {-Clear 4 mp_complex}
begin
  mpc_clear2(a,b);
  mpc_clear2(c,d);
end;


{---------------------------------------------------------------------------}
procedure mpc_clear5(var a,b,c,d,e: mp_complex);
  {-Clear 5 mp_complex}
begin
  mpc_clear2(a,b);
  mpc_clear2(c,d);
  mpc_clear(e);
end;


{---------------------------------------------------------------------------}
procedure mpc_conj(const a: mp_complex; var b: mp_complex);
  {-Return the complex conjugate b = a.re - i*a.im}
begin
  mpf_copyp(a.re,b.re);
  mpf_chs(a.im,b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_copy(const a: mp_complex; var b: mp_complex);
  {-Copy a to b including bitprecs}
begin
  mpf_copy(a.re, b.re);
  mpf_copy(a.im, b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_copyp(const a: mp_complex; var b: mp_complex);
  {-Copy a to b, preserve b's bitprecs}
begin
  mpf_copyp(a.re, b.re);
  mpf_copyp(a.im, b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_cos(const a: mp_complex; var b: mp_complex);
  {-Calculate b = cos(a)}
var
  st,ct,sh,ch: mp_float;
begin
  if mpc_is0(a) then mpc_set1(b)
  else begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
    mpf_initp4(st,ct,sh,ch,b.re.bitprec+8);
    if mp_error=MP_OKAY then begin
      mpf_sinhcosh(a.im, sh, ch);
      mpf_sincos(a.re, st, ct);
      mpf_mul(ct,ch,b.re);
      mpf_mul(st,sh,b.im);
      s_mpf_chs(b.im);
      mpf_clear4(st,ct,sh,ch);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_cosh(const a: mp_complex; var b: mp_complex);
  {-Calculate b = cosh(a)}
var
  st,ct,sh,ch: mp_float;
begin
  if mpc_is0(a) then mpc_set1(b)
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(b.re);
    {$endif}
    mpf_initp4(st,ct,sh,ch,b.re.bitprec+8);
    if mp_error=MP_OKAY then begin
      {HMF[41], 4.5.50: cosh(x + iy) = cos(y)*cosh(x) + i*sin(y)*sinh(x)}
      mpf_sinhcosh(a.re, sh, ch);
      mpf_sincos(a.im, st, ct);
      mpf_mul(ct,ch,b.re);
      mpf_mul(st,sh,b.im);
      mpf_clear4(st,ct,sh,ch);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_cot(const a: mp_complex; var b: mp_complex);
  {-Calculate b = cot(a)}
begin
  {cot(a) = i*coth(i*a)}
  {b=i*a}
  mpc_copyp(a,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
  mpc_coth(b,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
end;


{---------------------------------------------------------------------------}
procedure mpc_coth(const a: mp_complex; var b: mp_complex);
  {-Calculate b = coth(a)}
var
  st,ct,sh,ch: mp_float;
  prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) or mpc_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_coth');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a.im) then begin
    mpf_coth(a.re, b.re);
    mpf_set0(b.im);
  end
  else if s_mpf_is0(a.re) then begin
    {coth(iy) = i*cot(y)}
    mpf_cot(a.im, b.im);
    s_mpf_chs(b.im);
    mpf_set0(b.re);
  end
  else begin
    prec := b.re.bitprec+8;
    mpf_initp4(st,ct,sh,ch,prec);
    if mp_error=MP_OKAY then begin
      {a = x + iy}
      mpf_mul_2k(a.re,1,sh);
      mpf_mul_2k(a.im,1,st);
      if mpf_cmp_mag_dbl(sh,prec*0.6931472) > 0 then begin
        {|x| is large:  cosh(2x) ~ sinh(2|x|) ~ 0.5exp(2|x|)}
        {coth = sign(x) - 2sin(2y)exp(-2|x|)}
        mpf_set1(b.re);
        if s_mpf_is_neg(sh) then s_mpf_chs(b.re)
        else s_mpf_chs(sh);
        mpf_exp(sh,ch);
        mpf_sin(st,ct);
        s_mpf_incexp(ct,1);
        s_mpf_chs(ct);
        mpf_mul(ct,ch,b.im);
      end
      else begin
        if s_mpf_ldx(a.re) < 0 then begin
          {|x| < 0.5, see AMath, with q = coshm1(2x)+vers(2y)}
          {b.re = sinh(2x)/q, b.im := -sin(2y)/q}
          mpf_coshm1(sh,ch);
          mpf_sinh(sh,sh);
          {vers(2y)=2*sin(y)^2, here y<>0: avoid vers overhead}
          mpf_sin(a.im, ct);
          mpf_sqr(ct,ct);
          mpf_mul_2k(ct,1,ct);
          mpf_add(ch,ct,ch);
          mpf_sin(st,st);
          s_mpf_chs(st);
          mpf_div(sh,ch,b.re);
          mpf_div(st,ch,b.im);
        end
        else begin
          {HMF[1], 4.5.52}
          {sh = sinh(2x), ch = cosh(2x)}
          {st = sin(2y),  ct = cos(2y)}
          mpf_sinhcosh(sh, sh, ch);
          mpf_sincos(st, st, ct);
          s_mpf_chs(st);
          mpf_sub(ch,ct,ch);
          mpf_div(sh,ch,b.re);
          mpf_div(st,ch,b.im);
        end;
      end;
      mpf_clear4(st,ct,sh,ch);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_csc(const a: mp_complex; var b: mp_complex);
  {-Calculate the complex circular cosecant b = csc(a) = 1/sin(a)}
var
  z: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpc_initp(z, b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    mpc_sin(a,z);
    mpc_inv(z,b);
    mpc_clear(z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_csch(const a: mp_complex; var b: mp_complex);
  {-Calculate the complex hyperbolic cosecant b = csch(a) = 1/sinh(a)}
var
  z: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpc_initp(z, b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    mpc_sinh(a,z);
    mpc_inv(z,b);
    mpc_clear(z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_div(const a,b: mp_complex; var c: mp_complex);
  {-Calculate c = a/b}
var
  x: mp_float;
begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(c.re);
  {$endif}
  mpf_initp(x,c.re.bitprec+8);
  {c = a*conj(b)/abs(b)^2}
  if mp_error=MP_OKAY then begin
    {x = abs(b)^2}
    mpc_abs2(b,x);
    {c = a*conj(b)/x}
    s_mpc_mul(a,b,true,c);
    mpf_div(c.re, x, c.re);
    mpf_div(c.im, x, c.im);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_div_dbl(const a: mp_complex; b: double; var c: mp_complex);
  {-Calculate c = a/b}
begin
  mpf_div_dbl(a.re, b, c.re);
  mpf_div_dbl(a.im, b, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_div_mpf(const a: mp_complex; const b: mp_float; var c: mp_complex);
  {-Calculate c = a/b}
begin
  mpf_div(a.re, b, c.re);
  mpf_div(a.im, b, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_exch(var a,b: mp_complex);
  {-Exchange two mp_complexes (including bitprec)}
begin
  mpf_exch(a.re,b.re);
  mpf_exch(a.im,b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_exp(const a: mp_complex; var b: mp_complex);
  {-Calculate b = exp(a)}
var
  x: mp_float;
begin
  if mpc_is0(a) then mpc_set1(b)
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(b.re);
    {$endif}
    mpf_initp(x,b.re.bitprec+8);
    if mp_error=MP_OKAY then begin
      mpf_exp(a.re, x);
      mpf_sincos(a.im, b.im, b.re);
      mpf_mul(b.re,x,b.re);
      mpf_mul(b.im,x,b.im);
      mpf_clear(x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_exp2(const a: mp_complex; var b: mp_complex);
  {-Calculate b = 2^a = exp(a*ln(2))}
var
  x,y: mp_float;
begin
  if mpc_is0(a) then mpc_set1(b)
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(b.re);
    {$endif}
    mpf_initp2(x,y,b.re.bitprec+8);
    if mp_error=MP_OKAY then begin
      mpf_set_ln2(y);
      mpf_mul(y,a.im,y);
      mpf_exp2(a.re, x);
      mpf_sincos(y, b.im, b.re);
      mpf_mul(b.re,x,b.re);
      mpf_mul(b.im,x,b.im);
      mpf_clear2(x,y);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_exp10(const a: mp_complex; var b: mp_complex);
  {-Calculate b = 10^a = exp(a*ln(10))}
var
  x,y: mp_float;
begin
  if mpc_is0(a) then mpc_set1(b)
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(b.re);
    {$endif}
    mpf_initp2(x,y,b.re.bitprec+8);
    if mp_error=MP_OKAY then begin
      mpf_set_ln10(y);
      mpf_mul(y,a.im,y);
      mpf_exp10(a.re, x);
      mpf_sincos(y, b.im, b.re);
      mpf_mul(b.re,x,b.re);
      mpf_mul(b.im,x,b.im);
      mpf_clear2(x,y);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_expm1(const a: mp_complex; var b: mp_complex);
  {-Calculate b = exp(a)-1}
var
  lr, lm, bp: longint;
var
  w,x,y,z: mp_float;
  t: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) or mpc_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpc_expm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  lr := s_mpf_ldx(a.re);
  lm := s_mpf_ldx(a.im);
  bp := b.re.bitprec+8;
  if lr > lm then lm := lr;
  if mpf_is0(a.im) then begin
    {a is real}
    mpf_expm1(a.re, b.re);
    mpf_set0(b.im);
  end
  else if lm < -bp then begin
    if mpc_is0(a) then begin
      mpc_set0(b);
      exit;
    end;
    {very small a: expm1(a) = a + a^2/ + O(a^3)}
    {Note that using only one term would give the wrong}
    {result b.re = 0, if a.re = 0 and a.im <> 0}
    mpc_initp(t,bp);
    if mp_error=MP_OKAY then begin
      mpc_sqr(a,t);
      mpc_mul_2k(t,-1,t);
      mpc_add(a,t,b);
      mpc_clear(t);
    end;
    exit;
  end
  else if lm>0 then begin
    {|a| >= 1}
    mpc_exp(a,b);
    s_mpf_dec1(b.re);
  end
  else begin
    mpf_initp4(w,x,y,z,bp);
    if mp_error=MP_OKAY then begin
      {with z = expm1(a.re) we have (c.f. AMath manual)}
      {b.re = z - (1+z)*vers(a.im)}
      {b.im = (1+z)*sin(a.im)     }
      mpf_expm1(a.re, z);
      mpf_copy(z,w);
      s_mpf_inc1(w);
      mpf_sin(a.im, y);
      mpf_vers(a.im, x);
      mpf_mul(w,y, b.im);
      mpf_mul(w,x,x);
      mpf_sub(z,x,b.re);
      mpf_clear4(w,x,y,z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_init(var a: mp_complex);
  {-Initialize an mp_complex with default precision}
begin
  mpc_initp(a, mpf_get_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_init2(var a,b: mp_complex);
  {-Initialize two mp_complexes with default precision}
begin
  mpc_initp2(a,b,mpf_get_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_init3(var a,b,c: mp_complex);
  {-Initialize 3 mp_complexes with default precision}
begin
  mpc_initp3(a,b,c,mpf_get_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_init4(var a,b,c,d: mp_complex);
  {-Initialize 4 mp_complexes with default precision}
begin
  mpc_initp4(a,b,c,d,mpf_get_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_init5(var a,b,c,d,e: mp_complex);
  {-Initialize 5 mp_complexes with default precision}
begin
  mpc_initp5(a,b,c,d,e,mpf_get_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_initp(var a: mp_complex; prec: longint);
  {-Initialize an mp_complex with bit precision prec}
begin
  if mp_error<>MP_OKAY then exit;
  if prec<MPF_MIN_PREC then prec := MPF_MIN_PREC;
  if prec>MPF_MAX_PREC then prec := MPF_MAX_PREC;
  mpf_initp(a.re, prec);
  if mp_error=0 then begin
    mpf_initp(a.im, prec);
    if mp_error<>0 then mpf_clear(a.re);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_initp2(var a,b: mp_complex; prec: longint);
  {-Initialize two mp_complexes with bit precision prec}
var
  pa: array[0..1] of pmp_complex;
begin
  pa[0] := @a;
  pa[1] := @b;
  mpc_initp_multi_p(pa, prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_initp3(var a,b,c: mp_complex; prec: longint);
  {-Initialize 3 mp_complexes with bit precision prec}
var
  pa: array[0..2] of pmp_complex;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  mpc_initp_multi_p(pa, prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_initp4(var a,b,c,d: mp_complex; prec: longint);
  {-Initialize 4 mp_complexes with bit precision prec}
var
  pa: array[0..3] of pmp_complex;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  mpc_initp_multi_p(pa, prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_initp5(var a,b,c,d,e: mp_complex; prec: longint);
  {-Initialize 5 mp_complexes with bit precision prec}
var
  pa: array[0..4] of pmp_complex;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  mpc_initp_multi_p(pa, prec);
end;


{---------------------------------------------------------------------------}
procedure mpc_initp_multi_p(var pv: array of pmp_complex; prec: longint);
  {-Initialize with bit precision prec a list of mp_complexes given as a pointer}
  { vector; on error the already initialized mp_complexes will be cleared}
var
  i,k: integer;
begin
  if mp_error<>MP_OKAY then exit;
  for i:=low(pv) to high(pv) do begin
    mpc_initp(pv[i]^, prec);
    if mp_error<>MP_OKAY then begin
      {error, clear all previous mp_complexes}
      for k:=low(pv) to i-1 do mpc_clear(pv[k]^);
      break;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_inv(const a: mp_complex; var b: mp_complex);
  {-Calculate b = 1/a}
begin
  mpc_xdivc(1.0,a,b);
end;


{---------------------------------------------------------------------------}
function mpc_is0(const a: mp_complex): boolean;
  {-Return true if a=0}
begin
  mpc_is0 := mpf_is0(a.re) and mpf_is0(a.im);
end;


{---------------------------------------------------------------------------}
function mpc_is1(const a: mp_complex): boolean;
  {-Return true if a=1}
begin
  mpc_is1 := mpf_is0(a.im) and mpf_is1(a.re);
end;


{---------------------------------------------------------------------------}
function mpc_is1a(const a: mp_complex): boolean;
  {-Return true if a = +1 or -1}
begin
  mpc_is1a := mpf_is0(a.im) and mpf_is1a(a.re);
end;


{---------------------------------------------------------------------------}
function mpc_is_i(const a: mp_complex): boolean;
  {-Return true if a=I}
begin
  mpc_is_i := mpf_is0(a.re) and mpf_is1(a.im);
end;


{---------------------------------------------------------------------------}
function mpc_is_ia(const a: mp_complex): boolean;
  {-Return true if a=I or a= -I}
begin
  mpc_is_ia := mpf_is0(a.re) and mpf_is1a(a.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_ln(const a: mp_complex; var b: mp_complex);
  {-Calculate the natural logarithm b = ln(a); principal branch ln(|a|) + i*arg(a), accurate near |a|=1}
var
  x,y,z: mp_float;
  k: longint;
begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpf_initp3(x,y,z,b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    mpf_abs(a.re,x);
    mpf_abs(a.im,y);
    if mpf_is_lt(x,y) then mpf_exch(x,y);
    k := s_mpf_ldx(x);
    if (k<0) or (k>=2) then begin
      {x < 0.5 or x >= 2:  b.re = ln(abs(a))}
      mpf_hypot(x,y,x);
      mpf_ln(x,x);
    end
    else begin
      {Ref: Kahan[42], procedure CLOGS }
      {b.re = ln1p((x-1)*(x+1) + y*y)/2}
      mpf_copy(x,z);
      s_mpf_inc1(x);
      s_mpf_dec1(z);
      mpf_mul(x,z,x);
      mpf_sqr(y,y);
      mpf_add(x,y,x);
      mpf_ln1p(x,x);
      s_mpf_incexp(x, -1);
    end;
    mpc_arg(a,b.im);
    mpf_copyp(x,b.re);
    mpf_clear3(x,y,z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_ln1p(const a: mp_complex; var b: mp_complex);
  {-Calculate the natural logarithm b = ln(1+a)}
var
  lr, lm, bp: longint;
var
  x,y: mp_float;
  t: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) or mpc_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpc_ln1p');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  lr := s_mpf_ldx(a.re);
  lm := s_mpf_ldx(a.im);
  bp := b.re.bitprec+8;
  if lr > lm then lm := lr;
  if lm < -bp then begin
    if mpc_is0(a) then begin
      mpc_set0(b);
      exit;
    end;
    {very small a: ln(1+a) = a - a^2/ + O(a^3)}
    {Note that using only one term would give the wrong}
    {result b.re = 0, if a.re = 0 and a.im <> 0}
    mpc_initp(t,bp);
    if mp_error=MP_OKAY then begin
      mpc_sqr(a,t);
      mpc_mul_2k(t,-1,t);
      mpc_sub(a,t,b);
      mpc_clear(t);
    end;
    exit;
  end;
  mpf_initp2(x,y,bp);
  if mp_error=MP_OKAY then begin
    {compute b.re = y = ln(|1+a|)}
    if (lm>0) or ((lr >= 0) and (a.re.mantissa.sign=MP_NEG)) then begin
      {|a.im| >= 1  or  a.re <= -0.5; y = ln(sqrt(|1+a|^2))}
      s_mpf_add1(a.re, x);
      mpf_hypot(x,a.im,y);
      mpf_ln(y,y);
    end
    else begin
      {small a: y = 0.5*ln1p(2*a.re + |a|^2)}
      mpf_sqr(a.re, x);
      mpf_sqr(a.im, y);
      mpf_add(x,y,y);
      mpf_mul_2k(a.re, 1, x);
      mpf_add(x,y,y);
      mpf_ln1p(y,y);
      mpf_mul_2k(y, -1, y);
      s_mpf_add1(a.re, x);
    end;
    {b.im = arg(1+a)}
    mpf_arctan2(a.im, x, b.im);
    mpf_copyp(y, b.re);
    mpf_clear2(x,y);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_log10(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal branch of the base 10 logarithm of a, b=ln(z)/ln(10)}
var
  ln10: mp_float;
begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpc_ln(a,b);
  if mp_error<>MP_OKAY then exit;
  mpf_initp(ln10, b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    mpf_set_ln10(ln10);
    mpf_div(b.re, ln10, b.re);
    mpf_div(b.im, ln10, b.im);
    mpf_clear(ln10);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_mul(const a,b: mp_complex; var c: mp_complex);
  {-Calculate c = a*b}
begin
  s_mpc_mul(a,b,false,c);
end;


{---------------------------------------------------------------------------}
procedure mpc_mul_2k(const a: mp_complex; k: longint; var b: mp_complex);
  {-Calculate b = a*2^k}
begin
  mpf_mul_2k(a.re, k, b.re);
  mpf_mul_2k(a.im, k, b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_mul_dbl(const a: mp_complex; b: double; var c: mp_complex);
  {-Calculate c = a*b}
begin
  mpf_mul_dbl(a.re, b, c.re);
  mpf_mul_dbl(a.im, b, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_mul_mpf(const a: mp_complex; const b: mp_float; var c: mp_complex);
  {-Calculate c = a*b}
begin
  mpf_mul(a.re, b, c.re);
  mpf_mul(a.im, b, c.im);
end;


{---------------------------------------------------------------------------}
function mpc_not_init(const a: mp_complex): boolean;
  {-Sanity check if a is initialized, does not catch all cases!}
begin
  mpc_not_init := mpf_not_init(a.re) or mpf_not_init(a.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_pow(const a,b: mp_complex; var c: mp_complex);
  {-Calculate the principal value of the complex power c = a^b = exp(b*ln(a))}
var
  z: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;

  if s_mpf_is0(b.im) then begin
    if s_mpf_is0(b.re) then begin
      {b=0: c=1}
      mpc_set1(c);
      exit;
    end
    else if mpf_is1a(b.re) then begin
      {b=+-1: c=a or c=1/a}
      if s_mpf_is_neg(b.re) then mpc_inv(a,c)
      else mpc_copyp(a,c);
      exit;
    end;
  end;

  if mpc_is0(a) and s_mpf_is_gt0(b.re) then mpc_set0(c)
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(c.re);
    {$endif}
    mpc_initp(z, c.re.bitprec+8);
    if mp_error=MP_OKAY then begin
      {c = a^b = exp(b*ln(a))}
      mpc_ln(a,z);
      mpc_mul(b,z,z);
      mpc_exp(z,c);
      mpc_clear(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_nroot(const a: mp_complex; n: longint; var b: mp_complex);
  {-Calculate the nth principal root b = a^(1/n) = exp(ln(a)/n)}
begin
  if n=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpc_nroot: n=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mpc_ln(a,b);
  mpf_div_int(b.re,n,b.re);
  mpf_div_int(b.im,n,b.im);
  mpc_exp(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpc_nroot1(n: longint; var a: mp_complex);
  {-Calculate the principal nth root of unity a = exp(2*Pi*i/n)}
var
  x: mp_float;
  k: longint;
begin
  if mp_error<>MP_OKAY then exit;
  k := abs(n);
  if (k=1) or (k=2) then begin
    mpc_set1(a);
    if k=2 then s_mpf_chs(a.re);
    exit;
  end
  else if k=4 then begin
    {a = +- I}
    mpf_set0(a.re);
    if n>0 then mpf_set1(a.im)
    else mpf_set_dbl(a.im,-1.0);
    exit;
  end;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a.re);
  {$endif}
  mpf_initp(x, a.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    mpf_set_pi2k(x,1);
    mpf_div_int(x,n,x);
    mpc_cis(x,a);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_sec(const a: mp_complex; var b: mp_complex);
  {-Calculate the complex circular secant b = sec(a) = 1/cos(a)}
var
  z: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpc_initp(z, b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    mpc_cos(a,z);
    mpc_inv(z,b);
    mpc_clear(z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_sech(const a: mp_complex; var b: mp_complex);
  {-Calculate the complex hyperbolic secant b = sech(a) = 1/cosh(a)}
var
  z: mp_complex;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpc_initp(z, b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    mpc_cosh(a,z);
    mpc_inv(z,b);
    mpc_clear(z);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_set0(var a: mp_complex);
  {-Set a=0}
begin
  mpf_set0(a.re);
  mpf_set0(a.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_set1(var a: mp_complex);
  {-Set a=1}
begin
  mpf_set1(a.re);
  mpf_set0(a.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_seti(var a: mp_complex);
  {-Set a=i}
begin
  mpf_set0(a.re);
  mpf_set1(a.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_set_dbl(var a: mp_complex; x,y: double);
  {-Set a = x + iy with x,y double, error if x,y are NAN or INF}
begin
  mpf_set_dbl(a.re,x);
  mpf_set_dbl(a.im,y);
end;


{$ifndef EXT64}
{---------------------------------------------------------------------------}
procedure mpc_set_ext(var a: mp_complex; x,y: extended);
  {-Set a = x + iy with x,y extended, error if x,y are NAN or INF}
begin
  mpf_set_ext(a.re,x);
  mpf_set_ext(a.im,y);
end;
{$endif}


{---------------------------------------------------------------------------}
procedure mpc_set_mpf(var a: mp_complex; const x,y: mp_float);
  {-Set a = x + iy}
begin
  mpf_copy(x,a.re);
  mpf_copy(y,a.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_setp_mpf(var a: mp_complex; const x,y: mp_float);
  {-Set a = x + iy, bitprecs of a are preserved}
begin
  mpf_copyp(x,a.re);
  mpf_copyp(y,a.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_sin(const a: mp_complex; var b: mp_complex);
  {-Calculate b = sin(a)}
var
  st,ct,sh,ch: mp_float;
begin
  if mpc_is0(a) then mpc_set0(b)
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(b.re);
    {$endif}
    mpf_initp4(st,ct,sh,ch,b.re.bitprec+8);
    if mp_error=MP_OKAY then begin
      mpf_sinhcosh(a.im, sh, ch);
      mpf_sincos(a.re, st, ct);
      mpf_mul(st,ch,b.re);
      mpf_mul(ct,sh,b.im);
      mpf_clear4(st,ct,sh,ch);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_sincos(const a: mp_complex; var s,c: mp_complex);
  {-Calculate s = sin(a), c = cos(a)}
var
  st,ct,sh,ch: mp_float;
  p: longint;
begin
  if mpc_is0(a) then begin
    mpc_set0(s);
    mpc_set1(c);
  end
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(s.re);
      _CheckBitPrec(c.re);
    {$endif}
    p := s.re.bitprec;
    if p<c.re.bitprec then p := c.re.bitprec;
    mpf_initp4(st,ct,sh,ch,p+8);
    if mp_error=MP_OKAY then begin
      mpf_sinhcosh(a.im, sh, ch);
      mpf_sincos(a.re, st, ct);
      mpf_mul(st,ch,s.re);
      mpf_mul(ct,sh,s.im);
      mpf_mul(ct,ch,c.re);
      mpf_mul(st,sh,c.im);
      s_mpf_chs(c.im);
      mpf_clear4(st,ct,sh,ch);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_sinh(const a: mp_complex; var b: mp_complex);
  {-Calculate b = sinh(a)}
var
  st,ct,sh,ch: mp_float;
begin
  if mpc_is0(a) then mpc_set0(b)
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(b.re);
    {$endif}
    mpf_initp4(st,ct,sh,ch,b.re.bitprec+8);
    if mp_error=MP_OKAY then begin
      {HMF[41], 4.5.49: sinh(x + iy) = cos(y)*sinh(x) + i*sin(y)*cosh(x)}
      mpf_sinhcosh(a.re, sh, ch);
      mpf_sincos(a.im, st, ct);
      mpf_mul(ct,sh,b.re);
      mpf_mul(st,ch,b.im);
      mpf_clear4(st,ct,sh,ch);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_sinhcosh(const a: mp_complex; var s,c: mp_complex);
  {-Calculate s = sinh(a), c = cosh(a)}
var
  st,ct,sh,ch: mp_float;
  p: longint;
begin
  if mpc_is0(a) then begin
    mpc_set0(s);
    mpc_set1(c);
  end
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(s.re);
      _CheckBitPrec(c.re);
    {$endif}
    p := s.re.bitprec;
    if p<c.re.bitprec then p := c.re.bitprec;
    mpf_initp4(st,ct,sh,ch,p+8);
    if mp_error=MP_OKAY then begin
      mpf_sinhcosh(a.re, sh, ch);
      mpf_sincos(a.im, st, ct);
      mpf_mul(ct,sh,s.re);
      mpf_mul(st,ch,s.im);
      mpf_mul(ct,ch,c.re);
      mpf_mul(st,sh,c.im);
      mpf_clear4(st,ct,sh,ch);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_sqr(const a: mp_complex; var b: mp_complex);
  {-Calculate b = a*a}
var
  x,y: mp_float;
begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpf_initp2(x,y,b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    {b.re = (a.re - a.im)*(a.re + a.im)}
    {b.im = 2*z.re*z.im}
    mpf_sub(a.re, a.im, x);
    mpf_add(a.re, a.im, y);
    mpf_mul(a.re,a.im,b.im);
    s_mpf_incexp(b.im, 1);
    mpf_mul(x,y,b.re);
    mpf_clear2(x,y);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_sqrt(const a: mp_complex; var b: mp_complex);
  {-Calculate the principal square root b = sqrt(a)}
var
  x,y: mp_float;
begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b.re);
  {$endif}
  mpf_initp2(x,y,b.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    mpc_abs(a,x);
    if s_mpf_is0(x) then mpc_set0(b)
    else begin
      if s_mpf_is_ge0(a.re) then begin
        {a.re >=0:  b.re = sqrt((x+a.re)/2), b.im=a.im/(2*c)}
        mpf_add(x,a.re,x);
        s_mpf_incexp(x, -1);
        mpf_sqrt(x,x);
        mpf_copyp(x,b.re);
        mpf_div(a.im,x,x);
        mpf_mul_2k(x, -1, b.im)
      end
      else begin
        {a.re < 0:  b.im = sqrt((x-a.re)/2)*sign(a.im), b.re=a.im/(2*d)}
        mpf_sub(x,a.re,x);
        s_mpf_incexp(x, -1);
        mpf_sqrt(x,x);
        if s_mpf_is_neg(a.im) then s_mpf_chs(x);
        mpf_div(a.im,x,y);
        mpf_mul_2k(y, -1, b.re);
        mpf_copyp(x,b.im);
      end;
    end;
    mpf_clear2(x,y);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_sub(const a,b: mp_complex; var c: mp_complex);
  {-Calculate c = a-b}
begin
  mpf_sub(a.re, b.re, c.re);
  mpf_sub(a.im, b.im, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_sub_dbl(const a: mp_complex; b: double; var c: mp_complex);
  {-Calculate c = a-b}
begin
  mpf_sub_dbl(a.re, b, c.re);
  mpf_copy(a.im, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_sub_mpf(const a: mp_complex; const b: mp_float; var c: mp_complex);
  {-Calculate c = a-b}
begin
  mpf_sub(a.re, b, c.re);
  mpf_copy(a.im, c.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_tan(const a: mp_complex; var b: mp_complex);
  {-Calculate b = tan(a)}
begin
  {tan(a) = -i*tanh(i*a)}
  {b=i*a}
  mpc_copyp(a,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.re);
  mpc_tanh(b,b);
  mpf_exch(b.re,b.im);
  s_mpf_chs(b.im);
end;


{---------------------------------------------------------------------------}
procedure mpc_tanh(const a: mp_complex; var b: mp_complex);
  {-Calculate b = tanh(a)}
var
  st,ct,sh,ch: mp_float;
  prec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpc_not_init(a) or mpc_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_tanh');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a.im) then begin
    mpf_tanh(a.re, b.re);
    mpf_set0(b.im);
  end
  else if s_mpf_is0(a.re) then begin
    {tanh(iy) = i*tan(y)}
    mpf_tan(a.im, b.im);
    mpf_set0(b.re);
  end
  else begin
    prec := b.re.bitprec+8;
    mpf_initp4(st,ct,sh,ch,prec);
    if mp_error=MP_OKAY then begin
      {a = x + iy}
      if s_mpf_ldx(a.re) < 0 then begin
        {|x| < 0.5, avoid cancellation for cosh(2x) + cos(2y)}
        {with more accurate but slightly slower computation}
        mpf_sinhcosh(a.re, sh, ch);
        mpf_sincos(a.im, st, ct);
        {cosh(2x) + cos(2y) = 2*[sinh(x)^2 + cos(y)^2)]}
        {sinh(2x) = 2sinh(x)cosh(x), sin(2y)=2sin(y)cos(y)}
        {sh = sinh(x), ch = cosh(x)}
        {st = sin(y),  ct = cos(y) }
        {q  = ct^2 + sh^2, b.re = sh*ch/q, b.im = st*ct/q}
        mpf_mul(sh,ch,ch);
        mpf_mul(st,ct,st);
        mpf_sqr(sh,sh);
        mpf_sqr(ct,ct);
        mpf_add(ct,sh,ct);
        mpf_div(ch,ct,b.re);
        mpf_div(st,ct,b.im);
      end
      else if mpf_cmp_mag_dbl(a.re,prec*0.3465736) > 0 then begin
        mpf_mul_2k(a.re,1,sh);
        mpf_mul_2k(a.im,1,st);
        {|x| is large:  cosh(2x) ~ sinh(2|x|) ~ 0.5exp(2|x|)}
        {tanh = sign(x) + 2sin(2y)exp(-2|x|)}
        mpf_set1(b.re);
        if s_mpf_is_neg(sh) then s_mpf_chs(b.re)
        else s_mpf_chs(sh);
        mpf_exp(sh,ch);
        mpf_sin(st,ct);
        s_mpf_incexp(ct,1);
        mpf_mul(ct,ch,b.im);
      end
      else begin
        mpf_mul_2k(a.re,1,sh);
        mpf_mul_2k(a.im,1,st);
        {HMF[41], 4.5.51}
        {sh = sinh(2x), ch = cosh(2x)}
        {st = sin(2y),  ct = cos(2y)}
        mpf_sinhcosh(sh, sh, ch);
        mpf_sincos(st, st, ct);
        mpf_add(ct,ch,ch);
        mpf_div(sh,ch,b.re);
        mpf_div(st,ch,b.im);
      end;
      mpf_clear4(st,ct,sh,ch);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpc_xdivc(a: double; const b: mp_complex; var c: mp_complex);
  {-Calculate c = a/b}
var
  x: mp_float;
begin
  if a=0.0 then mpc_set0(c)
  else begin
    {$ifdef MPC_ArgCheck}
      _CheckBitPrec(c.re);
    {$endif}
    mpf_initp(x,c.re.bitprec+8);
    {c = a*conj(b)/abs(b)^2}
    if mp_error=MP_OKAY then begin
      {x = a/abs(b)^2}
      mpc_abs2(b,x);
      mpf_divr_dbl(a,x,x);
      {c = conj(b)*x}
      mpc_conj(b,c);
      mpf_mul(c.re, x, c.re);
      mpf_mul(c.im, x, c.im);
      mpf_clear(x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpc_mul(const a,b: mp_complex; bc: boolean; var c: mp_complex);
  {-Calculate c = a*b or c=a*conj(b) if bc}
var
  x,y,z: mp_float;
begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(c.re);
  {$endif}
  mpf_initp3(x,y,z,c.re.bitprec+8);
  if mp_error=MP_OKAY then begin
    {z = a.re*b.re -+ a.im*b.im}
    mpf_mul(a.re, b.re, x);
    mpf_mul(a.im, b.im, y);
    if bc then mpf_add(x,y,z)
    else mpf_sub(x,y,z);
    {c.im := a.im*b.re -+ a.re*b.im}
    mpf_mul(a.re, b.im, x);
    mpf_mul(a.im, b.re, y);
    if bc then mpf_sub(y,x,c.im)
    else mpf_add(x,y,c.im);
    mpf_copy(z,c.re);
    mpf_clear3(x,y,z);
  end;
end;


end.
