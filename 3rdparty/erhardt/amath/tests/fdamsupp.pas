unit FDAMSupp;

{MPArith support routines for D/AMath test programs}

interface


{$i STD.INC}

{$ifdef BIT16}
  {$N+}
{$endif}

uses
  BTypes, mp_types, mp_real;

{$i mp_conf.inc}


(*************************************************************************

 DESCRIPTION   :  MPArith support routines for D/AMath test programs}

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 REMARK        :  The routines should not be used as stand-alone functions!
                  The precision handling is appropriate only for checking
                  the AMath and DAMath functions! Many routines do not have
                  error handling, argument and/or range checks.


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------

 0.0.01   30.06.14  W.Ehrhardt  Separate unit forked from t_amathx
 0.0.02   01.07.14  we          mpf_rem_2pi from t_xdamat
 0.0.03   06.10.14  we          mpf_logistic
 0.0.04   31.05.15  we          mpf_versint
 0.0.05   19.06.15  we          mpf_sinhc
 0.0.06   25.06.15  we          mpf_sinhmx

**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2007-2015 Wolfgang Ehrhardt

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

procedure mpf_arccosd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccos(a), b in degrees}

procedure mpf_arccotcd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccotc(a), b in degrees}

procedure mpf_arccotd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccot(a), b in degrees}

procedure mpf_arcsind(const a: mp_float; var b: mp_float);
  {-Calculate b = arcsin(a), b in degrees}

procedure mpf_arctand(const a: mp_float; var b: mp_float);
  {-Calculate b = arctan(a), b in degrees}

procedure mpf_cbrt(const a: mp_float; var b: mp_float);
  {-Calculate b = a^(1/3) }

procedure mpf_cosd(const a: mp_float; var b: mp_float);
  {-Calculate b = cos(a), a in degrees}

procedure mpf_cotd(const a: mp_float; var b: mp_float);
  {-Calculate b = cot(a), a in degrees}

procedure mpf_covers(const a: mp_float; var b: mp_float);
  {-Calculate b = 1-sin(a)}

procedure mpf_deg2rad(const a: mp_float; var b: mp_float);
  {-Convert a in degrees to b in radians}

procedure mpf_exp3(const a: mp_float; var b: mp_float);
  {-Calculate b = 3^a}

procedure mpf_exp5(const a: mp_float; var b: mp_float);
  {-Calculate b = 5^a}

procedure mpf_exp7(const a: mp_float; var b: mp_float);
  {-Calculate b = 7^a}

procedure mpf_expmx2h(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(-0.5*x^2) }

procedure mpf_exprel(const a: mp_float; var b: mp_float);
  {-Calculate b = (exp(a)-1)/a, 1 for a=0}

procedure mpf_expt1p(const a,b: mp_float; var c: mp_float);
  {-Calculate c = (1+a)^b, a>-1}

procedure mpf_expx2(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(a*abs(a)) }

procedure mpf_langevin(const a: mp_float; var b: mp_float);
  {-Calculate the Langevin function b = L(a) = coth(a) - 1/a, L(0) = 0}

procedure mpf_ln1mexp(const a: mp_float; var b: mp_float);
  {-Calculate b := ln(1-exp(a)), a<0}

procedure mpf_ln1pmx(const x: mp_float; var y: mp_float);
  {-Calculate y = ln(1+x)-x}

procedure mpf_logistic(const a: mp_float; var b: mp_float);
  {-Calculate b = 1/(1+exp(-a))}

procedure mpf_logit(const a: mp_float; var b: mp_float);
  {-Calculate b := ln(a/(1-a)), 0<a<1)}

procedure mpf_rad2deg(const a: mp_float; var b: mp_float);
  {-Convert a in radians to b in degrees}

procedure mpf_rem_2pi(const a: mp_float; var b: mp_float);
  {-Calculate b = a mod (2*Pi)}

procedure mpf_sincPi(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(Pi*a)/(Pi*a)}

procedure mpf_sind(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(a), a in degrees}

procedure mpf_sinhc(const a: mp_float; var b: mp_float);
  {-Calculate b = sinh(a)/a}

procedure mpf_sinhmx(const a: mp_float; var b: mp_float);
  {-Calculate b = sinh(a)-a}

procedure mpf_tand(const a: mp_float; var b: mp_float);
  {-Calculate b = tan(a), a in degrees}

procedure mpf_versint(const a: mp_float; var b: mp_float);
  {-Return b := versint(a) = int(vers(a),a) = a - sin(a)}


implementation


{---------------------------------------------------------------------------}
procedure mpf_exp3(const a: mp_float; var b: mp_float);
  {-Calculate b = 3^a}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(t,3.0);
    mpf_expt(t,a,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_exp5(const a: mp_float; var b: mp_float);
  {-Calculate b = 5^a}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(t,5.0);
    mpf_expt(t,a,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_exp7(const a: mp_float; var b: mp_float);
  {-Calculate b = 7^a}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(t,7.0);
    mpf_expt(t,a,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sincPi(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(Pi*a)/(Pi*a)}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_trig_ex(a, true, nil, @t, nil);
    mpf_div(t,a,b);
    mpf_set_pi(t);
    mpf_div(b,t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_ln1pmx(const x: mp_float; var y: mp_float);
  {-Calculate y = ln(1+x)-x}
var
  c: mp_float;
begin
  mpf_initp(c,y.bitprec+32);
  if mp_error<>MP_OKAY then exit;
  mpf_ln1p(x,c);
  mpf_sub(c,x,y);
  mpf_clear(c);
end;


{---------------------------------------------------------------------------}
procedure mpf_expx2(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(a*abs(a)) }
begin
  mpf_abs(a,b);
  mpf_mul(a,b,b);
  mpf_exp(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_expmx2h(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(-0.5*x^2) }
begin
  mpf_mul(a,a,b);
  mpf_mul_ext(b,-0.5,b);
  mpf_exp(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_cbrt(const a: mp_float; var b: mp_float);
  {-Calculate b = a^(1/3) }
var
  an: boolean;
begin
  an := s_mpf_is_neg(a);
  mpf_abs(a,b);
  mpf_ln(b,b);
  mpf_div_ext(b,3.0,b);
  mpf_exp(b,b);
  if an then s_mpf_chs(b);
end;


{---------------------------------------------------------------------------}
procedure mpf_exprel(const a: mp_float; var b: mp_float);
  {-Calculate b = (exp(a)-1)/a, 1 for a=0}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_expm1(a,t);
    mpf_div(t,a,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_covers(const a: mp_float; var b: mp_float);
  {-Calculate b = 1-sin(a)}
var
  e,t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp2(e,t,b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_set1(e);
    mpf_sin(a,t);
    mpf_sub(e,t,t);
    mpf_copyp(t,b);
    mpf_clear2(e,t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_versint(const a: mp_float; var b: mp_float);
  {-Return b := versint(a) = int(vers(a),a) = a - sin(a)}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_versint');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  if s_mpf_ldx(a)>0 then mpf_initp(t, b.bitprec+32)
  else mpf_initp(t, 2*b.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_sin(a,t);
    mpf_sub(a,t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_deg2rad(const a: mp_float; var b: mp_float);
  {-Convert a in degrees to b in radians}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp(t,b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_set_pi(t);
    mpf_mul(t,a,t);
    mpf_div_d(t,180,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_rad2deg(const a: mp_float; var b: mp_float);
  {-Convert a in radians to b in degrees}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp(t,b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_set_pi(t);
    mpf_div(a,t,t);
    mpf_mul_d(t,180,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_rem_2pi(const a: mp_float; var b: mp_float);
  {-Calculate b = a mod (2*Pi)}
var
  oddm: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  s_mpf_mod_pi2k(a,1,b,oddm);
end;


{---------------------------------------------------------------------------}
procedure mpf_langevin(const a: mp_float; var b: mp_float);
  {-Calculate the Langevin function b = L(a) = coth(a) - 1/a, L(0) = 0}
var
  s,t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp2(s,t,b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_coth(a,s);
    mpf_inv(a,t);
    mpf_sub(s,t,t);
    mpf_copyp(t,b);
    mpf_clear2(s,t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_ln1mexp(const a: mp_float; var b: mp_float);
  {-Calculate b := ln(1-exp(a)), a<0}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  mpf_initp(t, b.bitprec+32);
  if mp_error=MP_OKAY then begin
    if mpf_cmp_ext(a,-0.693147180559945309)>0 then begin
      {range -ln(2) < a < 0}
      mpf_expm1(a,t);
      s_mpf_chs(t);
      mpf_ln(t,b);
    end
    else begin
      mpf_exp(a,t);
      s_mpf_chs(t);
      mpf_ln1p(t,b);
    end;
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_logistic(const a: mp_float; var b: mp_float);
  {-Calculate b = 1/(1+exp(-a))}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  mpf_initp(t,b.bitprec+8);
  if mp_error=MP_OKAY then begin
    {check exp(-a) < 2^(-t.bitprec)}
    if s_mpf_ldx(a) > 1.45*ln(t.bitprec+2.0) then begin
      if s_mpf_is_ge0(a) then mpf_set0(b)
      else begin
        {a negative, 1/(1+exp(-a)) ~ 1/exp(-a) = exp(a)}
        mpf_exp(a,b);
      end;
    end
    else begin
      mpf_chs(a,t);
      mpf_exp(t,t);
      s_mpf_inc1(t);
      mpf_inv(t,b);
    end;
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_logit(const a: mp_float; var b: mp_float);
  {-Calculate b := ln(a/(1-a)), 0<a<1)}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  mpf_initp(t, b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_set1(t);
    mpf_sub(t,a,t);
    mpf_div(a,t,t);
    mpf_ln(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_expt1p(const a,b: mp_float; var c: mp_float);
  {-Calculate c = (1+a)^b, a>-1}
var
  t: mp_float;
  p: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_expt1p');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is_le0(a) then p:=0
  else begin
    p := abs(s_mpf_ldx(a));
    if p<32 then p := 32;
  end;
  mpf_initp(t, c.bitprec+p);
  if mp_error=MP_OKAY then begin
    mpf_ln1p(a,t);
    mpf_mul(t,b,t);
    mpf_exp(t,c);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sind(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(a), a in degrees}
begin
  mpf_deg2rad(a,b);
  mpf_sin(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_sinhc(const a: mp_float; var b: mp_float);
  {-Calculate b = sinh(a)/a}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec+32);
  if mp_error=MP_OKAY then begin
    mpf_sinh(a,t);
    mpf_div(t,a,t);
    mpf_copyp(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sinhmx(const a: mp_float; var b: mp_float);
  {-Calculate b = sinh(a)-a}
var
  t: mp_float;
  p: longint;
begin
  if mp_error<>MP_OKAY then exit;
  p := b.bitprec+16;
  if s_mpf_ldx(a)<=1 then p := 2*p;
  mpf_initp(t, p);
  if mp_error=MP_OKAY then begin
    mpf_sinh(a,t);
    mpf_sub(t,a,t);
    mpf_copyp(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_cosd(const a: mp_float; var b: mp_float);
  {-Calculate b = cos(a), a in degrees}
begin
  mpf_deg2rad(a,b);
  mpf_cos(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_tand(const a: mp_float; var b: mp_float);
  {-Calculate b = tan(a), a in degrees}
begin
  mpf_deg2rad(a,b);
  mpf_tan(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_cotd(const a: mp_float; var b: mp_float);
  {-Calculate b = cot(a), a in degrees}
begin
  mpf_deg2rad(a,b);
  mpf_cot(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arctand(const a: mp_float; var b: mp_float);
  {-Calculate b = arctan(a), b in degrees}
begin
  mpf_arctan(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccotd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccot(a), b in degrees}
begin
  mpf_arccot(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccotcd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccotc(a), b in degrees}
begin
  mpf_arccotc(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arcsind(const a: mp_float; var b: mp_float);
  {-Calculate b = arcsin(a), b in degrees}
begin
  mpf_arcsin(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccosd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccos(a), b in degrees}
begin
  mpf_arccos(a,b);
  mpf_rad2deg(b,b);
end;


end.
