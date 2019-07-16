unit AMQuat;

{AMath based quaternion routines}

interface

{$i STD.INC}

{$ifdef BIT16}
{$N+,F+}
{$endif}

(*************************************************************************

 DESCRIPTION   :  AMath based quaternion routines


 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    : References used in this unit, main index in amath_info.txt/references


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.01     09.11.16  W.Ehrhardt  Initial BP7 version: qnorm, qabs, qnorm_im, qabs_im
                                qadd, qsub, qmul, qdiv, qexp, qcos, qln
 0.02     10.11.16  we          qcosh, qsinh, qsin, qtan, qtanh, qsqrt
 0.03     10.11.16  we          new qsqrt using cx_sqrt
 0.04     10.11.16  we          qarg, qinv, qaxis
 0.05     12.11.16  we          use AMCplx, c2qm, use_complex_ex, qarccos
 0.06     12.11.16  we          use_complex, new qtanh, qtan
 0.07     13.11.16  we          qnorm_l1, qneg, qconj
 0.08     13.11.16  we          qcot, qcoth
 0.09     14.11.16  we          qarcsin, qarctan, qarccosh, qarcsinh, qarctanh, qarccoth
 0.10     26.11.16  we          qarccot, qcbrt
 0.11     01.12.16  we          qsec, qsech, qarccot; fix qmul/qdiv to allow @c=@a or =@b
 0.12     22.08.17  we          references, rename norm/abs etc, comments
 0.13     18.10.17  we          qarccsch, qarccsc, qarcsech
 0.14     19.10.17  we          qcsc, qsech, qmulx, qpowx
 0.15     19.10.17  we          renamed some functions, test cases
 0.16     20.10.17  we          separate unit, qunit
 0.17     21.10.17  we          qdot, renamed unit
 0.18     21.10.17  we          use qabs in qdiv
 0.19     23.10.17  we          branch cut for arcsech
 0.20     21.06.18  we          Fix some erroneous comments

***************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2016-2018 Wolfgang Ehrhardt

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


type
  Quaternion = record
                 r:     extended; {real or scalar part     }
                 x,y,z: extended; {imaginary or vector part}
               end;

{#Z+}
{---------------------------------------------------------------------------}
{-----------------------  Operations and basics  ---------------------------}
{---------------------------------------------------------------------------}
{#Z-}

function  qabs(const a: quaternion): extended;
  {-Return the magnitude of a = sqrt(r^2 + x^2 + y^2 + z^2)}

function  qabs_im(const a: quaternion): extended;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the magnitude of the imaginary part of a}

procedure qadd(const a,b: quaternion; var c: quaternion);
  {-Compute the sum c = a + b}

function  qarg(const a: quaternion): extended;
  {-Return the argument of a, arg(a) = atan2(|x,y,z|, r)}

procedure qconj(const a: quaternion; var b: quaternion);
  {-Compute the conjugate of a, b=(r,-x,-y,-z)}

procedure qdiv(const a,b: quaternion; var c: quaternion);
  {-Compute the quotient c = a * (1/b)}

function  qdot(const a,b: quaternion): extended;
 {-Compute the scalar (dot) product a.b}

procedure qinv(const a: quaternion; var b: quaternion);
  {-Compute the inverse b = 1/a}

procedure qmul(const a,b: quaternion; var c: quaternion);
  {-Compute the product c = a * b}

procedure qmulx(const a: quaternion; b: extended; var c: quaternion);
  {-Compute the product c = a * b}

procedure qneg(const a: quaternion; var b: quaternion);
  {-Compute the negative of a, b = -a}

function  qnorm(const a: quaternion): extended;
  {-Return the norm of a = r^2 + x^2 + y^2 + z^2}

procedure qsub(const a,b: quaternion; var c: quaternion);
  {-Compute the difference c = a - b}

procedure qunit(const a: quaternion; var b: quaternion);
  {-Compute the unit b = a/abs(a)}

{#Z+}
{---------------------------------------------------------------------------}
{----------------------------  Functions  ----------------------------------}
{---------------------------------------------------------------------------}
{#Z-}
procedure qarccos(const a: quaternion; var b: quaternion);
  {-Return the inverse circular cosine of b := arccos(a)}

procedure qarccosh(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic cosine b := arccosh(a)}

procedure qarccot(const a: quaternion; var b: quaternion);
  {-Return the inverse circular cotangent b := arccot(a)}

procedure qarccoth(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic cotangent b := arccoth(a)}

procedure qarccsc(const a: quaternion; var b: quaternion);
  {-Return the inverse trigonometric cosecant b = arccsc(a)}

procedure qarccsch(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic cosecant b = arccsch(a)}

procedure qarcsec(const a: quaternion; var b: quaternion);
  {-Return the inverse circular secant b = arcsec(a)}

procedure qarcsech(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic secant b = arcsech(a)}

procedure qarcsin(const a: quaternion; var b: quaternion);
  {-Return the inverse circular sine of b := arcsin(a)}

procedure qarcsinh(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic sine b := arcsinh(a)}

procedure qarctan(const a: quaternion; var b: quaternion);
  {-Return the inverse circular tangent of b := arctan(a)}

procedure qarctanh(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic tangent b := arctanh(a)}

procedure qcbrt(const a: quaternion; var b: quaternion);
  {-Compute the cube root b = cbrt(a) = a^(1/3)}

procedure qcos(const a: quaternion; var b: quaternion);
  {-Compute the circular cosine b = cos(a)}

procedure qcosh(const a: quaternion; var b: quaternion);
  {-Compute the hyperbolic cosine b = cosh(a)}

procedure qcot(const a: quaternion; var b: quaternion);
  {-Return the circular cotangent b = cot(a)}

procedure qcoth(const a: quaternion; var b: quaternion);
  {-Return the hyperbolic cotangent b = coth(a)}

procedure qcsc(const a: quaternion; var b: quaternion);
  {-Return the circular secant b = csc(a)}

procedure qcsch(const a: quaternion; var b: quaternion);
  {-Return the hyperbolic secant  b = sech(a)}

procedure qexp(const a: quaternion; var b: quaternion);
  {-Compute the exponential function b = exp(a)}

procedure qln(const a: quaternion; var b: quaternion);
  {-Compute the logarithm b = ln(a)}

procedure qpowx(const a: quaternion; b: extended; var c: quaternion);
  {-Return the real power of w = a^b = exp(b*ln(a))}

procedure qsec(const a: quaternion; var b: quaternion);
  {-Return the circular secant  b = sec(a)}

procedure qsech(const a: quaternion; var b: quaternion);
  {-Return the hyperbolic secant b := sech(a)}

procedure qsin(const a: quaternion; var b: quaternion);
  {-Compute the circular sine b = cos(a)}

procedure qsinh(const a: quaternion; var b: quaternion);
  {-Compute the hyperbolic cosine b = cosh(a)}

procedure qsqrt(const a: quaternion; var b: quaternion);
  {-Compute the square root b = sqrt(a)}

procedure qtan(const a: quaternion; var b: quaternion);
  {-Return the circular tangent of b := tan(a)}

procedure qtanh(const a: quaternion; var b: quaternion);
  {-Return the hyperbolic tangent of b := tanh(a)}


implementation

uses
  AMath, AMCmplx;


{type for complex function call}
type
  tcproc = procedure (const z: complex; var w: complex);


{Most (inverse) trigonometric/hyperbolic quaternion function qf(a,b)}
{are basically computed with the corresponding function cf(z,w) with}
{z = Re(q) + i*abs(Im(q)), w = cf(z) and the mapping to quaternions }
{b.r = Re(b) = Re(w), [b.x,b.y,b.z] = Im(b) = Im(w)*Im(q)/abs(Im(q))}
{There are some technically additions from branch buts etc.         }

{See https://en.wikipedia.org/wiki/Quaternionic_analysis or more detailed in}
{https://de.wikipedia.org/wiki/Quaternion#Fortsetzungen_komplexer_Funktionen}
{See also D. W. Harder, C++ Complex Numbers, Quaternions, Octonions,        }
{Sedenions, etc., https://ece.uwaterloo.ca/~dwharder/C++/CQOST/             }


{---------------------------------------------------------------------------}
procedure c2qm(ai: extended; usei: boolean; const u: complex; const a: quaternion; var b: quaternion);
  {-Build b from u, use u.im if ai=0 and usei=true}
var
  f: extended;
begin
  f := u.im;
  b.r := u.re;
  if ai=0.0 then begin
    if usei then begin
      b.x := f;
      b.y := 0.0;
      b.z := 0.0;
      exit;
    end
  end
  else f := u.im/ai;
  b.x := f*a.x;
  b.y := f*a.y;
  b.z := f*a.z;
end;


{---------------------------------------------------------------------------}
procedure use_complex_ex(func: tcproc; ai: extended; usei: boolean; const a: quaternion; var b: quaternion);
  {-Compute b = func(a), ai=abs(im(a)), set usei if complex imaginary part}
  { should be use if ai=0, e.g. for branch cuts etc}
var
  u,v: complex;
begin
  u.re := a.r;
  u.im := ai;
  func(u,v);
  c2qm(ai, usei, v, a, b);
end;


{---------------------------------------------------------------------------}
procedure use_complex_std(func: tcproc; const a: quaternion; var b: quaternion);
  {-Compute b = func(a) using complex function}
var
  ai: extended;
  u,v: complex;
begin
  ai := qabs_im(a);
  u.re := a.r;
  u.im := ai;
  func(u,v);
  c2qm(ai, false, v, a, b);
end;


{---------------------------------------------------------------------------}
function qnorm(const a: quaternion): extended;
  {-Return the norm of a = r^2 + x^2 + y^2 + z^2}
begin
  with a do qnorm := r*r + x*x + y*y + z*z;
end;


{---------------------------------------------------------------------------}
function qabs(const a: quaternion): extended;
  {-Return the magnitude of a = sqrt(r^2 + x^2 + y^2 + z^2)}
begin
  with a do qabs := hypot(hypot(r,x), hypot(y,z));
end;


{---------------------------------------------------------------------------}
function qabs_im(const a: quaternion): extended;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the magnitude of the imaginary part of a}
begin
  with a do qabs_im := hypot3(x,y,z);
end;


{---------------------------------------------------------------------------}
function qarg(const a: quaternion): extended;
  {-Return the argument of a, arg(a) = atan2(|x,y,z|, r)}
begin
  with a do qarg := arctan2(hypot3(x,y,z),r);
end;


{---------------------------------------------------------------------------}
procedure qadd(const a,b: quaternion; var c: quaternion);
  {-Compute the sum c = a + b}
begin
  c.r := a.r + b.r;
  c.x := a.x + b.x;
  c.y := a.y + b.y;
  c.z := a.z + b.z;
end;


{---------------------------------------------------------------------------}
procedure qsub(const a,b: quaternion; var c: quaternion);
  {-Compute the difference c = a - b}
begin
  c.r := a.r - b.r;
  c.x := a.x - b.x;
  c.y := a.y - b.y;
  c.z := a.z - b.z;
end;


{---------------------------------------------------------------------------}
procedure qconj(const a: quaternion; var b: quaternion);
  {-Compute the conjugate of a, b=(r,-x,-y,-z)}
begin
  b.r := +a.r;
  b.x := -a.x;
  b.y := -a.y;
  b.z := -a.z;
end;


{---------------------------------------------------------------------------}
procedure qneg(const a: quaternion; var b: quaternion);
  {-Compute the negative of a, b = -a}
begin
  b.r := -a.r;
  b.x := -a.x;
  b.y := -a.y;
  b.z := -a.z;
end;


{---------------------------------------------------------------------------}
procedure qmul(const a,b: quaternion; var c: quaternion);
  {-Compute the product c = a * b}
var
  t: quaternion;
begin
  {use local t to allow  @c=@a  or  @c=@b}
  t.r := a.r*b.r - a.x*b.x - a.y*b.y - a.z*b.z;
  t.x := a.r*b.x + a.x*b.r + a.y*b.z - a.z*b.y;
  t.y := a.r*b.y - a.x*b.z + a.y*b.r + a.z*b.x;
  t.z := a.r*b.z + a.x*b.y - a.y*b.x + a.z*b.r;
  c.r := t.r;
  c.x := t.x;
  c.y := t.y;
  c.z := t.z;
end;


{---------------------------------------------------------------------------}
procedure qmulx(const a: quaternion; b: extended; var c: quaternion);
  {-Compute the product c = a * b}
begin
  c.r := b*a.r;
  c.x := b*a.x;
  c.y := b*a.y;
  c.z := b*a.z;
end;


{---------------------------------------------------------------------------}
procedure qdiv(const a,b: quaternion; var c: quaternion);
  {-Compute the quotient c = a * (1/b)}
var
  t:  quaternion;
  bn: extended;
begin
  bn := qabs(b);
  {use local t to allow  @c=@a  or  @c=@b}
  t.r := ( a.r*b.r + a.x*b.x + a.y*b.y + a.z*b.z) / bn;
  t.x := (-a.r*b.x + a.x*b.r - a.y*b.z + a.z*b.y) / bn;
  t.y := (-a.r*b.y + a.x*b.z + a.y*b.r - a.z*b.x) / bn;
  t.z := (-a.r*b.z - a.x*b.y + a.y*b.x + a.z*b.r) / bn;
  c.r := t.r / bn;
  c.x := t.x / bn;
  c.y := t.y / bn;
  c.z := t.z / bn;
end;


{---------------------------------------------------------------------------}
function qdot(const a,b: quaternion): extended;
 {-Compute the scalar (dot) product a.b}
begin
  qdot := a.r*b.r + a.x*b.x + a.y*b.y + a.z*b.z
end;


{---------------------------------------------------------------------------}
procedure qinv(const a: quaternion; var b: quaternion);
  {-Compute the inverse b = 1/a}
var
  w: extended;
begin
  w   := qnorm(a);
  b.r := a.r/w;
  b.x := -a.x/w;
  b.y := -a.y/w;
  b.z := -a.z/w;
end;


{---------------------------------------------------------------------------}
procedure qexp(const a: quaternion; var b: quaternion);
  {-Compute the exponential function b = exp(a)}
var
  u,v,w: extended;
begin
  u := exp(a.r);
  v := qabs_im(a);
  w := u*sinc(v);
  b.r := u*cos(v);
  b.x := w*a.x;
  b.y := w*a.y;
  b.z := w*a.z;
end;


{---------------------------------------------------------------------------}
procedure qcos(const a: quaternion; var b: quaternion);
  {-Compute the circular cosine b = cos(a)}
var
  v,w: extended;
  s,c: extended;
begin
  sincos(a.r,s,c);
  v := qabs_im(a);
  w := -s*sinhc(v);
  b.r := c*cosh(v);
  b.x := w*a.x;
  b.y := w*a.y;
  b.z := w*a.z;
end;


{---------------------------------------------------------------------------}
procedure qcosh(const a: quaternion; var b: quaternion);
  {-Compute the hyperbolic cosine b = cosh(a)}
var
  v,w: extended;
  s,c: extended;
begin
  sinhcosh(a.r,s,c);
  v := qabs_im(a);
  w := s*sinc(v);
  b.r := c*cos(v);
  b.x := w*a.x;
  b.y := w*a.y;
  b.z := w*a.z;
end;


{---------------------------------------------------------------------------}
procedure qsin(const a: quaternion; var b: quaternion);
  {-Compute the circular sine b = cos(a)}
var
  v,w: extended;
  s,c: extended;
begin
  sincos(a.r,s,c);
  v := qabs_im(a);
  w := c*sinhc(v);
  b.r := s*cosh(v);
  b.x := w*a.x;
  b.y := w*a.y;
  b.z := w*a.z;
end;


{---------------------------------------------------------------------------}
procedure qsinh(const a: quaternion; var b: quaternion);
  {-Compute the hyperbolic cosine b = cosh(a)}
var
  v,w: extended;
  s,c: extended;
begin
  sinhcosh(a.r,s,c);
  v := qabs_im(a);
  w := c*sinc(v);
  b.r := s*cos(v);
  b.x := w*a.x;
  b.y := w*a.y;
  b.z := w*a.z;
end;


{---------------------------------------------------------------------------}
procedure qln(const a: quaternion; var b: quaternion);
  {-Compute the logarithm b = ln(a)}
var
  v,w: extended;
begin
  v := qabs_im(a);
  w := arctan2(v,a.r);
  if v=0.0 then begin
    b.r := ln(abs(a.r));
    b.x := w;
    b.y := 0;
    b.z := 0;
  end
  else begin
    w   := w/v;
    b.r := ln(hypot(a.r, v));
    b.x := w*a.x;
    b.y := w*a.y;
    b.z := w*a.z;
  end;
end;


{---------------------------------------------------------------------------}
procedure qsqrt(const a: quaternion; var b: quaternion);
  {-Compute the square root b = sqrt(a)}
var
  u,v,w: extended;
begin
  w := qabs_im(a);
  if (w=0.0) and (a.r<0.0) then begin
    {pure negative real}
    b.r := 0.0;
    b.x := sqrt(abs(a.r));
    b.y := 0.0;
    b.z := 0.0;
  end
  else begin
    cx_sqrt(a.r, w, u, v);
    if w<>0.0 then v := v/w;
    b.r := u;
    b.x := v*a.x;
    b.y := v*a.y;
    b.z := v*a.z;
  end;
end;


{---------------------------------------------------------------------------}
procedure qpowx(const a: quaternion; b: extended; var c: quaternion);
  {-Return the real power of w = a^b = exp(b*ln(a))}
begin
  qln(a,c);
  qmulx(c,b,c);
  qexp(c,c);
end;


{---------------------------------------------------------------------------}
procedure qcbrt(const a: quaternion; var b: quaternion);
  {-Compute the cube root b = cbrt(a) = a^(1/3)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}ccbrt, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarccos(const a: quaternion; var b: quaternion);
  {-Return the inverse circular cosine of b := arccos(a)}
var
  ai: extended;
  ui: boolean;
begin
  ai := qabs_im(a);
  ui := (ai=0.0) and (abs(a.r)>1.0);
  use_complex_ex({$ifdef FPC_ProcVar}@{$endif}carccos, ai, ui, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarcsin(const a: quaternion; var b: quaternion);
  {-Return the inverse circular sine of b := arcsin(a)}
var
  ai: extended;
  ui: boolean;
begin
  ai := qabs_im(a);
  ui := (ai=0.0) and (abs(a.r)>1.0);
  use_complex_ex({$ifdef FPC_ProcVar}@{$endif}carcsin, ai, ui, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarctan(const a: quaternion; var b: quaternion);
  {-Return the inverse circular tangent of b := arctan(a)}
begin
  {complex branch cut on imag.axis}
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}carctan, a, b);
end;


{---------------------------------------------------------------------------}
procedure qtanh(const a: quaternion; var b: quaternion);
  {-Return the hyperbolic tangent of b := tanh(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}ctanh, a, b);
end;


{---------------------------------------------------------------------------}
procedure qtan(const a: quaternion; var b: quaternion);
  {-Return the circular tangent of b := tan(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}ctan, a, b);
end;


{---------------------------------------------------------------------------}
procedure qcot(const a: quaternion; var b: quaternion);
  {-Return the circular cotangent b = cot(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}ccot, a, b);
end;


{---------------------------------------------------------------------------}
procedure qcoth(const a: quaternion; var b: quaternion);
  {-Return the hyperbolic cotangent b = coth(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}ccoth, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarcsec(const a: quaternion; var b: quaternion);
  {-Return the inverse circular secant b = arcsec(a)}
var
  ai: extended;
  ui: boolean;
begin
  ai := qabs_im(a);
  ui := (ai=0.0) and (abs(a.r)<1.0);
  use_complex_ex({$ifdef FPC_ProcVar}@{$endif}carcsec, ai, ui, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarccosh(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic cosine b := arccosh(a)}
var
  ai: extended;
  ui: boolean;
begin
  ai := qabs_im(a);
  ui := (ai=0.0) and (a.r<1.0);
  use_complex_ex({$ifdef FPC_ProcVar}@{$endif}carccosh, ai, ui, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarcsinh(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic sine b := arcsinh(a)}
begin
  {complex branch cut on imag.axis}
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}carcsinh, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarctanh(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic tangent b := arctanh(a)}
var
  ai: extended;
  ui: boolean;
begin
  ai := qabs_im(a);
  ui := (ai=0.0) and (abs(a.r)>1.0);
  use_complex_ex({$ifdef FPC_ProcVar}@{$endif}carctanh, ai, ui, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarccoth(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic cotangent b := arccoth(a)}
var
  ai: extended;
  ui: boolean;
begin
  ai := qabs_im(a);
  ui := (ai=0.0) and (abs(a.r)<1.0);
  use_complex_ex({$ifdef FPC_ProcVar}@{$endif}carccoth, ai, ui, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarccot(const a: quaternion; var b: quaternion);
  {-Return the inverse circular cotangent b := arccot(a)}
begin
  {complex branch cut on imag.axis}
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}carccot, a, b);
end;


{---------------------------------------------------------------------------}
procedure qsec(const a: quaternion; var b: quaternion);
  {-Return the circular secant  b = sec(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}csec, a, b);
end;


{---------------------------------------------------------------------------}
procedure qcsch(const a: quaternion; var b: quaternion);
  {-Return the hyperbolic secant  b = sech(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}ccsch, a, b);
end;


{---------------------------------------------------------------------------}
procedure qcsc(const a: quaternion; var b: quaternion);
  {-Return the circular secant b = csc(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}ccsc, a, b);
end;



{---------------------------------------------------------------------------}
procedure qsech(const a: quaternion; var b: quaternion);
  {-Return the hyperbolic secant b := sech(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}csech, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarccsch(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic cosecant b = arccsch(a)}
begin
  use_complex_std({$ifdef FPC_ProcVar}@{$endif}carccsch, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarccsc(const a: quaternion; var b: quaternion);
  {-Return the inverse trigonometric cosecant b = arccsc(a)}
var
  ai: extended;
  ui: boolean;
begin
  ai := qabs_im(a);
  ui := (ai=0.0) and (abs(a.r)<1.0);
  use_complex_ex({$ifdef FPC_ProcVar}@{$endif}carccsc, ai, ui, a, b);
end;


{---------------------------------------------------------------------------}
procedure qarcsech(const a: quaternion; var b: quaternion);
  {-Return the inverse hyperbolic secant b = arcsech(a)}
var
  ai: extended;
  ui: boolean;
begin
  ai := qabs_im(a);
  ui := (ai=0.0) and ((a.r<0.0) or (a.r>1.0));
  use_complex_ex({$ifdef FPC_ProcVar}@{$endif}carcsech, ai, ui, a, b);
end;


{---------------------------------------------------------------------------}
procedure qunit(const a: quaternion; var b: quaternion);
  {-Compute the unit b = a/abs(a)}
begin
  qmulx(a, 1.0/qabs(a), b);
end;


end.
